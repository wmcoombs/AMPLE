%AMPLE 1.1: A Material Point Learning Environment
%--------------------------------------------------------------------------
% Author: William Coombs
% Date:   27/08/2020
% Description:
% Large deformation elasto-plastic (EP) material point method (MPM) code
% based on an updated Lagrangian (UL) descripition of motion with a 
% quadrilateral background mesh. 
%
%--------------------------------------------------------------------------
% See also:
% SETUPGRID             - analysis specific information
% ELEMMPINFO            - material point-element information
% DETEXTFORCE           - external forces
% DETFDOFS              - mesh unknown degrees of freedom
% LINSOLVE              - linear solver
% DETMPS                - material point stiffness and internal force
% UPDATEMPS             - update material points
% POSTPRO               - post processing function including vtk output
%--------------------------------------------------------------------------
clear;
addpath('constitutive','functions','plotting','setup','splitting');        
split_count = 1;
ghost_type = 0;%0 - no ghost; 1 - surface ghost; 2 - bulk ghost
crit = 0.5;
crit_stop = 0.6;
refine = true;
refine_count = 0;
max_refines = 10;
refine_mats = {};
[lstps,g,mpData,mesh] = setupGrid_collapse_ghost(1); 
plot_data_disp = zeros(max_refines,lstps);
plot_data_load = zeros(max_refines,lstps);
plot_data_error = zeros(max_refines,lstps);
[nodes,nD] = size(mesh.coord);                                              % number of nodes and dimensions
[nels,nen] = size(mesh.etpl);                                               % number of elements and nodes/element
nDoF = nodes*nD;   
cache_frct = zeros(lstps,nDoF);
cache_uvw  = zeros(lstps,nDoF);
plot_lstp = zeros(max_refines,1);
last_lstp = 1;
delete('cached_nodes\*.csv')
delete('outframes\*.png')
mkdir cached_nodes;
mkdir output;
mkdir load_disp_refine\;
mkdir outframes\;
total_timer = tic;
close all;
framecount = 1;
figure(1)
pbaspect([2 1 1])
while refine && refine_count < max_refines
refine_count = refine_count + 1;
stats_total_nr_steps = 0;
stats_total_nr_cache_steps = 0;
stats_total_lstp_cached = 0;
stats_total_lstp_uncached = 0;
stats_total_cache_updates = 0;
[lstps,g,mpData,mesh] = setupGrid_collapse_ghost(1);                                          % setup information
for i=1:split_count
    %mpData = split_mps(mesh,mpData,2*ones(length(mpData),1));
    mpData = split_mps(mesh,mpData,3*ones(length(mpData),1));
    %mpData = split_mps(mesh,mpData,to_split);
end
disp("Refining MPS");
for i=1:(refine_count-1)
    %mpData = split_mps(mesh,mpData,2*refine_mats{i});
    mpData = split_mps(mesh,mpData,refine_mats{i});
end
NRitMax = 50; tol = 1e-9;                                                   % Newton Raphson parameters
use_cache = true;
if refine_count == 1
    use_cache = false;
end
[nodes,nD] = size(mesh.coord);                                              % number of nodes and dimensions
[nels,nen] = size(mesh.etpl);                                               % number of elements and nodes/element
nDoF = nodes*nD;                                                            % total number of degrees of freedom
nmp  = length(mpData);                                                      % number of material points
lstp = 0;                                                                   % zero loadstep counter (for plotting function)
uvw  = zeros(nDoF,1);                                                       % zeros displacements (for plotting function)
run postPro;                                                                % plotting initial state & mesh
disp(nmp)
%pool = 1;
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool("threads");
end
%mpData_a = mpData
%mpData_b = mpData
data_load = zeros(lstps,1);
data_disp = zeros(lstps,1);
data_error = zeros(lstps,1);
data_split = zeros(lstps,nmp);
tic;
try
for lstp=1:lstps                                                            % loadstep loop
  fprintf(1,'\n%s %4i %s %4i\n','loadstep ',lstp,' of ',lstps);             % text output to screen (loadstep)
  [mesh,mpData] = elemMPinfo(mesh,mpData);                                  % material point - element information
  Ks   = mesh.ks*ghostPenalty(mesh);                                        % Ghost stabilisation matrix
  if ~enable_ghost
    Ks = Ks * 0;
  end

  fext = detExtForce(nodes,nD,g,mpData);                                    % external force calculation (total)
  fext = fext*lstp/lstps;                                                   % current external force value
  oobf = fext;                                                              % initial out-of-balance force
  fErr = 1;                                                                 % initial error
  frct = zeros(nDoF,1);                                                     % zero the reaction forces
  uvw  = zeros(nDoF,1);                                                     % zero the displacements
  fd   = detFDoFs(mesh);                                                    % free degrees of freedom
  NRit = 0;                                                                 % zero the iteration counter
  Kt   = 0;                                                                 % zero global stiffness matrix
    has_used_cache = false;
  while (fErr > tol) && (NRit < NRitMax) || (NRit < 2)                      % global equilibrium loop
    if refine_count > 1 && use_cache && NRit == 0 && lstp < last_lstp
      has_used_cache = true;
      %uvw_cache = readmatrix(sprintf("cached_nodes/uvw_%d.csv",lstp));
      %frct_cache = readmatrix(sprintf("cached_nodes/frct_%d.csv",lstp));

      uvw_cache=cache_uvw(lstp,:)';
      frct_cache=cache_frct(lstp,:)';
      %Kt = readmatrix(sprintf("cached_nodes/kt_%d.csv",lstp));
      uvw = uvw_cache;
      frct = frct_cache;
      fgst = Ks*uvw;
      [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
      oobf = (fext-fint+frct-fgst);                                               % out-of-balance force
      fErr = norm(oobf)/norm(fext+frct+eps); 
      fprintf(1,'Cached error norm %8.3e\n',fErr);                                 % normalised oobf error
      if (fErr < 1e-10)
        fprintf(1,'Cached value acceped\n');
        fprintf(1,'Cached error norm %8.3e\n',fErr);
        break;
      end
      %fprintf(1,'Cached uvw norm %8.3e\n',norm(duvw));   % text output to screen (NR error)#
      [duvw,drct] = linSolve(mesh.bc,Kt+Ks,oobf,NRit,fd);                        % linear solver
      uvw  = uvw+duvw;                                                        % update displacements
      frct = frct+drct;                                                       % update reaction forces
      fgst = Ks*uvw;                                                          % ghost stabilisation force
      [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
      oobf = (fext-fint+frct-fgst);                                           % out-of-balance force
      fErr = norm(oobf)/norm(fext+frct+eps);                                  % normalised oobf error
      NRit = NRit+1;                                                          % increment the NR counter
      %fprintf(1,'%s %2i %s %8.3e\n','  iteration ',NRit,' NR error ',fErr);   % text output to screen (NR error)
      fprintf(1,'Single step error norm %8.3e\n',fErr);   % text output to screen (NR error)#
      
      if (fErr > 1e-1)
        fprintf(1,'Cached value rejected\n');
        use_cache = false;
      end
    else
    [duvw,drct] = linSolve(mesh.bc,Kt+Ks,oobf,NRit,fd);                        % linear solver
    uvw  = uvw+duvw;                                                        % update displacements
    frct = frct+drct;                                                       % update reaction forces
    fgst = Ks*uvw;                                                          % ghost stabilisation force
    [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
    oobf = (fext-fint+frct-fgst);                                           % out-of-balance force
    fErr = norm(oobf)/norm(fext+frct+eps);                                  % normalised oobf error
    NRit = NRit+1;                                                          % increment the NR counter
    fprintf(1,'%s %2i %s %8.3e\n','  iteration ',NRit,' NR error ',fErr);   % text output to screen (NR error)
    end
  end
  if has_used_cache
    stats_total_nr_cache_steps = stats_total_nr_cache_steps + NRit;
    stats_total_lstp_cached = stats_total_lstp_cached + 1;
    %fprintf(1,'Duvw norm %8.3e\n',norm(uvw_cache-uvw));   % text output to screen (NR error)
  else
    stats_total_nr_steps = stats_total_nr_steps + NRit;
    stats_total_lstp_uncached = stats_total_lstp_uncached + 1;
  end

  %if NRit > 2 && use_cache
  %    fprintf(1,'Given up cached values\n');
  %    use_cache = false;
  %end
  error = 0;
  if has_used_cache
      x = mesh.coord(:,1);
      y = mesh.coord(:,2);
      derr = (uvw_cache - uvw);
      %derr = derr/max(derr);
      scale = 100;
      dx = derr(1:nD:end);
      dy = derr(2:nD:end);
      error = sqrt(sum(dx.*dx + dy.*dy));
  end
  pos = [mpData.mpC];
  lp = [mpData.lp];
  pos_x = pos(1:2:end);
  pos_y = pos(2:2:end);
  lp_x = lp(1:2:end)*2;
  lp_y = lp(2:2:end)*2;
  %pos_y = pos_y(pos_x==pos_x(1));
  %pos_x = pos_x(pos_x==pos_x(1));
  %data_split_m = split_critera(mesh,mpData,crit*0.25);
  data_split_l = split_critera(mesh,mpData,crit);
  data_split_l(data_split_l==1) = 1;
  data_split_l(data_split_l==2) = 1;
  data_split_l(data_split_l==3) = 1;
  data_split = data_split_l;
  %data_split = (0.5*(data_split_m.*(1-data_split_l))) + data_split_l;
  data_load(lstp) = lstp;
  data_disp(lstp) = max(pos_x+(0.5*lp_x));
  data_error(lstp) = error;
  plot_data_disp(refine_count,lstp) = data_disp(lstp);
  plot_data_load(refine_count,lstp) = data_load(lstp);
  plot_data_error(refine_count,lstp) = error;
  plot_lstp(refine_count) = lstp;
  %scatter(pos_x,pos_y,[],data_split);#
  if true
      if true
      colours = zeros(nmp,3);
      positions = [(pos_x-lp_x*0.5)', (pos_y-lp_y*0.5)', lp_x',lp_y'];
      colours(:,1) = data_split;
      colours(:,2) = data_split==0;
      cla;
      for i=1:nmp
        rectangle('Position', positions(i,:), 'FaceColor', colours(i,:));
      end
      colormap(colours);
      xlim([0,mesh.h(1)]);
      caxis([0,1]);
      xlim([0,32]);
      ylim([0,16]);
      %ylim([0, 50]);
    
      error = 0;
      if has_used_cache
          x = mesh.coord(:,1);
          y = mesh.coord(:,2);
          derr = (uvw_cache - uvw);
          %derr = derr/max(derr);
          scale = 100;
          dx = derr(1:nD:end);
          dy = derr(2:nD:end);
          error = sqrt(sum(dx.*dx + dy.*dy));
          fprintf(1,'Cache error norm: %8.3e\n',error);
    
          hold on;
          quiver(x,y,dx*scale,dy*scale,'AutoScale','off');
      end
      else
      %ylim([0, max(pos_y)* 1.1]);%figure(2)
      figure(2);
      title("Load-displacement graph")
      cla;
      hold on;
      for i = 1:(refine_count-1)
        plot(plot_data_load(i,1:plot_lstp(i)),plot_data_error(i,1:plot_lstp(i)));
      end
      plot(data_load(1:lstp),data_error(1:lstp));
      end
      drawnow;
      pause(0.5)
  end

  %if NRit > 2 && use_cache
  %    fprintf(1,'Given up cached values\n');
  %    use_cache = false;
  %end
  if ~use_cache || NRit > 3
    fprintf(1,"Cache updated\n");
    stats_total_cache_updates = stats_total_cache_updates+1;
    cache_uvw(lstp,:) = uvw;
    cache_frct(lstp,:) = frct;
    %writematrix(uvw,sprintf("cached_nodes/uvw_%d.csv",lstp))
    %writematrix(frct,sprintf("cached_nodes/frct_%d.csv",lstp))
  end
  %writematrix(duvw,sprintf("cached_nodes/duvw_%d.csv",lstp))
  %writematrix(drct,sprintf("cached_nodes/drct_%d.csv",lstp))
  mpData = updateMPs(uvw,mpData,mesh);                                           % update material points
  too_long_crit = split_critera(mesh,mpData,crit_stop);
  if(any(too_long_crit))
    max([mpData.lp])/mesh.h(2)
    throw(MException("AMPLE:domain_size","MP domain too long"));
  end
  drawnow;
  title(sprintf("Refine count: %d - Loadstep: %d - MP count: %d - Cache error: %8.3e",refine_count,lstp,nmp,error))
  %saveas(gcf,(sprintf("outframes/frame_05%d.png",framecount)));
  print(gcf,sprintf("outframes/frame_%05d.png",framecount),'-dpng','-r300');
  framecount =  framecount + 1;
  run postPro;                                                              % Plotting and post processing 
  refine = false;
end    
catch e %e is an MException struct
  fprintf(1,'The identifier was:\n%s',e.identifier);
  fprintf(1,'There was an error! The message was:\n%s',e.message);
  refine = true;
  new_split = split_critera(mesh,mpData,crit);
  refine_mats{refine_count} = new_split;  
  last_lstp = lstp;
  if all(new_split == 0)
      refine = false;
      fprintf(1,'Simulation stopped, but due to a non-length scale issue\n');
      %throw(e);
  end
end
run postPro; 
time = toc;
fprintf(1,'\n');
fprintf(1,'Average throughput: %8.3e\n',lstp/time);
fprintf(1,'Load steps: %i - cached : %i\n',lstp,stats_total_lstp_cached);
fprintf(1,'Total nrs - %i\n',stats_total_nr_steps + stats_total_nr_cache_steps);
fprintf(1,'Average nrs cached - %8.3e \n',stats_total_nr_cache_steps/stats_total_lstp_cached);
fprintf(1,'Average nrs uncached - %8.3e \n',stats_total_nr_steps/stats_total_lstp_uncached);
fprintf(1,'Cache updates - %i \n',stats_total_cache_updates);
data_load = data_load(1:lstp-1);
data_disp = data_disp(1:lstp-1);

writematrix([data_load,data_disp],sprintf("load_disp_refine/load_disp_%d.csv",refine_count));
%figure(2)
%title("Load-displacement graph")
%plot(data_load,data_disp);
end
total_time = toc(total_timer);
fprintf(1,'Average effective throughput: %8.3e\n',20/total_time);