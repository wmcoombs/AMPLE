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

refine = true;
refine_count = 1;
max_refines = 6;
refine_mats = {};
[lstps,g,mpData,mesh] = setupGrid(1); 
plot_data_disp = zeros(max_refines,lstps);
plot_data_load = zeros(max_refines,lstps);
plot_lstp = zeros(max_refines,1);
last_lstp = 1;
while refine && refine_count < max_refines
[lstps,g,mpData,mesh] = setupGrid(1);                                          % setup information
for i=1:split_count
    mpData = split_mps(mesh,mpData,2*ones(length(mpData),1));
    %mpData = split_mps(mesh,mpData,to_split);
end
disp("Refining MPS");
for i=1:(refine_count-1)
    mpData = split_mps(mesh,mpData,2*refine_mats{i});
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
%parpool("threads")
%mpData_a = mpData
%mpData_b = mpData
data_load = zeros(lstps,1);
data_disp = zeros(lstps,1);
data_split = zeros(lstps,nmp);
close all;
figure(1)
pbaspect([1 1 1])
tic;
try
for lstp=1:lstps                                                            % loadstep loop
  fprintf(1,'\n%s %4i %s %4i\n','loadstep ',lstp,' of ',lstps);             % text output to screen (loadstep)
  [mesh,mpData] = elemMPinfo(mesh,mpData);                                  % material point - element information
  fext = detExtForce(nodes,nD,g,mpData);                                    % external force calculation (total)
  fext = fext*lstp/lstps;                                                   % current external force value
  oobf = fext;                                                              % initial out-of-balance force
  fErr = 1;                                                                 % initial error
  frct = zeros(nDoF,1);                                                     % zero the reaction forces
  uvw  = zeros(nDoF,1);                                                     % zero the displacements
  fd   = detFDoFs(mesh);                                                    % free degrees of freedom
  NRit = 0;                                                                 % zero the iteration counter
  Kt   = 0;                                                                 % zero global stiffness matrix
pos = [mpData.mpC];
lp = [mpData.lp];
  pos_x = pos(1:2:end);
  pos_y = pos(2:2:end);
  lp_x = lp(1:2:end)*2;
  lp_y = lp(2:2:end)*2;
  %pos_y = pos_y(pos_x==pos_x(1));
  %pos_x = pos_x(pos_x==pos_x(1));
  crit = 0.9;
  %data_split_m = split_critera(mesh,mpData,crit*0.25);
  data_split_l = split_critera(mesh,mpData,crit*0.9);
  data_split = data_split_l;
  %data_split = (0.5*(data_split_m.*(1-data_split_l))) + data_split_l;
  data_load(lstp) = lstp;
  data_disp(lstp) = max(pos_y);
  plot_data_disp(refine_count,lstp) = data_disp(lstp);
  plot_data_load(refine_count,lstp) = data_load(lstp);
  plot_lstp(refine_count) = lstp;
  %scatter(pos_x,pos_y,[],data_split);
  %colours = zeros(nmp,3);
  %positions = [(pos_x-lp_x*0.5)', (pos_y-lp_y*0.5)', lp_x',lp_y'];
  %colours(:,1) = data_split;
  %colours(:,2) = data_split==0;
  %cla;
  %for i=1:nmp
  %  rectangle('Position', positions(i,:), 'FaceColor', colours(i,:));
  %end
  %colormap(colours);
  %xlim([0,mesh.h(1)]);
  %caxis([0,1]);
  %ylim([0, 50]);
  %drawnow;
  %pause(0.01)
  %ylim([0, max(pos_y)* 1.1]);%figure(2)
  title("Load-displacement graph")
  cla;
  hold on;
  for i = 1:(refine_count-1)
    plot(plot_data_load(i,1:plot_lstp(i)),plot_data_disp(i,1:plot_lstp(i)));
  end
  plot(data_load(1:lstp),data_disp(1:lstp));
  
  while (fErr > tol) && (NRit < NRitMax) || (NRit < 2)                      % global equilibrium loop
    if refine_count > 1 && use_cache && NRit == 0 && lstp < last_lstp
      uvw_cache = readmatrix(sprintf("cached_nodes/uvw_%d.csv",lstp));
      frct_cache = readmatrix(sprintf("cached_nodes/frct_%d.csv",lstp));
      uvw = uvw_cache;
      frct = frct_cache;
      [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
      oobf = (fext-fint+frct);                                                % out-of-balance force
      fErr = norm(oobf)/norm(fext+frct+eps);                                  % normalised oobf error
      if (fErr < 1e-10)
        fprintf(1,'Cached value acceped\n');
        fprintf(1,'Cached error norm %8.3e\n',fErr);
        break;
      end
      %fprintf(1,'Cached uvw norm %8.3e\n',norm(duvw));   % text output to screen (NR error)#

      [duvw,drct] = linSolve(mesh.bc,Kt,oobf,NRit,fd);
      uvw  = uvw+duvw;                                                        % update displacements
      frct = frct+drct;                                                       % update reaction forces
      [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
      oobf = (fext-fint+frct);                                                % out-of-balance force
      fErr = norm(oobf)/norm(fext+frct+eps);                                  % normalised oobf error
      NRit = NRit+1;                                                          % increment the NR counter
      fprintf(1,'Single step error norm %8.3e\n',fErr);   % text output to screen (NR error)#
      
      if (fErr > 1e-3)
        fprintf(1,'Cached value rejected\n');
        use_cache = false;
      end
    else
    [duvw,drct] = linSolve(mesh.bc,Kt,oobf,NRit,fd);                        % linear solver
    %end
    uvw  = uvw+duvw;                                                        % update displacements
    frct = frct+drct;                                                       % update reaction forces
    [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
    oobf = (fext-fint+frct);                                                % out-of-balance force
    fErr = norm(oobf)/norm(fext+frct+eps);                                  % normalised oobf error
    NRit = NRit+1;                                                          % increment the NR counter
    fprintf(1,'%s %2i %s %8.3e\n','  iteration ',NRit,' NR error ',fErr);   % text output to screen (NR error)
    end
  end
  if use_cache
    %fprintf(1,'Duvw norm %8.3e\n',norm(uvw_cache-uvw));   % text output to screen (NR error)
  end
  %if NRit > 2 && use_cache
  %    fprintf(1,'Given up cached values\n');
  %    use_cache = false;
  %end
  if ~use_cache || NRit > 3
  writematrix(uvw,sprintf("cached_nodes/uvw_%d.csv",lstp))
  writematrix(frct,sprintf("cached_nodes/frct_%d.csv",lstp))
  end
  %writematrix(duvw,sprintf("cached_nodes/duvw_%d.csv",lstp))
  %writematrix(drct,sprintf("cached_nodes/drct_%d.csv",lstp))
  mpData = updateMPs(uvw,mpData,mesh);                                           % update material points
  too_long_crit = split_critera(mesh,mpData,1);
  if(any(too_long_crit))
    max([mpData.lp])/mesh.h(2)
    throw(MException("AMPLE:domain_size","MP domain too long"));
  end
  drawnow;
  run postPro;                                                              % Plotting and post processing 
  refine = false;
end    
catch e %e is an MException struct
  fprintf(1,'The identifier was:\n%s',e.identifier);
  fprintf(1,'There was an error! The message was:\n%s',e.message);
  refine = true;
  new_split = split_critera(mesh,mpData,0.6);
  refine_mats{refine_count} = new_split;  
  refine_count = refine_count + 1;
  last_lstp = lstp;
  if all(new_split == 0)
      refine = false;
      fprintf(1,'Simulation stopped, but due to a non-length scale issue\n');
  end
end
time = toc;
fprintf(1,'Average throughput: %8.3e\n',time/lstp);
data_load = data_load(1:lstp-1);
data_disp = data_disp(1:lstp-1);

writematrix([data_load,data_disp],sprintf("load_disp_refine/load_disp_%d.csv",refine_count));
%figure(2)
%title("Load-displacement graph")
%plot(data_load,data_disp);
end
