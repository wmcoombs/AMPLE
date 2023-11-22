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
for split_count = 9:9
[lstps,g,mpData,mesh] = setupGrid(1);                                          % setup information
%mpData = split_mps(mesh,mpData,2*ones(length(mpData),1));
%mpData = split_mps(mesh,mpData,2*ones(length(mpData),1));
for i=1:split_count
    %to_split = zeros(length(mpData),1);
    %nmp = length(mpData);
    %for mp=1:nmp 
    %    if mpData(mp).mpC(2) < 1
    %        to_split(mp) = 2;
    %    end
    %end
    mpData = split_mps(mesh,mpData,2*ones(length(mpData),1));
    %mpData = split_mps(mesh,mpData,to_split);
end
NRitMax = 50; tol = 1e-9;                                                   % Newton Raphson parameters
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
figure()
pbaspect([1 1 1])
%try
for lstp=1:lstps                                                            % loadstep loop
  fprintf(1,'\n%s %4i %s %4i\n','loadstep ',lstp,' of ',lstps);             % text output to screen (loadstep)
  [mesh,mpData] = elemMPinfo(mesh,mpData);                                  % material point - element information
  fext = detExtForce(nodes,nD,g,mpData);                                    % external force calculation (total)
  
  fext = fext*lstp/lstps;                                                   % current external force value
  %fext = fext*1/(1+log(lstps/lstp));                                                   % current external force value

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
  crit = 1;
  data_split_m = split_critera(mesh,mpData,crit*0.25);
  data_split_l = split_critera(mesh,mpData,crit*0.5);
  data_split = (0.5*(data_split_m.*(1-data_split_l))) + data_split_l;
  data_load(lstp) = lstp;
  data_disp(lstp) = max(pos_y);
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
  %drawnow;
  %pause(0.1)
  %ylim([0, max(pos_y)* 1.1])
  title("Load-displacement graph")
  cla;
  hold on;
  plot(data_load(1:lstp),data_disp(1:lstp));
  while (fErr > tol) && (NRit < NRitMax) || (NRit < 2)                      % global equilibrium loop
    [duvw,drct] = linSolve(mesh.bc,Kt,oobf,NRit,fd);                        % linear solver
    uvw  = uvw+duvw;                                                        % update displacements
    frct = frct+drct;                                                       % update reaction forces
    [fint,Kt,mpData] = detMPs(uvw,mpData);                                  % global stiffness & internal force
    oobf = (fext-fint+frct);                                                % out-of-balance force
    fErr = norm(oobf)/norm(fext+frct+eps);                                  % normalised oobf error
    NRit = NRit+1;                                                          % increment the NR counter
    fprintf(1,'%s %2i %s %8.3e\n','  iteration ',NRit,' NR error ',fErr);   % text output to screen (NR error)
  end                           
  mpData = updateMPs(uvw,mpData,mesh);                                           % update material points
  too_long_crit = split_critera(mesh,mpData,1);
  if(any(too_long_crit))
    max([mpData.lp])/mesh.h(2)
    throw(MException("AMPLE:domain_size","MP domain too long\n"));
  end
  drawnow;
  run postPro;                                                              % Plotting and post processing 
end    
%catch e %e is an MException struct
%fprintf(1,'The identifier was:\n%s',e.identifier);
%fprintf(1,'There was an error! The message was:\n%s',e.message);
%end
data_load = data_load(1:lstp-1);
data_disp = data_disp(1:lstp-1);
writematrix([data_load,data_disp],sprintf("split_output/load_disp_%d.csv",split_count));
figure()
title("Load-displacement graph")
plot(data_load,data_disp);
end