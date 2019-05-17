# AMPLE
AMPLE - A Material Point Learning Environment

AMPLE is a quasi-static implicit implementation of the material point method in MATLAB.  
More informatio about AMPLE can be obtained from the project webapges:
https://wmcoombs.github.io/

AMPLE is an elasto-plastic large deformation material point code with a regular quadrilateral background mesh 
(the main code is ample.m).   The continuum framework is based on an updated Lagrangian formation and two 
different constitutive models are included: linear elasticity and a linear elastic-perfectly plastic model 
with a von Mises yield surface.  
 
It is suggested that users start by looking at the “index.html” file within the “documentation” folder.  
The html files were complied with M2HTML (https://www.artefact.tk/software/matlab/m2html/) and they should 
provide  an overview of the code and how it all fits together.  The scripts should also be compatible with 
MATLAB’s “help” function.  For example if you type
 
>> help ample
 
MATLAB should return:
 
 ample: A Material Point Learning Environment
 --------------------------------------------------------------------------
  Author: William Coombs
  Date:   29/01/2019
  Description:
  Large deformation elasto-plastic (EP) material point method (MPM) code
  based on an updated Lagrangian (UL) descripition of motion with a 
  quadrilateral background mesh. 
 
 --------------------------------------------------------------------------
  See also:
  setupGrid             - analysis specific information
  elemMPinfo            - material point-element information
  detExtForce           - external forces
  detFDoFs              - mesh unknown degrees of freedom
  linSolve              - linear solver
  detMPs                - material point stiffness and internal force
  updateMPs             - update material points
  postPro               - post processing function including vtk output
 -------------------------------------------------------------------------- 
 
The “See also:” part should contain hyperlinks to the script’s sub functions (if applicable). 
