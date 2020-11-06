%--------------------------------------------------------------------------
%                     Aorta+systemic Compartement
%                              6-Elements 
%                             Constant R1
%                           -------------
% Copyright: These codes can be used with author premission.
% Author:  Yasser Aboelkassem
% Department of Mechanical Engineering
% San Diego State University
% Year: November/2020
%---------------------------------------------------------------------
function dy =Func_6Elem_ConstR(t,y)
%--------------------------------------------------------------------------
% Load Exp. Data (Aortic valve flowrate ``Qav" and and Aortica pressure ``pa"
%--------------------------------------------------------------------------
%array=dlmread('Qav_measured.DAT');
%t_Exp=array(:,1);
%Qav_Exp=array(:,2); % it is assumed that Qav is in ml/s
% array=dlmread('pa_measured.DAT');
% tp=array(:,1);
% pa_measured=array(:,2);  % it is assumed that pa is in mmHg
%--------------------------------
% Step 2: State variables
%--------------------------------
Qa=y(1);
G =y(2);
Ps=y(3);
%-------------------------------
% ODEs Governing Equations
%--------------------------------
global Eta_o Co C1 R1 Rp L1 t_Qav_Exp Qav_Exp Tp;
tnew=rem(t,Tp);                                        % reminder so we always between (0,1) and can find Qav
Qav=interp1(t_Qav_Exp,Qav_Exp,tnew);                      % Get Qav from Exp. Data function 

Pa=G/Co+ Eta_o*(Qav-Qa);

dQa=(1/L1)*(Pa-Ps-R1*Qa);  % use Velocity profile formula for R(t)*QSc
dG =Qav-Qa;
dPs=(Qa-Ps/Rp)/C1;
dy=[dQa;dG;dPs];
