% Program for the identification of the model parameters based on measured
% Qav and pa. Program also gives solution for pressures and flow rates in
% the time domain.
clear all; close all; clc;
tic
global Phasor_Qav Phasor_pa Phasor_pa_model freq   % phasors and frequencies

%-----------------------------------------------------
% Load Exp data: (Exp1,Kass,segers1, segers2)
%-----------------------------------------------------
%------------------------------------
% Aortic Valve Flow Rate and Pressure
%-------------------------------------
Flag=7;       % 1 = Nichols Book Adolescent , 2= Nichols Book Middle, 3 = Nichols Book Elder, 4 = Milnor, 5-Kass, 6=Segers_Pig, 7= Segeres_Dog
if (Flag==1)
array1=dlmread('QAV_Nicols_Adolescent.POD');
array2=dlmread('PAV_Nicols_Adolescent.POD');
elseif (Flag==2)
array1=dlmread('QAV_Nichols_Middle.DAT');
array2=dlmread('PAV_Nichols_Middle.DAT');    
elseif (Flag==3)
array1=dlmread('QAV_Nichols_Elderly.DAT');
array2=dlmread('PAV_Nichols_Elderly.DAT');
elseif (Flag==4)
array1=dlmread('QAV_Milnor.DAT');
array2=dlmread('PAV_Milnor.DAT');
elseif (Flag==5)
array1=dlmread('QAV_Kass.POD');
array2=dlmread('PAV_Kass.POD');
elseif (Flag==6)
array1=dlmread('QAV_Segers_Pig.POD');
array2=dlmread('PAV_Segers_Pig.POD');
PSO_param=dlmread('PSO_Optm_Parm_WK6_Segers_Pig.DAT');
elseif (Flag==7)
array1=dlmread('QAV_Segers_Dog.POD');
array2=dlmread('PAV_Segers_Dog.POD');
PSO_param=dlmread('PSO_Optm_Parm_WK6_Segers_Dog.DAT');
elseif (Flag==8) % New Data
array1=dlmread('QAV_Nichols_28.POD');
array2=dlmread('PAV_Nichols_28.POD');
array3=dlmread('Zabs_Nichols_28.POD');
array4=dlmread('Zang_Nichols_28.POD');
elseif (Flag==9) % New Data
array1=dlmread('QAV_Nichols_52.POD');
array2=dlmread('PAV_Nichols_52.POD');
array3=dlmread('Zabs_Nichols_52.POD');
array4=dlmread('Zang_Nichols_52.POD');
elseif (Flag==10) % New Data
array1=dlmread('QAV_Nichols_68.POD');
array2=dlmread('PAV_Nichols_68.POD');
array3=dlmread('Zabs_Nichols_68.POD');
array4=dlmread('Zang_Nichols_68.POD');
elseif (Flag==11) % New Data
array1=dlmread('QAV_Nichols_Normotensive.POD');
array2=dlmread('PAV_Nichols_Normotensive.POD');
array3=dlmread('Zabs_Nichols_Normotensive.POD');
array4=dlmread('Zang_Nichols_Normotensive.POD');
PSO_param=dlmread('PSO_Optm_Parm_WK6_Nichols_Normotensive.DAT');
elseif (Flag==12) % New Data
array1=dlmread('QAV_Nichols_Mild_Hypertension.POD');
array2=dlmread('PAV_Nichols_Mild_Hypertension.POD');
array3=dlmread('Zabs_Nichols_Mild_Hypertension.POD');
array4=dlmread('Zang_Nichols_Mild_Hypertension.POD');
PSO_param=dlmread('PSO_Optm_Parm_WK6_Nichols_Mild_Hypertension.DAT');
elseif (Flag==13) % New Data
array1=dlmread('QAV_Nichols_Severe_Hypertension_corr.POD');
array2=dlmread('PAV_Nichols_Severe_Hypertension_corr.POD');
array3=dlmread('Zabs_Nichols_Severe_Hypertension.POD');
array4=dlmread('Zang_Nichols_Severe_Hypertension.POD');
PSO_param=dlmread('PSO_Optm_Parm_WK6_Nichols_Severe_Hypertension.DAT');
end
tQ=array1(:,1);
Qav_measured=array1(:,2); % it is assumed that Qav is in ml/s
% figure()
% plot(tQ,Qav_measured)
%----------------
tp=array2(:,1);
pa_measured=array2(:,2);  % it is assumed that pa is in mmHg
% figure()
% plot(tp,pa_measured)
%--------------------
% NOTE: sampling time step is not necessarily constant, but periods must be equal

Nh=60; %60; % number of harmonics in the frequency domain to be used in the least square method
N_points=2*(Nh+1); % number of points within one cardiac cycle
dt=(tQ(length(tQ))-tQ(1))/(N_points-1);  % sampling time for the selected number of harmonics    
time=linspace(tQ(1), tQ(length(tQ)), N_points); % time vector with a new sampling time
Qav=interp1(tQ, Qav_measured, time);  % Qav sampled with a constant frequency
pa =interp1(tp, pa_measured , time);  % pa sampled with a constant frequency
% Find the spectrum of measured signals
[Phasor_Qav, Nh1, freq] = ffourier_def (Qav,dt);  % Fourier decomposition of Qav
[Phasor_pa , Nh1, freq] = ffourier_def (pa, dt);  % Fourier decomposition of pa

R_tot=imag(Phasor_pa(1))/imag(Phasor_Qav(1));    % total peripheral resistance at omega=0

% model parameters;
%      eta=param(1)
%      C0=param(2)
%      L=param(3)
%      r=param(4)
%      C1=param(5)
%      Rp=param(6)
% Note that Rp+r=R_tot ******************

% param:    eta       C0       L           r           C1         Rp
par_min=[    0        0        0           0           0      0.7*R_tot];  % lower bounds of model parameters
par_max=[0.3*R_tot  3/R_tot  0.1*R_tot  0.3*R_tot    3/R_tot      R_tot ];  % upper bounds of model parameters

% define initatial values of model parameters
%param0=0.5*(par_min+par_max);
param0=[0.05*R_tot 2/R_tot .003*R_tot 0.05*R_tot 0.2/R_tot 0.95*R_tot];
options=optimset('MaxFunEvals',100000,'TolX',1.e-7,'MaxIter',100000,'TolFun',1.e-7);   
[param,RMS,residual,exitflag,output] = lsqnonlin(@Goal_function,param0,par_min,par_max,options);  % Least square method
  

% Print results
fprintf('%s %7.4f %s\n', 'eta=',param(1),' mmHg*s/ml'); 
fprintf('%s %7.4f %s\n', 'C0 =',param(2),' ml/mmHg');   
fprintf('%s %7.4f %s\n', 'L  =',param(3),' mmHg*s^2/ml');  
fprintf('%s %7.4f %s\n', 'R1 =',param(4),' mmHg*s/ml'); 
fprintf('%s %7.4f %s\n', 'C1 =',param(5),' ml/mmHg');  
fprintf('%s %7.4f %s\n\n', 'Rp =',param(6),' mmHg*s/ml');
fprintf('%s %7.2f %s\n', ' ******* RMS =',sqrt(RMS),'  mmHg'); 

% Plot measured and reconstructed pressure pa
Pfit=bfourier_def1(Phasor_pa_model, time);
figure(1)
plot(tp,pa_measured,'o',time,Pfit,'-'); 
grid on
legend('measured pressure','reconstructed pressure')

% plot pressure ps and Q-s
% model parameters;
%      eta=param(1)
%      C0=param(2)
%      L=param(3)
%      R1=param(4)
%      C1=param(5)
%      Rp=param(6)
omega=2*pi*freq;  % vector of circular frequencies
for i=1:length(omega)
Z1(i)=param(6)/(1+1i*omega(i)*param(5)*param(6));  % impedance of C1 and Rp in parallel
Z2(i)=param(4)+1i*omega(i)*param(3);               % impedance of r and L in series
Z3_rec(i)=1i*omega(i)*param(2)/(1+1i*omega(i)*param(2)*param(1)); % reciprocal of impedance of eta and C0 in series
end
% solutions for pressures and flow rates in the frequency domain
Phasor_Qeta=Z3_rec.*Phasor_pa_model;
Phasor_Qa=Phasor_Qav-Phasor_Qeta;
Phasor_ps=Phasor_pa_model-Z2.*Phasor_Qa;
Phasor_QR=Phasor_ps/param(6);
Phasor_QC1=Phasor_Qa-Phasor_QR;
% solutions for pressures and flow rates in the time domain
Qeta=bfourier_def1(Phasor_Qeta, time);
Qa=bfourier_def1(Phasor_Qa, time);
ps=bfourier_def1(Phasor_ps, time);
QR=bfourier_def1(Phasor_QR, time);
QC1=bfourier_def1(Phasor_QC1, time);
% plot pressures
figure (2)
plot(time,Pfit,'k',time,ps,'r'); 
grid on
legend('pa','ps')
figure (3)
plot(time,Qav,'k',time,Qa,'r',time,Qeta,'b',time,QR,'g'); 
grid on
legend('Qav','Qa','Qeta','QR')
%--------------------------
% Export data for Tecplot
%-------------------------
if (Flag==1)
   fid1=fopen('Pa_Spect_WK6_Pa_Nicols_Adolescent.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Nicols_Adolescent.dat','wt');   
elseif (Flag==2)
   fid1=fopen('Pa_Spect_WK6_Pa_Nichols_Middle.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Nichols_Middle.dat','wt');
elseif (Flag==3)
   fid1=fopen('Pa_Spect_WK6_Pa_Nichols_Elderly.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Nichols_Elderly.dat','wt');
elseif (Flag==4)
   fid1=fopen('Pa_Spect_WK6_Pa_Milnor.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Milnor.dat','wt');
elseif (Flag==5)
   fid1=fopen('Pa_Spect_WK6_Pa_Kass.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Kass.dat','wt');
elseif (Flag==6)
   fid1=fopen('Pa_Spect_WK6_Pa_Segers_Pig.dat','wt');
   fid2=fopen('Optm_Spect_Parm_WK6_Segers_Pig.dat','wt');
   fid3=fopen('Zabs_Spect_PSO_Imped_WK6_Segers_Pig.dat','wt');
   fid4=fopen('Zang_Spect_PSO_Imped_WK6_Segers_Pig.dat','wt');
elseif (Flag==7)
   fid1=fopen('Pa_Spect_WK6_Pa_Segers_Dog.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Segers_Dog.dat','wt');
   fid3=fopen('Zabs_Spect_PSO_Imped_WK6_Segers_Dog.dat','wt');
   fid4=fopen('Zang_Spect_PSO_Imped_WK6_Segers_Dog.dat','wt');
elseif (Flag==8)
   fid1=fopen('Pa_Spect_WK6_Pa_Nichols_28.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Nichols_28.dat','wt');
elseif (Flag==9)
   fid1=fopen('Pa_Spect_WK6_Pa_Nichols_52.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Nichols_52.dat','wt');
elseif (Flag==10)
   fid1=fopen('Pa_Spect_WK6_Pa_Nichols_68.dat','wt');
   fid2=fopen('Optm_Parm_WK6_Nichols_68.dat','wt'); 
elseif (Flag==11)
   fid1=fopen('Pa_Spect_WK6_Nichols_Normotensive.dat','wt');
   fid2=fopen('Optm_Spect_Parm_WK6_Nichols_Normotensive.dat','wt');
   fid3=fopen('Zabs_Spect_PSO_Imped_WK6_Nichols_Normotensive.dat','wt');
   fid4=fopen('Zang_Spect_PSO_Imped_WK6_Nichols_Normotensive.dat','wt');
elseif (Flag==12)
   fid1=fopen('Pa_Spect_WK6_Nichols_Mild_Hypertension.dat','wt');
   fid2=fopen('Optm_Spect_Parm_WK6_Nichols_Mild_Hypertension.dat','wt');
   fid3=fopen('Zabs_Spect_PSO_Imped_WK6_Nichols_Mild_Hypertension.dat','wt');
   fid4=fopen('Zang_Spect_PSO_Imped_WK6_Nichols_Mild_Hypertension.dat','wt');
elseif (Flag==13)
   fid1=fopen('Pa_Spect_WK6_Nichols_Severe_Hypertension.dat','wt');
   fid2=fopen('Optm_Spect_Parm_WK6_Nichols_Severe_Hypertension.dat','wt');
   fid3=fopen('Zabs_Spect_PSO_Imped_WK6_Nichols_Severe_Hypertension.dat','wt');
   fid4=fopen('Zang_Spect_PSO_Imped_WK6_Nichols_Severe_Hypertension.dat','wt');
end
%------- Write Data-------------------
% Time Frictional Pressure Drop (FPD)
FPD= param(4).*Qa;
% figure (6)
% plot(time,FPD,'k'); 
% grid on
% legend('FPD')

% fprintf(fid1,'Variables="t","Pa","Ps","Qav","Qa","Qeta","QR","FPD"');
%    for i=1:length(time)
%        fprintf(fid1,'%f %f %f %f %f %f %f %f\n',time(i),Pfit(i),ps(i),Qav(i),Qa(i),Qeta(i),QR(i),FPD(i));
%    end
% fclose(fid1);

fprintf(fid1,'Variables="t","Pa"');
   for i=1:length(time)
       fprintf(fid1,'%f %f\n',time(i),Pfit(i));
   end
fclose(fid1);
%------- Write Optimized Parameters--------------
fprintf(fid2,'%f %f %f %f %f %f %f\n',param(1),param(2),param(3),param(4),param(5),param(6),sqrt(RMS));
fclose(fid2);
toc
exitflag

%-------------------------------------------
% Find The numerical Total Input Impedence
%-------------------------------------------
Za=Z2;      % impedence of L and R1
Zs=Z1;      % impedence of Rs and C1
Z_Num= (Za+Zs)./(1+(Za+Zs).*Z3_rec); % total numerical input impedence
Z_Magn_Num = abs(Z_Num);    %magnitude
Z_Phase_Num = angle(Z_Num)*180/pi; %phase angle
%----------------------------
%-------------------------------------------
% Find The Experimental Total Input Impedence
%-------------------------------------------
Z_Exp=Phasor_pa./Phasor_Qav;
Z_Magn_Exp = abs(Z_Exp);    %magnitude
Z_Phase_Exp = angle(Z_Exp)*180/pi; %phase angle
%---------------------------------
% Adjust Data as a column vector
%--------------------------------
%- Numerical
Z_Magn_Num = Z_Magn_Num';   %magnitude
Z_Phase_Num = Z_Phase_Num'; %phase angle
%- Exp
Z_Magn_Exp=Z_Magn_Exp';
Z_Phase_Exp=Z_Phase_Exp';
freq=freq';
%-----
time=time';
FPD=FPD';
%--------------------------------------------------------------------------
%--------------------------------------------------
% 
%  Find Impedence based on PSO-Optimized Parameters
%
%--------------------------------------------------
for i=1:length(omega)
Z1_PSO(i)=PSO_param(6)/(1+1i*omega(i)*PSO_param(5)*PSO_param(6));  % impedance of C1 and Rp in parallel
Z2_PSO(i)=PSO_param(4)+1i*omega(i)*PSO_param(3);               % impedance of r and L in series
Z3_PSO_rec(i)=1i*omega(i)*PSO_param(2)/(1+1i*omega(i)*PSO_param(2)*PSO_param(1)); % reciprocal of impedance of eta and C0 in series
end

%-------------------------------------------
% Find The numerical Total Input Impedence
%-------------------------------------------
Za_PSO=Z2_PSO;      % impedence of L and R1
Zs_PSO=Z1_PSO;      % impedence of Rs and C1
Z_Num_PSO= (Za_PSO+Zs_PSO)./(1+(Za_PSO+Zs_PSO).*Z3_PSO_rec); % total numerical input impedence
Z_Magn_Num_PSO = abs(Z_Num_PSO);    %magnitude
Z_Phase_Num_PSO = angle(Z_Num_PSO)*180/pi; %phase angle

Z_Magn_Num_PSO = Z_Magn_Num_PSO';   %magnitude
Z_Phase_Num_PSO = Z_Phase_Num_PSO'; %phase angle


%---------------------------------------------
% Plot Magnitude and Phase of Input Impedence: Exp, Numerical and Scanned
%---------------------------------------------
if (Flag==8)||(Flag==9)||(Flag==10) ||(Flag==11) ||(Flag==12) ||(Flag==13)
%----------------------------------------------------
% Obtain measured and exact (scanned) Impedance-data
%----------------------------------------------------
freq_measured=array3(:,1);
Zabs_measured_scanned=array3(:,2)./1330; % Magnitude: we divided by 1330 to convert from (dyn*s/cm3) to (mmHg*s/ml)
Zang_measured_scanned=array4(:,2);       % Phase: degrees    
%-------------------------------------------
figure()
plot(freq(1:10),Z_Magn_Exp(1:10),'o',freq(1:10),Z_Magn_Num(1:10),'r',freq_measured,Zabs_measured_scanned,'-*'   , freq(1:10),Z_Magn_Num_PSO(1:10),'g')
figure()
plot(freq(1:10),Z_Phase_Exp(1:10),'o',freq(1:10),Z_Phase_Num(1:10),'r',freq_measured,Zang_measured_scanned,'-*' ,freq(1:10),Z_Phase_Num_PSO(1:10),'g')
else
figure()
plot(freq(1:10),Z_Magn_Exp(1:10),'o',freq(1:10),Z_Magn_Num(1:10),'r',freq(1:10),Z_Magn_Num_PSO(1:10),'g')
figure()
plot(freq(1:10),Z_Phase_Exp(1:10),'o',freq(1:10),Z_Phase_Num(1:10),'r',freq(1:10),Z_Phase_Num_PSO(1:10),'g')    
end
%-------------------------------------
% Write Impedence mangitude data
%-------------------------------------
fprintf(fid3,'Variables="freq","Z_Abs_Exp","Z_Abs_LSFIT","Z_Abs_PSO"');
   for i=1:length(omega)
       freq=i;
       fprintf(fid3,'%f %f %f %f\n',freq, Z_Magn_Exp(i), Z_Magn_Num(i), Z_Magn_Num_PSO(i));
   end
fclose(fid3);

%-------------------------------------
% Write Impedence phase data
%-------------------------------------
fprintf(fid4,'Variables="freq","Z_phse_Exp","Z_phse_LSFIT","Z_phse_PSO"');
   for i=1:length(omega)
       freq=i;
       fprintf(fid4,'%f %f %f %f\n',freq, Z_Phase_Exp(i), Z_Phase_Num(i), Z_Phase_Num_PSO(i));
   end
fclose(fid4);








