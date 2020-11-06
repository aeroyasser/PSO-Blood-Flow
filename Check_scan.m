% Program for the identification of the model parameters based on measured
% Qav and pa. Program also gives solution for pressures and flow rates in
% the time domain.
clear all; close all; clc;

%-----------------------------------------------------
% Load Exp data: (Exp1,Kass,segers1, segers2)
%-----------------------------------------------------
%------------------------------------
% Aortic Valve Flow Rate and Pressure
%-------------------------------------
Flag=13;       % 1 = Nichols Book Adolescent , 2= Nichols Book Middle, 3 = Nichols Book Elder, 4 = Milnor, 5-Kass, 6=Segers_Pig, 7= Segeres_Dog
if (Flag==8) % New Data
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
elseif (Flag==12) % New Data
array1=dlmread('QAV_Nichols_Mild_Hypertension.POD');
array2=dlmread('PAV_Nichols_Mild_Hypertension.POD');
array3=dlmread('Zabs_Nichols_Mild_Hypertension.POD');
array4=dlmread('Zang_Nichols_Mild_Hypertension.POD');
elseif (Flag==13) % New Data
array1=dlmread('QAV_Nichols_Severe_Hypertension.POD');
array2=dlmread('PAV_Nichols_Severe_Hypertension.POD');
array3=dlmread('Zabs_Nichols_Severe_Hypertension.POD');
array4=dlmread('Zang_Nichols_Severe_Hypertension.POD');
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

Nh=60;  % number of harmonics in the frequency domain to be used in the least square method
N_points=2*(Nh+1); % number of points within one cardiac cycle
dt=(tQ(length(tQ))-tQ(1))/(N_points-1);  % sampling time for the selected number of harmonics    
time=linspace(tQ(1), tQ(length(tQ)), N_points); % time vector with a new sampling time
Qav=interp1(tQ, Qav_measured, time);  % Qav sampled with a constant frequency
pa =interp1(tp, pa_measured , time);  % pa sampled with a constant frequency
% Find the spectrum of measured signals
[Phasor_Qav, Nh1, freq] = ffourier_def (Qav,dt);  % Fourier decomposition of Qav
[Phasor_pa , Nh1, freq] = ffourier_def (pa, dt);  % Fourier decomposition of pa

R_tot=imag(Phasor_pa(1))/imag(Phasor_Qav(1))*1330    % total peripheral resistance at omega=0

% Reduce the highest frequency in spectra to Nmax
Nmax=22;
Phasor_Qav(Nmax+1:Nh1)=0+0i;
Phasor_pa(Nmax+1:Nh1)=0+0i;
Qav_filter=bfourier_def1(Phasor_Qav, time);
pa_filter =bfourier_def1(Phasor_pa, time);

% Plot original and filtered pressure
figure(1)
plot(tp,pa_measured,'r-',time,pa_filter,'k-'); 
grid on
legend('original pressure','filtered pressure')

% Plot original and filtered pressure
figure(2)
plot(tQ,Qav_measured,'r-',time,Qav_filter,'k-'); 
grid on
legend('original flow rate','filtered flow rate')

 if (Flag==8)||(Flag==9)||(Flag==10) ||(Flag==11) ||(Flag==12) ||(Flag==13)
%----------------------------------------------------
% Obtain measured and exact (scanned) Impedance-data
%----------------------------------------------------
 Zabs=array3(:,2)./1330;    % Magnitude: we divided by 1330 to convert from (dyn*s/cm3) to (mmHg*s/ml)
 Zang=array4(:,2)*pi/180;   % Phase: degrees 
 Nhs=length(Zabs);
%-------------------------------------------
 end
 
 
% Correct pressure and flow rate to obtain meaured input impedance
Phasor_pa_corr=Phasor_pa;
Phasor_Qav_corr=Phasor_Qav;
for i=1:Nhs
    ih=i+1;
    Zin(ih)=Zabs(i)*cos(Zang(i))+1i*Zabs(i)*sin(Zang(i));
    Phasor_pa_corr(ih)=0.5*(Phasor_pa(ih)+Zin(ih)*Phasor_Qav(ih));
    Phasor_Qav_corr(ih)=0.5*(Phasor_Qav(ih)+Phasor_pa(ih)/Zin(ih));
end

Qav_corr=bfourier_def1(Phasor_Qav_corr, time);
pa_corr =bfourier_def1(Phasor_pa_corr, time);

% Plot original and corrected pressure
figure(3)
plot(tp,pa_measured,'r-',time,pa_corr,'b-'); 
grid on
legend('original pressure','corrected pressure')

% Plot original and corrected pressure
figure(4)
plot(tQ,Qav_measured,'r-',time,Qav_corr,'b-'); 
grid on
legend('original flow rate','corrected flow rate')

freq_m(2:Nhs+1)=freq(2:Nhs+1);
for i=2:Nhs+1
    Zin_corr(i)=Phasor_pa_corr(i)/Phasor_Qav_corr(i);
    Zin_m(i)=Phasor_pa(i)/Phasor_Qav(i);
end
figure (5)
plot(freq_m(2:Nhs+1),abs(Zin_corr(2:Nhs+1)),'o',freq_m(2:Nhs+1),Zabs,'-')
figure (6)
plot(freq_m(2:Nhs+1),atan2(imag(Zin_corr(2:Nhs+1)),real(Zin_corr(2:Nhs+1)))*180/pi,'o',freq_m(2:Nhs+1),Zang*180/pi,'-')

if (Flag==8) % New Data
   fid1=fopen('QAV_Nichols_28_corr.POD','wt');
   for i=1:length(time); fprintf(fid1,'%f %f\n',time(i),Qav_corr(i)); end;
   fid2=fopen('PAV_Nichols_28_corr.POD','wt');
   for i=1:length(time); fprintf(fid2,'%f %f\n',time(i),pa_corr(i)); end;
elseif (Flag==9) % New Data
   fid1=fopen('QAV_Nichols_52_corr.POD','wt');
   for i=1:length(time); fprintf(fid1,'%f %f\n',time(i),Qav_corr(i)); end;
   fid2=fopen('PAV_Nichols_52_corr.POD','wt');
   for i=1:length(time); fprintf(fid2,'%f %f\n',time(i),pa_corr(i)); end;
elseif (Flag==10) % New Data
   fid1=fopen('QAV_Nichols_68_corr.POD','wt');
   for i=1:length(time); fprintf(fid1,'%f %f\n',time(i),Qav_corr(i)); end;
   fid2=fopen('PAV_Nichols_68_corr.POD','wt');
   for i=1:length(time); fprintf(fid2,'%f %f\n',time(i),pa_corr(i)); end;
elseif (Flag==11) % New Data
   fid1=fopen('QAV_Nichols_Normotensive_corr.POD','wt');
   for i=1:length(time); fprintf(fid1,'%f %f\n',time(i),Qav_corr(i)); end;
   fid2=fopen('PAV_Nichols_Normotensive_corr.POD','wt');
   for i=1:length(time); fprintf(fid2,'%f %f\n',time(i),pa_corr(i)); end;
elseif (Flag==12) % New Data
   fid1=fopen('QAV_Nichols_Mild_Hypertension_corr.POD','wt');
   for i=1:length(time); fprintf(fid1,'%f %f\n',time(i),Qav_corr(i)); end;
   fid2=fopen('PAV_Nichols_Mild_Hypertension_corr.POD','wt');
   for i=1:length(time); fprintf(fid2,'%f %f\n',time(i),pa_corr(i)); end;
elseif (Flag==13) % New Data
   fid1=fopen('QAV_Nichols_Severe_Hypertension_corr.POD','wt');
   for i=1:length(time); fprintf(fid1,'%f %f\n',time(i),Qav_corr(i)); end;
   fid2=fopen('PAV_Nichols_Severe_Hypertension_corr.POD','wt');
   for i=1:length(time); fprintf(fid2,'%f %f\n',time(i),pa_corr(i)); end;
end
