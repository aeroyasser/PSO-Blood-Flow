% Goal function for the identification of parameters of 6-element
% Windkessel model defined by (eta + C0) // [(L + r) + (C1 // Rp)]
function [RMSp]=Goal_function(param)

global Phasor_Qav Phasor_pa Phasor_pa_model freq   % phasors and frequencies

% model parameters;
%      eta=param(1)
%      C0=param(2)
%      L=param(3)
%      R1=param(4)
%      C1=param(5)
%      Rp=param(6)

Nh=length(freq);   % number of harmonics
RMSp(1:Nh)=0;  
for i=1:Nh
    omega=2*pi*freq(i);  % circular frequency
% Input impedance of the 6-element Windkessel
    Z1=param(6)/(1+1i*omega*param(5)*param(6));  % impedance of C1 and Rp in parallel
    Z2=param(4)+1i*omega*param(3);               % impedance of r and L in series
    Z3_rec=1i*omega*param(2)/(1+1i*omega*param(2)*param(1)); % reciprocal of impedance of eta and C0 in series
    Zin=(Z3_rec+1/(Z1+Z2))^(-1);  % Input impedance
% Phasor of pa from the mathematical model
    Phasor_pa_model(i)=Zin*Phasor_Qav(i);
    RMSp(i)=sqrt(( real(Phasor_pa(i)-Phasor_pa_model(i))^2+ ...  % looking in time domain, this is the sum of sqares 
                   imag(Phasor_pa(i)-Phasor_pa_model(i))^2)/2);  % of pressure differences devided by the number of points
end