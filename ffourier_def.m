% This function computes the coefficient of the Fourier series (forward Fourier decomposition, by definition)
% This is les efficient than the Matlab fft function, but it is not problem since we deal with small vector y
function [Fazor, Nh, freq] = ffourier_def (y, dt)
%    Inlet:
% y = vector containing values of a considered function in time domain
% it is neccessary to give values at the beginning and at the end of time interval (for K divisions there will be K+1 points)
% dt = sampling time (period T=K*dt)
%    Outlet:
% Fazor=S+iC , S is the coefficient multiplying sin(k*omega*time), 
%              C is the coefficient multiplying cos(k*omega*time)
% Nh = is the number oh harmonics (k=1:Nh), including the harmonic at omega=0 
% freq = vector of frequencies

Kp1=length(y);      
K=Kp1-1;           
T=K*dt;             
f0=1/T ;             % basic frequency
Nh=fix((Kp1+1)/2);   
fmax=(Nh-1)*f0;      % maksimal frequency
freq=linspace(0,fmax,Nh);  
Fazor(1:Nh)=0+1i; 

% the first harmonic (at omega=0)
C0=0.5*(y(1)+y(Kp1)); 
for k=2:1:K
    C0=C0+y(k); 
end
C0=C0/K; 
Fazor(1)=0+1i*C0;

% all other harmonics
for i=2:1:Nh  
    n=i-1;    
    S=0; 
    C=0.5*(y(1)+y(Kp1)); 
  for k=2:1:K
    fi=n*2*pi*(k-1)/K;
    S=S+y(k)*sin(fi); 
    C=C+y(k)*cos(fi); 
  end
   S=2*S/K;
   C=2*C/K;

    Fazor(i)=S+1i*C; % Phasor of i-th harmonic
end




