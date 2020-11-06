%This function calculates values of y in given time instances ba using Fourier coefficients (phasor vector)
% by definition y(i)=sum(k=1_Nh)[S*sin(k*omega0*t(i))+C*cos(k*omega0*t(i))]
function [y]=bfourier_def1(Fazor, t)
% y = output values of y at time t
% t = vector of time instances at which the function y will be calculated
% Fazor = Phasor (S+iC) containing coefficients of Fourier series

Nh=length(Fazor);     % number of harmonics
Ktot=length(t);       % number of time points
T=t(Ktot)-t(1);       % period


y(1:1:Ktot) = imag(Fazor(1));   
for np = 2:Nh  
    n=np-1;  
    fi=n*2*pi*t./T;
    y = y + imag(Fazor(np))*cos(fi) + real(Fazor(np))*sin(fi);
end

