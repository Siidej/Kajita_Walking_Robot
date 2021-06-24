clear all
close all
clc
NL = 240;   % 1.6/0.005-> predicted steps
Nl = 160;
T = 0.005;
fs = 1/T;


A = [1, T, (T^2)/2; 0, 1, T; 0, 0, 1];
b = [(T^3)/6; (T^2)/2; T];
c = [1 0 -0.083];
Q = 1;
%Qx = [1 0 0; 0 1 0; 0 0 1]; 
R = 1e-06;
[K,P] = dlqr(A,b,c'*c, R);


% f = zeros(1,Nl);
for i = 1:NL
   f(i) = ((R + transpose(b)*P*b)^(-1))*transpose(b)*((transpose(A-b*K))^(i-1))*transpose(c)*1;  
end

% time = 0:T:1-T;
% plot(time, -f)

p =load('signalCarre.mat');
test=p.Scenario{1}.Values;
vec = 0:T:10;
test2 = resample(test,[0:T:10]);
p_ref_all = test2.Data;
%size(p_ref_all)

