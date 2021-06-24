clear 
close all 
clc
zc = 0.814;
g = 9.81; % gravity
dt = 0.005;
Qe = 1;
R = 1e-06;
%%%%%%%%%%%%%%%%%% Change Preview Time Here For Test %%%%%%%%%%%%%%%%%%
t_preview = 1.6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumEch = t_preview/dt;
simulation_time = 10;
% Matrix A, B, C, and D from cart table model
A = [0 1 0;
     0 0 1;
     0 0 0];
B = [0; 0; 1];
C = [1 0 -zc/g];
D = 0;

% Convert continuous system to discrete system
sys_c = ss(A, B, C, D);            
sys_d = c2d(sys_c, dt);
[A_d, B_d, C_d, D_d] = ssdata(sys_d);

% A, B, C matrix for LQI control
A_tilde = [1 C_d*A_d;
           zeros(3,1) A_d];
B_tilde = [C_d*B_d;
           B_d];
C_tilde = [1 0 0 0];

% Q matrix
Q = [Qe 0 0 0;
     0 0 0 0;
     0 0 0 0;
     0 0 0 0];
% Find optimal gain and ricatti equation with dlqr function
% K = optimal gain
% P = ricatti equation
[K, P] = dlqr(A_tilde, B_tilde, Q, R);

Gi = K(1);
Gx = K(2:end);

% Calculating preview gain
N = 0:dt:t_preview;
Gp = zeros(1, length(N));
%Gd(1,1) = -Gi; % First Gd = -Gi

% % Original formula in the paper
% Ac_tilde = A_tilde - B_tilde * (R+B_tilde'*P*B_tilde)^(-1)*B_tilde'*P*A_tilde;
% % Simplified formula
Ac_tilde = A_tilde - B_tilde * K;

I_tilde = [1;0;0;0];
% First X_tilde
X_tilde = -Ac_tilde'*P*I_tilde;

for i = 2:length(N)
    % % Gd(i) is calculated from previous X_tilde
    Gp(1,i) = (R+B_tilde'*P*B_tilde)^(-1)*B_tilde'*X_tilde;
    % % Calculating next X_tilde
    X_tilde = Ac_tilde'*X_tilde;
end
plot(N, -Gp)
xlabel('time [s]');
ylabel('previw gain')
axis([0, 2, 0, 1500]);
% hold on
p =load('signaux.mat');
signalEs=p.Scenario{1}.Values;
signalCa=p.Scenario{2}.Values;
vec = 0:dt:10;
signalEsTs = resample(signalEs,[0:dt:10]);
signalCaTs = resample(signalCa,[0:dt:10]);
p_ref_all_Es = signalEsTs.Data;
p_ref_all_Ca = signalCaTs.Data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_es(:,1) = [0; 0; 0];
x_ca(:,1) = [0; 0; 0];
u_es(1) = 0;
u_ca(1) = 0;
p_es(1) = 0;
p_ca(1) = 0;

for k = 2 : simulation_time/dt - NumEch
    x_es(:,k) = A_d*x_es(:,k-1)+B_d*u_es(k-1);
    x_ca(:,k) = A_d*x_ca(:,k-1)+B_d*u_ca(k-1);
    
    p_es(k) = C_d*x_es(:,k-1);
    p_ca(k) = C_d*x_ca(:,k-1);
    
    e_es(k) = p_es(k-1) - p_ref_all_Es(k-1);
    e_ca(k) = p_ca(k-1) - p_ref_all_Ca(k-1);
    SumGainXRefEs = 0;
    SumGainXRefCa = 0;
    for j = 2 : NumEch
        SumGainXRefEs = SumGainXRefEs + Gp(j)*p_ref_all_Es(k+j);
        SumGainXRefCa = SumGainXRefCa + Gp(j)*p_ref_all_Ca(k+j);
    end
    u_es(k) = -Gi*sum(e_es)-Gx*x_es(:,k) - SumGainXRefEs;
    u_ca(k) = -Gi*sum(e_ca)-Gx*x_ca(:,k) - SumGainXRefCa;
end
time = 0:dt:simulation_time - t_preview - dt;
figure(2)
plot(time,p_ref_all_Es(1:length(time)))
figure(3)
subplot(2,1,1)
plot(time,p_es, time,x_es(1,:), time,p_ref_all_Es(1:length(time)))
ylabel('x [m]')
axis([0, 7, -0.05, 1]);
grid on
legend("ZMP", "CoM", "ZMP ref")

subplot(2,1,2)
plot(time,p_ca, time,x_ca(1,:), time,p_ref_all_Ca(1:length(time)))
xlabel('time [s]');
ylabel('y [m]')
axis([0, 7, -0.15, 0.15]);
grid on
legend("ZMP", "CoM", "ZMP ref")
