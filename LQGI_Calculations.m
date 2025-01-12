%% Clear Workspace
clc; clear; close all

%% Parameters

[m , M, h, a, r, I, I_x, I_y, I_z, g] = deal(0.141, 2.313, 0.254, 0.165, 0.0615, 2.5e-4, 0.1986, 0.1942, 0.0056, 9.81);
% reference 
rho_ref = 2;
psi_ref = pi/4;
theta_ref = 0;

%% Linear State Space Realisation

A = [0, 0, 0, 1, 0, 0;...
       0, 0, 0, 0, 1, 0;... 
       0, 0, 0, 0, 0, 1;...
       0, 0, -1*(M^2*g*h^2*r^2)/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2), 0, 0, 0;...
       0, 0, 0, 0, 0, 0;... 
       0, 0, (2*g*h*M^2*r^2 + 4*g*h*m*M*r^2 + 4*I*g*h*M)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2), 0, 0, 0];

B = [0, 0;... 
       0, 0;... 
       0, 0;... 
       (r*(M*h^2 + M*r*h + I_y))/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2), (r*(M*h^2 + M*r*h + I_y))/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2);...
       -1*(a*r)/(2*m*a^2*r^2 + 2*I*a^2 + I_z*r^2), (a*r)/(2*m*a^2*r^2 + 2*I*a^2 + I_z*r^2);...
       -1*(4*I + 2*M*r^2 + 4*m*r^2 + 2*M*h*r)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2), -1*(4*I + 2*M*r^2 + 4*m*r^2 + 2*M*h*r)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2)];

C = eye(6);
C_I = C(1:2, :);
C_Kalman = C(1:3, :);
D = zeros(6, 2);
D_I = D(1:2, :);
D_Kalman = D(1:3, :);

states = {'rho' 'psi' 'theta' 'rho_dot' 'psi_dot' 'theta_dot'};
inputs = {'T_l' 'T_r'};
outputs = {'rho' 'psi' 'theta' 'rho_dot' 'psi_dot' 'theta_dot'};

sys_ol = ss(A, B, C, D, 'statename', states, 'inputname', inputs,'outputname', outputs);
G = tf(sys_ol);
display(G)

%% Compute Pole Placing Gain (LQI)
% Q chosen for faster response in rho by higher penalization
% and a gradual transition to linear motion by penalizing psi less
Q = diag([10, 0.01, 15, 1, 0.1, 1, 10, 10]); 
% R chosen to ensure 
R = 0.1*eye(2);
sys_k = ss(A, B, C_I, D_I);
K_I = lqi(sys_k,  Q, R, 0);
display(K_I)
