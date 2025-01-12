%% Clear Workspace
clc; clear; close all

%% Parameters

[m , M, h, a, r, I, I_x, I_y, I_z, g] = deal(0.141, 2.313, 0.254, 0.165, 0.0615, 2.5e-4, 0.1986, 0.1942, 0.0056, 9.81);

%% Linear State Space Realisation

A = [0, 0, 0, 1, 0, 0;...
       0, 0, 0, 0, 1, 0;...
       0, 0, 0, 0, 0, 1;...
       0, 0, -(M^2*g*h^2*r^2)/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2), 0, 0, 0;...
       0, 0, 0, 0, 0, 0;...
       0, 0, (2*g*h*M^2*r^2 + 4*g*h*m*M*r^2 + 4*I*g*h*M)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2), 0, 0, 0];

B = [0, 0;...
      0, 0;...
      0, 0;...
      (r*(M*h^2 + M*r*h + I_y))/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2), (r*(M*h^2 + M*r*h + I_y))/(2*I*I_y + 2*I*M*h^2 + I_y*M*r^2 + 2*I_y*m*r^2 + 2*M*h^2*m*r^2);...
    -(a*r)/(2*m*a^2*r^2 + 2*I*a^2 + I_z*r^2), (a*r)/(2*m*a^2*r^2 + 2*I*a^2 + I_z*r^2);...
    -(4*I + 2*M*r^2 + 4*m*r^2 + 2*M*h*r)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2), -(4*I + 2*M*r^2 + 4*m*r^2 + 2*M*h*r)/(4*I*I_y + 4*I*M*h^2 + 2*I_y*M*r^2 + 4*I_y*m*r^2 + 4*M*h^2*m*r^2)];

C = zeros(3, 6);
C(1:3, 1:3) = eye(3); % outputs are rho, psi and theta

D = zeros(3, 2);

%% Check for Observability

W_o = obsv(A, C);
r = rank(W_o);

if r == size(A, 1) % rank n(number of states)
    disp("(A, C) is observable.")
else
    disp("(A, C) is not observable.")
end