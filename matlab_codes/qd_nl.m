%% Quad-Rotor Dynamics for Non-Linear Model
% x = [phi theta psi p q r u v w x  y  z ]'
%      1   2     3   4 5 6 7 8 9 10 11 12
% X = [p q r phi theta psi z Vz y Vy x  Vx]'
%      1 2 3 4   5     6   7 8  9 10 11 12
% 1  ==  4
% 2  ==  5
% 3  ==  6
% 4  ==  1
% 5  ==  2
% 6  ==  3
% 7  ==  12
% 8  ==  10
% 9  ==  8
% 10 ==  11
% 11 ==  9
% 12 ==  7
% u = [netT Mx My Mz]'
%      1    2  3  4

%% 
% phi_dot    =  p + r[c(phi)t(theta)] + q[s(phi)t(theta)]
% theta_dot  =  q[c(phi)] - r[s(phi)]
% psi_dot    =  r(c(phi)/c(theta)) + q(s(phi)/c(theta))
% p_dot      =  (Iy - Iz)rq/Ix + (_tx + _twx)/Ix
% q_dot      =  (Iz - Ix)pr/Iy + (_ty + _twy)/Iy
% r_dot      =  (Ix - Iy)pq + (_tz + _twz)/Iz
% u_dot      =  rv - qw - g[s(theta)] + fwx/m
% v_dot      =  pw - ru + g[s(phi)c(theta)] + fwy/m
% w_dot      =  qu - pv + g[c(theta)c(phi)] + (fwz - ft)/m
% x_dot      =  w[s(phi)s(psi)] + c(phi)c(psi)s(theta)] - 
%               v[c(phi)s(psi) - c(psi)s(phi)s(theta)] + 
%               u[c(psi)c(theta)]
% y_dot      =  v[c(phi)c(psi) + s(phi)s(psi)s(theta)] -
%               w[c(psi)s(phi) - c(phi)s(psi)s(theta)] + 
%               u[c(theta)s(psi)]
% z_dot      =  w[c(phi)c(theta)] - u[s(theta)] + v[c(theta)s(phi)]

%% Code 
clc
clear all
% close all

%% Constt
global m g Ix Iy Iz tow_x tow_y tow_z tow_wx tow_wy tow_wz fwx fwy fwz ft
 m = 1.104; % mass
 g = 9.81;  % gravity
 d = 1.225; % density
 Ix  = 0.008562874765838073; 
 Iy  = 0.008788914621963906;
 Iz  = 0.015570395039175332;
 tow_x = 0;
 tow_y = 0;
 tow_z = 0;
 tow_wx = 0;
 tow_wy = 0;
 tow_wz = 0;
 fwx = 0;
 fwy = 0;
 fwz = 0;
%  ft = m*g;
%  ft = 0;
momentArm   = 0.225; %half of quadcopter diagonal

%% Initial Values
x_0 = zeros(12,1);
% x_0(4) = 1.0; % phi
% x_0(5) = 0.1; % theta
% x_0(6) = 0.1; % psi
% x_0(7) = 1;
% x_0(9)  = 1;
% x_0(11) = 1.0; % x
% x_0(9) = 1.0;  % y
% x_0(7) = 1.0;  % z
%% CGL Nodes
global ptspan 
ncgl = 40;
tspan = zeros(ncgl+1,1);
for i=1:(ncgl+1)
tspan(i) = cos((pi*(i-1))/(ncgl));
end
tspan;
t0 = 0.00;
tf = 15.00;

for i=1:(ncgl+1)
ptspan(i) = (tf/2.0)*(tspan(i)+1.0);
end 
ptspan = fliplr(ptspan);
ptspan = ptspan';

%% ODE Solver
%% Load Gain Matrix from Linear System
load K_matrix
global k
k = K;
%% Load Input Control
load inp_u.mat
global U
U = u;
% figure;
% plot(ptspan, U');
% title('Input Control from LQR');

% load U_unchanged.mat;
% U_ = U_';
% save param
% x = zeros(12,ncgl+1)
odeoptions = odeset('RelTol',1e-5,'AbsTol',1e-7);
[t, x] = ode45(@comp, ptspan, x_0, odeoptions);

figure;
plot(t,x)
legend('p', 'q', 'r', 'phi', 'theta', 'psi', 'z', 'Vz', 'y', 'Vy', 'x',  'Vx');
title('Non Linear States');
figure;
hold on
plot(t,x(:,8:8 ))
% plot(t,x(:,9:9))
% plot(t,x(:,7:7))
legend('Vz')
title('Non Linear Vz')
% hold off