%% Linear Model
% X = [p q r phi theta psi z Vz y Vy x  Vx]'
%      1 2 3 4   5     6   7 8  9 10 11 12
% u = [netT Mx My Mz]'
%      1    2  3  4
% quad1.mass = 1.104;
% quad1.gravity = 9.81;
% quad1.density = 1.225;
% quad1.Ixx  = 0.008562874765838073;
% quad1.Iyy  = 0.008788914621963906;
% quad1.Izz  = 0.015570395039175332;
% quad1.momentArm   = 0.225; %half of quadcopter diagonal

clc
clear
close all;
load('final_matrices.mat')
% load A_B_Q_R
% load A_B
A = A_trim;
B = B_trim;
clear A_trim;
clear B_trim;
clc
N = zeros(12,4);
Q = eye(12,12);
R = eye(4,4);
% x0 = zeros(12,1) + 1.0;
x0 = zeros(12,1);
% x0(4) = 1.0; % phi
% x0(5) = 0.1; % theta
% x0(6) = 0.1; % psi
% x0(11) = 1.0; % x
% x0(9) = 1.0;  % y
x0(7) = 1.0;  % z
% x0(4)=0.1;x0(5)=0.1;x0(6)=0.1;

ncgl = 40;
tspan = zeros(ncgl+1,1);
for i=1:(ncgl+1)
tspan(i) = cos((pi*(i-1))/(ncgl));
end
tspan;
t0 = 0.00;
tf = 12.00;

for i=1:(ncgl+1)
ptspan(i) = (tf/2.0)*(tspan(i)+1.0);
end
ptspan = fliplr(ptspan);
ptspan = ptspan'
[K,S,e] = lqr(A,B,Q,R,N);
% ptspan = [0 6];
% [t, x] = ode45(@(t,x)(A*x), ptspan, x0);
[t, x] = ode45(@(t,x)(A-B*K)*x, ptspan, x0);
% x(:,7:7);
% load guess_z.mat;
% x(:,7:7) = guess_z;
% plot(t,x(1:2,1:1))
% figure()
% plot(t,x(:,7:8))
% legend('z','w(Vz)')
% title('z and w(vz) in Linear System')
u = -K*x';
save('K_matrix','K');
save('inp_u','u');
figure;
plot(t,u)
legend('netT', 'Mx', 'My', 'Mz');
title('Linear LQR Controls')

figure;
plot(t,x)
legend('p', 'q', 'r', 'phi', 'theta', 'psi', 'z', 'Vz', 'y', 'Vy', 'x',  'Vx');
title('Linear LQR States')

figure;
hold on
plot(t,x(:,8:8))
% plot(t,x(:,9:9))
% plot(t,x(:,7:7))
legend('Vz')
title('Linear Vz')
%% rearrange x and u to give it as input to IPOPT

x = x'
x = fliplr(x);
x = x';
X = reshape(x,[1,(ncgl+1)*12]);

U = u
U = fliplr(U);
U = U'
% figure
% plot(t,U)

U = reshape(U,[1,(ncgl+1)*4]);
X_U = [X,U]'
%%

% [cgl, w] = clencurt(ncgl);
% w;

fileID = fopen('x.txt','w');
fprintf(fileID,'%3f\n',X_U);
fclose(fileID);
type('x.txt');

% Sir x has been generated. No need to run further line of code. Copy this x into IPOPT and play with it.
% put a breakpoint here
%%
save all;
[tsim1, xsim1] = ode45(@dxdt1, t, x0);
% figure()
% plot(t,flip(x))
% hold on
% figure
% plot(t,xsim1)
% % legend('ODE','ODE LQR')
% title('ODE LQR states')
% hold off

% figure()
% nth=9;
% plot(t,flip(x(:,nth:nth)))
% hold on
% plot(t,xsim1(:,nth:nth))
% legend('ODE','ODE LQR')
% title('y')
% hold off

%%
% figure()
% nth=11;
% plot(t,flip(x(:,nth:nth)))
% hold on
% plot(t,xsim1(:,nth:nth))
% legend('ODE','ODE LQR')
% title('x')
% hold off


% figure()
% nth=12;
% plot(t,flip(x(:,nth:nth)))
% hold on
% plot(t,xsim1(:,nth:nth))
% legend('ODE','ODE LQR')
% title('x')
% hold off

%%

% load U_changed.mat
load U_unchanged.mat
figure
plot(t,U_)
title('IPOPT cntrols')
U_ = U_';
save all;
[tsim2, xsim2] = ode45(@dxdt2, t, x0);
figure;
plot(t, xsim2)
title('ODE IPOPT')
%%

figure;
changed_u = [1.0000
1.0000
1.0000
0.9999
0.9997
0.9994
0.9987
0.9976
0.9961
1.7939
1.7910
1.7873
1.7828
1.7773
1.7709
1.7635
1.7552
1.7459
1.7359
1.7251
0.9136
0.9016
0.8893
0.8766
0.8638
0.8511
0.8385
0.8261
0.8142
0.3028
0.2920
0.2820
0.2728
0.2645
0.2571
0.2508
0.2456
0.2415
0.2385
0.2367
0.2361];
% plot(t,changed_u)
plot(t,xsim1(:,7:7))
hold on
plot(t,xsim2(:,7:7))
legend('ODE LQR z','ODE IPOPT z')
% title('diff in z in changed x')
title('z changed')
hold off
% saveas(f,'/Users/rishabh/Documents/MATLAB/pertubrations','png');
% figure
% stem(t,x);
% figure
% scatter(t,x)
% plot(t,x,'.')
