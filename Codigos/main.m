% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Nov 2, 2022
%
% file: main.m
%
% brief: this file holds the simulation of the "ball on the plate" system
% with the proposed controller using H infinity control and TS-Fuzzy model
% representation.

clear;clc;close all;

g = 9.8;
eps = 1e-6;

% Ball characteristics values
m = 0.13; % Mass in Kg
R = 0.03; % Radius in meters
Ib = 2/5 * m * R^2; % Inertial Momentum

Kv = m/(m+Ib/(R^2));

% extremal values definition
dtheta1_min = -0.1;     dtheta1_max = +0.1;
dtheta2_min = -0.1;     dtheta2_max = +0.1;

dtheta1_2_min = 0;                          dtheta1_2_max = dtheta1_max^2;
dtheta2_2_min = 0;                          dtheta2_2_max = dtheta2_max^2;
dth1th2_min = -dtheta1_max*dtheta2_max;     dth1th2_max = +dtheta1_max*dtheta2_max;
sincTheta1_min = 0.988;                     sincTheta1_max = 1;
sincTheta2_min = 0.988;                     sincTheta2_max = 1;

Hi_min = Kv*dtheta1_2_min;             Hi_max = Kv*dtheta1_2_max;
Hj_min = Kv*dtheta2_2_min;             Hj_max = Kv*dtheta2_2_max; 
Hk_min = Kv*dth1th2_min;               Hk_max = Kv*dth1th2_max; 
Hp_min = Kv*g*sincTheta1_min;          Hp_max = Kv*g*sincTheta1_max;
Hq_min = Kv*g*sincTheta2_min;          Hq_max = Kv*g*sincTheta2_max;

A{1} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{2} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{3} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{4} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{5} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{6} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{7} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{8} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{9} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{10} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{11} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{12} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{13} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{14} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{15} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{16} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{17} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{18} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{19} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{20} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{21} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{22} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{23} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{24} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_min   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{25} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{26} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{27} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{28} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_min   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{29} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{30} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{31} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            

A{32} = [0        0        0        0        1 0 0 0;
        Hi_max   Hk_max   Hp_max   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_max   Hj_max   0        Hq_max   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];

Bu = [1 0;%0 0
      1 0;%0 0
      1 0;%0 0
      1 0;
      0 1;%0 0
      0 1;%0 0
      0 1;%0 0
      0 1];

Bw = [0 0;
      0 0;%1 0;
      0 0;
      0 0;
      0 0;
      0 0;%0 1;
      0 0;
      0 0];

Cz = eye(8);
Dzu = zeros(8,2);
Dzw = zeros(8,2);

%% Temporal simulation

out = q_LPV_Hinf_optimization(eps, A, Bu, Bw, Cz, Dzu, Dzw);

if (out.sol.problem ~= 0)
    out.sol.info

else
    K = out.K;

% TODO: definir os valores da perturbação
w = [0;0]; %improviso

x0 = [0.01; 0.001; 0.001; 0.001; 0.01; 0.001; 0.001; 0.001];

[t,xs] = ode45(@(t,xs)TS_System(t,xs,w,A,Bu,Bw,K), [0,100], x0);

%% Graph of the evolution of the states

figure(1)

subplot(2,2,1)
plot(t, xs(:,1), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,5), 'r-', 'LineWidth', 1.5);
xlabel('$t(seconds)$','Interpreter','latex','FontSize',16);
ylabel('$x(m)$','Interpreter','latex','FontSize',16);
lgd1 = legend('$p_{x}$','$p_{y}$');
set(lgd1,'Interpreter','latex','FontSize',16);
legend('boxoff')
legend('Location','best')
grid minor

subplot(2,2,2)
plot(t, xs(:,2), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,6), 'r-', 'LineWidth', 1.5);
xlabel('$t(seconds)$','Interpreter','latex','FontSize',16);
ylabel('$x(m/s)$','Interpreter','latex','FontSize',16);
lgd2 = legend('$\dot{p}_x$','$\dot{p}_y$');
set(lgd2,'Interpreter','latex','FontSize',16);
legend('boxoff')
legend('Location','best')
grid minor

subplot(2,2,3)
plot(t, xs(:,3), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,7), 'r-', 'LineWidth', 1.5);
xlabel('$t(seconds)$','Interpreter','latex','FontSize',16);
ylabel('$x(rad)$','Interpreter','latex','FontSize',16);
lgd3 = legend('$\theta_x$','$\theta_y$');
set(lgd3,'Interpreter','latex','FontSize',16);
legend('boxoff')
legend('Location','best')
grid minor

subplot(2,2,4)
plot(t, xs(:,4), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,8), 'r-', 'LineWidth', 1.5);
xlabel('$t(seconds)$','Interpreter','latex','FontSize',16);
ylabel('$x(rad/s)$','Interpreter','latex','FontSize',16);
lgd4 = legend('$\dot{\theta}_x$','$\dot{\theta}_y$');
set(lgd4,'Interpreter','latex','FontSize',16);
legend('boxoff')
legend('Location','best')
grid minor
end

%% System Definition
% calculate and return the next states of the system
function dx = TS_System(t,x,w,A,Bu,Bw,K)

Ax = zeros(8);
Kx = zeros(1,8);


[h,verif] = pertinence(x);

if (verif~=1)
    disp("It is not a simplex")
else
    for i=1:32
        Ax = Ax + h(i)*A{i};
        Kx = Kx + h(i)*K{i};
    end

    u = Kx*x;

    dx = Ax*x + Bu*u + Bw*w;
end
end
