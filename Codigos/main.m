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

% Ball characteristics values
m = 0.13; % Mass in Kg
R = 0.03; % Radius in meters
Ib = 2/5 * m * R^2; % Inertial Momentum


Kv = m/(m+Ib/(R^2));

% extremal values definition
dtheta1_min = -10;
dtheta1_max = +10;
dtheta2_min = -10;
dtheta2_max = +10;
sincTheta1_min = 0.988;
sincTheta1_max = 1;
sincTheta2_min = 0.988;
sincTheta2_max = 1;

Hi_min = Kv*dtheta1_min^2;     Hi_max = Kv*dtheta1_max^2;
Hj_min = Kv*dtheta2_min^2;     Hj_max = Kv*dtheta2_max^2; 
Hk_min = Kv*dtheta1_min*dtheta2_min;   Hk_max = Kv*dtheta1_max*dtheta2_max; 
Hp_min = Kv*sincTheta1_min;     Hp_max = Kv*sincTheta1_max;
Hq_min = Kv*sincTheta2_min;     Hq_max = Kv*sincTheta2_max;

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

Bu = [0 0;
      0 0;
      0 0;
      1 0;
      0 0;
      0 0;
      0 0;
      0 1];

Bw = [0 0;
      1 0;
      0 0;
      0 0;
      0 0;
      0 1;
      0 0;
      0 0];

Cz = eye(8);
Dzu = zeros(8,2);
Dzw = zeros(8,2);

%% Temporal simulation

out = q_LPV_Hinf_optimization(A, Bu, Bw, Cz, Dzu, Dzw);

K = out.K;

% TODO: definir os valores da perturbação
w = [0;0]; %improviso

x0 = [0.5; 0; 0; 0; 0.5; 0; 0; 0];

[t,xs] = ode45(@(t,xs)TS_System(t,xs,w,A,Bu,Bw,K), [0,10], x0);

%% Graph of the evolution of the states

figure(1)

subplot(2,2,1)
plot(t, xs(:,1), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,5), 'r-', 'LineWidth', 1.5);
xlabel('$t (seconds)$','Interpreter','latex');
ylabel('$x (meters)$','Interpreter','latex');
lgd1 = legend('~P_1','~P_2','Interpreter','latex');
set(lgd1,'Interpreter','latex');

subplot(2,2,2)
plot(t, xs(:,2), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,6), 'r-', 'LineWidth', 1.5);
xlabel('$t (seconds)$','Interpreter','latex');
ylabel('$x (meters/second)$','Interpreter','latex');
lgd2 = legend('~\dot{P}_1','~\dot{P}_2','Interpreter','latex');
set(lgd2,'Interpreter','latex');

subplot(2,2,3)
plot(t, xs(:,3), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,7), 'r-', 'LineWidth', 1.5);
xlabel('$t (seconds)$','Interpreter','latex');
ylabel('$x (radians)$','Interpreter','latex');
lgd3 = legend('~\theta_1','~\theta_2','Interpreter','latex');
set(lgd3,'Interpreter','latex');

subplot(2,2,4)
plot(t, xs(:,4), 'b-', 'LineWidth', 1.5);
hold on
plot(t, xs(:,8), 'r-', 'LineWidth', 1.5);
xlabel('$t (seconds)$','Interpreter','latex');
ylabel('$x (radians/second)$','Interpreter','latex');
lgd4 = legend('~\dot{\theta}_1','~\dot{\theta}_2','Interpreter','latex');
set(lgd4,'Interpreter','latex');

% figure settings
%set(gca,'FontSize',14,'Fontname','Times New Roman');
%xlabel('$t(s)$','Interpreter','latex','FontSize',20);
%ylabel('$x(t)$','Interpreter','latex','FontSize',20);
%box on
%leg = legend('$~x_1$','$~x_2$','$~x_3$','$~x_4$','$~x_5$','$~x_6$','$~x_7$','$~x_8$');
%set(leg,'Interpreter','latex','FontSize',18);
%legend('boxoff')
%legend('Location','best')

%% System Definition
% calculate and return the next states of the system
function dx = TS_System(t,x,w,A,Bu,Bw,K)

Ax = zeros(8);
Kx = zeros(1,8);

for i=1:32
    Ax = Ax + pertinence(x)*A{i};
    Kx = Kx + pertinence(x)*K{i};
end

u = Kx*x;

dx = Ax*x + Bu*u + Bw*w;

end
