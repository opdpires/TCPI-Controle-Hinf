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

%% Initializations
% TODO: definir valores da bolinha
m = 0;
Ib = 0;
R = 0;

Kv = m/(m+Ib/(R^2));

% TODO: definir esses valores limite
Hi_min = Kv*0;     Hi_max = Kv*0;
Hj_min = Kv*0;     Hj_max = Kv*0; 
Hk_min = Kv*0*0;   Hk_max = Kv*0*0; 
Hp_min = Kv*0;     Hp_max = Kv*0;
Hq_min = Kv*0;     Hq_max = Kv*0; 

A{1} = [0        0        0        0        1 0 0 0;
        Hi_min   Hk_min   Hp_min   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        Hk_min   Hj_min   0        Hq_min   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];

% TODO: continuar como se fosse um num em binario

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

x0 = [0.5; 0; 0; 0; 0.5; 0; 0; 0];

[t,xs] = ode45(@(t,xs)TS_System(t,x,w,A,Bu,Bw,K), [0,10], x0);

%% Graph of the evolution of the states

figure(1)
hold on

plot(t, xs(:,1), 'b-', 'LineWidth', 1.5);
plot(t, xs(:,2), 'r-', 'LineWidth', 1.5);
% TODO: plotar outros estados de forma que a figura permaneca legivel

% figure settings
set(gca,'FontSize',14,'Fontname','Times New Roman');
xlabel('$t(s)$','Interpreter','latex','FontSize',20);
ylabel('$x(t)$','Interpreter','latex','FontSize',20);
box on
leg = legend('$~x_1$','$~x_2$','$~x_3$','$~x_4$','$~x_5$','$~x_6$','$~x_7$','$~x_8$');
set(leg,'Interpreter','latex','FontSize',18);
legend('boxoff')
legend('Location','best')

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