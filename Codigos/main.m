% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Nov 11, 2022
%
% file: main.m
%
% brief: this file holds the simulation of the "ball on the plate" system
% with the proposed controller using H infinity control and TS-Fuzzy model
% representation.

clear;clc;close all;

parameters; % Load parameters constants

eps = 1e-5;


% Building space state matrices
A = matrizA(Hi_min,Hi_max,Hj_min,Hj_max,Hk_min,Hk_max,Kv,g);

Bu = [0 0;
      0 0;
      0 0;
      1 0;
      0 0;
      0 0;
      0 0;
      0 1];

Bw = [0 0;
      1 0;%
      0 0;
      0 0;
      0 0;
      0 1;%
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

    
    [t,xs] = ode45(@(t,xs)TS_System(t,xs,A,Bu,Bw,K,ref_px(t),ref_py(t),w), tspan, x0);
    
    %% Graph of the evolution of the states
    
    figure(1)
    
    subplot(2,2,1)
    plot(t, xs(:,1), 'b-', 'LineWidth', 1.5);
    hold on
    plot(t, xs(:,5), 'r-', 'LineWidth', 1.5); hold on;
    plot(t,ref_px(t),'k'); hold on;
    plot(t,ref_py(t),'c','LineWidth',3);
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
    
    figure(2)
    plot3(xs(:,1),xs(:,5),t)
    %axis([-0.015,0.015,-0.015,0.015])
    xlabel('$p_x (m)$','Interpreter','latex','FontSize',16);
    ylabel('$p_y (m)$','Interpreter','latex','FontSize',16);
    zlabel('$t(s)$','Interpreter','latex','FontSize',16);
    grid minor
    xlim([-0.01 0.01])
    ylim([-0.01 0.01])

end

for i = 1:8
    fprintf("K(%d) = \n",i)
    disp(K{i})
end



