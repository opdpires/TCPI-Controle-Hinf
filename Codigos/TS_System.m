% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques <hugofm2@gmail.com>
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Dec 07, 2022
%
% file: TS_System.m
%
% brief: this function calculates and return the next states of the system,
% based on the space state Matrix, reference and distubance

function dx = TS_System(t,x,A,Bu,Bw,K,ref_px,ref_py,w)

Ax = zeros(8);
Kx = zeros(1,8);

% Definig the pertubation values

% if t < 10
%     w = [0;0];
% else
%     w = exp(-0.25*t)*[sin(t);cos(t)];
% end

% xref = [0.;0;0;0;0;0;0;0];
% xref = [0.01*sin(t);0;0;0;0.01*cos(t);0;0;0];
xref = [ref_px;0;0;0;ref_py;0;0;0];

[h,verif] = pertinence(x);

if (verif~=1)
    disp("It is not a simplex")
else
    for i=1:8
        Ax = Ax + h(i)*A{i};
        Kx = Kx + h(i)*K{i};
    end

    u = Kx*(x-xref);

    dx = Ax*(x-xref) + Bu*u + Bw*w;
end
end
