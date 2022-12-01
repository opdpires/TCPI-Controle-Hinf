% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Nov 11, 2022
%
% file: pertinence_8_vertices.m
%
% brief: this function calculates the pertinence values of the actual
% states of the system to each vertex of the politope.

function [h,verif] = pertinence_8_vertices(xt)
% h is the pertinence vector. Each vertex is associated with a pertinence value. 
%
% xt is the state vector. The states that are associated with a parameter
% are going to be used in the computation of the  pertinence values for
% both vertices related to those states.
%
% verif is 1 if the sum of the pertinence values is equal to 1. Otherwise,
% it is equal to 0.

%% Variables inicialization
h = zeros(8,1);
sum = 0;
g = 9.8;

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

Hi_min = Kv*dtheta1_2_min;             Hi_max = Kv*dtheta1_2_max;
Hj_min = Kv*dtheta2_2_min;             Hj_max = Kv*dtheta2_2_max; 
Hk_min = Kv*dth1th2_min;               Hk_max = Kv*dth1th2_max;

%% Calculation of the pertinence values associated with each parameter

Hi(1) = (Hi_max - Kv*xt(4)^2)/(Hi_max-Hi_min); % h1
Hi(2) = 1-Hi(1);

Hj(1) = (Hj_max - Kv*xt(8)^2)/(Hj_max-Hj_min); % h2
Hj(2) = 1-Hj(1);

Hk(1) = (Hk_max - Kv*xt(4)*xt(8))/(Hk_max-Hk_min); % h3
Hk(2) = 1-Hk(1);

%% Calculates the pertinence values for the equivalent parameter "rho"

h(1)  = Hi(1)*Hj(1)*Hk(1);
h(2)  = Hi(1)*Hj(1)*Hk(2);
h(3)  = Hi(1)*Hj(2)*Hk(1);
h(4)  = Hi(1)*Hj(2)*Hk(2);
h(5)  = Hi(2)*Hj(1)*Hk(1);
h(6)  = Hi(2)*Hj(1)*Hk(2);
h(7)  = Hi(2)*Hj(2)*Hk(1);
h(8)  = Hi(2)*Hj(2)*Hk(2);

%% Verification if it really IS a simplex

% "sum" must be equal to 1 in the end of the "for" loop
for i=1:8
    sum = sum + h(i);
end

fprintf("sum = %d\n",sum)

% verif == 1 -> it is a simplex
% verif == 0 -> there's something wrong (it's not a simplex)
verif = (sum < 1.01 && sum > 0.99);

end
