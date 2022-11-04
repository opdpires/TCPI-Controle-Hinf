% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Nov 2, 2022
%
% file: pertinence.m
%
% brief: this function calculates the pertinence values of the actual
% states of the system to each vertex of the politope.

function [h,verif] = pertinence(xt)
% h is the pertinence vector. Each vertex is associated with a pertinence value. 
%
% xt is the state vector. The states that are associated with a parameter
% are going to be used in the computation of the  pertinence values for
% both vertices related to those states.
%
% verif is 1 if the sum of the pertinence values is equal to 1. Otherwise,
% it is equal to 0.

%% Variables inicialization
h = zeros(32,1);
sum = 0;
g = 9.8;

% Ball characteristics values
m = 0.13; % Mass in Kg
R = 0.03; % Radius in meters
Ib = 2/5 * m * R^2; % Inertial Momentum


Kv = m/(m+Ib/(R^2));

% extremal values definition
dtheta1_min = -0.1;
dtheta1_max = +0.1;
dtheta2_min = -0.1;
dtheta2_max = +0.1;
sincTheta1_min = 0.988;
sincTheta1_max = 1;
sincTheta2_min = 0.988;
sincTheta2_max = 1;

Hi_min = Kv*dtheta1_min^2;     Hi_max = Kv*dtheta1_max^2;
Hj_min = Kv*dtheta2_min^2;     Hj_max = Kv*dtheta2_max^2; 
Hk_min = Kv*dtheta1_min*dtheta2_min;   Hk_max = Kv*dtheta1_max*dtheta2_max; 
Hp_min = Kv*g*sincTheta1_min;     Hp_max = Kv*g*sincTheta1_max;
Hq_min = Kv*g*sincTheta2_min;     Hq_max = Kv*g*sincTheta2_max;

%% Calculation of the pertinence values associated with each parameter

Hi(1) = (Hi_max - Kv*xt(4)^2)/(Hi_max-Hi_min); % h1
Hi(2) = 1-Hi(1);

Hj(1) = (Hj_max - Kv*xt(8)^2)/(Hj_max-Hj_min); % h2
Hj(2) = 1-Hj(1);

Hk(1) = (Hk_max - Kv*xt(4)*xt(8))/(Hk_max-Hk_min); % h3
Hk(2) = 1-Hk(1);

Hp(1) = (Hp_max - Kv*sin(xt(3))/xt(3))/(Hp_max-Hp_min); % h4
Hp(2) = 1-Hp(1);

Hq(1) = (Hq_max - Kv*sin(xt(7))/xt(7))/(Hq_max-Hq_min); % h5
Hq(2) = 1-Hq(1);

%% Calculates the pertinence values for the equivalent parameter "rho"

h(1)  = Hi(1)*Hj(1)*Hk(1)*Hp(1)*Hq(1);
h(2)  = Hi(1)*Hj(1)*Hk(1)*Hp(1)*Hq(2);
h(3)  = Hi(1)*Hj(1)*Hk(1)*Hp(2)*Hq(1);
h(4)  = Hi(1)*Hj(1)*Hk(1)*Hp(2)*Hq(2);
h(5)  = Hi(1)*Hj(1)*Hk(2)*Hp(1)*Hq(1);
h(6)  = Hi(1)*Hj(1)*Hk(2)*Hp(1)*Hq(2);
h(7)  = Hi(1)*Hj(1)*Hk(2)*Hp(2)*Hq(1);
h(8)  = Hi(1)*Hj(1)*Hk(2)*Hp(2)*Hq(2);
h(9)  = Hi(1)*Hj(2)*Hk(1)*Hp(1)*Hq(1);
h(10) = Hi(1)*Hj(2)*Hk(1)*Hp(1)*Hq(2);
h(11) = Hi(1)*Hj(2)*Hk(1)*Hp(2)*Hq(1);
h(12) = Hi(1)*Hj(2)*Hk(1)*Hp(2)*Hq(2);
h(13) = Hi(1)*Hj(2)*Hk(2)*Hp(1)*Hq(1);
h(14) = Hi(1)*Hj(2)*Hk(2)*Hp(1)*Hq(2);
h(15) = Hi(1)*Hj(2)*Hk(2)*Hp(2)*Hq(1);
h(16) = Hi(1)*Hj(2)*Hk(2)*Hp(2)*Hq(2);
h(17) = Hi(2)*Hj(1)*Hk(1)*Hp(1)*Hq(1);
h(18) = Hi(2)*Hj(1)*Hk(1)*Hp(1)*Hq(2);
h(19) = Hi(2)*Hj(1)*Hk(1)*Hp(2)*Hq(1);
h(20) = Hi(2)*Hj(1)*Hk(1)*Hp(2)*Hq(2);
h(21) = Hi(2)*Hj(1)*Hk(2)*Hp(1)*Hq(1);
h(22) = Hi(2)*Hj(1)*Hk(2)*Hp(1)*Hq(2);
h(23) = Hi(2)*Hj(1)*Hk(2)*Hp(2)*Hq(1);
h(24) = Hi(2)*Hj(1)*Hk(2)*Hp(2)*Hq(2);
h(25) = Hi(2)*Hj(2)*Hk(1)*Hp(1)*Hq(1);
h(26) = Hi(2)*Hj(2)*Hk(1)*Hp(1)*Hq(2);
h(27) = Hi(2)*Hj(2)*Hk(1)*Hp(2)*Hq(1);
h(28) = Hi(2)*Hj(2)*Hk(1)*Hp(2)*Hq(2);
h(29) = Hi(2)*Hj(2)*Hk(2)*Hp(1)*Hq(1);
h(30) = Hi(2)*Hj(2)*Hk(2)*Hp(1)*Hq(2);
h(31) = Hi(2)*Hj(2)*Hk(2)*Hp(2)*Hq(1);
h(32) = Hi(2)*Hj(2)*Hk(2)*Hp(2)*Hq(2);

%% Verification if it really IS a simplex

% "sum" must be equal to 1 in the end of the "for" loop
for i=1:32
    sum = sum + h(i);
end

% verif == 1 -> it is a simplex
% verif == 0 -> there's something wrong (it's not a simplex)
verif = (sum == 1);

end
