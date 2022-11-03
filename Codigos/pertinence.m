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

% Definindo valores da bolinha
m = 0.13; % Massa da bolinha, em Kg
R = 0.03; % Raio da bolinha, em m
Ib = 2/5 * m * r^2; % Momento de inércia


Kv = m/(m+Ib/(R^2));

% TODO: definir esses valores limite
Hi_min = Kv*0;     Hi_max = Kv*0;
Hj_min = Kv*0;     Hj_max = Kv*0; 
Hk_min = Kv*0*0;   Hk_max = Kv*0*0; 
Hp_min = Kv*0;     Hp_max = Kv*0;
Hq_min = Kv*0;     Hq_max = Kv*0;

%% Calculation of the pertinence values associated with each parameter

Hi(1) = (Hi_max - xt(1))/(Hi_max-Hi_min);
Hi(2) = 1-Hi(1);

Hj(1) = (Hj_max - xt(2))/(Hj_max-Hj_min);
Hj(2) = 1-Hj(1);

Hk(1) = (Hk_max - xt(3))/(Hk_max-Hk_min);
Hk(2) = 1-Hk(1);

Hp(1) = (Hp_max - xt(4))/(Hp_max-Hp_min);
Hp(2) = 1-Hp(1);

Hq(1) = (Hq_max - xt(5))/(Hq_max-Hq_min);
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
