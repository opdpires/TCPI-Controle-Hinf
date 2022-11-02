% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Nov 1, 2022
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
Hi = zeros(2,1);
Hi = zeros(2,1);
Hi = zeros(2,1);
Hi = zeros(2,1);
Hi = zeros(2,1);

h = zeros(32,1);
sum = 0;

% TODO: temos que encontrar os valores reais desses limites
Hi_min = 0; Hi_max = 0;
Hj_min = 0; Hj_max = 0;
Hk_min = 0; Hk_max = 0;
Hp_min = 0; Hp_max = 0;
Hq_min = 0; Hq_max = 0;

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

h(1) = Hi(1)*Hj(1)*Hk(1)*Hp(1)*Hq(1);
h(2) = Hi(1)*Hj(1)*Hk(1)*Hp(1)*Hq(2);
h(3) = Hi(1)*Hj(1)*Hk(1)*Hp(2)*Hq(1);
h(4) = Hi(1)*Hj(1)*Hk(1)*Hp(2)*Hq(2);
h(3) = Hi(1)*Hj(1)*Hk(2)*Hp(1)*Hq(1);
% TODO: Fazer para cada combinacao
% estou fazendo como se fosse um numero em binario
% .
% .
% .
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