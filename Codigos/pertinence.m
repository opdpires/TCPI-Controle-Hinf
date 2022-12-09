% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Nov 11, 2022
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
h = zeros(8,1);
parameters; % Load parameters constants
sum = 0;

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
