% authors: Pedro Otávio Fonseca Pires <piresopd@gmail.com>
%          Hugo Ferreira Marques
%          Gabriela Gonçalves Amorim de Oliveira
%
% created on: Nov 2, 2022
%
% file: q_LPV_Hinf_optimization.m
%
% brief: this function calculates the gain matrix for each vertex of the
% politope. 
% It uses H infinity control and Takagi-Sugeno Fuzzy model
% representation. An optimization problem is obtained and as we solve it,
% we obtain the gain matrices.

function out = q_LPV_Hinf_optimization(eps, A, Bu, Bw, Cz, Dzu, Dzw)
% out is the output of the function. It contains, mainly, the gain matrix,
% the value of the limit of the gain of the disturbance to the output and
% also contains the value of the P matrix of the Lyapunov function.
%
% A, Bu and Bw are the matrices of the state equation.
% Cz, Dzu and Dzw are the matrices of internal variable equation.

%% Initializations

nx = size(A{1},1) ; % number of states
nu = size(Bu,2) ; % number of inputs
[~,r] = size(A) ; % number of fuzzy rules

sigma = sdpvar(1,1,'full');
X = sdpvar(nx,nx,'symmetric');

for i=1:r
    Z{i} = sdpvar(nu,nx,'full');
end
%% Vertices Matrices

M = [sigma  0;
     0      sigma];

for i=1:r
    Gamma{i} = [(A{i}*X)'+A{i}*X+(Bu*Z{i})'+Bu*Z{i}           Bw                    (Cz*X+Dzu*Z{i})';
                Bw'                                           -M                    Dzw';
                Cz*X+Dzu*Z{i}                                 Dzw                   -eye(nx,nx)];
end

%% LMIs definition

LMIs = [];

% LMI condition to define a limit for the values of the gain matriz
xi = 1e-1;
for i = 1:r
    XI{i} = [xi*xi*eye(nu) Z{i};
            Z{i}'            eye(nx)];

    LMIs = [LMIs, XI{i} >= zeros(nu+nx)];
end

% Lyapunov function must be positive defined	
LMIs = [LMIs, X>=eps*eye(nx)];
	
% main H infinity LMI condition for each vertex of the politope
for i=1:r
    LMIs = [LMIs, Gamma{i} <= -eps*eye(size(Gamma{1}))];
end
    
%% Optimization
obj=[sigma];
options=sdpsettings('verbose',0,'warning',1,'solver','mosek','showprogress',0);
out.sol = optimize(LMIs,obj,options);
%warning('off','YALMIP: strict');
out.p=min(checkset(LMIs));
out.t=out.sol.solvertime;
out.variables=size(getvariables(LMIs),2);
	
out.LMIs=LMIs;

%% Output definition   
if out.p>0 || out.sol.problem==0

    P = inv(value(X));
    P = eye(nx)/(value(X));

    for i=1:r
        Zx{i} = value(Z{i});
        K{i} = value(Z{i})*P;
    end

    out.feas=1;
    out.P=P;
    out.Z = Zx;
	out.gamma=sqrt(value(sigma));
	out.K=K;
else
	out.feas=0;
end

end