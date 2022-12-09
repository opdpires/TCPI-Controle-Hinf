% Ball characteristics values
m = 0.13; % Mass in Kg
R = 0.03; % Radius in meters
Ib = 2/5 * m * R^2; % Inertial Momentum
g = 9.8; % Gravity

Kv = m/(m+Ib/(R^2)); % Constant

% Control values
x0 = [0.05; 0; 0; 0; -0.05; 0; 0; 0]; % Initial condition
tspan = 0:0.1:100;
ref_px = @(t)0.01*sin(t); % px reference
ref_py = @(t)0.01*cos(t); % py reference

w = [0;0]; % Disturbance

% extremal values definition
dtheta1_min = -1;     dtheta1_max = +1;
dtheta2_min = -1;     dtheta2_max = +1;

dtheta1_2_min = -1;                          dtheta1_2_max = dtheta1_max^2;
dtheta2_2_min = -1;                          dtheta2_2_max = dtheta2_max^2;
dth1th2_min = -dtheta1_max*dtheta2_max;     dth1th2_max = +dtheta1_max*dtheta2_max;

Hi_min = Kv*dtheta1_2_min;             Hi_max = Kv*dtheta1_2_max;
Hj_min = Kv*dtheta2_2_min;             Hj_max = Kv*dtheta2_2_max; 
Hk_min = Kv*dth1th2_min;               Hk_max = Kv*dth1th2_max;
