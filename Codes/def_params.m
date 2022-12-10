function [Kv,g,Hi_min,Hi_max,Hj_min,Hj_max,Hk_min,Hk_max] = def_params()

% Ball characteristics values
m = 0.13; % Mass in Kg
R = 0.03; % Radius in meters
Ib = 2/5 * m * R^2; % Inertial Momentum
g = 9.8; % Gravity

Kv = m/(m+Ib/(R^2)); % Constant

% extremal values definition
dtheta1_min = -1;     dtheta1_max = +1;
dtheta2_min = -1;     dtheta2_max = +1;

dtheta1_2_min = 0;                          dtheta1_2_max = dtheta1_max^2;
dtheta2_2_min = 0;                          dtheta2_2_max = dtheta2_max^2;
dth1th2_min = -dtheta1_max*dtheta2_max;     dth1th2_max = +dtheta1_max*dtheta2_max;

Hi_min = Kv*dtheta1_2_min;             Hi_max = Kv*dtheta1_2_max;
Hj_min = Kv*dtheta2_2_min;             Hj_max = Kv*dtheta2_2_max; 
Hk_min = Kv*dth1th2_min;               Hk_max = Kv*dth1th2_max;

end