function [A,Bu,Bw,Cz,Dzu,Dzw] = system_definition(Hi_min,Hi_max,Hj_min,Hj_max,Hk_min,Hk_max,Kv,g)

A{1} = [0        1        0        0        0      0   0    0;
        Hi_min   0        Kv*g     0        Hk_min 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_min   0        0        0        Hj_min 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];
            

A{2} = [0        1        0        0        0      0   0    0;
        Hi_min   0        Kv*g     0        Hk_max 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_max   0        0        0        Hj_min 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];
            

A{3} = [0        1        0        0        0      0   0    0;
        Hi_min   0        Kv*g     0        Hk_min 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_min   0        0        0        Hj_max 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];
            

A{4} = [0        1        0        0        0      0   0    0;
        Hi_min   0        Kv*g     0        Hk_max 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_max   0        0        0        Hj_max 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];
            

A{5} = [0        1        0        0        0      0   0    0;
        Hi_max   0        Kv*g     0        Hk_min 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_min   0        0        0        Hj_min 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];

A{6} = [0        1        0        0        0      0   0    0;
        Hi_max   0        Kv*g     0        Hk_max 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_max   0        0        0        Hj_min 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];
            

A{7} = [0        1        0        0        0      0   0    0;
        Hi_max   0        Kv*g     0        Hk_min 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_min   0        0        0        Hj_max 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];
            

A{8} = [0        1        0        0        0      0   0    0;
        Hi_max   0        Kv*g     0        Hk_max 0   0    0;
        0        0        0        1        0      0   0    0;
        0        0        0        0        0      0   0    0;
        0        0        0        0        0      1   0    0;
        Hk_min   0        0        0        Hj_max 0   Kv*g 0; 
        0        0        0        0        0      0   0    1;
        0        0        0        0        0      0   0    0];

Bu = [0 0;
      0 0;
      0 0;
      1 0;
      0 0;
      0 0;
      0 0;
      0 1];

Bw = [0 0;
      1 0;
      0 0;
      0 0;
      0 0;
      0 1;
      0 0;
      0 0];

Cz = eye(8);
Dzu = zeros(8,2);
Dzw = zeros(8,2);

end
