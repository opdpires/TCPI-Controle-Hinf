from itertools import combinations

import itertools
from itertools import permutations
 
# initialize lists
list_1 = ["Hi_min", "Hi_max"]
list_2 = ["Hj_min", "Hj_max"]
list_3 = ["Hk_min", "Hk_max"]
list_4 = ["Hp_min", "Hp_max"]
list_5 = ["Hq_min", "Hq_max"]
 
res = list(itertools.product(list_1, list_2,list_3,list_4,list_5))

# print(res)

# print(len(res))
# print

for i,lista in enumerate(res):
    estrutura = f"""
A{ {i+1} } = [0        0        0        0        1 0 0 0;
        {lista[0]}   {lista[2]}   {lista[3]}   0        0 0 0 0;
        0        0        0        0        0 0 1 0;
        0        0        0        0        0 0 0 0;
        0        0        0        0        0 1 0 0;
        {lista[2]}   {lista[1]}   0        {lista[4]}   0 0 0 0; 
        0        0        0        0        0 0 0 1;
        0        0        0        0        0 0 0 0];
            """
    print(estrutura)
