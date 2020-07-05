from sympy import *
import pickle
from multiprocessing import Pool

T_raw_file = open('T_raw.pkl', 'rb')
T_raw = pickle.load(T_raw_file)

T_sp = Matrix(T_raw)
T_trans = T_sp.T.tolist()
# print(T_trans.shape)
T = []
for i in range(len(T_trans)):
    print(f'simplifying {i+1}/{len(T_trans)}')
    T.append(simplify(sum(T_trans[i])))


T_final_raw = open('T_final.pkl', 'wb')
pickle.dump(T, T_final_raw)


