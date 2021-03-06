from sympy import *
from dot_product import *
from shape_gen import *
from torsion_shape_gen import torsion_shape_gen
from symbols_util import sym_suite, sym_str_gen, diff_wrapper, T_fusalage_gen, row_recduce_wrapper, U_wrapper
from multiprocessing import Pool
import pickle
THREAD_NUM = 20

x, alpha, alpha_dt, gamma, gamma_dt, c, x_f, cg, qb, qb_dt, M = symbols('x, alpha, alpha_dt, gamma, gamma_dt, c, x_f, cg, qb, qb_dt, M')
X, Y, X_dt, Y_dt, qt, qt_dt = symbols('X, Y, X_dt, Y_dt, qt, qt_dt')
y, L, t = symbols('y, L, t')
M_fusalage = symbols('M_fusalage')
E, I, G, J = symbols('E, I, G, J') # Young's modulus, second moment I_yy, shear modulus, polar moment I_p

parameters = {E:200*10**9, I:7.85*10**(-5), G:78*10**9, J:15.74*10**(-5), c:1, L:1, x_f:0.4, cg:0.5, M:251200, M_fusalage:251200} #assuming r=0.1, density=8*10**6
twoD_param = {X:0, Y:0, X_dt:0, Y_dt:0, gamma:0, gamma_dt:0}

replace_dict = {**parameters, **twoD_param}

mm_raw = open('matrix_from_letm.pkl', 'rb')
MM = pickle.load(mm_raw)
U_raw = open('U_final.pkl', 'rb')
U = Matrix(pickle.load(U_raw))

A, b = MM[0], MM[1]
A_sym = A.xreplace(replace_dict)
b_sym = b.xreplace(replace_dict)
U_sym = U.xreplace(replace_dict)
rhs = (b_sym - U_sym).tolist()

# print(A_sym)
# print(b_sym)
# print(U)
# print(rhs)

def rhs_simplify(i):
    print(f'simplifying {i+1}/{len(rhs)}')
    out = simplify(rhs[i][0])
    print(f'finished {i+1}/{len(rhs)}')
    return out

p = Pool(THREAD_NUM)
r = [r for r in range(len(rhs))]
rhs_s = p.map(rhs_simplify, r)

A_raw = open('A_b_U.pkl', 'wb')
pickle.dump([A_sym, b_sym, U_sym, rhs_s], A_raw)
