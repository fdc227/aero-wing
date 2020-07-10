from sympy import *
from dot_product import *
from shape_gen import *
from torsion_shape_gen import torsion_shape_gen
from symbols_util import sym_suite, sym_str_gen, diff_wrapper, T_fusalage_gen, row_recduce_wrapper, U_wrapper
from multiprocessing import Pool
import pickle
from sympy.solvers.solveset import *
THREAD_NUM = 20

x, alpha, alpha_dt, gamma, gamma_dt, c, x_f, cg, qb, qb_dt, M = symbols('x, alpha, alpha_dt, gamma, gamma_dt, c, x_f, cg, qb, qb_dt, M')
X, Y, X_dt, Y_dt, qt, qt_dt = symbols('X, Y, X_dt, Y_dt, qt, qt_dt')
y, L, t = symbols('y, L, t')
M_fusalage = symbols('M_fusalage')
E, I, G, J = symbols('E, I, G, J')

qb0_dt_dt, qb0_dot_dt_dt, qt0_dt_dt, qb1_dt_dt, qb1_dot_dt_dt, qt1_dt_dt, qb2_dt_dt, qb2_dot_dt_dt, qt2_dt_dt, qb3_dt_dt, qb3_dot_dt_dt, qt3_dt_dt, qb4_dt_dt, qb4_dot_dt_dt, qt4_dt_dt, qb5_dt_dt, qb5_dot_dt_dt, qt5_dt_dt, qb6_dt_dt, qb6_dot_dt_dt, qt6_dt_dt, qb7_dt_dt, qb7_dot_dt_dt, qt7_dt_dt, qb8_dt_dt, qb8_dot_dt_dt, qt8_dt_dt, qb9_dt_dt, qb9_dot_dt_dt, qt9_dt_dt, qb10_dt_dt, qb10_dot_dt_dt, qt10_dt_dt, qb11_dt_dt, qb11_dot_dt_dt, qt11_dt_dt, qb12_dt_dt, qb12_dot_dt_dt, qt12_dt_dt, qb13_dt_dt, qb13_dot_dt_dt, qt13_dt_dt, qb14_dt_dt, qb14_dot_dt_dt, qt14_dt_dt, qb15_dt_dt, qb15_dot_dt_dt, qt15_dt_dt, qb16_dt_dt, qb16_dot_dt_dt, qt16_dt_dt, qb17_dt_dt, qb17_dot_dt_dt, qt17_dt_dt, qb18_dt_dt, qb18_dot_dt_dt, qt18_dt_dt, qb19_dt_dt, qb19_dot_dt_dt, qt19_dt_dt, qb20_dt_dt, qb20_dot_dt_dt, qt20_dt_dt = symbols('qb0_dt_dt, qb0_dot_dt_dt, qt0_dt_dt, qb1_dt_dt, qb1_dot_dt_dt, qt1_dt_dt, qb2_dt_dt, qb2_dot_dt_dt, qt2_dt_dt, qb3_dt_dt, qb3_dot_dt_dt, qt3_dt_dt, qb4_dt_dt, qb4_dot_dt_dt, qt4_dt_dt, qb5_dt_dt, qb5_dot_dt_dt, qt5_dt_dt, qb6_dt_dt, qb6_dot_dt_dt, qt6_dt_dt, qb7_dt_dt, qb7_dot_dt_dt, qt7_dt_dt, qb8_dt_dt, qb8_dot_dt_dt, qt8_dt_dt, qb9_dt_dt, qb9_dot_dt_dt, qt9_dt_dt, qb10_dt_dt, qb10_dot_dt_dt, qt10_dt_dt, qb11_dt_dt, qb11_dot_dt_dt, qt11_dt_dt, qb12_dt_dt, qb12_dot_dt_dt, qt12_dt_dt, qb13_dt_dt, qb13_dot_dt_dt, qt13_dt_dt, qb14_dt_dt, qb14_dot_dt_dt, qt14_dt_dt, qb15_dt_dt, qb15_dot_dt_dt, qt15_dt_dt, qb16_dt_dt, qb16_dot_dt_dt, qt16_dt_dt, qb17_dt_dt, qb17_dot_dt_dt, qt17_dt_dt, qb18_dt_dt, qb18_dot_dt_dt, qt18_dt_dt, qb19_dt_dt, qb19_dot_dt_dt, qt19_dt_dt, qb20_dt_dt, qb20_dot_dt_dt, qt20_dt_dt')

var_list = [qb0_dt_dt, qb0_dot_dt_dt, qt0_dt_dt, qb1_dt_dt, qb1_dot_dt_dt, qt1_dt_dt, qb2_dt_dt, qb2_dot_dt_dt, qt2_dt_dt, qb3_dt_dt, qb3_dot_dt_dt, qt3_dt_dt, qb4_dt_dt, qb4_dot_dt_dt, qt4_dt_dt, qb5_dt_dt, qb5_dot_dt_dt, qt5_dt_dt, qb6_dt_dt, qb6_dot_dt_dt, qt6_dt_dt, qb7_dt_dt, qb7_dot_dt_dt, qt7_dt_dt, qb8_dt_dt, qb8_dot_dt_dt, qt8_dt_dt, qb9_dt_dt, qb9_dot_dt_dt, qt9_dt_dt, qb10_dt_dt, qb10_dot_dt_dt, qt10_dt_dt, qb11_dt_dt, qb11_dot_dt_dt, qt11_dt_dt, qb12_dt_dt, qb12_dot_dt_dt, qt12_dt_dt, qb13_dt_dt, qb13_dot_dt_dt, qt13_dt_dt, qb14_dt_dt, qb14_dot_dt_dt, qt14_dt_dt, qb15_dt_dt, qb15_dot_dt_dt, qt15_dt_dt, qb16_dt_dt, qb16_dot_dt_dt, qt16_dt_dt, qb17_dt_dt, qb17_dot_dt_dt, qt17_dt_dt, qb18_dt_dt, qb18_dot_dt_dt, qt18_dt_dt, qb19_dt_dt, qb19_dot_dt_dt, qt19_dt_dt, qb20_dt_dt, qb20_dot_dt_dt, qt20_dt_dt]


T_raw = open('T_final_v2.pkl','rb')
T_package= pickle.load(T_raw)
T_list, var = T_package[0], T_package[1]
# print(T_list[1])
# print(var)


def lc_mt(i):
    term = T_list[i]
    out = linear_coeffs(term, *var_list)
    print(f'{i+1}/{len(T_list)} finished')
    return out


# m_raw = open('matrix_from_lc_mt.pkl', 'wb')
# pickle.dump(m_final, m_raw)

# A, b = linear_eq_to_matrix(T_list, *var)
# mm_raw = open('matrix_from_letm.pkl', 'wb')
# m_package = [A, b]
# pickle.dump(m_package, mm_raw)

###  Windows specific ######

if __name__ == "__main__":
    p = Pool(THREAD_NUM)
    r = [r for r in range(len(T_list))]
    m_final = p.map(lc_mt, r)
        
