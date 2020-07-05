from sympy import *
from dot_product import *
from shape_gen import *
from torsion_shape_gen import torsion_shape_gen
from symbols_util import sym_suite, sym_str_gen, diff_wrapper, T_fusalage_gen, row_recduce_wrapper, U_wrapper, list_subs_wrapper
from taylor_expansion import Sin, Cos
from multiprocessing import Pool
import pickle
THREAD_NUM = 20

x, alpha, alpha_dt, gamma, gamma_dt, c, x_f, cg, qb, qb_dt, M = symbols('x, alpha, alpha_dt, gamma, gamma_dt, c, x_f, cg, qb, qb_dt, M')
X, Y, X_dt, Y_dt, qt, qt_dt = symbols('X, Y, X_dt, Y_dt, qt, qt_dt')
y, L, t = symbols('y, L, t')
M_fusalage = symbols('M_fusalage')
E, I, G, J = symbols('E, I, G, J')

T_rotation = Rational(1/2) * M/c * integrate(((x-cg)*(alpha_dt+gamma_dt))**2, (x, 0, c))

# print(simplify(T_rotation.expand()))

A = Matrix([[Cos(gamma), -Sin(gamma)],[Sin(gamma), Cos(gamma)]])

A_ = Matrix([[Cos(gamma), Sin(gamma)],[-Sin(gamma), Cos(gamma)]])

# A_ = A.inv()

print(simplify(A_.expand()))

V_local = Matrix([[0],[qb_dt]]) + alpha_dt * (x_f-cg) * Matrix([[-Sin(alpha)],[Cos(alpha)]])

V_global = A_ * Matrix([[X_dt],[Y_dt]]) + V_local + Matrix([[-qb_dt],[cg-x_f]])

T_trans = Rational(1/2) * M * Dot(V_global, V_global)

T = (T_rotation + T_trans).subs({alpha:qt, alpha_dt:qt_dt})

# print(simplify(T))

shape_func = shape_gen(4)
for i in range(4):
    shape_func[i] = shape_func[i].subs({x:y})
torsion_func = torsion_shape_gen()
# print(torsion_func)
# print(shape_func)

variables_str_list = sym_str_gen(['qb', 'qb_dot', 'qt'], [2,2,2], [0, 21])
print(variables_str_list)

var_list_func, var_list_func_sym, var_list_func_dt, var_list_func_dt_sym, var_list_func_dt_dt, var_list_func_dt_dt_sym, replace_dict = sym_suite(variables_str_list)
# print(replace_dict)

qb_section = []
qb_dt_section = []
qt_section = []
qt_dt_section = []

for i in range(20):
    qb_section.append(var_list_func[3*i]*shape_func[0] + var_list_func[3*i+1]*shape_func[1] + var_list_func[3*i+3]*shape_func[2] + var_list_func[3*i+4]*shape_func[3])
    qb_dt_section.append(var_list_func_dt[3*i]*shape_func[0] + var_list_func_dt[3*i+1]*shape_func[1] + var_list_func_dt[3*i+3]*shape_func[2] + var_list_func_dt[3*i+4]*shape_func[3])
    qt_section.append(var_list_func[3*i+2]*torsion_func[0] + var_list_func[3*i+5]*torsion_func[1])
    qt_dt_section.append(var_list_func_dt[3*i+2]*torsion_func[0] + var_list_func_dt[3*i+5]*torsion_func[1])

# print(qt_dt_section)
# T_list = []
# for i in range(10):

########## T and U list definitions ###########
T_list = []
# T_expanded = simplify(T.expand())
for i in range(20):
    print(f'T substituting {i+1}/20')
    expr = T.xreplace({qb:qb_section[i],qb_dt:qb_dt_section[i], qt:qt_section[i], qt_dt:qt_dt_section[i]})
    print(f'{i+1}/20 finishd')
    T_list.append(expr)

# T_substitute = list_subs_wrapper(T, [qb, qb_dt, qt, qt_dt], [qb_section, qb_dt_section, qt_section, qt_dt_section])
# Pt = Pool(THREAD_NUM)
# rt = [r for r in range(20)]
# T_list = Pt.map(T_substitute, rt)

# print(T_list[0])

T_ode_gen = diff_wrapper(T_list, var_list_func_dt, [y, 0, L], replace_dict)

U_list = []
for i in range(20):
    expr = Rational(1/2)*E*I*diff(diff(qb_section[i], y), y)**2 + Rational(1/2)*G*J*diff(qt_section[i], y)**2
    U_list.append(expr)

U_ode_gen = U_wrapper(U_list, var_list_func, [[y, 0, L]], replace_dict)

#//////////////////////////////////////////////#

################ T gen #########################

p = Pool(THREAD_NUM)
R = [r for r in range(20)]
T_ode_list = p.map(T_ode_gen, R)

T_fusalage = T.subs({M:M_fusalage, qb:var_list_func[30], qb_dt: var_list_func_dt[30], qt:var_list_func[32], qt_dt:var_list_func_dt[32]})
T_fusalage_list = T_fusalage_gen(T_fusalage, var_list_func_dt, replace_dict)
print(T_fusalage_list)

T_ode_list.append(T_fusalage_list)

T_ode_list_transpose = Matrix(T_ode_list).T.tolist()
T_fianl_gen = row_recduce_wrapper(T_ode_list_transpose)
Q = Pool(THREAD_NUM)
Rq = [r for r in range(len(T_ode_list_transpose))]
T_final = Q.map(T_fianl_gen, Rq)

#////////////////////////////////////////////////#

############### U gen ############################

Pp = Pool(THREAD_NUM)
U_ode_list = Pp.map(U_ode_gen, R)
U_ode_list_transpose = Matrix(U_ode_list).T.tolist()
U_final_gen = row_recduce_wrapper(U_ode_list_transpose)
Qu = Pool(THREAD_NUM)
Ru = [r for r in range(len(U_ode_list_transpose))]
U_final = Qu.map(U_final_gen, Ru)

#////////////////////////////////////////////////#

T_package = [T_final, var_list_func_dt_dt_sym]
T_final_raw = open('T_final_v2.pkl', 'wb')
pickle.dump(T_package, T_final_raw)

U_final_raw = open('U_final.pkl', 'wb')
pickle.dump(U_ode_list, U_final_raw)


