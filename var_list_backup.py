from sympy import *
from symbols_util import sym_suite, sym_str_gen
import pickle

variables_str_list = sym_str_gen(['qb', 'qb_dot', 'qt'], [2,2,2], [0, 21])

var_list_func, var_list_func_sym, var_list_func_dt, var_list_func_dt_sym, var_list_func_dt_dt, var_list_func_dt_dt_sym, replace_dict = sym_suite(variables_str_list)

out = [var_list_func_sym, var_list_func_dt_sym]

sym_raw = open('sym_standalone.pkl', 'wb')
pickle.dump(out, sym_raw)



