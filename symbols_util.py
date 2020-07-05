from sympy import *
t = symbols('t')

def sym_suite(var_list):
    var_list_func = []
    var_list_func_sym = []
    var_list_func_dt = []
    var_list_func_dt_sym = []
    var_list_func_dt_dt = []
    var_list_func_dt_dt_sym = []
    for sym in var_list:
        globals()[sym] = Function(sym)(t)
        var_list_func.append(Function(sym)(t))
        var_list_func_sym.append(symbols(sym))
        var_list_func_dt.append(diff(Function(sym)(t), t))
        var_list_func_dt_sym.append(symbols(sym+'_dt'))
        var_list_func_dt_dt.append(diff(diff(Function(sym)(t), t), t))
        var_list_func_dt_dt_sym.append(symbols(sym+'_dt_dt'))

    replace_dict = {}
    for i in range(len(var_list)):
        replace_dict[var_list_func[i]] = var_list_func_sym[i]
        replace_dict[var_list_func_dt[i]] = var_list_func_dt_sym[i]
        replace_dict[var_list_func_dt_dt[i]] = var_list_func_dt_dt_sym[i]
    
    return var_list_func, var_list_func_sym, var_list_func_dt, var_list_func_dt_sym, var_list_func_dt_dt, var_list_func_dt_dt_sym, replace_dict

def sym_str_gen(sym_list, insert_position_list, num_range_list):
    out = []
    for i in range(num_range_list[0],num_range_list[1]):
        for j in range(len(sym_list)):
            out.append(sym_list[j][0:insert_position_list[j]]+str(i)+sym_list[j][insert_position_list[j]:])
    return out

# def diff_wrapper(T_list, diff_var_list, integrate_var_list, replace_dict): #integrate_var_dict = [[y, 0, L]]
#     global T_ode_gen
#     def T_ode_gen(i):
#         out = []
#         n = len(diff_var_list)
#         for j in range(n):
#             local_diff = diff(diff(T_list[i], diff_var_list[j]), t).expand()
#             print(f'j:{j+1}/{n} i:{i} differentiated')
#             print(local_diff)
#             local_int = 0
#             for term in integrate_var_list:
#                 local_int += integrate(local_diff, (term[0], term[1], term[2]))
#             print(f'j:{j+1}/{n} i:{i} integrated')
#             local_int = local_int.xreplace(replace_dict)
#             out.append(local_int)
#         return out
#     return T_ode_gen

def diff_wrapper(T_list, diff_var_list, integrate_var_list, replace_dict): #integrate_var_list = [y, 0, L]
    global T_ode_gen
    def T_ode_gen(i):
        out = []
        n = len(diff_var_list)
        expr = T_list[i]
        local_int = integrate(expr, (integrate_var_list[0], integrate_var_list[1], integrate_var_list[2]))
        print(f'i:{i} integrated')
        for j in range(n):
            local_diff = diff(diff(local_int, diff_var_list[j]), t)
            print(f'j:{j+1}/{n} i:{i} differentiated')
            out.append(local_diff.xreplace(replace_dict))
        return out
    return T_ode_gen

def U_wrapper(U_list, diff_var_list, integrate_var_dict, replace_dict):
    global U_ode_gen
    def U_ode_gen(i):
        expr = U_list[i]
        n = len(diff_var_list)
        out = []
        local_int = 0
        for term in integrate_var_dict:
            local_int += integrate(expr, (term[0], term[1], term[2]))
            print(f'i:{i+1} of U integrated')
        for j in range(n):
            local_diff = diff(local_int, diff_var_list[j])
            out.append(local_diff.xreplace(replace_dict))
            print(f'j:{j+1}/{n} i:{i+1} of U differentiated')
        return out
    return U_ode_gen
        

def T_fusalage_gen(T, diff_var_list, replace_dict):
    out = []
    n = len(diff_var_list)
    for i in range(n):
        local_diff = diff(diff(T, diff_var_list[i]),t)
        out.append(local_diff.xreplace(replace_dict))
    return out

def row_recduce_wrapper(T_ode_list_transpose):
    global T_final_gen
    def T_final_gen(i):
        n = len(T_ode_list_transpose)
        out = simplify(sum(T_ode_list_transpose[i]))
        print(f'{i+1}/{n} row reduced')
        return out
    return T_final_gen

def list_subs_wrapper(T, subs_symbols, subs_lists):
    global T_substitute
    def T_substitute(i):
        out = T.xreplace({subs_symbols[0]:subs_lists[0], subs_symbols[1]:subs_lists[1], 
                          subs_symbols[2]:subs_lists[2], subs_symbols[3]:subs_lists[3]}).expand()
        print(f'T substitute {i+1}')
        return out
    return T_substitute
        
