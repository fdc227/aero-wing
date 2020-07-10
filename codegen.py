from sympy import *
from sympy.printing.ccode import C99CodePrinter
from sympy.printing.codeprinter import Assignment
import pickle
from symbols_util import row_recduce_wrapper
from multiprocessing import Pool
THREAD_NUM = 20

sym_list_raw = open('sym_standalone.pkl', 'rb')
out = pickle.load(sym_list_raw)
var_list_func_sym, var_list_func_dt_sym = out[0], out[1]

A_raw = open('A_b_U.pkl', 'rb')
A_list = pickle.load(A_raw)
A, b, U, RHS = A_list[0], A_list[1], A_list[2], A_list[3]
# print(U)

# X = MatrixSymbol('X', 63, 63)
# Y = MatrixSymbol('Y', 63, 1)

# print(Matrix(X).inv()[0,0])
# Z = X.inv()*Y
# print(Matrix(Z).shape)

# sub_exprs, simplified = cse(U)
# print(sub_exprs)
# print(simplified)

class CMatrixPrinter(C99CodePrinter):
    def _print_ImmutableDenseMatrix(self, expr):
        sub_exprs, simplified = cse(expr)
        lines = []
        for var, sub_expr in sub_exprs:
            lines.append('double ' + self._print(Assignment(var, sub_expr)))
        M = MatrixSymbol('M', *expr.shape)
        return '\n'.join(lines) + '\n' + self._print(Assignment(M, Matrix(simplified)))

p = CMatrixPrinter()
# print(p.doprint(A))
print(p.doprint(Matrix(RHS)))

# var_list = var_list_func_sym + var_list_func_dt_sym
# varstr = ''
# for i in range(len(var_list)):
#     varstr += f'double {var_list[i]} = var_list[{i}];\n' 
# print(varstr)

# print(A)