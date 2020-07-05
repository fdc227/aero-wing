from sympy import *

def Cos(x):
    out = 1 - Rational(1/2)*x**2 + Rational(1/24)*x**4 - Rational(1/720)*x**6
    return out

def Sin(x):
    out = x - Rational(1/6)*x**3 + Rational(1/120)*x**5 - Rational(1/5040)*x**7
    return out