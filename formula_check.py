from sympy import *

T, M, x, x_f, V, omega, c, x_G = symbols('T, M, x, x_f, V, omega, c, x_G')

dT = M/c * Rational(1/2) * (V+(x-x_f)*omega)**2

T = integrate(dT, (x, 0, c))

print(T.expand())

Tt_T = Rational(1/2)*M*(omega*(x_f-x_G))**2 + Rational(1/2)*M/c*integrate((omega*(x-x_G))**2, (x, 0, c)) + Rational(1/2)*M*V**2 + M * omega*(x_G-x_f)*V

print(Tt_T.subs({x_G:c/2}).expand())

Tt = Rational(1/2)*M*(omega*(x_f-x_G))**2 + Rational(1/2)*M/c*integrate((omega*(x-x_G))**2, (x, -c/2, c/2)) + Rational(1/2)*M*V**2 + M * omega*(x_G-x_f)*V

print(Tt.subs({x_G:0}).expand())

gamma = symbols('gamma')

A = Matrix([[cos(gamma), -sin(gamma)],[sin(gamma), cos(gamma)]])
A_ = simplify(A.inv())
print(A_)