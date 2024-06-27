from sympy import symbols, Poly

x = symbols('x')
f = (35/128)*x**9 - (180/128)*x**7 + (378/128)*x**5 - (420/128)*x**3 + (315/128)*x
g = (46623/1024)*x**9 - (113492/1024)*x**7 + (97015/1024)*x**5 - (34974/1024)*x**3 + (5850/1024)*x

h = f.subs(x, f.subs(x, f.subs(x, g.subs(x, g.subs(x, g)))))

# P = Poly(h, x)
print(h.subs(x, 0.001))