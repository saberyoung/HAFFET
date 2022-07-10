def poly6(x, a, b, c, d, e, f, g):
    return a*x**6+b*x**5+c*x**4+d*x**3+e*x**2+f*x+g

def poly5(x, a, b, c, d, e, f):
    return a*x**5+b*x**4+c*x**3+d*x**2+e*x+f

def poly4(x, a, b, c, d, e):
    return a*x**4+b*x**3+c*x**2+d*x+e

def poly3(x, a, b, c, d):
    return a*x**3+b*x**2+c*x+d

def poly2(x, a, b, c):
    return a*x**2+b*x+c

def linear(x, a, b):
    return a*x+b
