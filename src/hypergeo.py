from scipy import special as sp
from numpy import finfo
import math

EPS = finfo(float).eps *1000

def hyper2F1(a, b, c, z):
        if z> 1:
                raise ValueError("z<=1 z= {}".format(z))
        if z == 1 or math.isclose(z, 1.0, rel_tol=EPS):
            if c>a+b:
                return sp.gamma(c) * sp.gamma(c-(a+b))/ (sp.gamma(c-a) * sp.gamma(c-b))
            else:
                raise ValueError("if z=1, c must be larger than a+b, c:{}, a+b:{}".format(c, a+b))
        if -1<= z:
            return sp.hyp2f1(a,b,c,z)
        if z< -1:

            ab_abs = math.fabs(a-b)
            ab_int = int(ab_abs)

            if math.isclose(ab_abs, ab_int, rel_tol =  EPS):
                return (1-z)**(-a) * sp.hyp2f1(a, c-b, c ,z/(z-1))
                # same with '(1-z)**(-b)*sp.hyp2f1(c-a, b, c, z/(z-1))'
                # raise ValueError("a-b must not be in Z, |a-b|:{}".format(a-b))

            if c-a <0:
                    ca_abs = math.fabs(c-a)
                    ca_int = int(ca_abs)
                    if math.isclose(ca_abs, ca_int, rel_tol = EPS):
                        raise ValueError("c-a and c-b must not be negative integer c-a: {}, c-b:{}".format(c-a, c-b ))
            if c-b <0:
                cb_abs = math.fabs(c-b)
                cb_int = int(cb_abs)
                if math.isclose(cb_abs, cb_int, rel_tol= EPS):
                    raise ValueError("c-a and c-b must not be negative integer c-a: {}, c-b:{}".format(c-a, c-b ))

            w = 1/(1-z)

            c1 = (w**a) * (sp.gamma(c)*sp.gamma(b-a)/sp.gamma(b) * sp.gamma(c-a))
            c2 = (w**b) *(sp.gamma(c) * sp.gamma(a-b)/(sp.gamma(a) * sp.gamma(c-b)))

            return c1 *sp.hyp2f1(a, c-b, a-b+1, w) + c2* sp.hyp2f1(b, c-a, b-a+1, w) 

