from sympy import *
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt

g=0.1

E=g**(-1/3)*1e8
def chooseX1(x):
    angTmp=np.angle(x)

    ret=(angTmp<=-1/4*np.pi and angTmp>=-1/2*np.pi)
    return ret
def chooseX2(x):
    angTmp=np.angle(x)
    ret=(angTmp<-1/2*np.pi and angTmp>=-3/4*np.pi)
    return ret
coefs=[-1j*g,0,0,1,0,-E]

rts=np.roots(coefs)
rts=sorted(rts,key=np.angle,reverse=True)
x2=rts[4]
x1=rts[1]
for tmp in rts:
    if chooseX1(tmp):
        x1=tmp
    if chooseX2(tmp):
        x2=tmp

print(chooseX1(x1))
print(chooseX2(x2))

print(np.abs(x1)/g**(-1/3))
print(np.abs(x2)/g**(-1/3))