from sympy import *
import numpy as np
import cmath as cm
import matplotlib.pyplot as plt

g=0.01

E=g**(-1/3)/100
coefs=[-1j*g,0,0,1,0,-E]

rts=np.roots(coefs)
anglesInPi=[cm.polar(elem)[1]/np.pi for elem in rts]
normAndAngInPi=[(cm.polar(elem)[0]/g**(-1/3),cm.polar(elem)[1]/np.pi) for elem in rts]


print(normAndAngInPi)

