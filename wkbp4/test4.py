from sympy import *
import matplotlib.pyplot as plt
import numpy as np
import cmath as cm
from mpmath import mp
import mpmath
from scipy import  integrate
import scipy.optimize as sopt
import mpmath
from mpmath import mp
mp.dps=100
g=1e-4
E=2

coefs=[-1j*g,0,0,1,0,-E]
rtAll=np.roots(coefs)
rtAll=sorted(rtAll,key=np.angle,reverse=True)
angsAll=[np.angle(elem)/np.pi for elem in rtAll]
# print(rtAll)
# print(angsAll)
# print("g^{3/2}*E="+str(g**(2/3)*np.abs(E)))
def P(a):
    return 2*a-5*1j*a**4

def pert(a):
    return 1/P(a)*g**(1/3)*E

a01=np.exp(-1j*1/6*np.pi)
a02=np.exp(-1j*5/6*np.pi)

w1=a01*g**(-1/3)+pert(a01)
w2=a02*g**(-1/3)+pert(a02)
def f3(y):
    return 1j/6*y**6-1/3*y**3


y2=a02
y1=a01
I3=f3(y1)-f3(y2)
I2=y1-y2



sVal=1

def f4(y):
    return sVal*(1j*y**5-y**2)**(1/2)
def f5(y):
    return sVal*(1j*y**5-y**2)**(-1/2)
def I4(y2,y1):
    c1=np.real(y1)
    d1=np.imag(y1)

    c2=np.real(y2)
    d2=np.imag(y2)
    dx=1e-3
    N=int(np.abs(y2-y1)/dx)
    zAll=np.linspace(start=y2,stop=y1,num=N+1)
    oddVals4=[f4(zAll[j]) for j in range(1,N,2)]
    evenVals4=[f4(zAll[j]) for j in range(2,N,2)]
    return 1/3*(c1-c2)/N*(1+1j*(d1-d2)/(c1-c2))*(f4(zAll[0])+4*sum(oddVals4)+2*sum(evenVals4)+f4((zAll[N])))


def I5(y2,y1):
    lmd=0.1
    y2p=lmd*(y1-y2)+y2
    y1p=(1-lmd)*(y1-y2)+y2

    c1p=np.real(y1p)
    d1p=np.imag(y1p)

    c2p=np.real(y2p)
    d2p=np.imag(y2p)
    dx=1e-3
    N=int(np.abs(y1p-y2p)/dx)
    zAll=np.linspace(start=y2p,stop=y1p,num=N+1)
    oddVals5=[f5(zAll[j]) for j in range(1,N,2)]
    evenVals=[f5(zAll[j]) for j in range(2,N,2)]
    return 1/3*(c1p-c2p)/N*(1+1j*(d1p-d2p)/(c1p-c2p))*(f5(zAll[0])+4*sum(oddVals5)+2*sum(evenVals)+f5(zAll[N]))

n=100
Eest=2*((n+1/2)*np.pi-g**(-2/3)*I4(y2,y1))/I5(y2,y1)

def f(z, g, E):
    '''
    :param g: const
    :param z: point on x2x1
    :param E: trial eigenvalue
    :return: f value
    '''
    return (1j * g * z ** 5 - z ** 2 + E) ** (1 / 2)

def integralQuadrature(g,E,x1,x2):
    '''

    :param g: const
    :param E: trial eigenvalue
    :param x1: ending point
    :param x2: starting point
    :return:
    '''
    a1 = np.real(x1)
    b1 = np.imag(x1)

    a2 = np.real(x2)
    b2 = np.imag(x2)


    slope = (b1 - b2) / (a1 - a2)
    gFunc=lambda y:f(y+1j*(slope*(y-a2)+b2),g,E)
    return (1+1j*slope)*mpmath.quad(gFunc,[a2,a1])


def eqn(EIn,*data):
    '''

    :param EIn: trial eigenvalue, in the form of [re, im]
    :return:
    '''

    n,g=data

    E = EIn[0] + 1j * EIn[1]
    # print(E)

    x1, x2 = retX1X2(g,E)

    # dx = 1e-4
    # N = int(np.abs(x2 - x1) / dx)
    rst =integralQuadrature(g,E,x1,x2) - (n+1/2) * np.pi
    return np.real(rst), np.imag(rst)
def retX1X2(g, E):
    '''
    :param g: const
    :param E: trial eigenvalue
    :return: x1 and x2 of
    '''

    # E=float(mpmath.re(E))+1j*float(mpmath.im(E))
    # coefs = [-1j * g, 0, 0, 1, 0, -E]
    # print(E)

    # rootsAll = np.roots(coefs)
    # # print(rootsAll)
    # rootsAll = sorted(rootsAll, key=np.angle, reverse=True)
    # x1Val = rootsAll[3]
    # x2Val = rootsAll[4]
    a01=np.exp(-1j*np.pi/6)
    a02=np.exp(-1j*5/6*np.pi)
    y1=a01+1/P(a01)*g**(2/3)*E
    y2=a02+1/P(a02)*g**(2/3)*E
    x1Val=g**(-1/3)*y1
    x2Val=g**(-1/3)*y2

    # print(x1Val)
    return x1Val, x2Val

def computeOneSolution(inData):
    '''

    :param inData: [n,g]
    :return: [n,g, re(E), im(E)]
    '''
    n,g=inData


    eVecTmp=sopt.fsolve(eqn, [np.real(Eest), np.imag(Eest)],args=inData,maxfev=1000,xtol=1e-6)
    return [n,g,eVecTmp[0],eVecTmp[1]]


dataTmp=(n,g)
ret=computeOneSolution(dataTmp)
print("Eest = "+str(Eest))
print("|g^{2/3}*Eest| = "+str(np.abs(g**(2/3)*Eest)))
print("numerical solution = "+str(ret))
# lmdAll=np.linspace(start=0,stop=1,num=100)
# rts=[]
# intStart=a02
# intEnd=a01
# valsAll=[]
# for l in lmdAll:
#     xTmp=intStart+l*(intEnd-intStart)
#     yTmp=xTmp
#     valTmp=1j*yTmp**5-yTmp**2
#     prodTmp=g**(2/3)*E
#
#     rts.append(np.abs(valTmp/1))
#
#
#
# plt.figure()
# plt.plot(lmdAll,rts)
# plt.show()