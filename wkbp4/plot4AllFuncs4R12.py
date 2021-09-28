import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Pool
import scipy.special as sspecial
import mpmath
from mpmath import mp
mp.dps=100

# this script combines calculations from WKB
# and shooting as a verification for potential 4, V(x)=x^{2}-igx^{5}
# in region I-II
# WKB part

def P(a):
    return 2*a-5*1j*a**4
def retX1X2(g, E):
    '''
    :param g: const
    :param E: trial eigenvalue
    :return: x1 and x2 of
    '''
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll = np.roots(coefs)
    # print(rootsAll)
    rootsAll = sorted(rootsAll, key=np.angle, reverse=True)
    elemInMinusPi=[]
    for rtTmp in rootsAll:
        if np.angle(rtTmp)>-np.pi and np.angle(rtTmp)<0:
            elemInMinusPi.append(rtTmp)
    if len(elemInMinusPi)==2:
        x1Val = rootsAll[3]
        x2Val = rootsAll[4]
    elif len(elemInMinusPi)==3:
        x1Val=rootsAll[2]
        x2Val=rootsAll[3]
    ################################
    # a01 = np.exp(-1j * np.pi / 6)
    # a02 = np.exp(-1j * 5 / 6 * np.pi)
    # y1 = a01 + 1 / P(a01) * g ** (2 / 3) * E
    # y2 = a02 + 1 / P(a02) * g ** (2 / 3) * E
    # x1Val = g ** (-1 / 3) * y1
    # x2Val = g ** (-1 / 3) * y2
    # print(x1Val)
    return x1Val, x2Val

def retX1X2New(g, E):
    '''

    :param g: const
    :param E: trial eigenvalue
    :return: x1 and x2 of region I-II
    '''
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll = np.roots(coefs)
    # print(rootsAll)
    rtsSmallerThanMinusPiOver2=[]
    rtsLargerThanMinusPiOver2=[]
    for rtTmp in rootsAll:
        if np.angle(rtTmp)<-np.pi/2:
            rtsSmallerThanMinusPiOver2.append(rtTmp)
        else:
            rtsLargerThanMinusPiOver2.append(rtTmp)
    rtsSmallerThanMinusPiOver2=sorted(rtsSmallerThanMinusPiOver2,key=np.angle)
    rtsLargerThanMinusPiOver2=sorted(rtsLargerThanMinusPiOver2,key=np.angle)
    x2=rtsSmallerThanMinusPiOver2[-1]
    x1=rtsLargerThanMinusPiOver2[0]
    return x1,x2
def returnFivePairsOfRoots(g,E):
    '''

    :param g: const
    :param E: trial eigenvalue
    :return: an adjacent pair of roots, the first has smaller angle than the second root,
    the order is [x2, x1]
    '''
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll = np.roots(coefs)
    rootsSortedByAngle=sorted(rootsAll,key=np.angle)
    rst=[]
    for j in range(0,len(rootsSortedByAngle)):
        rst.append([rootsSortedByAngle[j],rootsSortedByAngle[(j+1)%len(rootsSortedByAngle)]])
    return rst
def f(z, g, E):
    '''
    :param g: const
    :param z: point on x2x1
    :param E: trial eigenvalue
    :return: f value
    '''
    return (1j * g * z ** 5 - z ** 2 + E) ** (1 / 2)
def fBranchOther(z,g,E):
    '''

    :param z: point on x2x1
    :param g: const
    :param E: trial eigenvalue
    :return: f value on another branch
    '''
    return -(1j * g * z ** 5 - z ** 2 + E) ** (1 / 2)

def simpInt(g, E, N, x1, x2):
    '''
    :param E: trial eigenvalue
    :param N: number of subintervals
    :param x1:
    :param x2:
    :return: simpson integral
    '''
    a1 = np.real(x1)
    b1 = np.imag(x1)

    a2 = np.real(x2)
    b2 = np.imag(x2)

    deltaA = (a1 - a2) / N
    slope = (b1 - b2) / (a1 - a2)

    zAll = []

    for j in range(0, N + 1):
        # j=0,1,...,N
        aj = deltaA * j + a2
        bj = slope * (aj - a2) + b2
        zAll.append(aj + 1j * bj)
    fValsOdd = [f(zAll[j], g, E) for j in range(1, N, 2)]
    fValsEven = [f(zAll[j], g, E) for j in range(2, N, 2)]

    return 1 / 3 * deltaA * (1 + 1j * slope) * (f(zAll[0], g, E)
                                                + 4 * sum(fValsOdd)
                                                + 2 * sum(fValsEven)
                                                + f(zAll[N], g, E)
                                                )

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

def integralQuadratureAnotherBranch(g,E,x1,x2):
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
    gFunc = lambda y: fBranchOther(y + 1j * (slope * (y - a2) + b2), g, E)
    return (1 + 1j * slope) * mpmath.quad(gFunc, [a2, a1])



def eqn(EIn,*data):
    '''

    :param EIn: trial eigenvalue, in the form of [re, im]
    :return:
    '''

    n,g=data

    E = EIn[0] + 1j * EIn[1]
    # print("current E="+str(E))
    # x1, x2 = retX1X2(g,E)
    x1, x2 = retX1X2New(g,E)
    # dx = 1e-4
    # N = int(np.abs(x2 - x1) / dx)
    intVal=integralQuadrature(g,E,x1,x2)
    # intEst=g**(-2/3)*2/5*np.cos(3/10*np.pi)*(g**(2/3)*E)**(7/10)*sspecial.beta(1/5,3/2)
    # print("intVal="+str(intVal))
    # print("intEst="+str(intEst))
    rst =intVal - (n+1/2) * np.pi
    return np.real(rst), np.imag(rst)

def eqnFivePairs(EIn,*data):
    """

    :param EIn: trial eigenvalue, in the form of [re, im]
    :param data:
    :return:
    """
    n,g =data
    E = EIn[0] + 1j * EIn[1]
    pairsAll=returnFivePairsOfRoots(g,E)
    rstValsAll=[]
    for onePair in pairsAll:
        x2Tmp,x1Tmp=onePair
        intValTmp=integralQuadrature(g,E,x1Tmp,x2Tmp)
        rstTmp=intValTmp- (n+1/2) * np.pi
        rstValsAll.append(rstTmp)
    rstValsAll=sorted(rstValsAll,key=np.abs)
    root0=rstValsAll[0]
    return np.real(root0), np.imag(root0)


def computeOneSolution(inData):
    '''

    :param inData: [n,g]
    :return: [n,g, re(E), im(E)]
    '''
    n,g=inData

    eVecTmp=sopt.fsolve(eqn, [np.abs(n + 0.5), 0],args=inData,maxfev=100,xtol=1e-2)

    return [n,g,eVecTmp[0],eVecTmp[1]]

#approximation of E
def f4(y):
    return (1j*y**5-y**2)**(1/2)
def f5(y):
    return (1j*y**5-y**2)**(-1/2)
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
def computeOneSolutionWithInit(inData):
    '''

    :param inData: [n,g,Eest]
    :return: [n,g,re(E), im(E)]
    '''
    n,g,Eest=inData

    eVecTmp=sopt.fsolve(eqn,[np.real(Eest),np.imag(Eest)],args=(n,g),maxfev=100,xtol=1e-3)
    # eVecTmp = sopt.fsolve(eqn, [np.abs(n+0.5), 0.01], args=(n, g), maxfev=100, xtol=1e-3)
    return [n,g,eVecTmp[0],eVecTmp[1]]


#-igx^5 dominant
def dom5E(n,g):
    rst=(
                5*(n+1/2)*np.pi*g**(1/5)/(
        2*np.cos(3/10*np.pi)*sspecial.beta(1/5,3/2)
    )
         )**(10/7)
    return rst
# a01=np.exp(-1j*1/6*np.pi)
# a02=np.exp(-1j*5/6*np.pi)
# g=0.1
# y2=a02
# y1=a01
# n=10
# Eest=2*((n+1/2)*np.pi-g**(-2/3)*I4(y2,y1))/I5(y2,y1)
# print(computeOneSolutionWithInit([n,g,Eest]))


