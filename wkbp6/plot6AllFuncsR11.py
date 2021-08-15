import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Pool
from scipy import integrate
import mpmath
# this script combines calculations from WKB
# and shooting as a verification for potential 6, V(x)=ix^{3}-gx^{4}
# in region I-I
# WKB part

def retX1X2(g,E):
    '''

    :param g: const
    :param E: trial eigenvalue
    :return: x1 and x2 of gx^{4}-ix^{3}+e=0
    '''
    coefs = [-g, 1j, 0, 0, -E]
    rootsAll = np.roots(coefs)
    # print(rootsAll)
    rootsAll = sorted(rootsAll, key=np.angle, reverse=True)
    x1=rootsAll[1]
    x2=rootsAll[0]
    return x1, x2


def f(z,g,E):
    '''

    :param z: point on x2x1
    :param g: const
    :param E: trial eigenvalue
    :return: f value
    '''
    return (g*z**4-1j*z**3+E)**(1/2)


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

def eqn(EIn,*data):
    '''

    :param EIn: trial eigenvalue, in the form of [re, im]
    :return:
    '''
    # print("entering eqn")
    n,g=data

    E = EIn[0] + 1j * EIn[1]
    x1, x2 = retX1X2(g,E)

    # dx = 1e-4
    # N = int(np.abs(x2 - x1) / dx)

    rst = integralQuadrature(g,E,x1,x2) - (n+1/2) * np.pi
    # print(rst)
    return np.real(rst), np.imag(rst)
def computeOneSolution(inData):
    '''

    :param inData: [n,g]
    :return: [n,g, re(E), im(E)]
    '''
    n,g=inData
    # print("entering")

    eVecTmp=sopt.fsolve(eqn, [np.abs(n + 0.5), np.abs(n+0.5)],args=inData,maxfev=1000,xtol=1e-1)

    return [n,g,eVecTmp[0],eVecTmp[1]]

#WKB part ends here