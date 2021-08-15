import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from datetime import datetime
from multiprocessing import Pool
import mpmath
from mpmath import mp
mp.dps=100
# this script combines calculations from WKB
# and shooting as a verification for potential 6, V(x)=ix^{3}-gx^{4}
# in region I-II
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
    x1=rootsAll[2]
    x2=rootsAll[3]
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

    n,g=data

    E = EIn[0] + 1j * EIn[1]
    x1, x2 = retX1X2(g,E)
    # dx = 1e-4
    # N = int(np.abs(x2 - x1) / dx)
    rst = integralQuadrature(g,E,x1,x2) - (n+1/2) * np.pi
    return np.real(rst), np.imag(rst)


def computeOneSolution(inData):
    '''

    :param inData: [n,g]
    :return: [n,g, re(E), im(E)]
    '''
    n,g=inData

    eVecTmp=sopt.fsolve(eqn, [np.abs(n + 0.5), 0],args=inData,maxfev=100,xtol=1e-7)

    return [n,g,eVecTmp[0],eVecTmp[1]]

#WKB part ends here

#shooting part
def Q(x,g,E):
    return 1j*x**3-g*x**4-E

def f1(x,g,E):
    QVal=Q(x,g,E)
    return -1/4*(3*1j*x**2-4*g*x**3)/QVal+QVal**(1/2)

def oneStepRk4(g, E, h, xCurr, yCurr, vCurr):
    '''

    :param E:
    :param xCurr:
    :param yCurr:
    :param vCurr:
    :return: yNext, vNext
    '''
    M0 = h * Q(xCurr, g, E) * yCurr

    M1 = h * Q(xCurr + 1 / 2 * h, g, E) * (yCurr + 1 / 2 * h * vCurr)

    M2 = h * Q(xCurr + 1 / 2 * h, g, E) * (yCurr + 1 / 2 * h * vCurr + 1 / 4 * h * M0)

    M3 = h * Q(xCurr + h, g, E) * (yCurr + h * vCurr + 1 / 2 * h * M1)

    yNext = yCurr + h * vCurr + 1 / 6 * h * (M0 + M1 + M2)
    vNext = vCurr + 1 / 6 * (M0 + 2 * M1 + 2 * M2 + M3)

    return yNext, vNext

def calculateBoundaryValue(E, *data):
    # L = 10
    # Ns = 1000
    g = data[0]
    E = E[0]


    # rootsAll = np.roots(coefs)

    # print(rootsAll)
    # rootsAll = sorted(rootsAll, key=np.angle, reverse=True)
    # x1Val = rootsAll[2]
    L = 10  # 4*np.abs(np.real(x1Val))
    # print(L)
    theta = -1 / 10 * np.pi
    hAbsEst = 1e-2
    Ns = int(L / hAbsEst)
    h = -L / Ns * np.exp(1j * (theta))
    x0 = L * np.exp(1j * theta)
    xAll = [x0]
    yAll = [1]

    vAll = [f1(x0, g, E)]
    for j in range(0, Ns):
        xCurr = xAll[-1]
        yCurr = yAll[-1]
        vCurr = vAll[-1]
        yNext, vNext = oneStepRk4(g, E, h, xCurr, yCurr, vCurr)
        xNext = xCurr + h
        xAll.append(xNext)
        yAll.append(yNext)
        vAll.append(vNext)
    # boundary condition: Re(vN/yN)=0
    return np.real(vAll[-1] / yAll[-1])
