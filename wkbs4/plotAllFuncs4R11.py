import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from datetime import datetime


# this script combines calculations from WKB
# and shooting as a verification for potential 4, V(x)=x^{2}-igx^{5}

# WKB part
def chooseX1(xTmp):
    '''

    :param xTmp: input complex x
    :return: true if arg(xTmp) is in (-1/14 pi, 3/14 pi); false otherwise
    '''
    angleXTmp = np.angle(xTmp)
    # print(angleXTmp)
    rst = (angleXTmp >= -1 / 14 * np.pi) and (angleXTmp <= 3 / 14 * np.pi)
    return rst


def chooseX2(xTmp):
    '''

    :param xTmp: input complex x
    :return: true if arg(xTmp) is in [-pi, -13/14 pi) U (11/14 pi, pi]
    '''
    angleTmp = np.angle(xTmp)
    ret = (angleTmp >= -np.pi and angleTmp <= -13 / 14 * np.pi) or (angleTmp >= 11 / 14 * np.pi and angleTmp < np.pi)
    return ret


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
    x1Val = rootsAll[2]
    x2Val = rootsAll[0]
    for rtTmp in rootsAll:
        if chooseX1(rtTmp):
            x1Val = rtTmp
        if chooseX2(rtTmp):
            x2Val = rtTmp
    # print(x1Val)
    return x1Val, x2Val


def f(z, g, E):
    '''
    :param g: const
    :param z: point on x2x1
    :param E: trial eigenvalue
    :return: f value
    '''
    return (1j * g * z ** 5 - z ** 2 + E) ** (1 / 2)


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

def testRoots(g,E):
    '''

    :param g: const
    :param E: eigenvalue
    :return: test if there are roots in left and right region R11
    '''
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll = np.roots(coefs)
    test1=False
    test2=False
    for rtTmp in rootsAll:
        test1=test1 or chooseX1(rtTmp)
    for rtTmp in rootsAll:
        test2=test2 or chooseX2(rtTmp)
    return (test1 and test2)



def eqn(EIn,*data):
    '''

    :param EIn: trial eigenvalue, in the form of [re, im]
    :return:
    '''

    n,g=data

    E = EIn[0] + 1j * EIn[1]
    x1, x2 = retX1X2(g,E)
    dx = 1e-4
    N = int(np.abs(x2 - x1) / dx)
    rst = simpInt(g,E, N, x1, x2) - (n+1/2) * np.pi
    return np.real(rst), np.imag(rst)

###functions for WKB end here

#shooting part
L = 10
Ns = 1000
h = -np.sign(L) * L / Ns

def Q(x, g,E):
    return x ** 2 - 1j * g * x ** 5 - E


def f1(x,g, E):
    return -1 / 4 * (2 * x - 5 * 1j * g * x ** 4) / Q(x,g, E) + Q(x,g, E) ** (1 / 2)

def oneStepRk4(g,E, xCurr, yCurr, vCurr):
    '''

    :param E:
    :param xCurr:
    :param yCurr:
    :param vCurr:
    :return: yNext, vNext
    '''
    M0 = h * Q(xCurr, g,E) * yCurr

    M1 = h * Q(xCurr + 1 / 2 * h, g,E) * (yCurr + 1 / 2 * h * vCurr)

    M2 = h * Q(xCurr + 1 / 2 * h,g, E) * (yCurr + 1 / 2 * h * vCurr + 1 / 4 * h * M0)

    M3 = h * Q(xCurr + h, g,E) * (yCurr + h * vCurr + 1 / 2 * h * M1)

    yNext = yCurr + h * vCurr + 1 / 6 * h * (M0 + M1 + M2)
    vNext = vCurr + 1 / 6 * (M0 + 2 * M1 + 2 * M2 + M3)

    return yNext, vNext

def calculateBoundaryValue(E,*data):
    xAll = [L]
    yAll = [1]
    g=data
    vAll = [f1(L,g, E)]
    for j in range(0, Ns):
        xCurr = xAll[-1]
        yCurr = yAll[-1]
        vCurr = vAll[-1]
        yNext, vNext = oneStepRk4(g,E, xCurr, yCurr, vCurr)
        xNext = xCurr + h
        xAll.append(xNext)
        yAll.append(yNext)
        vAll.append(vNext)
    # boundary condition: Re(vN/yN)=0
    return np.real(vAll[-1] / yAll[-1])


