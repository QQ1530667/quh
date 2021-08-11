import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from datetime import datetime

# this script calculates the eigenvalue for potential 4, V(x)=x^{2}-igx^{5}
#Region R11

g = 0.01


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


def retX1X2(E):
    '''

    :param E: trial eigenvalue
    :return: x1 and x2 of
    '''
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll = np.roots(coefs)
    # print(rootsAll)
    rootsAll=sorted(rootsAll,key=np.angle,reverse=True)
    x1Val = rootsAll[2]
    x2Val = rootsAll[0]
    for rtTmp in rootsAll:
        if chooseX1(rtTmp):
            x1Val = rtTmp
        if chooseX2(rtTmp):
            x2Val = rtTmp
    # print(x1Val)
    return x1Val, x2Val


def f(z, E):
    '''

    :param z: point on x2x1
    :param E: trial eigenvalue
    :return: f value
    '''
    return (1j * g * z ** 5 - z ** 2 + E) ** (1 / 2)


def simpInt(E, N, x1, x2):
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
    fValsOdd = [f(zAll[j], E) for j in range(1, N, 2)]
    fValsEven = [f(zAll[j], E) for j in range(2, N, 2)]

    return 1 / 3 * deltaA * (1 + 1j * slope) * (f(zAll[0], E)
                                                + 4 * sum(fValsOdd)
                                                + 2 * sum(fValsEven)
                                                + f(zAll[N], E)
                                                )


n = 6

def testRoots(E):
    coefs = [-1j * g, 0, 0, 1, 0, -E]
    rootsAll = np.roots(coefs)
    test1=False
    test2=False
    for rtTmp in rootsAll:
        test1=test1 or chooseX1(rtTmp)
    for rtTmp in rootsAll:
        test2=test2 or chooseX2(rtTmp)
    return (test1 and test2)

def eqn(EIn):
    '''

    :param EIn: trial eigenvalue, in the form of [re, im]
    :return:
    '''
    E = EIn[0] + 1j * EIn[1]
    x1, x2 = retX1X2(E)
    dx = 1e-4
    N = int(np.abs(x2 - x1) / dx)
    rst = simpInt(E, N, x1, x2) - (n+1/2) * np.pi
    return np.real(rst), np.imag(rst)

tStart=datetime.now()

eVec = sopt.fsolve(eqn, [10, 1], maxfev=100, xtol=1e-7)

tEnd=datetime.now()
print("time: ", tEnd-tStart)

eVal = eVec[0] + 1j * eVec[1]

print(eVal)
print("test :"+str(testRoots(eVal)))