import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from datetime import datetime

# this script combines calculations from WKB
# and shooting as a verification for potential 4, V(x)=x^{2}-igx^{5}
# in region I-II
# WKB part

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
    x1Val = rootsAll[3]
    x2Val = rootsAll[4]

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
