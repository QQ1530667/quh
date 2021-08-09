import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from datetime import datetime

# this script calculates eigenvalues using shooting in region R11,
# potential 4,  V(x)=x^{2}-igx^{5}

g = 0.01
L = 10
N = 1000
h = -np.sign(L) * L / N


def Q(x, E):
    return x ** 2 - 1j * g * x ** 5 - E


def f(x, E):
    return -1 / 4 * (2 * x - 5 * 1j * g * x ** 4) / Q(x, E) + Q(x, E) ** (1 / 2)


def oneStepRk4(E, xCurr, yCurr, vCurr):
    '''

    :param E:
    :param xCurr:
    :param yCurr:
    :param vCurr:
    :return: yNext, vNext
    '''
    M0 = h * Q(xCurr, E) * yCurr

    M1 = h * Q(xCurr + 1 / 2 * h, E) * (yCurr + 1 / 2 * h * vCurr)

    M2 = h * Q(xCurr + 1 / 2 * h, E) * (yCurr + 1 / 2 * h * vCurr + 1 / 4 * h * M0)

    M3 = h * Q(xCurr + h, E) * (yCurr + h * vCurr + 1 / 2 * h * M1)

    yNext = yCurr + h * vCurr + 1 / 6 * h * (M0 + M1 + M2)
    vNext = vCurr + 1 / 6 * (M0 + 2 * M1 + 2 * M2 + M3)

    return yNext, vNext


def calculateBoundaryValue(E):
    xAll = [L]
    yAll = [1]
    vAll = [f(L, E)]
    for j in range(0, N):
        xCurr = xAll[-1]
        yCurr = yAll[-1]
        vCurr = vAll[-1]
        yNext, vNext = oneStepRk4(E, xCurr, yCurr, vCurr)
        xNext = xCurr + h
        xAll.append(xNext)
        yAll.append(yNext)
        vAll.append(vNext)
    # boundary condition: Re(vN/yN)=0
    return np.real(vAll[-1] / yAll[-1])


tStart=datetime.now()
eVal=sopt.fsolve(calculateBoundaryValue,13,maxfev=100,xtol=1e-7)
tEnd=datetime.now()
print("time: ", tEnd-tStart)
print(eVal)