import numpy as np
import scipy.optimize as sopt
import matplotlib.pyplot as plt
from  multiprocessing import Pool
from datetime import datetime
# import mpmath
# from mpmath import mp
# mp.dps=100

def Q(x, g, E):
    return x ** 2 - 1j * g * x ** 5 - E

def f1(x,g,E):
   return -1 / 4 * (2 * x - 5 * 1j * g * x ** 4) / Q(x, g, E) + Q(x, g, E) ** (1 / 2)

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

def calculateBoundaryValue(E,*data):
    # L = 10
    # Ns = 1000
    E=E[0]
    g=data[0]




    # rootsAll = np.roots(coefs)

    # print(rootsAll)
    # rootsAll = sorted(rootsAll, key=np.angle, reverse=True)
    # x1Val = rootsAll[2]
    L = 15 # 4*np.abs(np.real(x1Val))
    # print(L)
    theta = - 4/ 21 * np.pi
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


def computeOneShootingSolution(inData):
    '''

    :param inData: [g,Eest]
    :return:

    '''
    g, Eest=inData


    eVal=sopt.fsolve(calculateBoundaryValue,Eest,args=(g),maxfev=100,xtol=1e-6)[0]
    return [g,eVal]