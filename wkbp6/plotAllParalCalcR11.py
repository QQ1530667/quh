from plotAllFuncsR11 import *

#this script calculates and plots eigenvalues, potential 6, V(x)=ix^{3}-gx^{4}
# in region I-I
# WKB part calculation
gnIndAll = np.linspace(start=-10, stop=-0, num=1)
gAll = [10 ** elem for elem in gnIndAll]

threadNum = 12
energyLevelMax = 3
levelsAll = range(0, energyLevelMax + 1)
# levelsAll=[2]
inDataAll = []
for nTmp in levelsAll:
    for gTmp in gAll:
        inDataAll.append((nTmp, gTmp))

tWKBParalStart = datetime.now()
pool1 = Pool(threadNum)
retAll = pool1.map(computeOneSolution, inDataAll)

tWKBParalEnd = datetime.now()
print("WKB time: ", tWKBParalEnd - tWKBParalStart)

#Data serialization
nWKBSctVals = []
gWKBSctVals = []
EWKBRealSctVals = []
EWKBImagSctVals = []
for itemTmp in retAll:
    nTmp, gTmp, ERe, EIm = itemTmp
    nWKBSctVals.append(nTmp)
    gWKBSctVals.append(gTmp)
    EWKBRealSctVals.append(ERe)
    EWKBImagSctVals.append(EIm)


print(EWKBRealSctVals)
print(EWKBImagSctVals)