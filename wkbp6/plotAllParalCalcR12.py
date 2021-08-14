from plotAllFuncsR12 import *

#this script calculates and plots eigenvalues, potential 6, V(x)=ix^{3}-gx^{4}
# in region I-II
# WKB part calculation

gnIndAll = np.linspace(start=-10, stop=-0, num=1)
gAll = [10 ** elem for elem in gnIndAll]

threadNum = 12
energyLevelMax = 30
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

# shooting part calculation
EShooting=[]
tShtStart=datetime.now()
for j in range(0,len(gAll)):
    gTmp=gAll[j]
    dataTmp=(gTmp)
    EMaxTmp=np.max(np.abs(EWKBRealSctVals))/1.5

    eVecTmp = []
    for eEstTmp in range(1, int(EMaxTmp) + 1):
        eValTmp = sopt.fsolve(calculateBoundaryValue, eEstTmp, args=dataTmp, maxfev=100, xtol=1e-7)[0]

        eVecTmp.append(eValTmp)
    EShooting.append(eVecTmp)


tShtEnd=datetime.now()
print("shooting time: ", tShtEnd-tShtStart)
#data serialization
gShtSctVals=[]
EShtSctVals=[]
for j in range(0,len(gAll)):
    gTmp=gAll[j]
    for eTmp in EShooting[j]:
        gShtSctVals.append(gTmp)
        EShtSctVals.append(eTmp)


tPltStart=datetime.now()
fig,ax=plt.subplots(figsize=(20,20))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=ix^{3}-gx^{4}$")
#plot WKB
sctWKB=ax.scatter(gWKBSctVals,EWKBRealSctVals,color="red",marker=".",s=50,label="WKB real part")

#plot shooting
scrSht=ax.scatter(gShtSctVals,EShtSctVals,color="green",marker="+",label="Shooting")

plt.legend()
plt.savefig("tmp12.png")
tPltEnd=datetime.now()
print("plotting time: ", tPltEnd-tPltStart)