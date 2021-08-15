from plotAllFuncsR11 import *

#this script calculates and plots eigenvalues, potential 6, V(x)=ix^{3}-gx^{4}
# in region I-I
# WKB part calculation
gnIndAll = np.linspace(start=-2, stop=np.log10(0.2), num=100)
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
    # if np.abs(ERe)>1e3:
    #     continue
    nWKBSctVals.append(nTmp)
    gWKBSctVals.append(gTmp)
    EWKBRealSctVals.append(ERe)
    EWKBImagSctVals.append(EIm)
#
# print(EWKBRealSctVals)
# print(EWKBImagSctVals)
tPltStart=datetime.now()
fig,ax=plt.subplots(figsize=(20,20))
# ax.set_xscale('log')
ax.set_yscale('symlog')
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=ix^{3}-gx^{4}$")
#plot WKB
sctWKB=ax.scatter(gWKBSctVals,EWKBRealSctVals,color="red",marker=".",s=20,label="WKB real part")


plt.legend()
plt.savefig("tmp11.png")
tPltEnd=datetime.now()
print("plotting time: ", tPltEnd-tPltStart)