from plotAllFuncs4R12 import *

# this script calculates and plots eigenvalues, potential 4, V(x)=x^{2}-igx^{5}
# Region I-II

gnIndAll = np.linspace(start=-7, stop=-0, num=1)
gAll = [10 ** elem for elem in gnIndAll]
EWKB = []
# calculate WKB eigenvalues

# calculate WKB eigenvalues
threadNum = 12
energyLevelMax = 3
levelsAll = range(0, energyLevelMax + 1)
inDataAll = []
for nTmp in levelsAll:
    for gTmp in gAll:
        inDataAll.append((nTmp, gTmp))

tWKBParalStart = datetime.now()
pool1 = Pool(threadNum)
retAll = pool1.map(computeOneSolution, inDataAll)

tWKBParalEnd = datetime.now()
print("WKB time: ", tWKBParalEnd - tWKBParalStart)
tPltStart = datetime.now()
# plot WKB
fig, ax = plt.subplots(figsize=(20, 20))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=x^{2}-igx^{5}$")

# data serialization
nSctVals = []
gSctVals = []
ERealSctVals = []
EImagSctVals = []
for itemTmp in retAll:
    nTmp, gTmp, ERe, EIm = itemTmp
    nSctVals.append(nTmp)
    gSctVals.append(gTmp)
    ERealSctVals.append(ERe)
    EImagSctVals.append(EIm)

sctRealPartWKB = ax.scatter(gSctVals, ERealSctVals, color="red", marker=".", s=50, label="WKB real part")
plt.legend()
plt.savefig("tmp2.png")

tPltEnd = datetime.now()
print("plotting time: ", tPltEnd - tPltStart)
print(retAll)