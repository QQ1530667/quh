from plotAllFuncs4R12 import *

# this script calculates and plots eigenvalues, potential 4, V(x)=x^{2}-igx^{5}
# Region I-II

gnIndAll = np.linspace(start=-7, stop=-0, num=1)
gAll = [10 ** elem for elem in gnIndAll]
EWKB = []
# calculate WKB eigenvalues

# calculate WKB eigenvalues

energyLevelMax = 30
tWKBStart = datetime.now()
for gTmp in gAll:
    nAndEPairs = []
    for nTmp in range(0, energyLevelMax + 1):
        dataTmp = (nTmp, gTmp)
        eVecTmp = sopt.fsolve(eqn, [np.abs(nTmp + 0.5), 0], args=dataTmp, maxfev=100, xtol=1e-7)
        eValTmp = eVecTmp[0] + 1j * eVecTmp[1]

        nAndEPairs.append([nTmp, eValTmp])
    EWKB.append(nAndEPairs)

tWKBEnd = datetime.now()

print("time for WKB: ", tWKBEnd - tWKBStart)

tPltStart = datetime.now()
fig, ax = plt.subplots(figsize=(20, 20))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=x^{2}-igx^{5}$")
# plot WKB
gSctVals = []
ERealSctVals = []
EImagSctVals = []
for j in range(0, len(gAll)):
    gTmp = gAll[j]
    nAndEPairsTmp = EWKB[j]
    if len(nAndEPairsTmp) > 0:
        for pairTmp in nAndEPairsTmp:
            gSctVals.append(gTmp)
            ERealSctVals.append(np.real(pairTmp[1]))
            # sctWKB = ax.scatter(gTmp, np.real(pairTmp[1]), color="red", marker=".", s=50)

sctRealPartWKB = ax.scatter(gSctVals, ERealSctVals, color="red", marker=".", s=50, label="WKB real part")

# plt.legend((sctRealPartWKB),("WKB real part"),loc="best",fontsize=15,scatterpoints=1)
plt.legend()
plt.savefig("tmp2.png")

tPltEnd = datetime.now()
print("plotting time: ", tPltEnd - tPltStart)
