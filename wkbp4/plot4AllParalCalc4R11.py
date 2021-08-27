from plot4AllFuncs4R11 import *

# this script calculates and plots eigenvalues, potential 4, V(x)=x^{2}-igx^{5}
# Region I-I

gnIndAll = np.linspace(start=np.log10(0.01), stop=np.log10(10), num=100)
gAll = [10 ** elem for elem in gnIndAll]
EWKB = []
# calculate WKB eigenvalues

# calculate WKB eigenvalues
threadNum = 12
energyLevelMax = 10
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

# ax.set_xscale('log')
# ax.set_yscale('symlog')
plt.yticks(np.arange(0,18,0.5))
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=x^{2}-igx^{5}$, region I-I")

# data serialization
nSctVals = []
gSctVals = []
ERealSctVals = []
EImagSctVals = []
#data serialization
for itemTmp in retAll:
    nTmp, gTmp, ERe, EIm = itemTmp
    nSctVals.append(nTmp)
    gSctVals.append(gTmp)
    ERealSctVals.append(ERe)
    EImagSctVals.append(EIm)

partitionByN=dict()
energyLevelSet=set(nSctVals)
for nTmp in energyLevelSet:
    partitionByN[nTmp]=[]

#data partition
for itemTmp in retAll:
    nTmp,gTmp,ERe,EIm=itemTmp
    partitionByN[nTmp].append([gTmp,ERe,EIm])

#plot for the same energy level
for nTmp in partitionByN.keys():
    gVecTmp=[]
    EReVecTmp=[]
    for vecTmp in partitionByN[nTmp]:
        gVecTmp.append(vecTmp[0])
        EReVecTmp.append(vecTmp[1])
    ax.plot(gVecTmp,EReVecTmp,color="red")
sctRealPartWKB = ax.scatter(gSctVals, ERealSctVals, color="red", marker=".", s=50, label="WKB real part")


#x^2 dominant
gDom2Vals=[0]*len(levelsAll)
EDom2Vals=[domE2(n) for n in levelsAll]

# -igx^5 dominant
ngEdomE5=[]
for nTmp in levelsAll:
    for gTmp in gAll:
        ngEdomE5.append([nTmp,gTmp,domE5(nTmp,gTmp)])
#data serialization for -igx^5
gDom5Vals=[]
nDom5Vals=[]
EDom5Vals=[]
for itemTmp in ngEdomE5:
    nDom5Vals.append(itemTmp[0])
    gDom5Vals.append(itemTmp[1])

    EDom5Vals.append(itemTmp[2])


#scatter plot of dom2
dom2Scatter=ax.scatter(gDom2Vals,EDom2Vals,color="green",marker="o",s=50,label="$x^{2}$ dominant")
#scatter plot of dom5
dom5Scatter=ax.scatter(gDom5Vals,EDom5Vals,color="blue",marker="+",s=50,label="$-igx^{5}$ dominant")




plt.legend()
plt.savefig("tmp11.png")

tPltEnd = datetime.now()
print("plotting time: ", tPltEnd - tPltStart)
