from shooting4ParalFuncsR12 import *
# this script calculates and plots shooting eigenvalues, potential 4, V(x)=x^{2}-igx^{5}
# Region I-II
num=100
startG=1e-3
stopG=1e-1

gnIndAll = np.linspace(start=np.log10(startG), stop=np.log10(stopG), num=num)
gAll = [10 ** elem for elem in gnIndAll]
EMax=30
inDataAll=[]
for g in gAll:
    for E in  range(1,EMax):
        inDataAll.append([g,E])

threadNum=12

pool1=Pool(threadNum)

tShootingStart=datetime.now()
retAll=pool1.map(computeOneShootingSolution,inDataAll)

tShootingEnd=datetime.now()
print("shooting time : ",tShootingEnd-tShootingStart)

# plot shooting
fig, ax = plt.subplots(figsize=(20, 20))

# ax.set_xscale('log')
ax.set_yscale('symlog')
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Shooting eigenvalues for potential $V(x)=x^{2}-igx^{5}$")

# data serialization
gShootingVals = []
EShootingVals=[]
for itemTmp in retAll:
    gTmp,ETmp=itemTmp
    gShootingVals.append(gTmp)
    EShootingVals.append(ETmp)


shootingScatter=ax.scatter(gShootingVals,EShootingVals,color="blue",marker=".",s=50,label="shooting")
plt.legend()

plt.savefig("shootingR12.png")