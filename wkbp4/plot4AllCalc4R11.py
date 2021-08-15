from plot4AllFuncs4R11 import *

#this script calculates and plots eigenvalues, potential 4, V(x)=x^{2}-igx^{5}
#Region I-I


gnIndAll=np.linspace(start=-3,stop=-2,num=2)
gAll=[10**elem for elem in gnIndAll]
EWKB=[]
EShooting=[]
#calculate WKB eigenvalues

energyLevelMax=30
tWKBStart=datetime.now()
for gTmp in gAll:
    nAndEPairs=[]
    for nTmp in range(0,energyLevelMax+1):
        dataTmp=(nTmp,gTmp)
        eVecTmp=sopt.fsolve(eqn,[(nTmp+0.5)**2,0],args=dataTmp,maxfev=100,xtol=1e-7)
        eValTmp=eVecTmp[0]+1j*eVecTmp[1]
        if testRoots(gTmp,eValTmp):
            nAndEPairs.append([nTmp,eValTmp])
    EWKB.append(nAndEPairs)

tWKBEnd=datetime.now()

print("time for WKB: ",tWKBEnd-tWKBStart)

#calculate shooting eigenvalue
# use same g values as in WKB
tShtStart=datetime.now()
for j in range(0,len(gAll)):
    gTmp=gAll[j]
    dataTmp = (gTmp)
    EMaxTmp=int(np.real(EWKB[j][-1][1]))

    eVecTmp=[]

    for eEstTmp in range(1,EMaxTmp+1):


        eValTmp=sopt.fsolve(calculateBoundaryValue,eEstTmp,args=dataTmp,maxfev=100,xtol=1e-7)[0]

        eVecTmp.append(eValTmp)
    EShooting.append(eVecTmp)


tShtEnd=datetime.now()
print("shooting time: ", tShtEnd-tShtStart)



tPltStart=datetime.now()
fig,ax=plt.subplots(figsize=(20,20))

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel("E")
ax.set_xlabel("g")
ax.set_title("Eigenvalues for potential $V(x)=x^{2}-igx^{5}$")
#plot WKB
for j in range(0,len(gAll)):
    gTmp=gAll[j]
    nAndEPairsTmp=EWKB[j]
    if len(nAndEPairsTmp)>0:
        for pairTmp in nAndEPairsTmp:
            sctWKB=ax.scatter(gTmp,np.real(pairTmp[1]),color="red",marker=".",s=50)

#plot shooting
for j in range(0,len(gAll)):
    gTmp=gAll[j]
    eShtVecTmp=EShooting[j]
    if len(eShtVecTmp)>0:
        for eTmp in eShtVecTmp:
            sctSht=ax.scatter(gTmp,eTmp,color="green",marker="+")



plt.legend((sctWKB,sctSht),("WKB","Shooting"),loc="best",fontsize=15)

plt.savefig("tmp.png")

tPltEnd=datetime.now()
print("plotting time: ", tPltEnd-tPltStart)