from plot4AllFuncs4R12 import *



g=0.1




nVals=range(0,11)
EVals=np.linspace(start=0,stop=100,num=100)

tCommence=datetime.now()

fig,ax=plt.subplots(figsize=(40,40))
for n in nVals:
    yTmps=[(n+1/2)*np.pi]*len(EVals)
    ax.plot(EVals,yTmps,color="black")

E2Vals=[]
int2Vals=[]
discarded=0
for ETmp in EVals:
    x1,x2=retX1X2(g,ETmp)
    intValTmp=integralQuadrature(g,ETmp,x1,x2)
    if np.abs(np.imag(intValTmp)/np.real(intValTmp))<1e-7:
        E2Vals.append(ETmp)
        int2Vals.append(np.real(intValTmp))
    else:
        discarded+=1

tEnd=datetime.now()
print("time: ",tEnd-tCommence)
print("discarded = "+str(discarded))
ax.plot(E2Vals,int2Vals,color="black")
ax.set_xlabel("E")
ax.set_ylabel("n")
plt.title("g="+str(g))
ytickVals=[(n+1/2)*np.pi for n in nVals]
yTickNames=[n for n in nVals]
plt.yticks(ytickVals,yTickNames)
plt.savefig("g"+str(g)+"energyIntersectionR12.png")