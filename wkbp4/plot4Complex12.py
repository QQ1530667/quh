from plot4AllFuncs4R12 import *



g=0.1

n=0
#plot horizontal line


fig,ax=plt.subplots(figsize=(20,20))
ax.scatter((n+1/2)*np.pi,0,color="red",marker="+",s=40)
r=3.25
#plot circle
thetaAll=np.linspace(start=0,stop=2*np.pi,num=100)
tStart=datetime.now()
intValsAll=[]
for thTmp in thetaAll:
    ETmp=r*np.exp(1j*thTmp)
    x1,x2=retX1X2New(g,ETmp)
    intValsAll.append(integralQuadrature(g,ETmp,x1,x2))

tEnd=datetime.now()
print("time : ",tEnd-tStart)
intReal=[]
intImag=[]

for intTmp in intValsAll:
    intReal.append(np.real(intTmp))
    intImag.append(np.imag(intTmp))
ax.scatter(intReal,intImag,color="black")





plt.savefig("g"+str(g)+"r"+str(r)+"energyIntersectionR12.png")
plt.close()

intDistAll=[np.abs(elem-(n+1/2)*np.pi) for elem in intValsAll]

intOrd=sorted(range(len(intDistAll)),key=lambda k: intDistAll[k])
th0Ind=intOrd[0]
th1Ind=intOrd[1]

e0=r*np.exp(1j*thetaAll[th0Ind])
e1=r*np.exp(1j*thetaAll[th1Ind])
print(e0)
print(e1)
