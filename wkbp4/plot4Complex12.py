from plot4AllFuncs4R12 import *



g=0.1

n=0
#plot horizontal line


fig,ax=plt.subplots(figsize=(20,20))
ax.scatter((n+1/2)*np.pi,0,color="red",marker="+",s=40)
r=12
#plot circle
thetaAll=np.linspace(start=-np.pi,stop=np.pi,num=100)
tStart=datetime.now()
int12345ValsAll=[]#contains 5 lists of integral values
for j in range(0,5):
    int12345ValsAll.append([])
for thTmp in thetaAll:
    ETmp=r*np.exp(1j*thTmp)
    # x1,x2=retX1X2New(g,ETmp)
    # int1ValsAll.append(integralQuadrature(g,ETmp,x1,x2))
    # int2ValsAll.append(integralQuadratureAnotherBranch(g,ETmp,x1,x2))
    rootPairs=returnFivePairsOfRoots(g,ETmp)
    for j in range(0,len(rootPairs)):
        x2Tmp,x1Tmp=rootPairs[j]
        int12345ValsAll[j].append(integralQuadrature(g,ETmp,x1Tmp,x2Tmp))


tEnd=datetime.now()
print("time : ",tEnd-tStart)
# int1Real=[]
# int1Imag=[]
# int2Real=[]
# int2Imag=[]
# for int1Tmp in int1ValsAll:
#     int1Real.append(np.real(int1Tmp))
#     int1Imag.append(np.imag(int1Tmp))
# for int2Tmp in int2ValsAll:
#     int2Real.append(np.real(int2Tmp))
#     int2Imag.append(np.imag(int2Tmp))
#
# ax.scatter(int1Real,int1Imag,color="black")
# ax.scatter(int2Real,int2Imag,color="red")
# plt.title("r="+str(r))
colorList=["b","g","r","m","k"]
for j in range(0,5):
    intTmpReal=[]
    intTmpImag=[]
    for intTmp in int12345ValsAll[j]:
        intTmpReal.append(np.real(intTmp))
        intTmpImag.append(np.imag(intTmp))

    ax.scatter(intTmpReal,intTmpImag,color=colorList[j])

plt.title("r="+str(r))


plt.savefig("./energyCircleR12/g"+str(g)+"r"+str(r)+"energyIntersectionR12.png")
plt.close()

