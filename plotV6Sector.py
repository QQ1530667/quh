import numpy as np

import matplotlib.pyplot as plt
from decimal import Decimal
#potential 6, V(x)=ix^{3}-gx^{4}
N=100


c=1e-1
#compute roots
g=0.01
E=g**(-1)*c
coefs=[-g,1j,0,0,-E]
rootsAll=np.roots(coefs)
rs=np.linspace(start=0,stop=np.max(np.abs(rootsAll))*1.1,num=N)

a0=1/3*np.pi

a1=0

a2=-1/3*np.pi

a3=-2/3*np.pi

a4=-np.pi

a5=-4/3*np.pi



p0x=[elem*np.sign(np.cos(a0)) for elem in rs]
p0y=[xtmp*np.tan(a0) for xtmp in p0x]

p1x=[elem*np.sign(np.cos(a1)) for elem in rs]
p1y=[xtmp*np.tan(a1) for xtmp in p1x]

p2x=[elem*np.sign(np.cos(a2)) for elem in rs]
p2y=[xtmp*np.tan(a2) for xtmp in p2x]

p3x=[elem*np.sign(np.cos(a3)) for elem in rs]
p3y=[xtmp*np.tan(a3) for xtmp in p3x]

p4x=[elem*np.sign(np.cos(a4)) for elem in rs]
p4y=[xtmp*np.tan(a4) for xtmp in p4x]

p5x=[elem*np.sign(np.cos(a5)) for elem in rs]
p5y=[xtmp*np.tan(a5) for xtmp in p5x]



b1=1/10*np.pi
b2=-3/10*np.pi
b3=-7/10*np.pi
b4=-11/10*np.pi

q0x=[0]*len(rs)
q0y=[elem for elem in rs]
q1x=[elem*np.sign(np.cos(b1)) for elem in rs]
q1y=[xtmp*np.tan(b1) for xtmp in q1x]

q2x=[elem*np.sign(np.cos(b2)) for elem in rs]
q2y=[xtmp*np.tan(b2) for xtmp in q2x]

q3x=[elem*np.sign(np.cos(b3)) for elem in rs]
q3y=[xtmp*np.tan(b3) for xtmp in q3x]

q4x=[elem*np.sign(np.cos(b4)) for elem in rs]
q4y=[xtmp*np.tan(b4) for xtmp in q4x]
fig,ax=plt.subplots(figsize=(20,20))

ax.spines['bottom'].set_color('grey')
ax.spines['bottom'].set_position('zero')
ax.spines['left'].set_color('grey')
ax.spines['left'].set_position('center')
# ax.set_yticks([])
# ax.set_xticks([])
ax.plot(p0x,p0y,p1x,p1y,p2x,p2y,p3x,p3y,p4x,p4y,p5x,p5y,color="blue")
ax.plot(q0x,q0y,q1x,q1y,q2x,q2y,q3x,q3y,q4x,q4y,color="red")
ax.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#color filling
ax.fill_between(q1x, q1y, p1y,color='gainsboro')
ax.fill_between(q4x,q4y,p4y,color="gainsboro")
ax.fill_between(p1x,p1y,q2y,color="aquamarine")
ax.fill_between(p4x,p4y,q3y,color="aquamarine")
#text
p=int(N/2)
#region right I-I

tx=q1x[p]
ty=1/3*(2*q1y[p]+p1y[p])
ax.text(tx,ty,"I-I",fontsize=16)
#region left I-I
tx=q4x[p]
ty=1/3*(2*q4y[p]+p4y[p])
ax.text(tx,ty,"I-I",fontsize=16)
#region right I-II
tx=p1x[p]
ty=1/3*(2*p1y[p]+q2y[p])
ax.text(tx,ty,"I-II",fontsize=16)
#region left I-II
tx=p4x[p]
ty=1/3*(2*p4y[p]+q3y[p])
ax.text(tx,ty,"I-II",fontsize=16)
#compute angles
rootsReversedAngle=sorted(rootsAll,key=np.angle,reverse=True)
anglesAllInPi=[np.angle(rt)/np.pi for rt in rootsReversedAngle]
print("angles in pi: ",anglesAllInPi)


sVal=30
#plot roots
for tmp in rootsAll:
    ax.scatter(np.real(tmp),np.imag(tmp),color="black",s=sVal)



ERe=np.real(E)
EIm=np.imag(E)
ax.set_title("$ix^{3}-gx^{4}-E=0$, $g=$"+str(g)+", $E=$"+'{:.2e}'.format(Decimal(str(ERe)))+"+i"+'{:.2e}'.format(Decimal(str(EIm))))



plt.savefig("./potential6/"+str(coefs)+".png")