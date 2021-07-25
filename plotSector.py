import numpy as np
import cmath as cm
import matplotlib.pyplot as plt


N=100
R=7400
c=1e18

rs=np.linspace(start=0,stop=R,num=N)

a1=7/14*np.pi
a2=3/14*np.pi
a3=-1/14*np.pi
a4=-5/14*np.pi
a5=-9/14*np.pi
a6=-13/14*np.pi
a7=-17/14*np.pi

p1x=[0]*len(rs)
p1y=rs

p2x=[elem*np.sign(np.cos(a2)) for elem in rs]
p2y=[xtmp*np.tan(a2) for xtmp in p2x]
p3x=[elem*np.sign(np.cos(a3)) for elem in rs]
p3y=[xtmp*np.tan(a3) for xtmp in p3x]
p4x=[elem*np.sign(np.cos(a4)) for elem in rs]
p4y=[xtmp*np.tan(a4) for xtmp in p4x]
p5x=[elem*np.sign(np.cos(a5)) for elem in rs]
p5y=[xtmp*np.tan(a5) for xtmp in p5x]

p6x=[elem*np.sign(np.cos(a6)) for elem in rs]
p6y=[xtmp*np.tan(a6) for xtmp in p6x]
p7x=[elem*np.sign(np.cos(a7)) for elem in rs]
p7y=[xtmp*np.tan(a7) for xtmp in p7x]
b1=1/4*np.pi
b2=-1/4*np.pi
b3=-3/4*np.pi
b4=-5/4*np.pi


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
ax.set_yticks([])
ax.set_xticks([])
ax.plot(p1x,p1y,p2x,p2y,p3x,p3y,p4x,p4y,p5x,p5y,p6x,p6y,p7x,p7y,color="blue")
ax.plot(q1x,q1y,q2x,q2y,q3x,q3y,q4x,q4y,color="red")
ax.fill_between(rs, p2y, p3y,color='gainsboro')
ax.fill_between(p3x,p3y,q2y,color="aquamarine")
ax.fill_between(p7x,p7y,p6y,color="gainsboro")
ax.fill_between(p6x,p6y,q3y,color="aquamarine")
p=int(N/2)
#region right I-I

tx=p2x[p]
ty=1/3*(2*p2y[p]+p3y[p])
ax.text(tx,ty,"I-I",fontsize=16)
#region left I-I
tx=p6x[p]
ty=1/3*(2*p7y[p]+p6y[p])
ax.text(tx,ty,"I-I",fontsize=16)
#region right I-II
tx=p3x[p]
ty=1/3*(2*q3y[p]+p3y[p])
ax.text(tx,ty,"I-II",fontsize=16)
#region left I-II
tx=p6x[p]
ty=1/3*(2*q3y[p]+p6y[p])
ax.text(tx,ty,"I-II",fontsize=16)

#values of g and E

g=0.01
E=g**(1/3)*c
ax.set_title("$x^{2}-igx^{5}-E=0$, $g=$"+str(g)+", $E=$"+str(E))

coefs=[-1j*g,0,0,1,0,-E]
rootsAll=np.roots(coefs)
print(np.abs(rootsAll))
sVal=30
for tmp in rootsAll:
    ax.scatter(np.real(tmp),np.imag(tmp),color="black",s=sVal)





plt.savefig("./roots1/"+str(coefs)+".png")