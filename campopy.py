import numpy as np
import matplotlib.pyplot as plt
import CampoLib as CL

print(CL.campo1.__doc__)

x,y,Bx1,By1 = CL.campo1(11,'conf0_10.xyz',0.0)

x,y,Bx2,By2 = CL.campo1(11,'parGG.xyz',0.0)
#Bx,By = CL.difcampo(100,Bx1,By1,Bx2,By2)

Bx = Bx2 - Bx1
By = By2 - By1

X,Y = np.meshgrid(x,y)
fig, ax = plt.subplots(figsize=(6,6))

ax.streamplot(X,Y,Bx,By)

ax.set_aspect('equal')

ax.set_title('Stream Plot of Reference State')

plt.show()

#x,y,Bx2,By2 = CL.campo(500,'linear1.xyz',0.2)


#Bx = Bx1-Bx2
#By = By1-By2

#ax.streamplot(X,Y,Bx,By)

#ax.set_aspect('equal')

#ax.set_title('Stream Plot of Reference State')

#plt.show()
