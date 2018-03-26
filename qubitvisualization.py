from sampling import *
from polyconfiregion import *
from parameters import *
from scipy.spatial import ConvexHull
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d as a3
import polytope.polytope as pc
from randomtest import *

#qubit example for visualization


Nqubit=1000
dqubit=2
mmts=sic2
#qubitdata=np.array([390, 211, 202, 197])/float(Nqubit)
qubitdata=simulate_measurements(rand_dm(dqubit), [mmts], Nqubit).Nm/float(Nqubit)
print qubitdata
eps=0.01

crqubit=polytopeCR(qubitdata,paulibasis,sic2,Nqubit,dqubit,eps,True,'CP')


fig = plt.figure()
ax=Axes3D(fig)
ax.set_aspect('equal')
res0=np.array(pc.extreme(crqubit)).transpose()
res0T=pc.extreme(crqubit)
ax.scatter(res0[0],res0[1],res0[2])
hull=ConvexHull(res0T)
for i in np.arange(len(hull.simplices)):
    square=[ tuple(res0T[hull.simplices[i,0]]), tuple(res0T[hull.simplices[i,1]]), tuple(res0T[hull.simplices[i, 2]])]
    face = a3.art3d.Poly3DCollection([square])
    face.set_alpha(.05)
    face.set_color('m')
    face.set_edgecolor('k')
    ax.add_collection3d(face)

'''

    ax.set_zlim(.25, .75)
    ax.set_xlim(-.25, .25)
    ax.set_ylim(-.25, .25)
    '''

#plot bloch sphere
# Make data
elev = 10.0
rot = 80.0 / 180 * np.pi
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 1. * np.outer(np.cos(u), np.sin(v))
y = 1. * np.outer(np.sin(u), np.sin(v))
z = 1. * np.outer(np.ones(np.size(u)), np.cos(v))
#plot the grid lines
a = np.array([-np.sin(elev / 180 * np.pi), 0, np.cos(elev / 180 * np.pi)])
b = np.array([0, 1, 0])
b =  + np.cross(a, b)  + a * np.dot(a, b)
ax.plot(np.sin(u),np.cos(u),0,color='k', linestyle = 'dashed', alpha=0.05)
horiz_front = np.linspace(0, np.pi, 100)
ax.plot(np.sin(horiz_front),np.cos(horiz_front),0,color='k', linestyle = 'dashed', alpha=0.05)
vert_front = np.linspace(np.pi / 2, 3 * np.pi / 2, 100)
ax.plot(a[0] * np.sin(u) + b[0] * np.cos(u), b[1] * np.cos(u), a[2] * np.sin(u) + b[2] * np.cos(u),color='k', linestyle = 'dashed', alpha=0.05)
ax.plot(a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front), b[1] * np.cos(vert_front), a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front),color='k', linestyle = 'dashed', alpha=0.05)
ax.plot(b[1] * np.cos(u),a[0] * np.sin(u) + b[0] * np.cos(u),  a[2] * np.sin(u) + b[2] * np.cos(u),color='k', linestyle = 'dashed', alpha=0.05)
ax.plot(b[1] * np.cos(vert_front), a[0] * np.sin(vert_front) + b[0] * np.cos(vert_front), a[2] * np.sin(vert_front) + b[2] * np.cos(vert_front),color='k', linestyle = 'dashed', alpha=0.05)


# Plot the surface.
ax.plot_surface(x, y, z, color='g', alpha=0.05)

plt.xlabel('x')
plt.ylabel('y')
plt.show()
