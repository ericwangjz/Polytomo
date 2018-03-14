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


#qubit example for visualization


Nqubit=1000
qubitdata=np.array([390., 211, 202, 197])/Nqubit
dqubit=2
eps=0.01

crqubit=polytopeCR(qubitdata,paulibasis,sic2,Nqubit,dqubit,eps,True,'CP')


fig = plt.figure()
ax=Axes3D(fig)

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
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 1. * np.outer(np.cos(u), np.sin(v))
y = 1. * np.outer(np.sin(u), np.sin(v))
z = 1. * np.outer(np.ones(np.size(u)), np.cos(v))
X = np.arange(-1, 1, 0.025)
Y = np.arange(-1, 1, 0.025)
X, Y = np.meshgrid(X, Y)



# Plot the surface.
ax.plot_surface(x, y, z, color='g', alpha=0.05)

plt.xlabel('x')
plt.ylabel('y')
plt.show()
