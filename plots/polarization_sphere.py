import sys
import numpy as np
import matplotlib.pyplot as plt

# 3D arrow copied from https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

chi_deg = 50.0
psi_deg = 70.0

chi = chi_deg*np.pi/180.0
psi = psi_deg*np.pi/180.0

class Arrow3D(FancyArrowPatch):

    def __init__(self, x, y, z, dx, dy, dz, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._xyz = (x, y, z)
        self._dxdydz = (dx, dy, dz)

    def draw(self, renderer):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        super().draw(renderer)
        
    def do_3d_projection(self, renderer=None):
        x1, y1, z1 = self._xyz
        dx, dy, dz = self._dxdydz
        x2, y2, z2 = (x1 + dx, y1 + dy, z1 + dz)

        xs, ys, zs = proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs) 
    
def _arrow3D(ax, x, y, z, dx, dy, dz, *args, **kwargs):
    '''Add a 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
    return arrow

setattr(Axes3D, 'arrow3D', _arrow3D)

# https://www.datylon.com/blog/data-visualization-for-colorblind-readers
psi_color = "#e57a77"
chi_color = "#7ca1cc"

major_color = "#f05039"
minor_color = "#1f449c"

plt.rc('font', family='serif', serif='cm10', size=21)
plt.rc('text', usetex=True)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
ax.set_box_aspect([1,1,1])

limit = 0.85
ax.set_xlim3d(-limit, limit)
ax.set_zlim3d(-limit, limit)
ax.set_ylim3d(-limit, limit)

plt.axis('off')

elev = 25
azim = 35
roll = 0
ax.view_init(elev, azim, roll)

# draw the x, y, and z axes
len=1.5
astyle="-|>,head_width=0.2,head_length=0.6"
ax.arrow3D(0,0,0, len,0,0, arrowstyle=astyle, mutation_scale=30, fc='black', ec='black')
ax.arrow3D(0,0,0, 0,len,0, arrowstyle=astyle, mutation_scale=30, fc='black', ec='black')
ax.arrow3D(0,0,0, 0,0,len, arrowstyle=astyle, mutation_scale=30, fc='black', ec='black')

zx=[0,0]
zy=[0,0]
zz=[-1.1,-0.9]
ax.plot(zx,zy,zz,color='k')

stretch=1.03

# label the x, y, and z axes
ax.text(len*stretch,0,0,"$S_1$", ha='right')
ax.text(0,len*stretch,0,"$S_2$")
ax.text(0,0,len*stretch,"$S_3$")

# Plot the equator
circle = np.linspace(0, 2*np.pi, 100)

x_lat = np.cos(circle)
y_lat = np.sin(circle)
z_lat = 0.0
ax.plot(x_lat, y_lat, z_lat, color='k', linestyle='dotted')

# create the arc through psi at the origin
N = 10
radius=1.0
phi = np.linspace(0, psi, N)

# reserve an extra point for the apex = (arcx[0], arcy[0])
psix = np.zeros(N+1)
psiy = np.zeros(N+1)
psiz = np.zeros(N+1)

psix[1:] = radius * np.cos(phi)
psiy[1:] = radius * np.sin(phi)

verts = [list(zip(psix, psiy, psiz))]
poly = Poly3DCollection(verts, alpha=0.6, color=psi_color)
ax.add_collection3d(poly)
ax.plot(psix[1:], psiy[1:], psiz[1:], color=major_color)

stretch=1.37
idx=N//2
ax.text(stretch*psix[idx], stretch*psiy[idx], stretch*psiz[idx], "$2\psi$", ha='left')

# Plot the meridian
half_circle = np.linspace(-np.pi/2, np.pi/2, 100)
x_long = np.cos(psi) * np.cos(half_circle)
y_long = np.sin(psi) * np.cos(half_circle)
z_long = np.sin(half_circle)
ax.plot(x_long, y_long, z_long, color='k', linestyle='dotted')

# create the arc through chi at the origin
phi = np.linspace(0, chi, N)

# reserve an extra point for the apex = (arcx[0], arcy[0])
chix = np.zeros(N+1)
chiy = np.zeros(N+1)
chiz = np.zeros(N+1)

chix[1:] = radius * np.cos(psi) * np.cos(phi)
chiy[1:] = radius * np.sin(psi) * np.cos(phi)
chiz[1:] = radius * np.sin(phi)

verts = [list(zip(chix, chiy, chiz))]
poly = Poly3DCollection(verts, alpha=0.6, color=chi_color)
ax.add_collection3d(poly)
ax.plot(chix[1:], chiy[1:], chiz[1:], color=minor_color)

stretch=1.24
ax.text(stretch*chix[idx], stretch*chiy[idx], stretch*chiz[idx], "$2\chi$", ha='center')

# finally, draw the polarization vector
r=1.8
S1 = r * np.cos(psi) * np.cos(chi)
S2 = r * np.sin(psi) * np.cos(chi)
S3 = r * np.sin(chi)

ax.arrow3D(0,0,0, S1,S2,S3, linewidth=2.0, arrowstyle=astyle, mutation_scale=30, fc='black', ec='black')
ax.text(S1,S2,S3,r'\textbf{\textit{S}}')

if np.size(sys.argv) > 1:
    plt.savefig(sys.argv[1], pad_inches=0)
else:
    plt.show()
