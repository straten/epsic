import math
import numpy as np
import matplotlib.pyplot as plt

# 3D arrow copied from https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

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
    '''Add an 3d arrow to an `Axes3D` instance.'''

    arrow = Arrow3D(x, y, z, dx, dy, dz, *args, **kwargs)
    ax.add_artist(arrow)
    return arrow

setattr(Axes3D, 'arrow3D', _arrow3D)

# https://www.datylon.com/blog/data-visualization-for-colorblind-readers
psi_color = "#e57a77"
chi_color = "#7ca1cc"

major_color = "#f05039"
minor_color = "#1f449c"

plt.rcParams['font.family'] = 'serif'
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.rc('font', size=32)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(111, projection='3d')
ax.set_box_aspect([1,1,1])

# Set up the figure axes, etc.
#ax.set_xlim(-1.0, 1.0)
#ax.set_ylim(-1.0, 1.0)
#ax.set_zlim(-1.0, 1.0)

plt.axis('off')

elev = 30
azim = 30
roll = 0
ax.view_init(elev, azim, roll)

# draw the x, y, and z axes
len=1.5
ax.arrow3D(-len,0,0, 2*len,0,0, arrowstyle="-|>", mutation_scale=30, fc='black', ec='black')
ax.arrow3D(0,-len,0, 0,2*len,0, arrowstyle="-|>", mutation_scale=30, fc='black', ec='black')
ax.arrow3D(0,0,-len, 0,0,2*len, arrowstyle="-|>", mutation_scale=30, fc='black', ec='black')

# label the x, y, and z axes
ax.text(1.01,-0.3,-0.1,"$S_1$")
ax.text(0.03,1.01,0,"$S_2$")
ax.text(0.03,0,1.01,"$S_3$")

chi_deg = 60.0
psi_deg = 60.0
r = 1.0

chi = chi_deg*math.pi/180.0
psi = psi_deg*math.pi/180.0

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
arcx = np.zeros(N+1)
arcy = np.zeros(N+1)
arcz = np.zeros(N+1)

arcx[1:] = radius * np.cos(phi)
arcy[1:] = radius * np.sin(phi)

verts = [list(zip(arcx, arcy, arcz))]
poly = Poly3DCollection(verts, alpha=0.6, color=psi_color)
ax.add_collection3d(poly)
#ax.plot(arcx[1:], arcy[1:], arcz[1:], color='k')

# Plot the meridian
half_circle = np.linspace(-np.pi/2, np.pi/2, 100)
x_long = np.cos(psi) * np.cos(half_circle)
y_long = np.sin(psi) * np.cos(half_circle)
z_long = np.sin(half_circle)
ax.plot(x_long, y_long, z_long, color='k', linestyle='dotted')

# create the arc through chi at the origin
phi = np.linspace(0, chi, N)
arcx[1:] = radius * np.cos(psi) * np.cos(phi)
arcy[1:] = radius * np.sin(psi) * np.cos(phi)
arcz[1:] = radius * np.sin(phi)

verts = [list(zip(arcx, arcy, arcz))]
poly = Poly3DCollection(verts, alpha=0.6, color=chi_color)
ax.add_collection3d(poly)
ax.plot(arcx[1:], arcy[1:], arcz[1:], color='k')

plt.show()
