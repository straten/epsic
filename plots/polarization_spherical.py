import math
import numpy as np
import matplotlib.pyplot as plt

# 3D arrow copied from https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D

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
ax = fig.add_subplot(projection='3d')

# Set up the figure axes, etc.
ax.set_xlim(-1.0, 1.0)
ax.set_ylim(-1.0, 1.0)
ax.set_zlim(-1.0, 1.0)

plt.axis('off')

elev = 30
azim = 25
roll = 0
ax.view_init(elev, azim, roll)

# draw the x, y, and z axes
ax.arrow3D(-1,0,0, 2,0,0, arrowstyle="-|>", mutation_scale=30, fc='black', ec='black')
ax.arrow3D(0,-1,0, 0,2,0, arrowstyle="-|>", mutation_scale=30, fc='black', ec='black')
ax.arrow3D(0,0,-1, 0,0,2, arrowstyle="-|>", mutation_scale=30, fc='black', ec='black')

# label the x, y, and z axes
ax.text(1.01,-0.3,-0.1,"$S_1$")
ax.text(0.03,1.01,0,"$S_2$")
ax.text(0.03,0,1.01,"$S_3$")

chi_deg = 30.0
psi_deg = 40.0
r = 1.0

chi = chi_deg*math.pi/180.0
psi = psi_deg*math.pi/180.0

plt.show()
