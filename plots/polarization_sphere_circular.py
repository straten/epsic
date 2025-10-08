import sys
import numpy as np
import matplotlib.pyplot as plt

# 3D arrow copied from https://gist.github.com/WetHat/1d6cd0f7309535311a539b42cccca89c

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

elev = 17
azim = 25
roll = 0

chi_deg = 25.0
psi_deg = 60.0

chi = chi_deg*np.pi/180.0
psi = psi_deg*np.pi/180.0

xi = np.arctan2( np.cos(psi) * np.cos(chi), np.sin(chi) )
zeta = np.arcsin( np.sin(psi) * np.cos(chi) )

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
zoff=0.1
ax.set_xlim3d(-limit, limit)
ax.set_zlim3d(-limit+zoff, limit+zoff)
ax.set_ylim3d(-limit, limit)

plt.axis('off')

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
ax.text(len*stretch,0,0,"$S_2$", ha='right')
ax.text(0,len*stretch,0,"$S_3$")
ax.text(0,0,len*stretch,"$S_1$")

# Plot the equator
circle = np.linspace(0, 2*np.pi, 100)

x_lat = np.cos(circle)
y_lat = np.sin(circle)
z_lat = 0.0
ax.plot(x_lat, y_lat, z_lat, color='k', linestyle='dotted')

# Plot from zero to psi solid
circle = np.linspace(0, psi, 20)

x_lat = np.cos(circle) 
y_lat = np.sin(circle)
z_lat = 0.0
ax.plot(x_lat, y_lat, z_lat, color='k')

stretch=1.05
idx=8
ax.text(stretch*x_lat[idx], stretch*y_lat[idx], -0.05, "$2\\psi$", ha='left', va='top')


#
# draw the arc through xi in the S1-S2 plane at the origin
#
N = 10
radius=1.0
phi = np.linspace(0, xi, N)

# reserve an extra point for the apex = (arcx[0], arcy[0])
xi_x = np.zeros(N+1)
xi_y = np.zeros(N+1)
xi_z = np.zeros(N+1)

xi_z[1:] = radius * np.cos(phi)
xi_x[1:] = radius * np.sin(phi)

verts = [list(zip(xi_x, xi_y, xi_z))]
poly = Poly3DCollection(verts, alpha=0.6, color=psi_color)
ax.add_collection3d(poly)
ax.plot(xi_x[1:], xi_y[1:], xi_z[1:], color=major_color)

stretch=1.05
idx=7
ax.text(stretch*xi_x[idx], stretch*xi_y[idx], stretch*xi_z[idx], "$2\\xi$", ha='right')

#
# Plot the meridian that passes through the Stokes vector
#
half_circle = np.linspace(-np.pi/2, np.pi/2, 100)
x_long = np.cos(psi) * np.cos(half_circle)
y_long = np.sin(psi) * np.cos(half_circle)
z_long = np.sin(half_circle)
ax.plot(x_long, y_long, z_long, color='k', linestyle='dotted')

# Plot the meridian from chi to LCP solid
half_circle = np.linspace(chi, np.pi/2, 20)
x_long = np.cos(psi) * np.cos(half_circle)
y_long = np.sin(psi) * np.cos(half_circle) 
z_long = np.sin(half_circle)
ax.plot(x_long, y_long, z_long, color='k')

idx = 8
stretch=1.1
ax.text(stretch*x_long[idx], stretch*y_long[idx], stretch*z_long[idx], "$\\theta$", ha='left')

#
# draw the arc through zeta at the origin
#
phi = np.linspace(0, zeta, N)

# reserve an extra point for the apex = (arcx[0], arcy[0])
zeta_x = np.zeros(N+1)
zeta_y = np.zeros(N+1)
zeta_z = np.zeros(N+1)

zeta_z[1:] = radius * np.cos(xi) * np.cos(phi)
zeta_x[1:] = radius * np.sin(xi) * np.cos(phi)
zeta_y[1:] = radius * np.sin(phi)

verts = [list(zip(zeta_x, zeta_y, zeta_z))]
poly = Poly3DCollection(verts, alpha=0.6, color=chi_color)
ax.add_collection3d(poly)
ax.plot(zeta_x[1:], zeta_y[1:], zeta_z[1:], color=minor_color)

#
# mark the right angle
#
dright = 0.07
# Plot the line segment along the "latitude" parallel to meridian from the pole to xi offset by dright
arc = np.linspace(xi-dright, xi, 10)
x_long = np.cos(dright) * np.sin(arc)
y_long = np.sin(dright)
z_long = np.cos(arc)
ax.plot(x_long, y_long, z_long, color='k', zorder=10)
# Plot the line segment along the great circle from the S1-S2 plane toward zeta, offset by dright
arc = np.linspace(0, dright, 10)
x_long = np.sin(xi-dright) * np.cos(arc)
y_long = np.sin(arc)
z_long = np.cos(xi-dright) * np.cos(arc)
ax.plot(x_long, y_long, z_long, color='k', zorder=10)

stretch=1.1
idx=5
ax.text(stretch*zeta_x[idx], stretch*zeta_y[idx], stretch*zeta_z[idx], "$2\\zeta$", ha='left')

#
# draw the polarization vector
#
r=1.8
S1 = r * np.cos(xi) * np.cos(zeta)
S2 = r * np.sin(xi) * np.cos(zeta)
S3 = r * np.sin(zeta)

ax.arrow3D(0,0,0, S2,S3,S1, linewidth=2.0, arrowstyle=astyle, mutation_scale=30, fc='black', ec='black')
ax.text(S2,S3,S1,r'\textbf{\textit{S}}')

if np.size(sys.argv) > 1:
    plt.savefig(sys.argv[1], pad_inches=0)
else:
    plt.show()
