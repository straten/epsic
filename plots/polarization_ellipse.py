import sys
import numpy as np
import matplotlib.pyplot as plt

# https://www.datylon.com/blog/data-visualization-for-colorblind-readers
psi_color = "#e57a77"
chi_color = "#7ca1cc"

major_color = "#f05039"
minor_color = "#1f449c"

plt.rc('font', family='serif', serif='cm10', size=32)
plt.rc('text', usetex=True)

fig = plt.figure(figsize=(8,8))

# Set up the figure axes, etc.
plt.xlim(-1.0, 1.0)
plt.ylim(-1.0, 1.0)
    
plt.axis('off')

# draw the x and y axes
hw, hl = 0.06, 0.10
plt.arrow(-1,0,2,0, head_width=hw, head_length=hl, fc='black', ec='black', length_includes_head=True)
plt.arrow(0,-1,0,2, head_width=hw, head_length=hl, fc='black', ec='black', length_includes_head=True)

# label the x and y axes
plt.text(1.0,0.05,"$x$")
plt.text(0.03,1.0,"$y$")

# draw the ellipse
chi_deg = 30.0
psi_deg = 30.0
r = 1.0

chi = chi_deg*np.pi/180.0
psi = psi_deg*np.pi/180.0

# 60 steps around ellipse (+1 because 0 and 2pi overlap)
N = 61
phi = np.linspace(0, 2*np.pi, N)
        
## Equation 4
xprime = r*np.cos(chi)*np.sin(phi)
yprime = r*np.sin(chi)*np.cos(phi)
        
## Inverse of Equations 2 and 3
x = xprime*np.cos(psi) - yprime*np.sin(psi)
y = xprime*np.sin(psi) + yprime*np.cos(psi)

# Set up the figure and the objects to hold plot data
plt.plot(x, y, color='k', linewidth=2)

# draw tangential arrow at pi/4 (one eighth of the way) around the ellipse
q=N//8
plt.arrow(x[q],y[q],x[q+1]-x[q],y[q+1]-y[q], head_width=hw, head_length=hl, fc='black', ec='black', length_includes_head=True)

# draw line through the major axis
q=N//4
u=(3*N)//4
linx=[x[q],x[u]]
liny=[y[q],y[u]]
plt.plot(linx, liny, color=major_color, linewidth=2)

# draw line through the minor axis
q=0
u=N//2
linx=[x[q],x[u]]
liny=[y[q],y[u]]
plt.plot(linx, liny, color=minor_color, linewidth=2)

# draw line connecting the ends of the minor and major axes in the third quadrant
q=(3*N)//4
u=N//2
linx=[x[q],x[u]]
liny=[y[q],y[u]]
plt.plot(linx, liny, color='k', linewidth=1)

# create the arc through psi at the origin
N = 10
arc_radius=0.35
phi = np.linspace(0, psi, N)

# reserve an extra point for the apex = (arcx[0], arcy[0])
arcx = np.linspace(0,1,N+1)
arcy = np.linspace(0,1,N+1)

arcx[1:] = arc_radius * np.cos(phi)
arcy[1:] = arc_radius * np.sin(phi)

plt.fill(arcx,arcy, color=psi_color, alpha=0.6)
plt.plot(arcx[1:], arcy[1:], color=major_color, linewidth=1)

# label the position angle
label_offset = 0.05
phimid=phi[2*N//5] # chosen by eye
xl = (arc_radius + label_offset) * np.cos(phimid)
yl = (arc_radius + label_offset) * np.sin(phimid)
plt.text(xl,yl,"$\psi$")

# create the arc through chi at the origin
N = 10
phi = np.linspace(psi-chi, psi, N)

arcx[1:] = x[q] + arc_radius * np.cos(phi)
arcy[1:] = y[q] + arc_radius * np.sin(phi)
arcx[0] = x[q]
arcy[0] = y[q]

plt.fill(arcx,arcy, color=chi_color, alpha=0.6)
plt.plot(arcx[1:], arcy[1:], color=minor_color, linewidth=1)

# label the ellipticity angle
phimid=phi[2*N//5] # chosen by eye
xl = x[q] + (arc_radius + label_offset) * np.cos(phimid)
yl = y[q] + (arc_radius + label_offset) * np.sin(phimid)
plt.text(xl,yl,"$\chi$")

if np.size(sys.argv) > 1:
    plt.savefig(sys.argv[1], bbox_inches='tight')
else:
    plt.show()