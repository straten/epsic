import sys
import numpy as np
import matplotlib.pyplot as plt

psi_deg = 35.0
chi_deg = 25.0
r = 0.9

chi = chi_deg*np.pi/180.0
psi = psi_deg*np.pi/180.0

# https://www.datylon.com/blog/data-visualization-for-colorblind-readers
psi_color = "#e57a77"
chi_color = "#7ca1cc"

major_color = "#f05039"
minor_color = "#1f449c"

# axis label offset
label_offset=0.05

plt.rc('font', family='serif', serif='cm10', size=28)
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
plt.text(1.0,label_offset,"$x$")
plt.text(label_offset,1.0,"$y$")

# draw the x' and y' axes with empty heads

# equation 2
xpx = np.cos(psi) 
xpy = np.sin(psi)

# equation 3
ypx = -np.sin(psi)
ypy = np.cos(psi)

# dashed lines
plt.plot([-xpx, xpx], [-xpy, xpy], color='black', dashes=[8, 4])
plt.plot([-ypx, ypx], [-ypy, ypy], color='black', dashes=[8, 4])

# emtpy arrow heads
plt.arrow(xpx,xpy,hl*xpx,hl*xpy, fill=False,
        head_width=hw, head_length=hl, fc='black', ec='black', length_includes_head=True)
plt.arrow(ypx,ypy,hl*ypx,hl*ypy, fill=False,
        head_width=hw, head_length=hl, fc='black', ec='black', length_includes_head=True)

# label the x' and y' axes
xpx *= 1+hl
xpy *= 1+hl
ypx *= 1+hl
ypy *= 1+hl

plt.text(xpx-xpy*label_offset,xpy+xpx*label_offset,"$x'$")
plt.text(ypx+ypy*label_offset,ypy-ypx*label_offset,"$y'$")

# draw the ellipse

# 60 steps around ellipse (+1 because 0 and 2pi overlap)
N = 61
phi = np.linspace(0, 2*np.pi, N)
        
## Equation 4
xprime = r*np.cos(chi)*np.sin(phi)
yprime = r*np.sin(chi)*np.cos(phi)
        
## Inverse of Equations 2 and 3
x = xprime*np.cos(psi) - yprime*np.sin(psi)
y = xprime*np.sin(psi) + yprime*np.cos(psi)

# Draw the ellipse
plt.plot(x, y, color='k', linewidth=2)

# draw tangential arrow at pi/4 (one eighth of the way) around the ellipse
q=N//8
plt.arrow(x[q],y[q],x[q+1]-x[q],y[q+1]-y[q], head_width=1.5*hw, head_length=1.5*hl, shape='full', overhang=0.25, fc='black', ec='black', length_includes_head=True)

# draw red line through the major axis
q=N//4
u=(3*N)//4
linx=[x[q],x[u]]
liny=[y[q],y[u]]
plt.plot(linx, liny, color=major_color, linewidth=2)

# draw blue line through the minor axis
q=0
u=N//2
linx=[x[q],x[u]]
liny=[y[q],y[u]]
plt.plot(linx, liny, color=minor_color, linewidth=2)

# draw line of length r that connects the ends of the minor and major axes in the third quadrant
# this line marks the ellipticity angle, chi
q=(3*N)//4
u=N//2
linx=[x[q],x[u]]
liny=[y[q],y[u]]
plt.plot(linx, liny, color='k', linewidth=1)

# draw two dotted lines perpendicular to the above line, out to where the length r can be labelled
# perpendicular offset vector
dx = (liny[1] - liny[0])
dy = (linx[0] - linx[1])
offset = 0.3 
dx_offset = dx * offset
dy_offset = dy * offset

# draw the first dotted line
plt.plot(
    [linx[0], linx[0] + dx_offset],
    [liny[0], liny[0] + dy_offset],
    linestyle=':',
    color='k'
)

# draw the second dotted line
plt.plot(
    [linx[1], linx[1] + dx_offset],
    [liny[1], liny[1] + dy_offset],
    linestyle=':',
    color='k'
)

# add the label "r" at the midpoint between the ends of the dotted lines
midx = 0.5 * (linx[1]+linx[0]) + dx_offset
midy = 0.5 * (liny[1]+liny[0]) + dy_offset

plt.text(midx, midy, "$r$", ha='center', va='center')

# draw a pair of arrows to mark the length

# parallel offset vector
dy = (liny[1] - liny[0])
dx = (linx[1] - linx[0])

space = 0.08
llen = 0.5 - space

plt.arrow(midx+space*dx, midy+space*dy, llen*dx, llen*dy, 
          head_width=hw*.6, head_length=hl*.6, fc='black', ec='black', length_includes_head=True)

plt.arrow(midx-space*dx, midy-space*dy, -llen*dx, -llen*dy,
          head_width=hw*.6, head_length=hl*.6, fc='black', ec='black', length_includes_head=True)

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
