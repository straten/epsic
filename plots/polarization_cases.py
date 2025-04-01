import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

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

def setup (ax, label, psi, chi):
    ax.axis('off')
    ax.set_aspect('equal')

    # draw the x and y axes
    hw, hl = 0.06, 0.10
    ax.arrow(-1,0,2,0, head_width=hw, head_length=hl, fc='black', ec='black', length_includes_head=True)
    ax.arrow(0,-1,0,2, head_width=hw, head_length=hl, fc='black', ec='black', length_includes_head=True)

    # label the x and y axes
    ax.text(1.0,label_offset,"$x$")
    ax.text(label_offset,1.0,"$y$")
    ax.text(1.0,1.0,label)

    r=0.8
    hw, hl = 0.1, 0.15

    x = r*np.cos(np.radians(psi))
    y = r*np.sin(np.radians(psi))
    
    if chi == 0:
        print(f"{label} x={x} y={y}")
        ax.arrow(0,0,x,y, lw=3, head_width=hw, head_length=hl, shape='full', overhang=0.25, fc=chi_color, ec=chi_color, length_includes_head=True)
        ax.arrow(0,0,-x,-y, lw=3, head_width=hw, head_length=hl, shape='full', overhang=0.25, fc=chi_color, ec=chi_color, length_includes_head=True)
    else:
        circle = patches.Circle((0, 0), radius=r, fill=False, ec=chi_color, lw=3)
        ax.add_patch(circle)
        v_x = 0.5 * hl * np.sin(np.radians(psi)) * np.sign(chi)
        v_y = -0.5 * hl * np.cos(np.radians(psi)) * np.sign(chi)
        ax.arrow(x,y,v_x,v_y, lw=3, head_width=hw, head_length=hl, shape='full', overhang=0.25, fc=chi_color, ec=chi_color, length_includes_head=True)


plt.rc('font', family='serif', serif='cm10', size=28)
plt.rc('text', usetex=True)

# Create a figure and a 2x3 grid of subplots.
fig, axes = plt.subplots(3, 2, figsize=(8, 12)) 

# Set up the figure axes, etc.
plt.xlim(-1.1, 1.1)
plt.ylim(-1.1, 1.1)

setup(axes[0,0],"$+Q$",0,0)
setup(axes[0,1],"$-Q$",90,0)
setup(axes[1,0],"$+U$",45,0)
setup(axes[1,1],"$-U$",-45,0)
setup(axes[2,0],"$+V$",45,45)
setup(axes[2,1],"$-V$",-45,-45)

if np.size(sys.argv) > 1:
    plt.savefig(sys.argv[1], bbox_inches='tight')
else:
    plt.show()
