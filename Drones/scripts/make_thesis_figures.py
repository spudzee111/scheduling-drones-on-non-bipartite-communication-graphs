"""Vector thesis illustrations from the exact example and benchmark graphs."""
from pathlib import Path
from fractions import Fraction as F
import sys
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from scs.model import Graph
from check_structural_examples import penny_graph,square_family
OUT=ROOT.parent/"figures"
COLORS=["#087f8c","#c94f2e","#8463a0","#608c3d","#d59827","#3978ad"]
plt.rcParams.update({"font.family":"DejaVu Serif","font.size":10,"pdf.fonttype":42})
RADIUS=5/12

def draw(ax,g,labels=False,drones=False):
    centers=g.centers
    for i,j,a in g.edges:
        x,y=centers[i];u,v=centers[j]
        dx,dy=u-x,v-y;length=math.hypot(dx,dy)
        ax.plot([x+RADIUS*dx/length,u-RADIUS*dx/length],
                [y+RADIUS*dy/length,v-RADIUS*dy/length],
                color="#a1aab1",lw=1,zorder=0)
    for aid,a in enumerate(g.arcs):
        x,y=centers[a.circle]
        start=360*a.start/g.Q;end=start+360*a.length/g.Q
        ax.add_patch(Arc((x,y),2*RADIUS,2*RADIUS,theta1=start,theta2=end,
                         lw=2.4,color=COLORS[g.arc_ring[aid]%len(COLORS)]))
    if labels:
        for i,(x,y)in enumerate(centers):
            ax.text(x,y,str(i+1),ha="center",va="center",fontsize=12)
    if drones:
        for c,a,color in[(0,0,"#152536"),(2,F(1,2),"#152536")]:
            x,y=centers[c];t=2*math.pi*float(a)
            ax.plot(x+RADIUS*math.cos(t),y+RADIUS*math.sin(t),"o",
                    ms=4,color=color,zorder=4)
            t1=t-.18;t2=t-.55
            ax.annotate("",xy=(x+RADIUS*math.cos(t2),y+RADIUS*math.sin(t2)),
                        xytext=(x+RADIUS*math.cos(t1),y+RADIUS*math.sin(t1)),
                        arrowprops=dict(arrowstyle="-|>",lw=1.2,color=color))
    xs,ys=zip(*centers)
    ax.set_xlim(min(xs)-.6,max(xs)+.6);ax.set_ylim(min(ys)-.6,max(ys)+.6)
    ax.set_aspect("equal");ax.axis("off")

def save(fig,name):
    fig.savefig(OUT/(name+".pdf"),bbox_inches="tight",pad_inches=.03)
    fig.savefig(OUT/(name+".png"),dpi=180,bbox_inches="tight",pad_inches=.03)
    plt.close(fig)

def main():
    OUT.mkdir(exist_ok=True)
    fig,ax=plt.subplots(figsize=(4.6,4.3));draw(ax,penny_graph(),labels=True)
    from matplotlib.lines import Line2D
    ax.legend([Line2D([0],[0],color=COLORS[0],lw=3),
               Line2D([0],[0],color=COLORS[1],lw=3)],
              [r"$R_A:\ L_A=23/6$",r"$R_B:\ L_B=13/6$"],
              loc="lower center",bbox_to_anchor=(.5,-.08),ncol=2,frameon=False)
    save(fig,"penny_no")
    fig=plt.figure(figsize=(7,4.8));grid=fig.add_gridspec(2,2,height_ratios=[1,1.2],hspace=.3)
    for ax,r in[(fig.add_subplot(grid[0,0]),1),(fig.add_subplot(grid[0,1]),2),
                (fig.add_subplot(grid[1,:]),5)]:
        draw(ax,square_family(r),drones=True)
        ax.set_title(rf"$n={4*r},\quad k={r+1}$",pad=5)
    save(fig,"two_drone_family")
    fig,axes=plt.subplots(1,3,figsize=(7,2.4))
    coords={"triangle":[(0,0),(1,0),(.5,math.sqrt(3)/2)],
            "square":[(0,0),(1,0),(1,1),(0,1)],
            "diamond":[(0,0),(1,0),(.5,math.sqrt(3)/2),(1.5,math.sqrt(3)/2)]}
    for ax,name in zip(axes,["triangle","square","diamond"]):
        g=Graph.read(ROOT/"benchmarks"/"vse"/(name+".json"))
        # Input circles use the same index order; only display scale changes.
        g.centers=[(x/(12/5),y/(12/5))for x,y in g.centers]
        draw(ax,g);ax.set_title(name.capitalize())
    save(fig,"enumeration_instances")
    print("Wrote three vector figures and PNG previews.")
if __name__=="__main__":main()
