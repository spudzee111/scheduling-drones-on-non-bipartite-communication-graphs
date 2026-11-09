"""Exact checks for the two structural examples in the thesis."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scs.model import Graph
from scs.dynamics import Engine
from run_vse_experiments import grid_properties, reference_properties

PENNY_EDGES = [(0,1,F(0)),(0,3,F(1,4)),(1,2,F(1,4)),
               (2,3,F(1,2)),(0,4,F(5,6)),(1,4,F(2,3)),
               (0,5,F(5,12)),(3,5,F(7,12))]
# Each coordinate a+b*sqrt(3) is a pair (a,b); these centers are unscaled.
PENNY_CENTERS = [((0,0),(0,0)),((1,0),(0,0)),((1,0),(1,0)),
                 ((0,0),(1,0)),((F(1,2),0),(0,F(-1,2))),
                 ((0,F(-1,2)),(F(1,2),0))]

def sign(a,b):
    if b == 0: return (a>0)-(a<0)
    if a == 0: return (b>0)-(b<0)
    if a>0 and b>0: return 1
    if a<0 and b<0: return -1
    t=a*a-3*b*b
    return ((t>0)-(t<0)) if a>0 else -((t>0)-(t<0))

def geometric_check():
    intended={(i,j):a for i,j,a in PENNY_EDGES}
    cos=[(1,0),(0,F(1,2)),(F(1,2),0),(0,0),
         (F(-1,2),0),(0,F(-1,2)),(-1,0),(0,F(-1,2)),
         (F(-1,2),0),(0,0),(F(1,2),0),(0,F(1,2))]
    for i,p in enumerate(PENNY_CENTERS):
        for j,q in enumerate(PENNY_CENTERS):
            if j<=i: continue
            vec=[(F(c)-F(a),F(d)-F(b))for(a,b),(c,d)in zip(p,q)]
            a=sum(x*x+3*y*y for x,y in vec)
            b=sum(2*x*y for x,y in vec)
            assert sign(a-F(25,36),b)>0  # physical distance >2, scale 12/5
            linked=sign(a-1,b)<=0  # distance <=12/5 = 2+rho, rho=2/5
            assert linked==((i,j)in intended)
            if linked:
                assert (a,b)==(1,0)
                h=int(intended[i,j]*12)
                assert vec==[cos[h],cos[(3-h)%12]]
            else:
                assert sign(a-2,b)>=0  # nonedges >=sqrt(2)

def penny_graph():
    import math
    centers=[(float(x[0])+float(x[1])*math.sqrt(3),
              float(y[0])+float(y[1])*math.sqrt(3))for x,y in PENNY_CENTERS]
    return Graph(6,PENNY_EDGES,name="six_circle_no",centers=centers)

def square_family(r):
    centers=[(2*j+x,(-1)**j*y)for j in range(r)for x in(0,1)for y in(0,1)]
    bearings={(1,0):F(0),(-1,0):F(1,2),(0,1):F(1,4),(0,-1):F(3,4)}
    edges=[]
    for i,(x,y)in enumerate(centers):
        for j,(u,v)in enumerate(centers):
            if j>i and (u-x)**2+(v-y)**2==1:
                edges.append((i,j,bearings[u-x,v-y]))
    return Graph(4*r,edges,name=f"alternating_squares_{r}",centers=centers)

def link_coordinate_table(graph):
    rec={}
    for r,arcs in enumerate(graph.rings):
        q=F(0)
        for aid in arcs:
            ar=graph.arcs[aid];s=graph.arc_sign[aid]
            q+=F(ar.length,graph.Q)
            port=ar.end if s==1 else ar.start
            e=graph.ports[ar.circle][port]
            rec.setdefault(e,[]).append((r,q%F(graph.ring_lengths[r],graph.Q),ar.circle))
    return rec

def main():
    geometric_check()
    g=penny_graph()
    lengths=[F(x,g.Q)for x in g.ring_lengths]
    assert lengths==[F(23,6),F(13,6)]
    div=F(gcd(*g.ring_lengths),g.Q)
    rows=[]
    C={1:set(),-1:set()};E={1:set(),-1:set()}
    for edge,coordinates in link_coordinate_table(g).items():
        if coordinates[0][0]==coordinates[1][0]:continue
        (r,a,u),(s,b,v)=sorted(coordinates)
        ce=1 if u==v else -1
        residues={1:(a-b)%div,-1:(a+b)%div}
        for direction in (1,-1):
            (C if direction==ce else E)[direction].add(residues[direction])
        i,j,_=g.edges[edge]
        rows.append(dict(link=[i+1,j+1],z_A=str(a),z_B=str(b),c=ce))
    assert all(C[s]==E[s]=={F(1,12)}for s in (1,-1))
    family=[]
    for r in range(1,17):
        h=square_family(r)
        assert len(h.rings)==r+1
        assert sorted(h.ring_lengths)==sorted([8,8]+[16]*(r-1))
        initial=(h.encode(0,0,-1),h.encode(2,F(1,2),-1))
        result=Engine(h).simulate(initial)
        assert result.connected
        missing,same=grid_properties(h,result)
        assert not missing and same
        assert F(result.period_ticks,h.Q)==4*r-2
        independent=reference_properties(h,initial)
        assert independent["period"]==4*r-2
        assert independent["covered_lengths"]==[F(1)]*(4*r)
        family.append(dict(squares=r,n=4*r,rings=r+1,period=4*r-2,
                           coverage=True,collision_free=True,connected=True))
    out=dict(penny=dict(lengths=list(map(str,lengths)),gcd=str(div),
                       link_coordinates=rows,safe_residues_equal=[],safe_residues_opposite=[],
                       exact_geometry_verified=True),
             family=family)
    path=ROOT/"results"/"vse"/"structural_examples.json"
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(out,indent=2)+"\n")
    print("PASS: exact six-circle geometry/residues; 16 family members independently replayed.")
if __name__=="__main__":main()
