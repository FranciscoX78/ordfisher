"""Reproducible Cayley graph for S4 (adjacent rank swaps), with a Tutte planar embedding."""
from itertools import permutations
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D

G = nx.Graph()
for p in permutations(range(4)):
    G.add_node(p)
    for j in range(3):
        q = list(p)
        q[j],q[j+1] = q[j+1],q[j]
        G.add_edge(p,tuple(q),kind=j+1)
assert (G.number_of_nodes(),G.number_of_edges()) == (24,36)
planar, embedding = nx.check_planarity(G)
assert planar
seen=set();faces=[]
for u,v in embedding.edges():
    if (u,v) not in seen:
        face=embedding.traverse_face(u,v)
        faces.append(face)
        for a,b in zip(face,face[1:]+face[:1]):seen.add((a,b))
outer = next(f for f in faces if len(f)==6)
# Barycentric (Tutte) embedding: place a hexagonal face on the boundary
pos={p:np.array([1.8*np.cos(2*np.pi*i/6+.18),1.8*np.sin(2*np.pi*i/6+.18)]) for i,p in enumerate(outer)}
inner=[p for p in G if p not in pos]
M=np.zeros((len(inner),len(inner)));Y=np.zeros((len(inner),2));idx={p:i for i,p in enumerate(inner)}
for p in inner:
    i=idx[p];M[i,i]=len(list(G.neighbors(p)))
    for n in G.neighbors(p):
        if n in idx:M[i,idx[n]]-=1
        else:Y[i]+=pos[n]
sol=np.linalg.solve(M,Y)
for p in inner:pos[p]=sol[idx[p]]
# Ensure no close label centers
from itertools import combinations
D=sorted((np.linalg.norm(pos[a]-pos[b]),a,b) for a,b in combinations(G.nodes,2))
print('Min pairwise vertex separation:',D[0]);print('faces:',len(faces),[len(f) for f in faces]);print('outer:',outer)

fig,ax=plt.subplots(figsize=(10.2,8.0))
colors={1:'#8296AF',2:'#237C86',3:'#D28B50'}
for j in (1,2,3):
    nx.draw_networkx_edges(G,pos,edgelist=[(u,v) for u,v,d in G.edges(data=True) if d['kind']==j],edge_color=colors[j],width=1.8,ax=ax)
node_colors=['#EEF3F8' if sum(p[i]>p[j] for i in range(4) for j in range(i+1,4))%2==0 else '#E5F3F0' for p in G]
nx.draw_networkx_nodes(G,pos,node_size=750,node_color=node_colors,edgecolors='#344D61',linewidths=1.3,ax=ax)
nx.draw_networkx_labels(G,pos,labels={p:''.join(map(str,p)) for p in G},font_size=7.8,font_weight='medium',ax=ax)
ax.legend(handles=[Line2D([0],[0],color=colors[j],lw=2.8,label=f'Swap ranks {j}–{j+1}') for j in (1,2,3)],frameon=False,loc='lower right',fontsize=10)
ax.set_aspect('equal');ax.axis('off');fig.subplots_adjust(left=.04,right=.96,top=.98,bottom=.04)
d=Path(__file__).parent/'figures';d.mkdir(exist_ok=True)
fig.savefig(d/'figD1_Cayley_m4.pdf',bbox_inches='tight',pad_inches=.08)
fig.savefig(d/'figD1_Cayley_m4.png',bbox_inches='tight',pad_inches=.08,dpi=200)
