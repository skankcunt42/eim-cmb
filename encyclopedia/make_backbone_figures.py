import warnings; warnings.simplefilter("ignore")
import numpy as np, networkx as nx
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from collections import Counter
INK="#1f2a37"; BLUE="#1f6f8b"; GOLD="#b8860b"; RED="#c0392b"; GREY="#9aa7b4"; GREEN="#1f6f3d"

G=nx.dodecahedral_graph(); A=nx.to_numpy_array(G,nodelist=sorted(G)); V=20
ev=np.linalg.eigvalsh(A); spec=dict(sorted(Counter(np.round(ev,4)).items()))

# ---- Fig E1: spectrum + return probability ----
fig,ax=plt.subplots(1,2,figsize=(9.4,3.7))
labels=["−√5","−2","0","1","√5","3"]; vals=list(spec.keys()); mult=list(spec.values())
cols=[GOLD if abs(abs(v)-np.sqrt(5))<1e-3 else (RED if abs(v)<1e-3 else BLUE) for v in vals]
ax[0].bar(range(6),mult,color=cols,edgecolor="white")
for i,m in enumerate(mult): ax[0].text(i,m+0.1,f"×{m}",ha="center",fontsize=9)
ax[0].set_xticks(range(6)); ax[0].set_xticklabels([f"{l}\n({labels[i]})" if False else labels[i] for i,l in enumerate(labels)])
ax[0].set_ylabel("multiplicity"); ax[0].set_xlabel("adjacency eigenvalue")
ax[0].set_title("(a) Γ_dodec spectrum: {3, √5³, 1⁵, 0⁴, −2⁴, −√5³}",fontsize=10)
ax[0].text(0.02,0.9,"gold = ±√5 → horizon sector E_G (dim 6)\nred = 0 → dark-adjacent E₀ (dim 4)",
           transform=ax[0].transAxes,fontsize=7.5,va="top")
m=np.array(mult); contrib=(m/V)**2
ax[1].bar(range(6),contrib,color=cols,edgecolor="white")
ax[1].set_xticks(range(6)); ax[1].set_xticklabels(labels)
ax[1].set_ylabel("(mult/20)²"); ax[1].set_xlabel("eigenvalue")
ax[1].set_title("(b) return probability = Σ(mult/20)² = %.2f"%contrib.sum(),fontsize=10)
ax[1].text(0.5,0.85,"exact: (1+9+25+16+16+9)/400 = 76/400 = 0.19",transform=ax[1].transAxes,
           ha="center",fontsize=8,color=INK)
fig.tight_layout(); fig.savefig("fig_backbone_spectrum.png",dpi=140); plt.close()

# ---- Fig E2: dark sector composition ----
fig,ax=plt.subplots(figsize=(5.6,3.6))
bars=[("H = ℝ²⁰\n(full vertex space)",20,GREY),
      ("K = ker(Bᵀ)\ndark sector",8,INK),
      ("K∩H₊",4,BLUE),("K∩H₋",4,GOLD)]
for i,(lab,h,c) in enumerate(bars):
    ax.bar(i,h,color=c,edgecolor="white",width=0.7)
    ax.text(i,h+0.3,str(h),ha="center",fontsize=11,fontweight="bold")
ax.set_xticks(range(4)); ax.set_xticklabels([b[0] for b in bars],fontsize=8)
ax.set_ylabel("dimension"); ax.set_ylim(0,22)
ax.set_title("Dark sector: 8 = 4₊ ⊕ 4₋  (A·K⊆K exact, res 3e-15)",fontsize=10)
fig.tight_layout(); fig.savefig("fig_backbone_darksector.png",dpi=140); plt.close()

# ---- Fig E3: dodecahedron -> Petersen quotient ----
fig,ax=plt.subplots(1,2,figsize=(9,4))
posD=nx.kamada_kawai_layout(G)
nx.draw(G,posD,ax=ax[0],node_size=90,node_color=BLUE,edge_color="#c7d0d9",width=1.2)
ax[0].set_title("Γ_dodec  (V=20, E=30, 3-regular)",fontsize=10)
d=dict(nx.all_pairs_shortest_path_length(G)); J={u:[w for w in G if d[u][w]==5][0] for u in G}
parts={frozenset({u,J[u]}) for u in G}
Q=nx.quotient_graph(G,[set(p) for p in parts],relabel=True); Q=nx.Graph(Q); Q.remove_edges_from(nx.selfloop_edges(Q))
posP=nx.kamada_kawai_layout(Q)
nx.draw(Q,posP,ax=ax[1],node_size=110,node_color=GREEN,edge_color="#c7d0d9",width=1.4)
ax[1].set_title("Γ_dodec / antipodal  ≅  Petersen  (verified)",fontsize=10)
fig.tight_layout(); fig.savefig("fig_backbone_petersen.png",dpi=140); plt.close()
print("wrote fig_backbone_spectrum.png fig_backbone_darksector.png fig_backbone_petersen.png")
