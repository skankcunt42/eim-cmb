"""Verify the basis-preserving decoherence result on the real 120-flag dodecahedral system,
and demonstrate the open contrast: a localized flag can never develop coherence (immediate/
permanent decoherence), while a genuine superposition is NOT decohered by the permutation
dynamics (coherence conserved) — the exact open question flagged in the Encyclopedia."""
import warnings; warnings.simplefilter("ignore")
import numpy as np, networkx as nx
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

G=nx.dodecahedral_graph(); ok,emb=nx.check_planarity(G)
# ordered pentagonal faces
faces=[]; seen=set()
for u in G:
    for v in G[u]:
        f=emb.traverse_face(u,v)
        if len(f)==5 and frozenset(f) not in seen:
            seen.add(frozenset(f)); faces.append(f)
assert len(faces)==12, len(faces)

# build 120 flags (v, edge, face_index)
def edges_of(face):
    return [frozenset((face[i],face[(i+1)%5])) for i in range(5)]
flags=[]
for fi,F in enumerate(faces):
    for e in edges_of(F):
        for v in e: flags.append((v,e,fi))
assert len(flags)==120, len(flags)
idx={fl:i for i,fl in enumerate(flags)}

edge_faces={}
for fi,F in enumerate(faces):
    for e in edges_of(F): edge_faces.setdefault(e,[]).append(fi)

def sigma1(fl):  # other edge of the face at v
    v,e,fi=fl; es=[x for x in edges_of(faces[fi]) if v in x]; e2=[x for x in es if x!=e][0]; return (v,e2,fi)
def sigma2(fl):  # other face sharing edge e
    v,e,fi=fl; f2=[x for x in edge_faces[e] if x!=fi][0]; return (v,e,f2)
def perm(func):
    P=np.zeros((120,120))
    for fl in flags: P[idx[func(fl)],idx[fl]]=1.0
    return P
P1,P2=perm(sigma1),perm(sigma2)
print("[flags] |flags|=%d  σ1²=I:%s  σ2²=I:%s  permutation(monomial):%s"%(
    len(flags), np.allclose(P1@P1,np.eye(120)), np.allclose(P2@P2,np.eye(120)),
    np.all((P1==0)|(P1==1))))

def channel(rho,eps=0.5):
    return (1-eps)*rho + (eps/2)*(P1@rho@P1) + (eps/2)*(P2@rho@P2)
def l1coh(rho):
    return float(np.sum(np.abs(rho))-np.sum(np.abs(np.diag(rho))))

# localized initial flag
loc=np.zeros((120,120),complex); loc[idx[flags[0]],idx[flags[0]]]=1.0
# superposition of two flags
a,b=idx[flags[0]],idx[flags[37]]
psi=np.zeros(120,complex); psi[a]=psi[b]=1/np.sqrt(2); sup=np.outer(psi,psi.conj())

T=25; cl=[l1coh(loc)]; cs=[l1coh(sup)]; rl,rs=loc.copy(),sup.copy()
for t in range(T):
    rl=channel(rl); rs=channel(rs); cl.append(l1coh(rl)); cs.append(l1coh(rs))
print("[localized]     coherence: start %.3f  end %.3f  (theorem: stays 0 — cannot form)"%(cl[0],cl[-1]))
print("[superposition] coherence: start %.3f  end %.3f  (OPEN: permutation dynamics does NOT decohere it)"%(cs[0],cs[-1]))

fig,ax=plt.subplots(figsize=(6.4,3.8))
ax.plot(range(T+1),cl,"-o",ms=3,color="#1f6f3d",label="localized flag  (interior → classical)")
ax.plot(range(T+1),cs,"-s",ms=3,color="#b8860b",label="flag superposition  (exterior → preserved)")
ax.set_xlabel("coordination step"); ax.set_ylabel("ℓ₁ coherence (Σ|off-diagonal|)")
ax.set_ylim(-0.05,1.15); ax.set_title("Basis-preserving decoherence on the 120-flag system",fontsize=11)
ax.legend(fontsize=8.5); ax.grid(alpha=0.25)
ax.text(0.98,0.5,"permutation (monomial) dynamics:\nlocalized state can never develop coherence;\na superposition's coherence is conserved,\nnot destroyed — the open exterior sector",
        transform=ax.transAxes,ha="right",va="center",fontsize=7.2,color="#334")
fig.tight_layout(); fig.savefig("fig_decoherence.png",dpi=140)
print("wrote fig_decoherence.png")
