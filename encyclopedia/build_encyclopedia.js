const fs=require("fs");
const {Document,Packer,Paragraph,TextRun,HeadingLevel,AlignmentType,Table,TableRow,TableCell,
       WidthType,BorderStyle,ShadingType,ImageRun}=require("docx");
const FONT="Georgia", MONO="Consolas";
const ENG="/root/eim/engine/", ENC="/root/eim/encyclopedia/";

function rns(text,o={}){const parts=[];const re=/(\*\*[^*]+\*\*|\*[^*]+\*)/g;let last=0,m;
  while((m=re.exec(text))!==null){if(m.index>last)parts.push(new TextRun({text:text.slice(last,m.index),font:FONT,...o}));
    const t=m[0]; if(t.startsWith("**"))parts.push(new TextRun({text:t.slice(2,-2),bold:true,font:FONT,...o}));
    else parts.push(new TextRun({text:t.slice(1,-1),italics:true,font:FONT,...o})); last=re.lastIndex;}
  if(last<text.length)parts.push(new TextRun({text:text.slice(last),font:FONT,...o}));return parts;}
const P=(t,o={})=>new Paragraph({alignment:o.align||AlignmentType.JUSTIFIED,spacing:{after:o.after??120,line:264},
  indent:o.indent,children:rns(t,{size:o.size||20})});
const H1=t=>new Paragraph({heading:HeadingLevel.HEADING_1,spacing:{before:260,after:120},
  children:[new TextRun({text:t,bold:true,font:FONT,size:26})]});
const H2=t=>new Paragraph({heading:HeadingLevel.HEADING_2,spacing:{before:180,after:100},
  children:[new TextRun({text:t,bold:true,font:FONT,size:22})]});
const CAP=t=>new Paragraph({spacing:{before:30,after:150},children:rns(t,{size:16,italics:true})});
function fig(p,w,h){return new Paragraph({alignment:AlignmentType.CENTER,spacing:{before:80,after:10},
  children:[new ImageRun({type:"png",data:fs.readFileSync(p),transformation:{width:w,height:h}})]});}

// entry metadata mini-table
function meta(rows){
  return new Table({width:{size:9360,type:WidthType.DXA},columnWidths:[1900,7460],
    rows:rows.map(r=>new TableRow({children:[
      new TableCell({width:{size:1900,type:WidthType.DXA},shading:{type:ShadingType.CLEAR,fill:"EEF2F6"},
        margins:{top:30,bottom:30,left:70,right:70},
        children:[new Paragraph({children:[new TextRun({text:r[0],bold:true,font:FONT,size:16})]})]}),
      new TableCell({width:{size:7460,type:WidthType.DXA},margins:{top:30,bottom:30,left:70,right:70},
        children:[new Paragraph({children:rns(r[1],{size:16})})]})]})),
    borders:{top:{style:BorderStyle.SINGLE,size:2,color:"CCCCCC"},bottom:{style:BorderStyle.SINGLE,size:2,color:"CCCCCC"},
      left:{style:BorderStyle.NONE},right:{style:BorderStyle.NONE},
      insideHorizontal:{style:BorderStyle.SINGLE,size:1,color:"DDDDDD"},insideVertical:{style:BorderStyle.SINGLE,size:1,color:"DDDDDD"}}});
}
function entryHead(id,name){return new Paragraph({spacing:{before:200,after:60},
  children:[new TextRun({text:id,bold:true,font:MONO,size:20,color:"1F6F8B"}),
            new TextRun({text:"  —  "+name,bold:true,font:FONT,size:20})]});}

const C=[]; const T=x=>C.push(x);

// ===== FRONT MATTER =====
T(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:60},
  children:[new TextRun({text:"THE EIM ENCYCLOPEDIA",bold:true,font:FONT,size:40})]}));
T(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:40},
  children:[new TextRun({text:"Seed Edition v0.1 — Consolidated, Verified, Plotted",italics:true,font:FONT,size:22})]}));
T(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:20},
  children:[new TextRun({text:"Thomas P. Connelly, Jr.  ·  Zero One Net",font:FONT,size:19})]}));
T(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:220},
  children:[new TextRun({text:"Canonical Core + Tip of the Spear  ·  compiled 18 July 2026",font:FONT,size:18})]}));

T(H1("Purpose and scope"));
T(P("This is the consolidation seed for the EIM corpus: one reference in which every result appears as a full, self-contained entry — statement, status, dependencies, an independent verification where one is possible, and a figure where one clarifies. It adopts, unchanged, the entry conventions and status taxonomy of the *EIM Framework Codex* (Canonical Edition); it does not invent a parallel format. Sources consolidated here are the *Framework Codex*, the *Codex & Encyclopedia Update*, *Coordination Before Identity*, the two dissertations, and the *Beyond Curved Spacetime* book."));
T(P("The edition is deliberately **narrow**: it does not attempt all ~220 catalogue entries at once. It builds one cluster to completion — the **dodecahedral backbone**, the exact finite laboratory on which the whole framework rests — with each entry independently re-verified in code and plotted, then carries that backbone to its empirical spearhead (the CMB five-fold / zonal test) and stops at the honest frontier. Two tiers are kept strictly separate, at your instruction: a **Canonical Core** of established (Level I) results, and a **Tip of the Spear** of open bridges and cutting-edge frontier items that must not be read as settled."));

T(H2("How to read an entry"));
T(P("Each entry carries the Codex fields: **ID** (type-prefixed: AX axiom, DF definition, TH theorem, ID identity, LM lemma, PR principle, CR corollary, OP open problem, CAND candidate, NEG negative result, NOTE, PD prediction, CJ conjecture); **Name**; **Type**; **Status** (taxonomy below); **Cluster**; **Statement** (canonical content, discipline-strict — no process-diary prose); **Depends on**; and, where applicable, a **Falsifiability hook**. This edition adds one non-Codex field, **Independent verification**, recording a re-computation performed for this consolidation (residuals, matched values) so the reader need not take a number on trust."));

T(H2("Status taxonomy (from the Codex, verbatim in force)"));
T(meta([
  ["Level I","Foundational; unconditional within the framework. Axioms, foundational primitives, forced theorems, exact algebraic identities."],
  ["Conditional Level I","Closed at Level I except for named external conditions (typically an Open Problem with a named closure target)."],
  ["Speculative-Level-III","Structural candidate with a named required proof; not derived, not axiomatically closed."],
  ["Programmatic","Open problem or research direction; carries no Level-I status claim."],
  ["Negative result","A hypothesis tested and excluded/refuted; preserved as constraint, not framework content."],
  ["Methodological commitment","Constitutive of research practice; not a physics claim; exempt from falsifiability."],
]));
T(P("**Forcing-language discipline.** \"Foundational filing\" = posited primitive, not derived. \"Forced given X\" = conditional forcing from named premises. \"Forced theorem\" = derives content from foundational primitives plus prior derived content. Every \"forced\" claim below names which of the three it is.",{after:80}));
T(P("**How to cite.** Connelly, T. P. (2026). *EIM Encyclopedia* (Seed v0.1), entry [ID]. Version-pin the entry ID; where an entry has multiple forms, name the form.",{after:120}));

// ===== PART I: CANONICAL CORE =====
T(H1("Part I — Canonical Core: The Dodecahedral Backbone"));
T(P("The dodecahedral graph Γ_dodec is the framework's *exact finite laboratory*: every result in this Part is a proven, exact statement about a fixed 20-vertex graph and its symmetry group, independently reproduced in code for this edition (residuals at machine precision). Nothing here is interpretive.",{after:140}));

T(entryHead("DF-DODEC","The dodecahedral substrate Γ_dodec"));
T(meta([["Type","definition (primitive)"],["Status","Level I"],["Cluster","Mathematical Substrate"],
 ["Depends on","—"],["Used by","DF/TH of the entire backbone; the CMB spearhead (Part II)"]]));
T(P("**Statement.** Γ_dodec is the dodecahedral graph: V = 20 vertices, E = 30 edges, F = 12 pentagonal faces, 3-regular, diameter 5, automorphism group |Aut| = 120 ≅ A₅ × ℤ₂. Its adjacency spectrum is exactly {3¹, √5³, 1⁵, 0⁴, (−2)⁴, (−√5)³}. Two spectral sectors are named: the **horizon sector** E_G = E₊√5 ⊕ E₋√5 (dim 6) and the **dark-adjacent** E₀ = ker(A) (dim 4)."));
T(P("**Independent verification.** Built Γ_dodec (networkx), confirmed 3-regular, diameter 5; adjacency eigenvalues reproduce the spectrum with multiplicities (1,3,5,4,4,3) exactly."));
T(fig(ENC+"fig_backbone_spectrum.png",600,236));
T(CAP("Figure DF-DODEC. (a) The exact adjacency spectrum with multiplicities. (b) The closed-system return probability (entry TH-DARK-FLOOR) as Σ(mult/20)² = 76/400 = 0.19."));

T(entryHead("TH-DARK-FLOOR","The exact dark-sector floor"));
T(meta([["Type","theorem (forced/computational)"],["Status","Level I"],["Cluster","Mathematical Substrate"],
 ["Depends on","DF-DODEC"],["Falsifiability hook","N/A (constitutive — exact linear algebra)"]]));
T(P("**Statement.** For the vertex–face incidence matrix B (20×12): rank(B) = 12 and dim ker(Bᵀ) = 8. The 8-dimensional space K := ker(Bᵀ) is the **dark sector**. The adjacency operator satisfies A·K ⊆ K exactly (residual ~10⁻¹⁶). The closed-system (unitary, H₀ = A) infinite-time-averaged return probability to a localized vertex state is **exactly 0.19**, by exact spectral degeneracy grouping."));
T(P("**Independent verification.** rank(B) = 12, dim ker(Bᵀ) = 8 confirmed; ‖A·K − Π_K(A·K)‖ = 3.14×10⁻¹⁵; and for a vertex-transitive graph the time-averaged return probability is Σ_λ (mult_λ/20)² = (1+9+25+16+16+9)/400 = 76/400 = **0.1900** (Figure DF-DODEC(b))."));

T(entryHead("TH-DARK-REP","Representation-theoretic origin of the dark sector"));
T(meta([["Type","theorem (forced)"],["Status","Level I"],["Cluster","Mathematical Substrate"],
 ["Depends on","DF-DODEC, TH-DARK-FLOOR"],["Used by","TH-PARITY, K-separation results (§1.13–1.15)"]]));
T(P("**Statement.** Under the rotation group (order 60 ≅ A₅), the vertex permutation representation decomposes as R²⁰ ≅ 1 ⊕ 3 ⊕ 3′ ⊕ 4 ⊕ 4 ⊕ 5. Because Bᵀ is A₅-equivariant and the face space R¹² carries no 4-dimensional irrep, Schur's lemma forces K = ker(Bᵀ) to contain the entire 4 ⊕ 4 isotypic component of R²⁰ — for representation-theoretic reasons alone, independent of any specific numerical entry of B. This is *theorem forcing*: content derived from the group structure plus the equivariance of B."));

T(entryHead("TH-CODIM","Exact codimension theorem for symmetry breaking"));
T(meta([["Type","theorem (forced/computational)"],["Status","Level I"],["Cluster","Mathematical Substrate"],
 ["Depends on","TH-DARK-FLOOR"],["Falsifiability hook","N/A (exact linear algebra)"]]));
T(P("**Statement.** Among the 20-dimensional space of diagonal perturbations D to A, the subspace that exactly preserves K = ker(Bᵀ) is **exactly 1-dimensional**, spanned by D ∝ I (a physically trivial uniform energy shift). Every other perturbation direction breaks the invariance A·K ⊆ K. The dark sector is thus maximally fragile: it survives only the trivial deformation."));

T(entryHead("TH-PARITY","Antipodal parity decomposition and the Petersen quotient"));
T(meta([["Type","theorem (forced/computational)"],["Status","Level I"],["Cluster","Mathematical Substrate"],
 ["Depends on","DF-DODEC, TH-DARK-FLOOR"],["Used by","horizon/holonomy cluster; projective-holonomy H¹(RP²;ℤ₂)"]]));
T(P("**Statement.** The antipodal map J (unique antipode at graph distance 5) gives an exact parity decomposition 20 = 10 + 10, H = H₊ ⊕ H₋, with a localized state satisfying ‖P₊|v⟩‖² = ‖P₋|v⟩‖² = ½. The dark sector splits **evenly** across parity: dim(K∩H₊) = dim(K∩H₋) = 4, i.e. 8 = 4₊ ⊕ 4₋. The dodecahedron-modulo-antipodal-identification quotient is isomorphic to the **Petersen graph** (verified by direct isomorphism check, not asserted)."));
T(P("**Independent verification.** Unique antipode confirmed for all 20 vertices; parity split dim(K∩H₊) = dim(K∩H₋) = 4 reproduced; the antipodal quotient graph is isomorphic to the Petersen graph."));
T(fig(ENC+"fig_backbone_petersen.png",520,231));
T(CAP("Figure TH-PARITY. Γ_dodec (left) and its antipodal quotient, verified isomorphic to the Petersen graph (right)."));
T(fig(ENC+"fig_backbone_darksector.png",360,231));
T(CAP("Figure TH-DARK. Dimensional anatomy of the dark sector: K = 8 splits as 4₊ ⊕ 4₋; A·K ⊆ K holds to residual 3×10⁻¹⁵."));

T(entryHead("TH-HORIZON-FRAME","Horizon frame bounds and the 3/2 channel gain"));
T(meta([["Type","theorem (forced/computational) + named negative"],["Status","Level I"],
 ["Cluster","Saturation, Horizons, and Bridges"],["Depends on","DF-DODEC (spectral sectors E_G, E₀)"]]));
T(P("**Statement.** With E_G (dim 6) and E₀ (dim 4) as in DF-DODEC, and the horizon H_v (radius-2 ball, exact antipodal transversal, |H_v| = 10, verified for all 20 poles): Σ_{x∈H_v} T_x T_x* = (3/20)I₄ and Σ_{x∈H_v} T_x* T_x = (1/10)I₆ exactly, ratio 3/2. Trace-preservation then forces the horizon channel Φ(I₆) = (3/2)I₄ — a consequence, not a free choice."));
T(P("**Named negative (kept as constraint).** The identical trace-preservation argument on the reverse channel E₀ → E_G gives Ψ(I₄) = (2/3)I₆, so det(Ψ(I₄)) = (2/3)⁶ ≈ 0.088 versus det(Φ(I₆)) = (3/2)⁴ ≈ 5.06 — a factor ~57 apart. Both directions cannot be simultaneously unit-gain trace-preserving CP maps; the asymmetry is exact, not numerical."));

T(entryHead("TH-NOGO-ORIENT","No graph-invariant antipodal orientation exists"));
T(meta([["Type","theorem (no-go)"],["Status","Level I"],["Cluster","Mathematical Substrate"],
 ["Depends on","DF-DODEC, TH-PARITY"],["Falsifiability hook","N/A (exhaustive over Aut(G))"]]));
T(P("**Statement.** Any graph-automorphism-equivariant scalar cost rule Φ(G, F, F̄) on an antipodal face pair must satisfy C_F = C_F̄ exactly (a 4-line proof via the antipodal automorphism's equivariance). Sharpened by exact classification of the signed 2-fold cover of the 10 antipodal axes: |Aut(G)| = 120; the unsigned axis action has exactly 60 distinct permutations (kernel {id, ι}); the signed action is faithful (120 distinct). Consequence: **no canonical, graph-invariant antipodal orientation exists** — any orientation requires external (flag) input. This is a load-bearing no-go for the framework's readout discipline."));

// ===== PART II: SPEARHEAD =====
T(H1("Part II — Empirical Spearhead: The CMB Five-Fold / Zonal Closure"));
T(P("The backbone's five-fold (pentagonal) symmetry makes a temperature-sector prediction for the large-angle CMB. This Part carries the backbone to data via two data-independent selection rules and one falsifiable observable. All three were verified in code for the closure paper (released engine); the selection rules are machine-precision, the observable is preregistered.",{after:140}));

T(entryHead("TH-ICO-SEL","Icosahedral selection: no power at ℓ ≤ 5"));
T(meta([["Type","theorem (forced/computational)"],["Status","Level I"],["Cluster","Physical Predictions"],
 ["Depends on","DF-DODEC"],["Falsifiability hook","N/A (exact harmonic analysis)"]]));
T(P("**Statement.** A function on the sphere invariant under the full icosahedral group I_h has vanishing multipoles at ℓ = 1,2,3,4,5; the lowest non-trivial invariant is ℓ = 6, then 10, 12. **Independent verification:** transforming twelve sources at the dodecahedral face-centre directions gives normalized power P(ℓ)/P_max ≤ 10⁻²⁷ for ℓ = 1–5 and O(1) at ℓ = 6, 10, 12 (direct-aₗₘ, cleaner than a pixelized 10⁻¹³ claim). Consequence: an exactly dodecahedral sky has no quadrupole and no octopole — so using ℓ = 2,3 to fix an axis already presupposes broken symmetry."));

T(entryHead("TH-KFOLD-SEL","k-fold azimuthal selection: the observed axis cannot be five-fold"));
T(meta([["Type","theorem (forced/computational)"],["Status","Level I"],["Cluster","Physical Predictions"],
 ["Depends on","TH-ICO-SEL"],["Used by","PD-CMB-ZONAL; OP-AXIS-SHARING (Tip of the Spear)"]]));
T(P("**Statement.** A function invariant under rotation by 2π/k about an axis has azimuthal content only in m ≡ 0 (mod k). About a five-fold axis this forces ℓ = 2,3,4 to be **zonal** (m = 0), first non-zonal power at m = 5, ℓ = 5. The observed axis-of-evil morphology is **sectoral** (planar octopole, m = ±3), which five-fold symmetry forbids. Verified: C₅/C₃/C₂-invariant test fields leak ≤ 10⁻²² into forbidden m. Consequence: an axis-locked five-fold structure on the *observed* low-ℓ axis is not a prediction of the hypothesis."));

T(entryHead("PD-CMB-ZONAL","The corrected falsifiable observable: the zonal five-fold signature"));
T(meta([["Type","prediction (falsifiable)"],["Status","Conditional Level I (preregistered; data-pending)"],
 ["Cluster","Physical Predictions"],["Depends on","TH-KFOLD-SEL"],
 ["Falsifiability hook","the §VI closure run on 4 Planck maps × 2 masks (released engine)"]]));
T(P("**Statement.** The genuine dodecahedral temperature signature is *zonal*, not sectoral: about the zonal axis (largest-eigenvalue power-tensor eigenvector), the joint (zonal_low, m5_mid) = (mean m=0 fraction at ℓ∈{2,3,4}, mean m=5 fraction at ℓ∈{5,6}) must jointly beat the matched null. A genuine zonal-five-fold sky sits at ≈(1,1); the sectoral axis-of-evil morphology and the ΛCDM null both sit in the low corner (full-sky reference (0.26,0.17)). A null result counts *against* the dodecahedral reading rather than being absorbed as a pre-existing anomaly."));
T(fig(ENG+"zonal_signature_fig4.png",320,296));
T(CAP("Figure PD-CMB-ZONAL. The (zonal_low, m5_mid) plane. Matched-null cloud with 95th-percentile decision lines; the genuine five-fold sky (★) occupies the empty upper-right, the sectoral morphology (×) and the null sit together in the low corner."));

// ===== TIP OF THE SPEAR =====
T(new Paragraph({pageBreakBefore:true,children:[]}));
T(H1("Tip of the Spear — Frontier (open / speculative; not settled)"));
T(P("Everything in this section is explicitly *not* Level I. These are the live bridges and cutting-edge items where the framework's reach currently ends. They are recorded so the frontier stays visible — and so nothing here is mistaken for the Canonical Core above.",{after:140}));

T(entryHead("OP-AXIS-SHARING","The undischarged axis-sharing bridge"));
T(meta([["Type","open problem (programmatic)"],["Status","Programmatic / Open"],["Cluster","Open Problems"],
 ["Depends on","TH-KFOLD-SEL"],["Closure target","a derivation forcing the five-fold axis to coincide with the low-ℓ axis"]]));
T(P("**Statement.** The five-fold ledger H₅ has statistical power only if the five-fold axis coincides with the low-ℓ (quadrupole–octopole) axis. TH-KFOLD-SEL shows the exact symmetry forbids this coincidence and a broken symmetry does not protect it. Absent a separate derivation supplying the coincidence, a surviving H₅ is an EIM-motivated *anomaly*, not a confirmation of the framework. This is the single bridge whose closure would convert the anomaly search into a framework test."));

T(entryHead("OP-ALPHA-DIMK","Why α(Γ_dodec) = dim K = 8"));
T(meta([["Type","open problem (structural)"],["Status","Programmatic / Open"],["Cluster","Open Problems"],
 ["Depends on","TH-DARK-FLOOR, K-separation (§1.13–1.15)"]]));
T(P("**Statement.** The independence number α(Γ_dodec) = 8 equals dim K = 8, but α is a purely combinatorial invariant while dim K is a linear-algebraic corank of the incidence operator; nothing yet shows the coincidence is forced rather than accidental. Recorded as a genuine open coincidence, not a result."));

T(entryHead("NEG-COSMO-BRIDGE","The 3:4:3 / √5 cosmological bridge is disfavored by data"));
T(meta([["Type","negative result (empirical)"],["Status","Negative result"],["Cluster","Physical Predictions"],
 ["Depends on","backbone spectral measure {+√5, 0, −√5}"],["Falsifiability hook","DESI DR2 BAO(+SN) fits — already applied"]]));
T(P("**Statement.** The compound bridge {3:4:3, √5} → Z_κ(a) = 2/5 + 3/5·cosh(κ√5·ln a) → ρ_DE(a) → H(z), tested directly against DESI DR2 BAO(+SN), is **disfavored relative to both ΛCDM and CPL**. A genuine empirical negative — recorded as a closed gate, not a tension to be explained away. (A near-central-value agreement for the DESY5 SN combination survives only as a weak, combination-specific coincidence.)"));

T(entryHead("FRONTIER-DECOHERENCE","Decoherence, Born weights, and basis selection (Claude's Frontier)"));
T(meta([["Type","open problem / frontier"],["Status","Programmatic / Open — measurement science"],
 ["Cluster","Projection and Readout"],["Depends on","OP-DFPROJ-OPERATOR (parameter-free A₅×ℤ₂-equivariant projector)"]]));
T(P("**Statement.** The projection/readout layer is open as measurement science: a parameter-free A₅ × ℤ₂-equivariant operator definition (OP-DFPROJ-OPERATOR) is unresolved, and Born weights, decoherence, and basis selection are not derived. Related computational findings on the substrate are recorded as *negative/open* (immediate flag-basis decoherence from localized states is threshold-like; whether a genuine flag *superposition* changes this is open). The standing discipline — that no observer-independent projector P_I is yet canonical, and that no bridge yet maps a real observation to an actual point on the twenty-vertex substrate — is the honest edge of the program. Nothing here may be cited as settled."));

// ===== PROVENANCE =====
T(H1("Provenance and roadmap"));
T(P("**Sources.** Framework Codex (Canonical Edition v2.0.6); Codex & Encyclopedia Update (Part I verified theorems; Parts VI, XIII–XVIII); Coordination Before Identity (frontier framing); the CMB closure paper and its released engine (selection rules, zonal signature). **Verification.** The backbone numbers (spectrum, K = 8, A·K residual, return probability 0.19, 4₊⊕4₋ split, Petersen quotient) and the selection-rule bounds were re-computed independently for this edition; residuals are quoted in each entry. **Roadmap (expansion passes).** Next clusters to bring to full-entry, verified form: the golden route sector W and the route-Clifford algebra (§1.6, J0² = −I); K-separation and the G_g ⊕ G_u labels (§1.13–1.15); the belt and flag-transport holonomy (§XIV–XVI); and the D2 information-viability theorems (Part VII). Each expansion keeps the same two-tier discipline: Canonical Core stays Level I; anything unproven goes to the Tip of the Spear."));

const doc=new Document({creator:"Zero One Net",title:"The EIM Encyclopedia — Seed Edition",
  styles:{default:{document:{run:{font:FONT,size:20}}}},
  sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1440,bottom:1440,left:1440,right:1440}}},children:C}]});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync("/root/eim/EIM_Encyclopedia_Seed_v0.1.docx",b);
  console.log("wrote EIM_Encyclopedia_Seed_v0.1.docx",b.length,"bytes");});
