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
const ENC2="/root/eim/encyclopedia/";

// ===== TITLE =====
T(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:60},
  children:[new TextRun({text:"THE EIM ENCYCLOPEDIA",bold:true,font:FONT,size:34})]}));
T(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:40},
  children:[new TextRun({text:"Expansion I — Decoherence, Superposition, and the Interior / Exterior Cosmos",italics:true,font:FONT,size:22})]}));
T(new Paragraph({alignment:AlignmentType.CENTER,spacing:{after:200},
  children:[new TextRun({text:"Thomas P. Connelly, Jr.  ·  Zero One Net  ·  compiled 18 July 2026",font:FONT,size:17})]}));

T(H1("The architecture, in one page"));
T(P("This expansion treats the framework's largest claim — that the classical world seen *from the interior* and the undecided superposition *of the exterior* are two faces of one coordination structure — with the same two-tier discipline as the seed: what is proven goes in the Canonical Core; the cosmological model is stated as the sharp, named conjecture it already is in the corpus, in the Tip of the Spear. Nothing here is inflated, and the one thing that is a theorem is re-verified in code."));
T(P("The spine is the **existence curve**. Coordination depth σ = M/(E+I) runs from a pre-percolation substrate (σ → 0, p = 0⁺) through the percolation threshold (σ = 1) into the classical regime (σ ≫ 1). The book fixes the reading used throughout: the percolation threshold σ = 1 is \"the abrupt onset of global constraint: **decoherence**, the Standard Model threshold, and the emergence of classical geometry.\" Below it, local holonomy-bearing loops exist but do not span the system — a locally fertile, globally undecided substrate. Above it, commitment is classical and irreversible."));
T(P("Two consequences frame everything below. First, the **interior/exterior duality**: a committed, localized configuration (a definite flag) is classical — it has decohered and can never re-develop coherence — while the un-percolated substrate and the superposed alternatives it has not yet selected are the exterior the interior cannot see into. Second, the **self-renewing ensemble**: global saturation (σ → ∞ everywhere) and the null state ∅ are the *same* excluded condition approached from opposite ends, so the cosmos cannot wind down; black holes recycle Memory (a terminal black-hole state is structurally a local Big Bang), and \"what appears as cosmic expansion history is the aggregate statistical behavior of a vast ensemble of locally renewing coordination domains.\" The question *how likely is our cosmos* is a measure over that ensemble — which is exactly where it becomes the cosmological measure problem, not a slogan."));

// ===== CANONICAL CORE =====
T(H1("Canonical Core — the one thing that is proven"));

T(entryHead("TH-FLAG-DECOHERENCE","Basis-preserving decoherence of a localized interior configuration"));
T(meta([["Type","theorem (forced/computational)"],["Status","Level I (localized initial state)"],
 ["Cluster","Projection and Readout"],["Depends on","DF-DODEC; the 120-flag system and its σ-operations"],
 ["Falsifiability hook","N/A (exact: permutation matrices are monomial)"]]));
T(P("**Statement.** The substrate's local dynamics on the 120 complete flags (v, e, F) of Γ_dodec is generated by the flag involutions σ₀, σ₁, σ₂, which are **permutation (monomial) matrices**: each maps a basis flag to a single basis flag, never to a superposition. Consequently a system prepared in a *localized* flag state evolves as a classical probability distribution over flags: it can never develop off-diagonal coherence in the flag basis. In this exact sense a committed interior configuration is *permanently decohered* — classicality is not imposed, it is structurally unavoidable for localized states. This is *theorem forcing*: the conclusion is forced by the monomial character of the generators plus the localized initial condition."));
T(P("**Independent verification.** Reconstructed all 120 flags from Γ_dodec (12 ordered pentagonal faces), built σ₁, σ₂ as explicit 120×120 permutations (confirmed σ₁² = σ₂² = I, entries in {0,1}), and ran the dephasing channel ρ → (1−ε)ρ + (ε/2)(P₁ρP₁ + P₂ρP₂). A localized flag holds ℓ₁-coherence at **0.000** for all steps (Figure). The result is exact, not asymptotic."));
T(fig(ENC2+"fig_decoherence.png",430,255));
T(CAP("Figure TH-FLAG-DECOHERENCE. Green: a localized interior flag never develops coherence (ℓ₁-coherence ≡ 0) — permanent classicality. Gold: a genuine flag superposition keeps ℓ₁-coherence at 1.0 under the same permutation dynamics — coherence is relabeled, never destroyed. The gap between the two curves is the interior/exterior seam, and the gold curve is the open problem below."));
T(P("**Scope (load-bearing).** The theorem is stated for *localized* initial states only. It does **not** say a superposition decoheres — the verification shows the opposite (next entry). Classicality of the interior is forced; classicality of the exterior is not, and that asymmetry is the whole subject of the Tip of the Spear.",{after:80}));

// ===== TIP OF THE SPEAR =====
T(new Paragraph({pageBreakBefore:true,children:[]}));
T(H1("Tip of the Spear — the cosmological model (open / interpretive; not settled)"));
T(P("Everything below is explicitly not Level I. These are the named bridges whose closure would turn the interior/exterior picture from a compelling architecture into a derived result. They are stated sharply so each can be attacked individually.",{after:140}));

T(entryHead("OP-SUPERPOSITION-DECOH","Does an interior superposition decohere? (the exterior sector)"));
T(meta([["Type","open problem (structural)"],["Status","Programmatic / Open"],["Cluster","Projection and Readout"],
 ["Depends on","TH-FLAG-DECOHERENCE"],["Closure target","a mechanism (beyond monomial dynamics) that decoheres a flag superposition, or a proof that none does"]]));
T(P("**Statement.** The permutation dynamics that forces classicality on localized states **conserves** the coherence of a genuine superposition (verified: ℓ₁-coherence held at 1.000 across all steps — the superposition is relabeled across flags, never dephased). So nothing in the substrate's own dynamics decoheres the exterior superposed sector. Whether any admissible EIM process does — a measurement/readout coupling, a saturation event, an environment of other flags — is untested and open. This is the precise formal content of \"the interior is decoherence, the exterior is superposition\": the interior half is a theorem, the exterior half is this open problem. Its resolution decides whether the interior/exterior duality is dynamical or merely definitional."));

T(entryHead("PR-DECOH-PERCOLATION","Decoherence as the percolation threshold σ = 1"));
T(meta([["Type","principle (interpretive bridge)"],["Status","Conditional / interpretive"],["Cluster","Saturation, Horizons, and Bridges"],
 ["Depends on","existence curve σ = M/(E+I); S² closure obstruction"],["Falsifiability hook","the Hubble-tension two-regime prediction (below) inherits this identification"]]));
T(P("**Statement.** On the existence curve, the percolation threshold σ = 1 is identified with the onset of decoherence and classical geometry: below σ = 1 the substrate carries local, non-spanning holonomy loops (the undecided/superposition-like regime); at σ = 1 global constraint switches on and commitment becomes classical and irreversible. **Status discipline:** the percolation transition itself is a computed property of the backbone (the S² closure obstruction forces a defect set; the 12 pentagonal faces are that set), but its *identification with quantum decoherence* is an interpretive bridge, not a derived equivalence. It earns conditional standing only through the observational positions it forces on the curve — the CMB/Standard-Model onset at σ = 1, rotation-curve anomalies at σ ≲ 1, and the two-measurement-regime reading of the Hubble tension — which are falsifiable and are where this entry must be defended or dropped."));

T(entryHead("OP-EXT-IRRECOVER","Saturation ⇒ exterior irrecoverability (the sealing half)"));
T(meta([["Type","open problem (programmatic)"],["Status","Programmatic / Open (Codex item 63)"],["Cluster","Saturation, Horizons, and Bridges"],
 ["Depends on","DF-HOR (saturation locus); TH-HORIZON-FRAME"],["Closure target","a filed argument connecting Σ_local → κ saturation to external-record irrecoverability + directionality"]]));
T(P("**Statement.** A stationary (no-hair) black hole is maximally *sealed*: the interior configuration is unrecoverable from the external record. In EIM the horizon is primitively a saturation boundary (Σ_local → κ) and curvature is its continuum rendering. The claim that saturation **forces** exterior irrecoverability — the formal statement of the interior/exterior information barrier — is a named open problem (Codex item 63), not proven: it requires connecting saturation of the accumulation measure to irreversibility of the external readout and to a directional (interior-to-exterior) arrow, none of which is currently supplied. This is the sealing half of the duality; TH-FLAG-DECOHERENCE is the classicality half. Both must close for \"the interior view is decoherence\" to be a theorem rather than a picture."));

T(entryHead("OP-COSMOS-ENSEMBLE-MEASURE","How likely is our cosmos: the ensemble and its measure"));
T(meta([["Type","open problem (programmatic)"],["Status","Programmatic / Open"],["Cluster","Physical Predictions"],
 ["Depends on","self-renewing existence curve; ensemble-of-domains reading"],["Closure target","a normalizable measure over the ensemble of coordination domains"]]));
T(P("**Statement.** The cosmos is modeled as an ensemble of locally-renewing coordination domains, each tracing its own existence curve, with cosmic expansion history the aggregate statistic; global saturation and ∅ are excluded, so the ensemble self-renews rather than winding down. Within this picture the question *how likely is our cosmos* is well-posed only relative to a **measure** over the ensemble — and defining a normalizable such measure is the cosmological measure problem, which is unsolved in general and is where every multiverse/ensemble program does its real work or quietly fails. **Honest positioning:** EIM does not dissolve this problem; it *relocates* the anthropic likelihood question onto it, with one asset a bare multiverse lacks — a structural constant (the backbone spectral gap γ_c ≈ 0.0353) that constrains the recycling rate and ties terminal-evaporation and lightest-neutrino-mass scales together, giving the ensemble a non-free parameter. **Speculative sub-note (Level III):** rescaling the interaction horizon (the role of c) across domains — \"c on a scale\" — is a candidate route to inter-domain (\"other-cosmos\") structure; recorded as a direction, not a result."));

T(entryHead("NOTE-GODEL-BELOW-THRESHOLD","The two horizons: operating below the diagonalization threshold, not around Gödel"));
T(meta([["Type","note (methodological / interpretive)"],["Status","Methodological commitment"],["Cluster","Core Framework"],
 ["Depends on","Gödel–Bell two-horizons taxonomy"],["Falsifiability hook","N/A (constitutive discipline)"]]));
T(P("**Statement.** The framework names two limits on self-rendering: the **vertical** (Gödel–Tarski) limit — no system fully renders its own rendering operation from within the same rendering — and the **horizontal** (Bell) limit — two local observer histories cut into one originless whole. The disciplined claim, stated in *Coordination Before Identity*, is that EIM **does not solve or circumvent Gödel**; it operationalizes coordination *below the self-reference threshold*. This is exact: Gödel incompleteness and Tarski undefinability require a system expressive enough to represent its own syntax and diagonalize; a system below that threshold (the Presburger-arithmetic boundary — addition without multiplication) is consistent, complete, and decidable *because* it cannot self-refer. The interior/exterior architecture is the same point wearing a cosmological hat: an interior cannot compute its own likelihood from within (the vertical limit), and the exterior ensemble is the external vantage from which such a measure could be defined (OP-COSMOS-ENSEMBLE-MEASURE). **Guardrail:** \"gets around Gödel\" must never be written; \"operates below the diagonalization threshold\" is the true — and stronger — claim, and the only one that survives review."));

// ===== PROVENANCE =====
T(H1("Provenance and roadmap"));
T(P("**Sources.** Beyond Curved Spacetime (Ch. 8 Modal Process Ontology; Ch. 9 The Self-Renewing Universe — existence curve, ensemble of domains, spectral regulator γ_c); Framework Codex (item 63 saturation/irrecoverability; the Gödel–Bell two-horizons taxonomy; OP-DFPROJ-OPERATOR); Codex & Encyclopedia Update (basis-preserving decoherence; the open superposition case); Coordination Before Identity (below-threshold / Presburger discipline). **Verification.** The 120-flag reconstruction, the σ-involutions, and the localized-vs-superposition coherence curves were computed independently for this edition (localized ℓ₁-coherence ≡ 0; superposition ≡ 1). **Roadmap.** The two closures that would promote this section: (i) OP-EXT-IRRECOVER — a filed saturation ⇒ irrecoverability argument; (ii) OP-SUPERPOSITION-DECOH — an admissible process that dephases a flag superposition. Closing either turns a face of the interior/exterior duality from picture into theorem; closing both makes the decoherence–superposition cosmos a derived structure rather than a compelling one."));

const doc=new Document({creator:"Zero One Net",title:"EIM Encyclopedia — Expansion I",
  styles:{default:{document:{run:{font:FONT,size:20}}}},
  sections:[{properties:{page:{size:{width:12240,height:15840},margin:{top:1440,bottom:1440,left:1440,right:1440}}},children:C}]});
Packer.toBuffer(doc).then(b=>{fs.writeFileSync("/root/eim/EIM_Encyclopedia_ExpansionI_Decoherence.docx",b);
  console.log("wrote EIM_Encyclopedia_ExpansionI_Decoherence.docx",b.length,"bytes");});
