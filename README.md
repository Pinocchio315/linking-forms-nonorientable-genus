# Smooth and topological nonorientable four-genera of low-crossing knots: code and data

This repository contains the programs and the recorded results of the computations of the paper [C] *Smooth and
Topological Nonorientable Four-Genera of Low-Crossing Knots* by Seong-Jin Lee.  [C] uses the results of the theory paper
[T] *Linking-Form Splittings and Four-Dimensional Genera of Knots*, arXiv:2610.40183.  Theorem, section and table
numbers below refer to the version of [C] of November 2026 and to arXiv:2610.40183v1.

- `scripts/` holds the programs.  Each script states its purpose, its usage and its output in its docstring.
- `data/` holds their outputs.  These are the records on which the counts, tables and theorems of [C] rest.
- `scripts/bands.py` (band moves on PD codes and the Smith normal form), `scripts/goeritz.py` (Goeritz matrices) and
  `scripts/open_sets.py` (the sets of knots of each stage) are modules that other scripts import.
- `vendor/` holds eight modules from the author's earlier code that some scripts import (see `vendor/README.md`).

## Checking [C] against the data

```
python scripts/paper_counts.py
```

This recomputes every count quoted in [C], and those of the section Algorithms and checks below, from `data/`.  Each
count is checked against the quoted value by an assertion, so the script stops at the first difference.  It prints
the counts under the names of the macros of [C].

The searches `mobius_open_search.py`, `donaldson_surfaces.py`, `surgery_drill.py` and `mobius_J_search.py` take their
knots from `scripts/open_sets.py`, which runs `paper_counts.py` only up to the set that each of them needs; so each of
them needs only the outputs of the stages before it, in the order given in the docstring of `open_sets.py`.

## Results, scripts and data

### In the order of [C]

| Result | Script | Data |
|---|---|---|
| [C] §2.1 (for [T] Examples 3.5 and 3.6 and Figure 2): the knots K and K′ built from the pretzel knot P(3,−3,1) in the chirality of [T], their Goeritz matrices (with the η of Gordon and Litherland, sign included), their identification as 8₁₀ and 8₁₁ up to mirror image, and their PD codes in the convention of KnotInfo, as quoted in [C] §2.1 | `small_examples.py` | `small_examples.txt` |
| [C] §3.1: KnotInfo values and ranges of γ₄ | `knotinfo_extract.py` | `gamma4_knotinfo_2026.10.5.json`, `open_gamma4_knotinfo_2026.10.5.json` |
| [C] Proposition 3.1(a) and §3.5: h_λ and h^⊕ for all prime knots with at most 13 crossings, and the congruences of [T] Lemma 2.5 | `linking_genus.py` | `linking_genus.json` |
| [C] Proposition 3.1(a) (h_λ = 3), rank two, and the binary forms of §3.6 | `rank2.py` | `rank2_targets.json`, `rank2_controls.json`, `rank2_checks.json` |
| [C] Proposition 3.1(b), rank one: the linking-form test | `rank1_linking.py`, `rank1_goeritz.py` | `rank1_targets.json`, `rank1_controls.json`, `rank1_validate.json`, `rank1_goeritz.json` |
| [C] Proposition 3.1(c) with Proposition 2.6: correction terms | `rank1_dinv.py` | `rank1_dinv_targets.json`, `rank1_dinv_mtargets.json`, `rank1_dinv_controls.json`, `rank1_dinv_mcontrols.json` |
| [C] Proposition 3.2(a): the double branched covers of the 316 open non-alternating knots that are not Montesinos, classified with the pipeline of the follow-up paper arXiv:2610.25430 (its modules in `vendor/`), and the correction terms of the 59 negative-definite plumbing trees with at most two bad vertices (57 knots excluded); controls; band searches for γ₄^top = 1 | `plumbing_dinv.py` | `plumbing_dinv.json`, `plumbing_dinv_gap.json`; input `plumbing_dinv_routeP/` (the verdicts of an independent implementation, compared in mode `summary`) |
| [C] Proposition 3.2(b) and §3.8 (11n_176): surgery descriptions Σ_K = ±S³_{D/2}(C) of the 50 open knots with unknotting number one (Montesinos trick), calibrated on 32 alternating knots, and the test of Proposition 2.6 for every V allowed by the bounds of Hom–Wu and Rasmussen (13 knots excluded) | `surgery_dinv.py`, `surgery_box.py` | `surgery_dinv.json`, `surgery_dinv_gap.json`, `surgery_box.json` |
| [C] Proposition 3.2(b), exact correction terms: for the other 33 open knots with unknotting number one, the V_j(C) determined from the UV = 0 knot Floer complex of C (Szabó's program), by the gradings of the hat complexes and over every lift to F₂[U,V], with the chirality of C fixed by an isometry, Regina's names, the Khovanov rank or finite covers; the test of Proposition 2.6 excludes 13n_3844 and 13n_3855 and passes for the other 31; controls (calibration, the 17 knots already excluded, 11n_176, nine connected sums of torus knots); mode gap: bands to knots J with Δ_J = 1 for the two excluded knots (Theorem B) | `surgery_cfk.py` | `surgery_cfk.json`, `surgery_cfk_gap.json` |
| [C] Proposition 3.2(c): Donaldson's theorem with definite checkerboard surfaces (Proposition 2.7) for the 231 knots left open by Proposition 3.2(a), (b) and the second search that are neither alternating nor Montesinos and have σ + 4 Arf ≡ ±2 (mod 8): the KnotInfo diagram, the stored diagrams of the two Möbius searches, the stored pilot diagrams and new diagrams from about 1,500 random sequences of Reidemeister moves per knot (a few hundred distinct diagrams); γ₄ = 2 for 64 knots, each certificate diagram stored and checked by a second program (diagram of K with the right chirality, Goeritz matrix recomputed, linking form against V + Vᵀ, independent embedding search); controls: 229 knots with Möbius bands, all 216 definite lattices embed | `donaldson_surfaces.py`, `donaldson_check.py`, `open_sets.py` | `donaldson_surfaces.json`, `donaldson_check.json`, `donaldson_pilot_diagrams.json` |
| [C] §3.2: surgery descriptions Σ_K = S³_{p/q}(C) by drilling short closed curves, for the 143 knots of Proposition 3.2(c) whose Σ_K is hyperbolic and bounds no negative-definite plumbing tree and whose KnotInfo unknotting number is not one: Σ_K drilled along at most 30 dual curves, an S³ filling recognised by Regina; descriptions for 21 knots (q = 1, 3, 5), V_j(C) for 8, and the test of Proposition 2.6 excludes 13n_2626, 13n_2632 and 13n_3616, which Proposition 3.2(c) excludes as well | `surgery_drill.py` | `surgery_drill.json` |
| [C] §3.2, the other 13 drilling descriptions: V_j(C) from the UV = 0 knot Floer complex over every lift (as in `surgery_cfk.py`); the test of Proposition 2.6 passes for 11, and for 2 the knot Floer homology of C was not computed within the time limit, so none is excluded; the 8 drilling knots with a verdict, recomputed as controls, keep it | `drill_cfk.py` | `drill_cfk.json` |
| [C] Proposition 3.3 and §3.8: surfaces built from band moves (punctured Klein bottles for 1,998 + 122 knots, surfaces with b₁ = 3 for 13a_2890 and 11n_133, Möbius bands for 4 + 905 knots in the first search; mode `mobius` searches the 1,280 knots of range [1,2] or [1,3] that pass both tests), and their independent rebuild | `klein_bottle_search.py`, `klein_bottle_rebuild.py` | `klein_bottle_search.json`, `klein_bottle_search_ranges.json`, `klein_bottle_search_mobius.json`, `klein_bottle_rebuild.json`, `klein_bottle_rebuild_mobius.json` |
| [C] §3.3: replay of the first search, mode `mobius`, on its stored diagrams (all stored results reproduced; see Algorithms and checks) | `klein_bottle_search.py mobius replay` | `klein_bottle_search_mobius_replay.txt` (printed output; the replay writes no data file) |
| [C] Proposition 3.3(a), the second search (§3.3): on the 305 knots left open by the first search and Propositions 3.1 and 3.2(a), (b); it also looks for the ribbon disc of J in other diagrams of J when SnapPy identifies J with a knot that KnotInfo lists as slice (the identification only guides the search); Möbius bands for 12 knots, all in diagrams of mode `mobius`, rebuilt and related to the KnotInfo diagrams by Reidemeister moves alone; 32 controls without a band; the search was stopped after 165 to 279 diagrams per knot left open; the Fox–Milnor condition for the 4 of the 12 with h_λ = 0 and signature 0 (it fails for 13a_4729) | `mobius_open_search.py`, `mobius_open_rebuild.py`, `fox_milnor.py wide` | `mobius_open_search.json.gz`, `mobius_open_rebuild.json`, `fox_milnor_wide.json` |
| [C] Proposition 3.3(a), the second construction continued after Proposition 3.2(c) on the 229 knots then left open: a ribbon disc looked for for every knot J of the bands of the second search in its stored diagrams (the J grouped by determinant, Alexander polynomial and volume, the groups ordered and pruned by necessary conditions for sliceness, then one band in up to 300 diagrams of J reached by Reidemeister moves); Möbius bands for 13a_1412, 13a_2898, 13a_3125 and 13n_3901, each rebuilt by `mobius_open_rebuild.py` and diagrammatic; no band for the 10 controls; 16 knots stopped at the time limit of 600 s; the Fox–Milnor condition holds for 13n_3901 (`fox_milnor.py wide` reads this search as well) | `mobius_J_search.py`, `mobius_open_rebuild.py` (REBUILD_IN), `fox_milnor.py wide` | `mobius_J_search.json.gz`, `mobius_J_search.groups.json.gz` (stages 1 and 2), `mobius_J_rebuild.json`, `fox_milnor_wide.json` |
| [C] §3.4: the Fox–Milnor condition for the 21 knots with a Möbius band of the first search, h_λ = 0 and signature 0 (Alexander polynomials by Fox calculus, compared with KnotInfo and factored in PARI and SymPy); it fails for 5 of them, which are therefore not topologically slice | `fox_milnor.py` | `fox_milnor.json` |
| [C] §2.3 and §3.5: comparison with the lattice criterion of Jabuka and Kelly | `jk_lattice.py` | `jk_targets.json`, `jk_controls.json`, `jk_ls14.json` |
| [C] Table 3: certificates of h_λ = 3 | `hthree_certificates.py` | `hthree_certificates.json` |
| [C] §3.6 and Table 3 (for [T] Theorem 4.6): Jabuka's rational Witt span ws, from the rational Witt class alone, for the knots with at most 13 crossings and the 196 knots with 14 crossings and h_λ = 3 | `witt_span.py` | `witt_span.json` |
| [C] Theorem D (§3.7): knots with 14 crossings | `rank1_ls14.py`, `rank2_ls14.py`, `rank1_dinv.py` | `rank1_ls14.json`, `rank2_ls14.json`, `rank1_dinv_ls14.json`, `rank1_dinv_ls14controls.json` |
| [C] §3.8: 11n_57, whose published value 1 the correction terms exclude, rechecked with the plumbing of the Seifert invariants that Regina recognises | `montesinos_dinv_check.py literature` | `montesinos_dinv_check_literature.json` |
| [C] Theorem B: the twelve knots of Proposition 3.1(c) with γ₄^top = 1 < γ₄ = 2 | `smooth_top_gap.py`, `smooth_top_gap_check.py`, `montesinos_dinv_check.py` | `smooth_top_gap.json`, `smooth_top_gap_check.json`, `montesinos_dinv_check.json` |
| [C] Theorem B: the 29 knots of Proposition 3.2(a), (b) with γ₄^top = 1, the knot Floer homology of the knot J of each band used (τ(J) = 0 for four of them) | `gap_new.py` | `gap_new.json` |
| Every count quoted in [C] | `paper_counts.py` | all of the above |

The usage lines of the scripts are:

```
# [C] Section 2.1: the knots K and K' of [T] Examples 3.5 and 3.6
python scripts/small_examples.py

# [C] Section 3: data, lower bounds, upper bounds, values and controls
python -I scripts/knotinfo_extract.py path/to/knotinfo_data_complete.csv [out_dir]
python scripts/linking_genus.py [workers]
python scripts/rank2.py controls | targets | checks [workers]
python scripts/rank1_linking.py validate | controls | targets [workers]
python scripts/rank1_goeritz.py [workers]
python scripts/rank1_dinv.py controls | targets | mcontrols | mtargets [workers]
python scripts/plumbing_dinv.py records path/to/negative-definite-plumbings | recognise | jsj | normal | algebraic | identify | dtest | controls | gap [workers] | summary
sage -python scripts/plumbing_dinv.py hyp | piece [workers]
python scripts/surgery_dinv.py calibration | targets | pilot | gap [workers] | summary
sage -python scripts/surgery_box.py
python scripts/surgery_cfk.py [workers] | gap [workers]
python scripts/donaldson_surfaces.py targets | controls [workers]   # DS_ONLY=knot,knot,... and DS_NEW=n restrict a run
python scripts/donaldson_check.py [workers]
python scripts/surgery_drill.py [workers] [redo]
python scripts/drill_cfk.py [workers]
python scripts/klein_bottle_search.py [ranges | mobius] [workers] [earlier_output.json ...] | [ranges | mobius] replay [workers]
python scripts/klein_bottle_rebuild.py [workers] [search_output.json ...]
python scripts/mobius_open_search.py [fresh] [workers] | status [workers] | finalize
python scripts/mobius_open_rebuild.py [workers]
python scripts/mobius_J_search.py [workers] | status
REBUILD_IN=data/mobius_J_search.json.gz REBUILD_OUT=data/mobius_J_rebuild.json python scripts/mobius_open_rebuild.py [workers]
python scripts/fox_milnor.py [wide]
python scripts/jk_lattice.py controls | targets [workers]
python scripts/hthree_certificates.py
python scripts/witt_span.py [workers]
python scripts/rank1_ls14.py path/to/knot_3-14data.csv [workers]
python scripts/rank2_ls14.py [workers]
python scripts/rank1_dinv.py ls14 | ls14controls path/to/knot_3-14data.csv [workers]
python scripts/jk_lattice.py ls14 path/to/knot_3-14data.csv [workers]
python scripts/montesinos_dinv_check.py [literature]

# [C] Section 4: locally flat but not smooth Möbius bands
python scripts/smooth_top_gap.py [workers] [earlier_output.json ...] | replay [workers]
python scripts/smooth_top_gap_check.py
python scripts/gap_new.py

# every count quoted in [C] and in this README
python scripts/paper_counts.py

```

Run them from the repository root.  They write to `data/` unless stated otherwise.

## Algorithms and checks

This section describes the computations of [C] and their checks in more detail than the paper does.  `paper_counts.py`
recomputes the counts quoted here as well.

### Linking forms

Jordan decompositions are computed by diagonalising V + Vᵀ, or a Goeritz matrix, over Z₍ₚ₎ in exact rational arithmetic,
as in [T] §2.1: at each step an entry of least p-adic valuation is chosen, and if no diagonal entry attains it, a basis
vector eᵢ is replaced by eᵢ + eⱼ.  The invariant h_λ is then given by [C] Theorem 2.2, with (R_h) evaluated through the
list of oddities of [T] Lemma 2.4.

The rank-one program (`rank1_linking.py`) decides, level by level, whether λ_K ≅ ⟨τ/ℓ⟩ ⊕ μ with μ metabolic; it refuses
inputs with an exceptional prime, where the orthogonal-summand condition is stronger than the Witt-level one ([T]
Proposition 4.9).  For the 7,893 knots of determinant below 200 it was compared with a search over all elements x ∈ H_K
with λ_K(x, x) = ±1/ℓ and, for each, over the subgroups of x^⊥ generated by at most two isotropic elements, which
contain a metaboliser whenever one exists for these determinants.

The rank-two program (`rank2.py`) enumerates a representative of every class of binary forms of determinant ±ℓ: reduced
forms in the definite case, and in the indefinite case the forms ax² + 2bxy + cy² with |2b| ≤ |a| ≤ √ℓ, together with
the isotropic forms when ℓ is a square.  This bound holds if |a| is the least nonzero absolute value of a form that
does not represent 0, after a shear that makes |2b| ≤ |a|: then |c| ≥ |a|, so ac < 0, since otherwise ℓ = b² − ac < 0,
and ℓ = b² + |ac| ≥ a².  For each odd form it decides prime by prime whether λ_Q is an orthogonal summand of λ_K with
metabolic complement.

### Correction terms

For an alternating knot, the correction terms of ±Σ_K = ∂X_G are the maxima of (c² + rk G)/4 over the characteristic
covectors c of the negative-definite Goeritz lattice G in each class, computed by enumeration (`rank1_dinv.py`).  For a
Montesinos knot the same is done for the star-shaped plumbing lattice; the plumbing is built from the Montesinos
notation with the code of arXiv:2609.09861 (in `vendor/`).  The independent check of the three Montesinos knots of [C]
Theorem B, and of 11n_57 ([C] §3.8), builds the plumbing from the Seifert invariants that Regina recognises and
enumerates characteristic covectors with the Fincke–Pohst algorithm in exact rational arithmetic
(`montesinos_dinv_check.py`).

### Plumbing trees and surgery descriptions

`plumbing_dinv.py` classifies the double branched covers of the knots of [C] Proposition 3.2 with the pipeline of
arXiv:2610.25430 (five of its modules are in `vendor/`, and the other functions that it needs are copied into the
script): the cover is recognised with Regina on random retriangulations, hyperbolicity is verified with SnapPy's
interval arithmetic in SageMath, and the normal forms of both orientations are computed from the recognised Seifert
pieces or, where Regina names no decomposition, from an algebraic diagram of K.  Every tree is identified with ±Σ_K by
triangulating its plumbed manifold, and the correction terms are computed and tested as in [C] §2.2.

`surgery_dinv.py` carries out the steps of [C] Proposition 3.2(b) for the knots with unknotting number one, including
the check of the labelling on lens spaces; the 32 calibration knots of [C] §3.5 are a fixed list in it.
`surgery_box.py` runs the test for all V between the bounds of [C] §2.2 on the surgery descriptions of slope ±D/2 (a
drilling route of another slope is left to `surgery_drill.py`), and `surgery_cfk.py` computes the exact V_j(C),
runs the test with them and searches bands for the two knots that they exclude.  `surgery_drill.py` runs the drilling
routes of `surgery_dinv.py`, with their checks, and the test of [C] Proposition 2.6 with both implementations.  The band
searches of [C] §4 for the knots of Proposition 3.2 use the code of `smooth_top_gap.py`, and `gap_new.py` recomputes the
knot Floer homology of the knot J of each band used in the proof of [C] Theorem B.

### Definite surfaces

`donaldson_surfaces.py` carries out the search of [C] Proposition 3.2(c) and its controls.  It uses the Goeritz matrices
of `rank1_goeritz.py`, whose sign is fixed on the KnotInfo diagram of 3_1: its signature is −2, and the matrix of its
Seifert surface is positive definite, so the program's matrices are −𝒢_V.  The embedding search is that of
`jk_lattice.py`, which places the rows of the lattice one at a time and keeps only partial embeddings that are canonical
under the symmetries of the standard lattice.  The new diagrams are made as in `mobius_open_search.py`.  The random
moves of spherogram cannot be repeated from a seed, so every certificate diagram is stored, and the diagrams of
`donaldson_pilot_diagrams.json` are scanned as stored diagrams.

`donaldson_check.py` checks the certificates and uses no code of the search: it traces the faces and forms the Goeritz
matrices again, computes LDLᵀ decompositions over Q, compares the matrices with the stored lattice by PARI's `qfisom`,
checks the chirality of the diagram by an isometry of exteriors in SnapPy and the sign by the linking form of V + Vᵀ
(with the Jordan data of `rank1_linking.py`), and repeats the embedding search over vectors in canonical form with its
own code.

### Bands and Alexander polynomials

Band moves are carried out on PD codes: a band in a face joining two of its edges, with t crossings between its long
sides, is inserted as in [C] Figure 1.  In the searches the determinant of a knot J is computed from the smaller of its
two Goeritz matrices; for every J with Δ_J = 1 both Goeritz matrices were checked to present the trivial group.
Alexander polynomials are computed by Fox calculus, with exact determinants in PARI; for the Fox–Milnor condition of [C]
§3.4 they are factored over Z in PARI.  This routine agrees with the Seifert matrices for all 12,965 prime knots with at
most 13 crossings and with knot Floer homology for 663 of them, and |Δ(−1)| agrees with the Goeritz determinant for all
79,352 knots on which it was compared.  Knots are identified with SnapPy and compared with the KnotInfo diagrams by knot
Floer homology; diagrams are compared with KnotInfo knots by isometry of their complements.

### Surfaces built from band moves

The first search of [C] §3.3 (`klein_bottle_search.py`) tries the bands in the faces of the KnotInfo diagram with |t|
≤ 6.  For the 49 knots with h_λ = 3, which serve as controls, and for the knots of [C] Proposition 3.3(a) and (b)
without a surface in the KnotInfo diagram, it also tries further diagrams made as in [C] §4; the controls were scanned
in 9,312 diagrams in all.  A knot J produced by a band is compared, through its determinant, its Alexander polynomial
and SnapPy's identification or its simplified diagram, with the knots that KnotInfo lists with value 1, and for
13a_2890 and 11n_133 also with those that it lists with value 2 and that have h_λ ≤ 2.  For the chosen bands the
Alexander polynomial of J was checked with a Seifert matrix and its knot Floer homology with that of the named knot,
the dual band was checked to give back K, and the band was rebuilt with the band operation of SnapPy's link module.
The bands of J and of J′ are searched in the diagrams that the earlier bands produce, so every surface is a sequence
of band moves on stored PD codes ending in a diagram that Reidemeister moves reduce to the unknot or to the
two-component unlink.

`klein_bottle_rebuild.py` shares no code with the search.  It rebuilds from the KnotInfo PD codes all 7,236 band records
that the search stored, and the surfaces of 13a_2890, 11n_133 and the four knots of [C] Proposition 3.3(a) of range
[1,4] or [1,5].  Of the records, 7,038 are complete surfaces, ending in the unknot or the two-component unlink; each of
the other 198 stores a band move to a knot J together with a surface of J in another diagram of J (93 records) or with
no surface of J (105 records).  Of these 198, only the surface for 13a_3126 is used in a proof, so every surface that a
proof uses is rebuilt completely.  The program performs each band move in two independent ways and checks that the
resulting diagrams are isomorphic, recognises the final unknot or unlink by Reidemeister moves alone, and relates the
further diagrams of K, and the diagram of J used for 13a_3126, to the diagrams that they come from by Reidemeister
moves; so these surfaces depend on the Reidemeister moves of spherogram, but on no identification of knots by SnapPy.

For the Möbius bands of [C] Proposition 3.3(a) and the last group of (b), the same search (mode `mobius`) looks, in the
knots of range [1,2] or [1,3] that pass both tests, for a band to a knot J′ followed by a band from J′ to the
two-component unlink, or for a band to the unknot, first in the KnotInfo diagram and then in further diagrams.  It
predicts the determinant of each band move from det K as in [C] §4, builds only the bands that can give a knot of square
determinant, and checks two properties of the unlink, a vanishing minor of the Fox matrix and the determinant 1 of each
component, before the Reidemeister moves; the 484 knots that a test excludes are searched in the same way.  The knots of
range [1,3] without a Möbius band are then searched for punctured Klein bottles as in (b), and `klein_bottle_rebuild.py`
rebuilds all these Möbius bands and punctured Klein bottles as well.  A replay of this search on its stored diagrams
reproduced the results for all 1,764 knots, the Möbius bands, the punctured Klein bottles and the rebuilt surfaces of
the knots J alike; only the number of hits identified through a simplified diagram, which no result uses, changed for
two knots.

### The second search for Möbius bands

`mobius_open_search.py` carries out the second search of [C] §3.3.  It scans each knot in rounds of forty diagrams,
first the stored diagrams of the first search and then six rounds of new diagrams (one round of each for the control
knots), and it stores the results of every round, and every new diagram together with the diagram that it was made from,
in a `.partial` file from which a stopped run resumes.  The search was stopped when the knots without a
Möbius band had been scanned in five to seven rounds, at most two of them of new diagrams, and its mode `finalize` wrote
`mobius_open_search.json.gz` from the stored rounds.  `mobius_open_rebuild.py` uses no code of the searches: it rebuilds
each of these Möbius bands with the functions of `klein_bottle_rebuild.py` and relates every diagram that they use to
the KnotInfo diagram of K, or to the diagram of J made by the first band, by Reidemeister moves.  `fox_milnor.py wide`
tests the Fox–Milnor condition for the knots with h_λ = 0 and signature 0 among them.

### Witt span

The values of ws in [C] §3.6 were recomputed without linking forms (`witt_span.py`).  The matrix V + Vᵀ, or for 14
crossings a Goeritz matrix, which has the same residues because the moves between spanning surfaces change the rational
Witt class by [⟨±1⟩] at most (S. Jabuka, *The concordance crosscap number and rational Witt span of a knot*,
arXiv:2012.14801v1, proof of Prop. 3.3), is diagonalised over Q in exact arithmetic, and the residues are read off from
the diagonal entries.  Forms of rank 0 and 1 are tested directly.  For rank 2 the program runs over the squarefree odd
determinants d made of a sign and primes of D, and at ∞, at 2 and at each prime of D over the local forms ⟨x, dx⟩, x a
square class, with the required residue; at the other primes such a form has Hasse invariant 1.  A rational form with
these localisations exists if and only if the product of the Hasse invariants is 1 (J. H. Conway and N. J. A. Sloane,
*Sphere Packings, Lattices and Groups*, 3rd ed., Ch. 15, (15), Thm. 5 and §5.3); the Hilbert symbols at 2 are given by
(7) of Jabuka's paper.  The values agree with [T] Theorem 4.6 for all 12,965 knots with at most 13 crossings and for the
196 knots with 14 crossings and h_λ = 3, and with ws(7_5) = ws(9_40) = 2 and ws(9_49) = 3 from Jabuka's Example 3.5.
For each of the 17 knots with ws = 2 and h_λ = 3, the program also finds an even binary form as in the proof of [T]
Theorem 4.6, for instance the form with rows (6, 3) and (3, −6) for 8_18 ([T] Example 4.7).  As a check, the rank-two
program, which enumerates the classes of binary forms, finds among the 49 knots with h_λ = 3 a binary form whose
discriminant form is an orthogonal summand of λ_K with metabolic complement exactly for these 17 knots, and every such
form is even; for the 6 knots that satisfy only (iii) of [C] Corollary 2.3(c), with σ(K) ≡ 4 (mod 8), it finds none.

### Figures

The knot diagrams of [C] Figure 2 are drawn from the KnotInfo PD codes and the band records in `data/` by the figure
program of the manuscripts, which is not part of this repository: the orthogonal layout of spherogram places crossings
and bends on a grid, and the smoothing of SnapPy's link editor turns the polygonal arcs into Bézier curves.  The program
checks that the four strands at every crossing appear counterclockwise in the order of the PD code, and labels every
edge by its PD label.  Each band is drawn inside the face named in the data, between the two named edges, with its t
crossings.

## Environment

The computations were run with Python 3.13.7 and the following packages:

- `snappy` 3.3.2 with `snappy_manifolds` 1.4, and its link module `spherogram` 2.4.1;
- `database_knotinfo` 2026.8.1;
- `cypari` 2.5.6 with PARI 2.15.4;
- `numpy` 2.5.2 and `sympy` 1.14.0;
- `knot_floer_homology` 1.2.2 (Szabó's program);
- `regina` 7.4, which recognises the double branched covers, the Seifert fibred spaces and the S³ fillings of
  `plumbing_dinv.py`, `surgery_dinv.py`, `surgery_drill.py` and `montesinos_dinv_check.py`;
- SageMath 10.7 with `snappy` 3.3.2 installed in it (`sage -pip install --user snappy==3.3.2`), for the scripts run
  with `sage -python`, among them the verified computations of SnapPy.

## Input data

Two input files are not redistributed here.  The scripts that need them take their path as an argument.

| Data | Source | SHA-256 |
|---|---|---|
| KnotInfo, version 2026.10.5: `knotinfo_data_complete.csv` | https://knotinfo.org, downloaded 2026-10-15 | `eb511ebc61204bc1c6d700d9da42257dcbb0369cb6756e89ed978a7821867c7d` |
| Lee and Sabloff, table of knots with at most 14 crossings: `knot_3-14data.csv` | https://github.com/LeeJoPing/Mobius-Band-Search, downloaded 2026-10-17 | `357c855fc76561723fd5bc35e0c1aefc0ee64d144f37d8ce4db331afff79db84` |

Seifert matrices, PD codes and Montesinos notations come from the package `database_knotinfo` 2026.8.1.  That
package has values of γ₄ only up to 10 crossings and for the non-alternating knots with 11 crossings, so all values
and ranges of γ₄ are taken from the data file of version 2026.10.5.  `knotinfo_extract.py` turns that file into `gamma4_knotinfo_2026.10.5.json` (SHA-256
`5b4388359b62336972a7bf54813cca59e9719b75d3be3e1e777a04d24073cfb8`) and `open_gamma4_knotinfo_2026.10.5.json`
(SHA-256 `c4816943218e4dd0bbb25b278009c7589087088e5887f5ba8fd516173b543b80`).

[C] cites KnotInfo version 2026.11.2: the package `database_knotinfo` 2026.11.2 of 2 November 2026, whose file
`csv_data/knotinfo_data_complete.csv` has SHA-256 `ed4447ff81dc50e2561764760d631356186c34e2dd35f89cc84e0a4d8fedd3d0`.
On the 12,965 prime knots with at most 13 crossings it agrees with the data file of version 2026.10.5 in every column
but the smooth concordance order of `11n_34`.  It agrees with the package 2026.8.1 in every column but these: the
values of γ₄, which are taken from the data file; 328 unknotting numbers, none of them 1 in either version; the
double slice genus of 84 knots; the Kauffman polynomial of `10_125`; and the smooth concordance order of `11n_34`.
`knotinfo_compare.py` makes this comparison:

```
python -I scripts/knotinfo_compare.py path/to/2026.10.5/knotinfo_data_complete.csv path/to/2026.8.1/knotinfo_data_complete.csv path/to/2026.11.2/knotinfo_data_complete.csv
```

## Conventions

- Knots are named as in KnotInfo: `8_10`, `11n_133`, `13a_4578`.
- KnotInfo records the least first Betti number of a nonorientable surface in B⁴, which is 1 for slice knots.  The
  [C] sets γ₄ = 0 for slice knots ([C] §3.1).
- λ_K is the linking form presented by V + Vᵀ for a Seifert matrix V.  It is minus the linking form of the double
  branched cover ([T] §2.3).
- The rational tangle 1/n of spherogram is the mirror image of a twist region of |n| crossings with incidence number η =
  sign n, as in [T] Figure 2.  So `small_examples.py` builds the pretzel knot P(3,−3,1) of [T] as the numerator closure
  of spherogram's tangles −1/3, 1/3 and −1.

## Reproducibility

- The linking-form and lattice computations use exact integer and rational arithmetic.
- The random Reidemeister moves of spherogram cannot be seeded.  So every diagram used in the searches for [C] Theorem
  4.1 and for the punctured Klein bottles is stored in `smooth_top_gap.json` and `klein_bottle_search.json`, and the
  mode `replay` of `smooth_top_gap.py` and of `klein_bottle_search.py` scans the stored diagrams again; on the final
  output the replay of `smooth_top_gap.py` reproduced the results for all 81 candidates of [C] §4.  The other searches
  store their diagrams as well (see Algorithms and checks).
- `open_sets.py` runs `paper_counts.py` only up to the set of knots that each search needs (see Checking [C] against the
  data).
- `montesinos_dinv_check.py` seeds the randomised triangulations of SnapPy and Regina, so it reproduces its output
  file exactly.
