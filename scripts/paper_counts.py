"""Recompute every count quoted in the computational paper [C] from the data files, and check each against the value
quoted there.

[C] is "Smooth and topological nonorientable four-genera of low-crossing knots", and [T] the theory paper "Linking-form
splittings and four-dimensional genera of knots" (arXiv:2610.40183), whose results [C] uses; comments name the result of
[C] that a count belongs to.  The counts of the README section Algorithms and checks are recomputed too.
Usage: python scripts/paper_counts.py   (from the repository root; prints the counts)

Each count is computed from data/*.json and data/*.txt.  An assertion records the value quoted in [C], so a
change in the data stops the script instead of silently changing a number.  The names are those of the LaTeX macros
of [C].
"""
import gzip
import json
import os
import re
from collections import Counter

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(name):
    with open(os.path.join(ROOT, "data", name)) as f:
        return json.load(f)


numbers = {}

# Knots with at most 13 crossings (KnotInfo 2026.10.5) and the signature lemma.
lg = load("linking_genus.json")
numbers["nKnotsThirteen"] = len(lg)
assert len(lg) == 12965
assert all(v["sigma_check"] for v in lg.values())                       # sigma = -sum of excesses (mod 8)
assert all((v["det"] - (-1) ** (v["sigma"] // 2)) % 4 == 0 for v in lg.values())  # D = (-1)^(sigma/2) (mod 4)

# Distribution of h_lambda (field h_witt); h^oplus (field h_summand) agrees, as no determinant reaches 3^9.
# [C] Proposition 3.1(a) quotes only h_lambda <= 3 and the number of knots with h_lambda = 3.
hdist = {h: sum(1 for v in lg.values() if v["h_witt"] == h) for h in range(4)}
assert sum(hdist.values()) == len(lg)                                    # h_lambda <= 3 for every knot
numbers["nHthree"] = hdist[3]
assert (hdist[0], hdist[1], hdist[2], hdist[3]) == (725, 7275, 4916, 49)
assert all(v["h_summand"] == v["h_witt"] for v in lg.values())
numbers["nMaxDetThirteen"] = max(v["det"] for v in lg.values())
assert numbers["nMaxDetThirteen"] == 663 < 3 ** 9
assert lg["8_10"]["h_witt"] == lg["8_10"]["h_summand"] == 1

# [C] Section 3: the census of the knots with at most 13 crossings whose value is open in KnotInfo 2026.10.5.
op = load("open_gamma4_knotinfo_2026.10.5.json")
rng = {k: v["smooth_4d_crosscap_number"].replace(" ", "") for k, v in op.items()}
count = lambda keys, r: sum(1 for k in keys if rng.get(k) == r)
numbers["nOpen"] = len(op)
for key, r in [("nOpenOneTwo", "[1,2]"), ("nOpenOneThree", "[1,3]"), ("nOpenTwoThree", "[2,3]"), ("nOpenTwoFour", "[2,4]"),
               ("nOpenOneFour", "[1,4]"), ("nOpenOneFive", "[1,5]")]:
    numbers[key] = count(op, r)
assert (numbers["nOpen"], numbers["nOpenOneTwo"], numbers["nOpenOneThree"], numbers["nOpenTwoThree"]) == (3525, 991, 773, 1751)

r1 = load("rank1_targets.json")                           # all open knots whose range starts at 1
ob1 = {k for k, v in r1.items() if v["obstructed"]}
assert len(r1) == 1768 and not any(v["obstructed_signed"] and not v["obstructed"] for v in r1.values())
numbers["nRankOneTargets"], numbers["nRankOneObstr"] = len(r1), len(ob1)
numbers["nRankOneObstrOneTwo"], numbers["nRankOneObstrOneThree"] = count(ob1, "[1,2]"), count(ob1, "[1,3]")
assert (len(ob1), numbers["nRankOneObstrOneTwo"], numbers["nRankOneObstrOneThree"]) == (251, 116, 135)

da = load("rank1_dinv_targets.json")                      # alternating, not obstructed by the linking form
oba = {k for k, v in da.items() if v["obstructed"]}
dm = load("rank1_dinv_mtargets.json")                     # Montesinos, not alternating, same condition; plus 11n_38
obm = {k for k, v in dm.items() if v["obstructed"]}
assert "11n_38" in obm and "11n_38" not in rng
numbers["nDAltTargets"], numbers["nDAltObstr"] = len(da), len(oba)
numbers["nDAltObstrOneTwo"], numbers["nDAltObstrOneThree"] = count(oba, "[1,2]"), count(oba, "[1,3]")
numbers["nDMontTargets"] = sum(1 for k in dm if k in rng)
numbers["nDMontObstr"] = sum(1 for k in obm if k in rng)
numbers["nDMontObstrOneTwo"], numbers["nDMontObstrOneThree"] = count(obm, "[1,2]"), count(obm, "[1,3]")
assert (len(da), len(oba), numbers["nDMontTargets"], numbers["nDMontObstr"]) == (887, 219, 20, 14)
nonalt = {k for k in r1 if not r1[k]["obstructed"] and "n" in k.split("_")[0]}
assert set(da) == {k for k in r1 if not r1[k]["obstructed"] and "a" in k.split("_")[0]}
assert set(dm) - {"11n_38"} <= nonalt
numbers["nNonAltPass"], numbers["nNonAltNoTest"] = len(nonalt), len(nonalt) - numbers["nDMontTargets"]
assert (len(nonalt), numbers["nNonAltNoTest"]) == (630, 610)
new2 = {k for k in ob1 | oba | obm if rng.get(k) == "[1,2]"}
new23 = {k for k in ob1 | oba | obm if rng.get(k) == "[1,3]"}
numbers["nNewTwo"], numbers["nNewTwoThree"] = len(new2), len(new23)
numbers["nNewAll"] = len(new2) + len(new23)
assert (len(new2), len(new23)) == (197, 287) and not (ob1 & oba) and not (ob1 & obm)

r2 = load("rank2_targets.json")                           # lower bound 2 (after the above) and upper bound at least 3
ob2 = {k for k, v in r2.items() if v["obstructed"]}
assert all(v["obstructed"] == v["obstructed_linking_only"] == v["obstructed_odd"] for v in r2.values())
numbers["nRankTwoTargets"], numbers["nRankTwoObstr"] = len(r2), len(ob2)
numbers["nThreeFromTwoThree"], numbers["nThreeFromOneThree"] = count(ob2, "[2,3]"), count(ob2, "[1,3]")
numbers["nNewThree"] = numbers["nThreeFromTwoThree"] + numbers["nThreeFromOneThree"]
assert (len(r2), len(ob2), numbers["nNewThree"], count(ob2, "[2,4]")) == (2044, 46, 45, 1)
h3 = {k for k, v in lg.items() if v["h_witt"] == 3}
assert h3 == ob2 | {"8_18", "9_49", "11n_133"}
numbers["nDObstrOneTwo"] = numbers["nDAltObstrOneTwo"] + numbers["nDMontObstrOneTwo"]
numbers["nDObstrOneThree"] = numbers["nDAltObstrOneThree"] + numbers["nDMontObstrOneThree"]
numbers["nTwoThreeFromOneThree"] = len(new23) - numbers["nThreeFromOneThree"]
numbers["nOpenOneFourFive"] = numbers["nOpenOneFour"] + numbers["nOpenOneFive"]
assert numbers["nOpenOneFive"] == 1                     # [C] Section 3.1 says "one knot"
assert (numbers["nDObstrOneTwo"], numbers["nDObstrOneThree"], numbers["nTwoThreeFromOneThree"]) == (81, 152, 283)
# The census table: the columns h_lambda >= 2 and [C] Proposition 3.1(c) are disjoint and add up to the knots with gamma_4 >= 2;
# the knots of range [1,3] with h_lambda = 3 are among those with h_lambda >= 2.
assert numbers["nRankOneObstrOneTwo"] + numbers["nDObstrOneTwo"] == numbers["nNewTwo"]
assert numbers["nRankOneObstrOneThree"] + numbers["nDObstrOneThree"] == numbers["nNewTwoThree"]

# [C] Proposition 3.3(b), (c): surfaces from band moves (scripts/klein_bottle_search.py).  Each of the 283 knots of range
# [1,3] with h_lambda <= 2 bounds a punctured Klein bottle K -> J -> unknot or K -> J -> J' -> two-component unlink,
# in the KnotInfo diagram; 13a_2890 bounds a surface with b_1 = 3 (K -> 7_5 -> J' -> unknot).  The 49 knots with
# h_lambda = 3 are controls: no band to a knot with KnotInfo value 1.
kb = load("klein_bottle_search.json")
ks, kk = kb["summary"], kb["knots"]
kt = {k for k, v in kk.items() if v["proposition_5_1b"] and v["h_witt"] <= 2}
assert kt == {k for k in new23 if lg[k]["h_witt"] <= 2} and len(kt) == numbers["nTwoThreeFromOneThree"]
assert ks["targets"] == ks["certified"] == ks["self_contained"] == len(kt) and not ks["not_certified"]
assert ks["contradictions"] == 0 and ks["certificates_failing_a_recheck"] == 0
assert all(kk[k]["representative"]["diagram_source"] == "knotinfo" for k in kt)
mob = Counter(kk[k]["representative"]["moebius_band_of_J"]["kind"] for k in kt)
numbers["nKleinCertified"] = len(kt)
numbers["nKleinUnknotBand"], numbers["nKleinRibbonBand"] = mob["band to the unknot"], mob["band to a ribbon knot"]
numbers["nKleinTZero"] = sum(kk[k]["representative"]["t"] == 0 for k in kt)
assert (len(kt), numbers["nKleinUnknotBand"], numbers["nKleinRibbonBand"], numbers["nKleinTZero"]) == (283, 194, 89, 277)
assert {abs(kk[k]["representative"]["t"]) for k in kt} == {0, 1}
ctrl = ks["controls"]
assert set(ctrl) == h3 and all(v["contradictions"] == 0 for v in ctrl.values())
assert all(v["hits"] == 0 for k, v in ctrl.items() if k != "13a_2890")
numbers["nKleinControlDiagrams"] = sum(v["diagrams"] for v in ctrl.values())
assert numbers["nKleinControlDiagrams"] == 9312
k2890 = ks["knot_13a_2890"]["representative"]
assert ks["knot_13a_2890"]["certified"] and k2890["self_contained"] and k2890["J"] == "7_5" and k2890["t"] == 1
assert count(ob2, "[2,4]") == 1 and "13a_2890" in ob2

# [C] Proposition 3.3(b) and [C] Section 3.8, mode ranges (data/klein_bottle_search_ranges.json).  Every knot of range
# [2,3] or [2,4] with h_lambda = 2 bounds a punctured Klein bottle built from band moves, so gamma_4 = gamma_top = 2;
# 11n_133 bounds a surface with b_1 = 3 (K -> 4_1 -> J' -> unknot); the 4 knots of range [1,4] or [1,5] bound Mobius
# bands and are not slice (signature), so gamma_4 = gamma_top = 1.  The rebuild (data/klein_bottle_rebuild.json)
# links every certificate to the KnotInfo PD code by band moves and Reidemeister moves, without SnapPy.
kr = load("klein_bottle_search_ranges.json")
rs, rk = kr["summary"], kr["knots"]
rt = {k for k, v in rk.items() if v["group"] == "target"}
upper = {k for k, v in rk.items() if v["group"] == "upper"}
assert rt == {k for k, r in rng.items() if r in ("[2,3]", "[2,4]") and lg[k]["h_witt"] == 2}
assert {k for k, r in rng.items() if r in ("[2,3]", "[2,4]")} == rt | {k for k in h3 if rng.get(k) in ("[2,3]", "[2,4]")}
assert upper == {k for k, r in rng.items() if r in ("[1,4]", "[1,5]")} and all(lg[k]["sigma"] != 0 for k in upper)
assert rs["targets"] == rs["certified"] == len(rt) and not rs["not_certified"] and rs["not_self_contained"] == ["13a_3126"]
assert rs["contradictions"] == 0 and rs["certificates_failing_a_recheck"] == 0 and rs["diagrams_not_isometric"] == 0
assert all(rs["upper"][k]["certified"] and rs["upper"][k]["upper_bound"] == 1 for k in upper)
numbers["nRangesTargets"] = len(rt)
numbers["nRangesTwoThree"], numbers["nRangesTwoFour"] = count(rt, "[2,3]"), count(rt, "[2,4]")
numbers["nRangesKI"], numbers["nRangesStageB"] = rs["representative_from_knotinfo_diagram"], len(rs["stage_B"])
assert numbers["nRangesKI"] + numbers["nRangesStageB"] == len(rt) and "13a_3126" in rs["stage_B"]
numbers["nRangesTZero"] = rs["representative_t"]["0"]
assert set(rs["representative_t"]) == {"0", "1"}
numbers["nRangesUnknotBand"] = rs["representative_surface_of_J"]["band to the unknot"]
numbers["nRangesRibbonBand"] = rs["representative_surface_of_J"]["band to a ribbon knot"]
assert (len(rt), numbers["nRangesTwoThree"], numbers["nRangesKI"], numbers["nRangesTZero"], numbers["nRangesUnknotBand"],
        numbers["nRangesRibbonBand"]) == (1715, 1710, 1698, 1675, 965, 750)
k133 = rs["knot_11n_133"]
assert k133["certified"] and k133["self_contained"] and k133["b1"] == 3 and k133["representative"]["J"] == "4_1"
assert k133["representative"]["diagram_source"] == "knotinfo" and "11n_133" in h3
assert rk["11n_133"]["contradictions"] == 0 and "1" not in rk["11n_133"]["J_gamma4_seen"]
numbers["nElevenControlDiagrams"] = len(rk["11n_133"]["diagrams"])
rb = load("klein_bottle_rebuild.json")["summary"]
b1, b0 = rb["klein_bottle_search_ranges.json"], rb["klein_bottle_search.json"]
assert b1["targets_representative_diagrammatic"] == len(rt) and b1["targets_representative_b1"] == {"2": len(rt)}
assert b0["targets_representative_diagrammatic"] == len(kt) and b0["targets_representative_b1"] == {"2": len(kt)}
assert not b1["certificates_failing"] and not b0["certificates_failing"]
assert b1["certificates_rebuilt"] == b1["certificates"] and b0["certificates_rebuilt"] == b0["certificates"]
assert b1["others"]["11n_133"]["representative_ok"] and b1["others"]["11n_133"]["b1"] == 3
assert b0["others"]["13a_2890"]["representative_ok"] and b0["others"]["13a_2890"]["b1"] == 3
assert all(b1["others"][k]["moebius_band_of_K_ok"] and b1["others"][k]["moebius_band_of_K_b1"] == 1 for k in upper)
# [C] Proposition 3.3(a), (b), the first search, mode mobius (data/klein_bottle_search_mobius.json): Moebius bands for the knots of range
# [1,2] or [1,3] that neither test obstructs, by our own bands (K -> J' -> two-component unlink, or K -> unknot), first in
# the KnotInfo diagram and then in diagrams of stage B; the knots of range [1,3] without one bound punctured Klein bottles,
# so their range becomes [1,2].  The knots that a test obstructs are controls and give no Moebius band.  The rebuild
# (data/klein_bottle_rebuild_mobius.json) redoes every band move.  A Moebius band gives gamma_4 <= 1, and gamma_4 = 1
# once K is not smoothly slice: for h_lambda = 1 (not metabolic) or a nonzero signature K is not even topologically
# slice, so gamma_top = 1 too; for the other knots KnotInfo lists the smooth 4-genus 1, and those whose Alexander
# polynomial is not of the form f(t) f(1/t) are not topologically slice either (Fox-Milnor, scripts/fox_milnor.py).
km = load("klein_bottle_search_mobius.json")
ms, mk = km["summary"], km["knots"]
mt = {k for k, v in mk.items() if v["group"] == "mobius"}
mc = {k for k, v in mk.items() if v["group"] == "control"}
assert mt == {k for k, r in rng.items() if r in ("[1,2]", "[1,3]")} - new2 - new23
assert mc == (new2 | new23) and ms["contradictions"] == 0 and ms["fit_failures"] == 0
assert ms["diagrams_not_isometric"] == 0 and ms["moebius_band_far"] == []
mb = {k for k in mt if mk[k]["moebius_band_of_K"]}
numbers["nMobiusOneTwo"], numbers["nMobiusOneThree"] = count(mb, "[1,2]"), count(mb, "[1,3]")
kinds_mobius_with_surface_of_J = sum(bool(c.get("moebius_band_of_J") or c.get("klein_bottle_of_J"))
                                     for v in mk.values() for c in v.get("certificates") or [])
numbers["nMobiusAll"] = len(mb)
assert ms["moebius_band_of_K"] == {"[1,2]": [numbers["nMobiusOneTwo"], count(mt, "[1,2]")],
                                   "[1,3]": [numbers["nMobiusOneThree"], count(mt, "[1,3]")]}
numbers["nMobiusKI"] = ms["moebius_band_diagram_source"]["knotinfo"]
numbers["nMobiusStageB"] = len(ms["moebius_band_stage_B"])
numbers["nMobiusUnknot"] = ms["moebius_band_kind"]["band to the unknot"]
numbers["nMobiusRibbon"] = ms["moebius_band_kind"]["band to a ribbon knot"]
assert numbers["nMobiusKI"] + numbers["nMobiusStageB"] == len(mb) == numbers["nMobiusUnknot"] + numbers["nMobiusRibbon"]
mkb = set(ms["range_one_three_to_one_two"])
assert mkb == {k for k in mt if rng[k] == "[1,3]"} - mb
mkbs = ms["klein_bottle"]
assert mkbs["knots"] == mkbs["certified"] == len(mkb) and mkbs["self_contained"] == [len(mkb), len(mkb)]
assert mkbs["contradictions"] == 0 and mkbs["certificates_failing_a_recheck"] == 0 and mkbs["diagrams_not_isometric"] == 0
numbers["nMobiusKlein"] = len(mkb)
numbers["nMobiusKleinKI"] = mkbs["representative_from_knotinfo_diagram"]
numbers["nMobiusControls"] = sum(v["knots"] for v in ms["controls"].values())
numbers["nMobiusControlDiagrams"] = sum(v["diagrams"] for v in ms["controls"].values())
assert numbers["nMobiusControls"] == len(mc) and all(v["moebius_bands"] == 0 for v in ms["controls"].values())
rbm = load("klein_bottle_rebuild_mobius.json")["summary"]["klein_bottle_search_mobius.json"]
assert rbm["mobius"]["moebius_bands_of_K_rebuilt"] == {r: [count(mb, r)] * 2 for r in ("[1,2]", "[1,3]")}
assert rbm["mobius"]["moebius_bands_of_K_b1"] == {"1": len(mb)} and rbm["mobius"]["moebius_bands_of_K_diagrammatic"] == len(mb)
assert not rbm["mobius"]["moebius_bands_of_K_failing"] and not rbm["certificates_failing"]
assert rbm["mobius"]["klein_bottles"]["rebuilt"] == rbm["mobius"]["klein_bottles"]["diagrammatic"] == len(mkb)
assert rbm["mobius"]["klein_bottles"]["b1"] == {"2": len(mkb)} and rbm["certificates_rebuilt"] == rbm["certificates"]
# not slice: h_lambda = 1, else a nonzero signature, else the smooth 4-genus 1 of KnotInfo
mh1 = {k for k in mb if lg[k]["h_witt"] == 1}
msig = {k for k in mb - mh1 if lg[k]["sigma"] != 0}
mki = mb - mh1 - msig
assert all(lg[k]["h_witt"] in (0, 1) for k in mb) and all(op[k]["smooth_four_genus"].strip() == "1" for k in mki)
numbers["nMobiusHOne"], numbers["nMobiusSig"], numbers["nMobiusKIgfour"] = len(mh1), len(msig), len(mki)
fm = load("fox_milnor.json")
mfm = set(fm["summary"]["fox_milnor_fails"])
assert set(fm["knots"]) == mki and fm["summary"]["agree_with_knotinfo"] == len(mki)
assert mfm == {k for k, v in fm["knots"].items() if not v["fox_milnor"]}
assert all(f["reciprocal"] and f["multiplicity"] == 1 for k in mfm for f in fm["knots"][k]["factors"])
assert mfm == {"12n_837", "13a_1431", "13a_1506", "13n_1967", "13n_2985"}       # named in [C]
numbers["nMobiusFoxMilnor"], numbers["nMobiusFoxMilnorHolds"] = len(mfm), len(mki - mfm)
assert (len(mb), numbers["nMobiusOneTwo"], numbers["nMobiusOneThree"], len(mkb), len(mh1), len(msig), len(mki)) == \
    (905, 541, 364, 122, 879, 5, 21)
# README, Algorithms and checks: the replay of the search of mode mobius on its stored diagrams
# (data/klein_bottle_search_mobius_replay.txt, the printed output of `klein_bottle_search.py mobius replay`) reproduced
# every stored result.
with open(os.path.join(ROOT, "data", "klein_bottle_search_mobius_replay.txt")) as f:
    rp = [line for line in f.read().splitlines() if not line.startswith("#")]
tail = rp[-7:]
assert tail[0] == f"replay: knots {len(mt) + len(mc)} Mobius band search differing from the stored output: []"
assert tail[1] == f"punctured Klein bottle searches {len(mkb)} differing from the stored output: []"
assert tail[2] == "knots whose number of hits identified by the simplified diagram changed: [['13n_1380', 3, 2], ['13n_1385', 4, 3]]"
assert tail[3] == f"certificates rebuilt: {rbm['certificates']} with a different J or isometry check: []"
assert tail[4] == "surfaces of J rebuilt: 363 failing: []" and kinds_mobius_with_surface_of_J == 363
assert tail[5] == f"Mobius bands of K rebuilt: {len(mb)} failing: []"
assert tail[6] == f"moebius_band in the KnotInfo diagrams: {len(mt) + len(mc)} bands {numbers['nMobiusKI']} differing from moebius_band_fast: []"
assert not any("Traceback" in line or "Error" in line for line in rp)
numbers["nMobiusSearched"] = len(mt) + len(mc)
assert numbers["nMobiusSearched"] == 1764
numbers["nKleinRebuilt"] = (b1["certificates_rebuilt"] + b0["certificates_rebuilt"] + rbm["certificates_rebuilt"]
                            + len(mb))
# README, Algorithms and checks: of the rebuilt band records, the complete surfaces (the Moebius bands of K and the
# records whose surface of J lies in the diagram of J that the first band produces), and the others, which store a band
# move to J with a surface of J in another diagram of J (among them the representative of 13a_3126) or with no surface
# of J.
kinds = Counter()
for f in ("klein_bottle_search.json", "klein_bottle_search_ranges.json", "klein_bottle_search_mobius.json"):
    for v in load(f)["knots"].values():
        for c in v.get("certificates") or []:
            sj = c.get("moebius_band_of_J") or c.get("klein_bottle_of_J")
            kinds["none" if not sj else "complete" if sj["diagram"] == "J made by the band" else "other"] += 1
assert kinds["complete"] == (b1["self_contained_certificates_rebuilt"] + b0["self_contained_certificates_rebuilt"]
                             + rbm["self_contained_certificates_rebuilt"])
assert list(b1["representatives_not_self_contained"]) == ["13a_3126"] and not b0["representatives_not_self_contained"]
numbers["nKleinComplete"] = kinds["complete"] + len(mb)
numbers["nKleinPartial"] = kinds["other"] + kinds["none"]
numbers["nKleinJOther"], numbers["nKleinJNone"] = kinds["other"], kinds["none"]
assert numbers["nKleinComplete"] + numbers["nKleinPartial"] == numbers["nKleinRebuilt"]
assert (numbers["nKleinComplete"], numbers["nKleinJOther"], numbers["nKleinJNone"]) == (7038, 93, 105)
numbers["nKleinAll"] = len(kt) + len(rt)
numbers["nGammaOne"] = len(upper) + len(mb)
numbers["nGammaTopOne"] = len(upper) + len(mh1) + len(msig) + len(mfm)
numbers["nGammaTopTwo"] = sum(1 for k in new2 | kt | rt if lg[k]["h_witt"] == 2)
numbers["nGammaThreeAll"] = numbers["nNewThree"] + 1
# [C] Propositions 3.1 and 3.3 (first search): the values they determine and the knots it leaves with gamma_4 in {1, 2}.
first_gamma_two = len(new2) + len(kt) + len(rt)
first_determined = numbers["nGammaOne"] + first_gamma_two + numbers["nGammaThreeAll"]
first_open_one_two = numbers["nOpenOneTwo"] - len(new2) - numbers["nMobiusOneTwo"]
numbers["nFirstOpenFromOneThree"] = len(mkb)
first_open = first_open_one_two + numbers["nFirstOpenFromOneThree"]
# no knot keeps the range [1,3] ([C] Theorem C: every other range [1,3] becomes [1,2])
numbers["nOpenLeftOneThree"] = numbers["nOpenOneThree"] - len(new23) - numbers["nMobiusOneThree"] - len(mkb)
assert numbers["nOpenLeftOneThree"] == 0
numbers["nBandsOneThree"] = len(kt) + numbers["nMobiusOneThree"] + len(mkb)
assert first_determined + first_open == numbers["nOpen"]
assert numbers["nOpenOneThree"] == len(new23) + numbers["nMobiusOneThree"] + len(mkb)
assert (numbers["nGammaOne"], numbers["nGammaTopOne"], first_gamma_two, numbers["nGammaThreeAll"], first_determined,
        first_open) == (909, 893, 2195, 46, 3150, 375)
left0 = {k for k in mt if not mk[k]["moebius_band_of_K"]}
assert len(left0) == first_open and mkb <= left0

# [C] Proposition 3.2(a), (b): correction terms for the non-alternating knots that are not Montesinos among the
# knots that Proposition 3.1 and the first search leave open, from negative-definite plumbing trees (scripts/plumbing_dinv.py, with the
# normal forms of the follow-up paper arXiv:2610.25430 and of the same pipeline) and, for unknotting number one, from
# surgery descriptions (scripts/surgery_dinv.py, and scripts/surgery_box.py for every V allowed by published bounds).
alt = lambda k: op[k]["alternating"] == "Y"
mont = lambda k: op[k]["montesinos_notation"].startswith("K(")
nd = {k for k in left0 if not alt(k) and not mont(k)}
first_open_alt = sum(1 for k in left0 if alt(k))
assert {k for k in left0 if mont(k) and not alt(k)} == {"12n_219"}
assert (len(nd), first_open_alt, len(left0 - nd) - first_open_alt) == (316, 58, 1)
pdv = load("plumbing_dinv.json")
ps = pdv["summary"]
assert set(pdv["targets"]) == nd and ps["targets"] == len(nd) and not ps["undecided"] and not ps["problems"]
cls = Counter(v["cls"] for v in pdv["targets"].values())
assert cls == Counter(ps["classes"]["all"])
numbers["nNDTargets"], numbers["nNDFollowup"], numbers["nNDNew"] = len(nd), ps["followup_records"], ps["new"]
(numbers["nNDHyp"], numbers["nNDMixed"], numbers["nNDNoDef"], numbers["nNDNoDefAlg"],
 numbers["nNDNegDef"]) = (cls[c] for c in ("hyperbolic", "mixed", "nodef", "nodef_alg", "negdef"))
numbers["nNDNoDefAll"] = numbers["nNDNoDef"] + numbers["nNDNoDefAlg"]
negd = {k for k, v in pdv["targets"].items() if v["cls"] == "negdef"}
assert negd == set(ps["negdef"]["knots"]) and not ps["negdef"]["not_eligible"]
assert all(pdv["targets"][k]["dtest"]["bad"] <= 2 and pdv["targets"][k]["dtest"]["det_is_D"]
           and pdv["targets"][k]["dtest"]["negative_definite"] for k in negd)
numbers["nNDBadTwo"] = ps["negdef"]["bad_vertices"]["2"]
pob = set(ps["obstructed"])
assert pob == {k for k in negd if pdv["targets"][k]["verdict"] == "obstructed"}
assert set(ps["passes"]) == negd - pob == {"13n_2413", "13n_3505"}
numbers["nNDPlumbObstr"], numbers["nNDPlumbPass"] = len(pob), len(negd - pob)
assert (numbers["nNDHyp"], numbers["nNDMixed"], numbers["nNDNoDef"], numbers["nNDNoDefAlg"], numbers["nNDNegDef"],
        numbers["nNDBadTwo"], len(pob)) == (176, 30, 47, 4, 59, 9, 57)
pck = ps["checks"]
assert pck["independent_dtest"] == {"agrees": 57, "H_1 not cyclic": 2}
assert pck["star_plumbing_route"]["agrees"] == 38 and pck["knotinfo_algebraic_agrees"] == "51 of 51"
numbers["nNDIndependent"], numbers["nNDStarRoute"] = pck["independent_dtest"]["agrees"], pck["star_plumbing_route"]["agrees"]
pcs = ps["controls"]
ff, npl, kfl = pcs["followup_forms"], pcs["new_pipeline"], pcs["known_fillings"]
assert ff["forms"] == 180 and ff["gamma4_one"] == 114 and not ff["gamma4_one_obstructed"] and not ff["independent_disagrees"]
assert npl["knots"] == 2874 and not npl["errors"] and not npl["obstructed"] and npl["dtests"] == 569
assert npl["regina_trees_identified"] == npl["regina_trees"] and not npl["routes_disagree"] and not npl["independent_disagrees"]
assert not npl["casson_walker_d_fails"] and not npl["tree_checks_fail"] and not npl["not_eligible_bad"]
assert kfl["knots"] == 1619 and kfl["knots_with_dtest"] == 1424 and not kfl["disagree"] and not kfl["errors"]
numbers["nNDControlsFollowup"], numbers["nNDControlsOne"] = ff["forms"], ff["gamma4_one"]
numbers["nNDControlsPipeline"], numbers["nNDControlsTrees"] = npl["knots"], npl["dtests"]
numbers["nNDControlsTreeKnots"] = sum(1 for v in pdv["controls"]["new_pipeline"].values() if v.get("trees"))
numbers["nNDControlsKnown"], numbers["nNDControlsKnownTested"] = kfl["knots"], kfl["knots_with_dtest"]
assert numbers["nNDControlsTreeKnots"] == 335
sdv = load("surgery_dinv.json")
sbx = load("surgery_box.json")
st, sb = sdv["summary"]["targets"], sbx["summary"]
su = {k for k in nd if op[k]["unknotting_number"].strip() == "1"}
assert set(sdv["targets"]) == su and st["knots"] == len(su) == 50 and not any(sbx["summary"]["checks"].values())
sob = set(sb["box_obstructed"])
assert sob <= su and set(sb["exact_obstructed"]) == set(st["obstructed"]) and (
    not sb["box_not_exact"] or set(sb["box_not_exact"]) == set(st["obstructed_by_bounds"]))
assert set(sb["exact_not_box"]) <= pob                       # these four are excluded by their plumbing trees
numbers["nSurgTargets"], numbers["nSurgBox"] = len(su), len(sob)
numbers["nSurgOnly"], numbers["nSurgAndPlumb"] = len(sob - pob), len(sob & pob)
cal = sdv["summary"]["calibration"]
assert cal["knots"] == 32 and not cal["disagree"] and list(cal["not_completed"]) == ["13a_1085"]
assert not cal["obstructed_with_knotinfo_value_1"] and not sb["checks"]["calibration_gamma4_one_box_obstructed"]
g4cal = load("gamma4_knotinfo_2026.10.5.json")
assert not [k for k in sb["calibration"]["box_obstructed"] if g4cal[k] == "1"]     # "none with KnotInfo value 1"
assert Counter(g4cal[k] for k in sdv["calibration"]) == Counter({"1": 14, "2": 18})
numbers["nSurgCalib"], numbers["nSurgCalibAgree"] = cal["knots"], len(cal["agree_with_goeritz"])
numbers["nSurgCalibBox"] = len(sb["calibration"]["box_obstructed"])
assert set(sb["calibration"]["box_obstructed"]) == {k for k, v in sdv["calibration"].items() if v["verdict"] == "obstructed"}
assert sb["pilot"] == {"11n_176": False} and sdv["summary"]["pilot"]["11n_176"]["passing"] == {"sign": 1, "l": 7, "M": 3}
# Exact V_j(C) for the other knots with unknotting number one (scripts/surgery_cfk.py): the UV = 0 knot Floer complex of C
# (Szabo's program) determines V_j(C), by the gradings of H_*(A^_s) and over every lift to F_2[U,V]; the chirality of C
# is fixed by an isometry, Regina's names, the Khovanov rank (L-space) or finite covers.  Controls: calibration,
# the knots already excluded, 11n_176 and connected sums of torus knots.
cfk = load("surgery_cfk.json")
cs = cfk["summary"]
cob = set(cs["excluded"])
assert cs["targets"] == len(su - pob - sob) == 33 and set(cs["excluded"]) | set(cs["passes"]) == su - pob - sob
assert not cs["undetermined"] and not any(cs["problems"].values()) and set(cs["excluded_by_box_with_chirality"]) == cob
assert cob == {"13n_3844", "13n_3855"}                          # named in [C] Section 3.2
ccal, coth, cten = cs["controls"]["calibration"], cs["controls"]["others"], cs["controls"]["tensor"]
assert ccal["knots"] == 32 and ccal["with_route"] == ccal["goeritz_agrees"] == numbers["nSurgCalibAgree"]
assert set(ccal["excluded"]) == set(sb["calibration"]["box_obstructed"])
assert coth["knots"] == coth["excluded"] == len(su & (pob | sob)) == 17
assert cs["controls"]["pilot"]["11n_176"]["verdict"] == "passes" and \
    cs["controls"]["pilot"]["11n_176"]["V"] == sdv["summary"]["pilot"]["11n_176"]["V"]
assert all(v["agrees"] and v["same_hat"] for v in cten.values())
numbers["nSurgExactTargets"], numbers["nSurgExact"], numbers["nSurgTensor"] = cs["targets"], len(cob), len(cten)
numbers["nSurgExactPass"] = len(cs["passes"])
numbers["nSurgExactControls"] = coth["knots"]
# a surgery description excludes the knots excluded by the bounds, by the exact V_j(C), and, again, the knots with
# unknotting number one that (a) excludes (the controls of surgery_cfk.py)
surg_excl = sob | cob | (su & pob)
numbers["nSurgExcl"], numbers["nSurgExclOnly"] = len(surg_excl), len(surg_excl - pob)
assert (numbers["nSurgExcl"], numbers["nSurgExclOnly"]) == (19, 13)
newd = pob | sob | cob
numbers["nNewD"] = len(newd)
numbers["nNDOneTwo"], numbers["nNDOneThree"] = count(newd, "[1,2]"), count(newd, "[1,3]")
assert newd <= nd and newd & mkb == {k for k in newd if rng[k] == "[1,3]"}
assert (len(sob), numbers["nSurgOnly"], numbers["nSurgAndPlumb"], len(newd), numbers["nNDOneTwo"],
        numbers["nNDOneThree"], numbers["nSurgCalibAgree"], numbers["nSurgCalibBox"]) == (13, 11, 2, 70, 40, 30, 31, 17)
assert all(lg[k]["h_witt"] == 1 for k in newd)                 # so gamma_top >= 1 for each of them ([C] Theorem 2.1)
if os.environ.get("PAPER_COUNTS_UPTO") == "newd":                 # scripts/open_sets.py: mobius_open_search.py
    raise SystemExit("stopped after newd")
# the knots of Proposition 3.2(a), (b) that stay open: neither method applies, or the test passes, or the bounds do not suffice
ndleft = nd - newd
numbers["nNDLeft"] = len(ndleft)
numbers["nNDLeftNoMethod"] = sum(1 for k in ndleft if pdv["targets"][k]["cls"] != "negdef" and k not in su)
numbers["nSurgLeft"] = len(ndleft & su)
assert (numbers["nNDLeft"], numbers["nNDLeftNoMethod"], numbers["nSurgLeft"]) == (246, 213, 31)
assert ndleft & su == set(cs["passes"])                       # the exact correction terms pass the test for these
assert numbers["nNDLeftNoMethod"] + numbers["nSurgLeft"] + numbers["nNDPlumbPass"] == numbers["nNDLeft"]
# [C] Section 3.8, 11n_176: Sigma_K = +-S^3_{63/2}(C), the reduced Khovanov homology over Z/2 has rank D = 63, so Sigma_K
# is an L-space (Ozsvath-Szabo, branched double covers, Cor. 1.2) and C an L-space knot; V(C) = torsion coefficients.
p176 = sdv["pilot"]["11n_176"]
assert p176["khovanov"]["rank"] == p176["khovanov"]["det"] == 63 and p176["khovanov"]["sigma_is_L_space"]
assert p176["V"] == [4, 3, 3, 2, 2, 2, 1, 1, 1, 0] and p176["V_source"] == "L-space knot" and not p176["obstructed"]
r176 = [r for r in p176["routes"] if r.get("kappa_hat")][0]
assert "K6_29(0,0)" in r176["kappa"]["identify"] and r176["kappa_hat"]["hfk"]["seifert_genus"] == 9
assert r176["kappa_hat"]["hfk"]["tau"] == 9 and r176["kappa_hat"]["hfk"]["fibered"]
# [C] Section 3.3, the second search for Moebius bands, on the knots left open after Proposition 3.2(a), (b)
# (scripts/mobius_open_search.py, stopped when the knots left open had been scanned in five to seven rounds of 40
# diagrams, at most two of them of new diagrams), every band rebuilt and related to the KnotInfo diagram of
# K by Reidemeister moves (scripts/mobius_open_rebuild.py).  All twelve bands lie in diagrams of mode mobius, and the
# ribbon disc of J in another diagram of J.  Bands to knots that KnotInfo lists as slice are not used.
mos = json.load(gzip.open(os.path.join(ROOT, "data", "mobius_open_search.json.gz"), "rt"))
moss = mos["summary"]
wide = set(moss["self_contained"]) | set(moss["self_contained_up_to_reidemeister_moves_of_J"])
assert {k for k, v in mos["knots"].items() if v["group"] == "target"} == left0 - newd and moss["targets"] == 305
reb = load("mobius_open_rebuild.json")["summary"]["moebius_bands"]
assert set(reb["rebuilt"]) == set(reb["diagrammatic"]) == wide and not reb["failing"] and reb["b1"] == {"1": len(wide)}
assert not moss["errors"] and not moss["certificates_failing_a_recheck"] and not moss["controls"]["contradictions"]
assert moss["diagrams_not_isometric"] == 0 and set(moss["knotinfo_slice_only"]) == {"13a_1412"}   # not used
numbers["nWideTargets"], numbers["nWideMobius"] = moss["targets"], len(wide)
numbers["nWideMobiusAlt"] = moss["by_kind"]["alternating"]["moebius_band"]
rep = {k: mos["knots"][k]["moebius_band"] for k in wide}
assert all(mos["knots"][k]["rounds"] == 1 and c["diagram"] < len(mos["knots"][k]["diagrams"]) and
           "parent" not in mos["knots"][k]["diagrams"][c["diagram"]] for k, c in rep.items())   # stored diagrams
numbers["nWideKI"] = sum(c["diagram_source"] == "knotinfo" for c in rep.values())
numbers["nWideReached"] = sum(c["surface_of_J"]["diagram"] == "diagram of J reached from J made by the band"
                              for c in rep.values())
numbers["nWideViaX"] = len(wide) - numbers["nWideReached"]
assert {k for k, c in rep.items() if len(c["surface_of_J"]["ribbon_bands"]) != 1} == {"13n_1374"}
assert len(rep["13n_1374"]["surface_of_J"]["ribbon_bands"]) == 2
opn = [v["diagrams"] for v in mos["knots"].values() if v["group"] == "target" and not v["moebius_band"]]
numbers["nWideOpenMin"], numbers["nWideOpenMax"] = min(map(len, opn)), max(map(len, opn))
numbers["nWideOpenNew"] = sum("parent" in d for ds in opn for d in ds)
ctl = [v for v in mos["knots"].values() if v["group"] == "control"]
numbers["nWideControls"], numbers["nWideControlDiagrams"] = len(ctl), sum(len(v["diagrams"]) for v in ctl)
assert moss["controls"]["knots"] == len(ctl) and moss["controls"]["diagrams"] == numbers["nWideControlDiagrams"]
assert all(not v["certificates"] and not v["knotinfo_slice_count"] for v in ctl)
assert (numbers["nWideKI"], numbers["nWideReached"], numbers["nWideOpenMin"], numbers["nWideOpenMax"],
        numbers["nWideOpenNew"], numbers["nWideControls"], numbers["nWideControlDiagrams"]) == (6, 10, 165, 279, 5960, 32, 2560)
# [C] Proposition 3.2(c): Donaldson's theorem with definite checkerboard surfaces of non-minimal diagrams
# (scripts/donaldson_surfaces.py; every certificate checked by scripts/donaldson_check.py, which shares no code with the
# search), and surgery descriptions by drilling short closed curves (scripts/surgery_drill.py), whose obstructions are all among
# the former.
oleft3 = left0 - newd - wide                                      # open before Proposition 3.2(c)
if os.environ.get("PAPER_COUNTS_UPTO") == "oleft3":               # scripts/open_sets.py: the searches that start here
    raise SystemExit("stopped after oleft3")
dsf, dck, sdr = load("donaldson_surfaces.json"), load("donaldson_check.json"), load("surgery_drill.json")
dts, dcs = dsf["targets"]["summary"], dsf["controls"]["summary"]
dkn = {k for k in oleft3 if not alt(k) and not mont(k)}
assert set(dsf["targets"]["knots"]) == {k for k in dkn if r1[k]["s"] in (2, 6)} and dts["knots"] == len(dkn)
assert {r1[k]["s"] for k in dkn} == {0, 2, 6}
assert dsf["meta"]["targets"]["convention_check"]["ok"] and dck["meta"]["convention_check"]["ok"]
dcert = set(dts["certified"])
assert set(dck["summary"]["passes"]) == dcert and not dck["summary"]["failing"]
assert not dcs["contradictions"] and dsf["meta"]["targets"]["new"] == dsf["meta"]["controls"]["new"]
assert dcs["searched"] == len(dsf["controls"]["knots"])
drl = set(sdr["summary"]["obstructed"])
assert set(sdr["knots"]) == {k for k in dkn if k not in su and pdv["targets"][k]["cls"] == "hyperbolic"}
assert drl <= dcert and not sdr["summary"]["errors"]
numbers["nOpenLeftThree"] = len(oleft3)
numbers["nOpenThreeAlt"] = sum(1 for k in oleft3 if alt(k))
numbers["nDefKnots"], numbers["nDefSearched"] = len(dkn), len(dsf["targets"]["knots"])
numbers["nDefSZero"] = numbers["nDefKnots"] - numbers["nDefSearched"]
numbers["nDefNew"], numbers["nDefDiagrams"] = dsf["meta"]["targets"]["new"], dts["diagrams"]
numbers["nDefWithSurface"], numbers["nDefCert"] = dts["with_definite_surface"], len(dcert)
numbers["nDefOneTwo"], numbers["nDefOneThree"] = count(dcert, "[1,2]"), count(dcert, "[1,3]")
assert numbers["nDefOneTwo"] + numbers["nDefOneThree"] == len(dcert)
numbers["nDefSignDecides"] = dck["summary"]["linking_form_sign_decides"]
assert dck["summary"]["searched"] == numbers["nDefSearched"] and not dck["summary"]["knotinfo_diagram_definite"]
assert numbers["nDefSignDecides"] == len(dcert)                   # the linking form decides the sign for every certificate
cc = [r["certificate"]["crossings"] for r in dsf["targets"]["knots"].values() if r["certificate"]]
numbers["nDefCertMinC"], numbers["nDefCertMaxC"] = min(cc), max(cc)
numbers["nDefControls"], numbers["nDefControlsSurface"] = dcs["searched"], dcs["with_definite_surface"]
numbers["nDefControlLattices"] = dcs["lattices"]
numbers["nDrillTargets"] = len(sdr["knots"])
numbers["nDrillDescr"] = len(sdr["summary"]["with_surgery_description"])
assert drl == {"13n_2626", "13n_2632", "13n_3616"}
# The descriptions without V_j in surgery_drill.py: the test with the exact V_j of [C] Proposition 3.2(b)
# (scripts/drill_cfk.py); the drilling knots with a verdict, recomputed there, keep it.
dcf = load("drill_cfk.json")
dcfs = dcf["summary"]
assert set(dcf["targets"]) == {k for k, v in sdr["knots"].items() if v["verdict"] == "no d-function"}
assert not dcfs["excluded"] and not any(v for v in dcfs["problems"].values() if v)
assert all(v["verdict"] == v["expected"] for v in dcfs["controls"].values())
dcf_none = set(dcfs["undetermined"])
assert all(any(r.get("flag", "").startswith("knot Floer homology of kappa not computed") for r in sdr["knots"][k]["routes"])
           for k in dcf_none)
numbers["nDrillPass"] = len(sdr["summary"]["passes"]) + len(dcfs["passes"])
numbers["nDrillDetermined"] = len(drl) + numbers["nDrillPass"]
numbers["nDrillNoHFK"] = len(dcf_none)
assert numbers["nDrillDetermined"] + numbers["nDrillNoHFK"] == numbers["nDrillDescr"]
assert (numbers["nDrillPass"], numbers["nDrillNoHFK"]) == (16, 2)
# The second construction continued after Proposition 3.2(c) on the knots then left open (the set oleft below): a ribbon
# disc looked for for every knot J of its bands (scripts/mobius_J_search.py), every Moebius band rebuilt by
# scripts/mobius_open_rebuild.py (data/mobius_J_rebuild.json).  wide_all: all Moebius bands of the second construction.
mjs = json.load(gzip.open(os.path.join(ROOT, "data", "mobius_J_search.json.gz"), "rt"))
wj = set(mjs["summary"]["moebius_band"])
wide_all = wide | wj
assert wj == {"13a_1412", "13a_2898", "13a_3125", "13n_3901"} and not wide & wj
numbers["nWideAll"] = len(wide_all)
# [C] Table 1: of the knots of range [1,3] left with punctured Klein bottles by the first search, these get Moebius
# bands in the second construction and count under gamma_4 = 1
numbers["nWideFromOneThree"], numbers["nKleinFromOneThreeLeft"] = len(wide_all & mkb), len(mkb - wide_all)
assert (numbers["nWideFromOneThree"], numbers["nKleinFromOneThreeLeft"]) == (8, 114)
# gamma_top = 1 for those with h_lambda = 1, or a nonzero signature, or Delta failing Fox-Milnor (fox_milnor.py wide)
fmw = load("fox_milnor_wide.json")["summary"]
assert set(fmw["fox_milnor_fails"]) | set(fmw["fox_milnor_holds"]) == \
    {k for k in wide_all if lg[k]["h_witt"] == 0 and lg[k]["sigma"] == 0}
wide_top = {k for k in wide_all if lg[k]["h_witt"] >= 1 or lg[k]["sigma"] != 0} | set(fmw["fox_milnor_fails"])
numbers["nWideTopOne"], numbers["nWideFMHolds"] = len(wide_top), len(fmw["fox_milnor_holds"])
assert (len(wide), numbers["nWideMobiusAlt"], len(wide_top), numbers["nWideFMHolds"]) == (12, 8, 12, 4)
# Final counts ([C] Theorem C, Section 3.4): Propositions 3.1 to 3.3.
numbers["nGammaOneAll"] = numbers["nGammaOne"] + len(wide_all)
numbers["nGammaTopOneAll"] = numbers["nGammaTopOne"] + len(wide_top)
numbers["nMobiusOneTwoAll"] = numbers["nMobiusOneTwo"] + count(wide_all, "[1,2]")
numbers["nMobiusOneThreeAll"] = numbers["nMobiusOneThree"] + count(wide_all, "[1,3]")
numbers["nGammaTwoAll"] = first_gamma_two + len(newd) + len(dcert)
numbers["nDetermined"] = first_determined + len(newd) + len(wide_all) + len(dcert)
numbers["nOpenLeftOneTwo"] = (first_open_one_two - numbers["nNDOneTwo"] - count(wide_all, "[1,2]")
                              - numbers["nDefOneTwo"])
numbers["nOpenLeftFromOneThree"] = (numbers["nFirstOpenFromOneThree"] - numbers["nNDOneThree"]
                                    - count(wide_all, "[1,3]") - numbers["nDefOneThree"])
numbers["nOpenLeft"] = numbers["nOpenLeftOneTwo"] + numbers["nOpenLeftFromOneThree"]
numbers["nGammaTwoOneTwo"] = numbers["nNewTwo"] + numbers["nNDOneTwo"] + numbers["nDefOneTwo"]
numbers["nNewDAll"] = numbers["nNewD"] + numbers["nDefCert"]            # [C] Proposition 3.2, (a)-(c)
numbers["nKleinTotal"] = numbers["nKleinAll"] + numbers["nMobiusKlein"]  # [C] Proposition 3.3(b)
numbers["nGammaTwoOneThree"] = numbers["nKleinCertified"] + numbers["nNDOneThree"] + numbers["nDefOneThree"]
assert numbers["nDetermined"] + numbers["nOpenLeft"] == numbers["nOpen"]
assert (numbers["nGammaTwoAll"], numbers["nDetermined"], numbers["nOpenLeft"], numbers["nOpenLeftFromOneThree"]) == \
    (2329, 3300, 225, 55)
# [C] Section 3.3: the continuation of the second construction ran on the knots left open after Proposition 3.2(c);
# its Moebius bands are rebuilt, and its controls (of mobius_open_search.py) get none.
oleft = oleft3 - dcert
mjss = mjs["summary"]
assert {k for k, v in mjs["knots"].items() if v["group"] == "target"} == oleft and wj <= oleft
assert not mjss["contradictions"] and not mjss["recheck_failing"]
jreb = load("mobius_J_rebuild.json")["summary"]["moebius_bands"]
assert set(jreb["rebuilt"]) == set(jreb["diagrammatic"]) == wj and not jreb["failing"] and jreb["b1"] == {"1": len(wj)}
jctl = {k for k, v in mjs["knots"].items() if v["group"] == "control"}
assert jctl <= {k for k, v in mos["knots"].items() if v["group"] == "control"}
assert all(not mjs["knots"][k]["certificates"] for k in jctl)
numbers["nJControls"] = len(jctl)
numbers["nJStopped"] = len(set(mjss["stopped_by_time"]) & oleft)
assert (numbers["nJControls"], numbers["nJStopped"]) == (10, 16)
# [C] Section 3.4: the knots that remain open are not smoothly slice, so their range is [1,2]: h_lambda = 1, else a
# nonzero signature, else the smooth 4-genus 1 of KnotInfo (as for the knots of [C] Proposition 3.3(a)).
ofinal = oleft - wj
assert len(ofinal) == numbers["nOpenLeft"] and all(lg[k]["h_witt"] in (0, 1) for k in ofinal)
oh1 = {k for k in ofinal if lg[k]["h_witt"] == 1}
osig = {k for k in ofinal - oh1 if lg[k]["sigma"] != 0}
oki = ofinal - oh1 - osig
assert all(op[k]["smooth_four_genus"].strip() == "1" for k in oki)
numbers["nOpenLeftHOne"], numbers["nOpenLeftSig"], numbers["nOpenLeftKIgfour"] = len(oh1), len(osig), len(oki)
assert (len(oh1), osig, oki) == (224, set(), {"13n_3505"})
# [C] Section 3.4: the Moebius bands of the second construction are not slice either, which gives gamma_4 = 1 for them.
wh1 = {k for k in wide_all if lg[k]["h_witt"] == 1}
wsig = {k for k in wide_all - wh1 if lg[k]["sigma"] != 0}
wki = wide_all - wh1 - wsig
assert all(op[k]["smooth_four_genus"].strip() == "1" for k in wki)
assert wki == set(fmw["fox_milnor_fails"]) | set(fmw["fox_milnor_holds"]) and fmw["fox_milnor_fails"] == ["13a_4729"]
numbers["nWideHOne"], numbers["nWideKIgfour"] = len(wh1), len(wki)
assert (len(wh1), wsig, len(wki)) == (9, {"13a_3125", "13a_4799"}, 5)
# Sharpness of h_lambda in the census ([C] Section 3.4): every prime knot with at most 13 crossings and
# h_lambda >= 2 has gamma_4 = gamma_top = h_lambda.  Values: our determinations (h_lambda = 3: the 46 knots of
# [C] Section 3.4, 9_49 and 11n_133 by [C] Section 3.8, 8_18 from KnotInfo; h_lambda = 2: [C] Section 3.4), and
# otherwise the exact KnotInfo value, which is then 2.
g4exact = load("gamma4_knotinfo_2026.10.5.json")
sharp = Counter()
for k, v in lg.items():
    h = v["h_witt"]
    if h < 2:
        continue
    if k in h3:
        assert k in ob2 | {"8_18", "9_49", "11n_133", "13a_2890"}
        sharp[(h, 3)] += 1
    elif k in new2 | kt | rt:
        sharp[(h, 2)] += 1
    else:
        assert str(g4exact[k]).strip() == "2" and k not in rng, k
        sharp[(h, 2)] += 1
assert set(sharp) == {(2, 2), (3, 3)}
numbers["nSharpTwo"], numbers["nSharpThree"] = sharp[(2, 2)], sharp[(3, 3)]
numbers["nSharpAll"] = numbers["nSharpTwo"] + numbers["nSharpThree"]
assert (numbers["nSharpTwo"], numbers["nSharpThree"]) == (4916, 49)
assert {k for k in ob2 if rng.get(k) == "[1,3]"} <= ob1
assert not any(rng.get(k) in ("[1,4]", "[1,5]") for k in ob1 | oba | obm)

# Controls and independent checks (13 crossings).
c1 = load("rank1_controls.json"); numbers["nControlsOne"] = len(c1); assert not any(v["obstructed"] for v in c1.values())
cd = load("rank1_dinv_controls.json"); numbers["nDAltControls"] = len(cd); assert not any(v["obstructed"] for v in cd.values())
cm = load("rank1_dinv_mcontrols.json"); numbers["nDMontControls"] = len(cm); assert not any(v["obstructed"] for v in cm.values())
val = load("rank1_validate.json"); numbers["nValidate"] = len(val)
assert all(v["agree"] for v in val.values()) and max(v["det"] for v in val.values()) < 200
gz = load("rank1_goeritz.json"); numbers["nGoeritzChecked"] = len(gz); assert all(v["agree"] for v in gz.values())
c2 = load("rank2_controls.json"); numbers["nRankTwoControls"] = len(c2)
assert {k for k, v in c2.items() if v["obstructed"]} == {"9_49", "11n_133"}
jt = load("jk_targets.json"); objk = {k for k, v in jt.items() if v["obstructed"]}
numbers["nJKTargets"], numbers["nJKObstr"] = len(jt), len(objk)
numbers["nJKLinking"] = len({k for k in jt if r1[k]["obstructed"]})
assert objk == {k for k in jt if r1[k]["obstructed"]} | {k for k in jt if k in da and da[k]["obstructed"]}
assert (len(jt), len(objk), numbers["nJKLinking"]) == (1028, 360, 141)
jc = load("jk_controls.json"); assert len(jc) == 3030 and not any(v["obstructed"] for v in jc.values())
# [C] Section 3.5: the signed criterion (Proposition 2.7 for the Goeritz forms of alternating diagrams) excludes none of
# the alternating controls; with the opposite sign it would exclude many of them.
assert all(v["definite"] for v in jc.values()) and not any(v["obstructed_signed"] for v in jc.values())
assert set(jc) == set(cd)                                          # the alternating controls of Proposition 2.6
assert {k for k, v in jt.items() if v["obstructed_signed"]} == objk     # signed and unsigned agree on the targets
numbers["nJKControls"] = len(jc)
numbers["nJKReversed"] = sum(1 for v in jc.values() if (v["s"] == 2 and not v["embeds"][0]) or (v["s"] == 6 and not v["embeds"][1]))
assert numbers["nJKReversed"] == 1769

# Fourteen crossings (table of Lee and Sabloff).
ls = load("rank1_ls14.json")
numbers["nLS"] = len(ls)
lsd = {h: sum(1 for v in ls.values() if v["h_witt"] == h) for h in range(4)}
assert sum(lsd.values()) == len(ls)                                      # h_lambda <= 3 for every knot ([C] Theorem 2.2)
numbers["nLSHthree"] = lsd[3]
assert (len(ls), lsd[0], lsd[1], lsd[2], lsd[3]) == (46972, 2195, 25996, 18585, 196)
assert all(v["consistent"] for v in ls.values()) and all(v["obstructed"] == v["obstructed_summand"] for v in ls.values())
band = {k for k, v in ls.items() if v["band"]}
lsarf = {k for k, v in ls.items() if v["arf"]}
lslink = {k for k, v in ls.items() if v["linking"] == "Obstructed"}
lsopen = {k for k in ls if k not in band and k not in lsarf and k not in lslink}
numbers["nLSBand"], numbers["nLSArf"], numbers["nLSLinking"], numbers["nLSOpen"] = len(band), len(lsarf), len(lslink), len(lsopen)
assert not any(ls[k]["obstructed"] for k in band)
assert all(ls[k]["obstructed"] for k in lsarf | lslink)
lsopen_lf = {k for k in lsopen if ls[k]["obstructed"]}
numbers["nLSOpenLinking"] = len(lsopen_lf)
dl = load("rank1_dinv_ls14.json"); obdl = {k for k, v in dl.items() if v["obstructed"]}
assert obdl <= lsopen and not (obdl & lsopen_lf)
numbers["nLSDTargets"], numbers["nLSDObstr"] = len(dl), len(obdl)
numbers["nLSNew"] = len(lsopen_lf | obdl)
assert (len(lsopen), len(lsopen_lf), len(dl), len(obdl), numbers["nLSNew"]) == (11247, 1082, 5189, 1020, 2102)
dlc = load("rank1_dinv_ls14controls.json"); numbers["nLSDControls"] = len(dlc)
assert not any(v["obstructed"] for v in dlc.values())
h3ls = {k for k, v in ls.items() if v["h_witt"] == 3}
numbers["nLSHthreeAlt"] = sum(1 for k in h3ls if k.startswith("K14a"))
assert numbers["nLSHthreeAlt"] == 83
r2ls = load("rank2_ls14.json"); assert {k for k, v in r2ls.items() if v["obstructed_odd"]} == h3ls
numbers["nMaxDetFourteen"] = max(v["det"] for v in ls.values()); assert numbers["nMaxDetFourteen"] == 1145
jl = load("jk_ls14.json"); objl = {k for k, v in jl.items() if v["obstructed"]}
assert objl == {k for k in jl if ls[k]["obstructed"]} | obdl
numbers["nLSJKTargets"], numbers["nLSJKObstr"] = len(jl), len(objl)
numbers["nLSJKLinking"] = len({k for k in jl if ls[k]["obstructed"]})
assert (len(jl), len(objl), numbers["nLSJKLinking"]) == (5754, 1585, 565)

# [C] Section 4: locally flat but not smooth Moebius bands.
gap = load("smooth_top_gap.json")
gs = gap["summary"]
numbers["nGapCandidates"], numbers["nGapAlt"], numbers["nGapMont"] = (gs["candidates"][k] for k in
                                                                    ("count", "alternating", "montesinos"))
assert (numbers["nGapCandidates"], numbers["nGapAlt"], numbers["nGapMont"]) == (81, 72, 9)
assert gs["candidates"]["h_witt_all_one"] and numbers["nGapCandidates"] == numbers["nDObstrOneTwo"]
assert set(gap["candidates"]) == {k for k in oba | obm if rng.get(k) == "[1,2]"}
numbers["nGapDiagrams"], numbers["nGapBands"] = gs["diagrams"], gs["bands_tried"]
numbers["nGapKnots"], numbers["nGapDetOne"], numbers["nGapHits"] = gs["knot_results"], gs["det1_knots"], gs["hits"]
numbers["nGapCertified"] = len(gs["certified"])
numbers["nGapNoForm"] = len(gs["not_certified_without_rank_one_form"])
numbers["nGapRest"] = numbers["nGapCandidates"] - numbers["nGapCertified"] - numbers["nGapNoForm"]
rest = set(gap["candidates"]) - set(gs["certified"]) - set(gs["not_certified_without_rank_one_form"])
numbers["nGapRestNoDetOne"] = sum(1 for k in rest if gap["candidates"][k]["det1_knots"] == 0)
numbers["nGapRestDetOne"] = len(rest) - numbers["nGapRestNoDetOne"]
assert (len(rest), numbers["nGapRest"], numbers["nGapRestDetOne"], numbers["nGapRestNoDetOne"]) == (61, 61, 55, 6)
assert (gs["diagrams"], gs["bands_tried"], gs["knot_results"], gs["det1_knots"], gs["hits"]) == \
    (24970, 39914728, 19761630, 22453, 6309)
assert gs["diagrams_not_isometric"] == gs["contradictions"] == gs["hits_failing_a_recheck"] == 0
assert gs["det1_knots_without_rank_one_form"] == 0 and gs["fit_failures"] == 0
numbers["nFoxDetChecks"] = gs["fox_det_checks"][0]; assert gs["fox_det_checks"] == [79352, 79352]
gv = gap["validation"]
numbers["nAlexValidated"], numbers["nAlexHFK"] = gv["fox_ok"], gv["hfk_ok"]
assert gv["fox_ok"] == gv["seifert_ok"] == gv["knots"] == 12965 and gv["hfk_ok"] == gv["hfk_checked"] == 663
gc = load("smooth_top_gap_check.json")
assert gc["all_checks_pass"] and gc["hits_with_tau_zero"] == 0 and len(gc["knots"]) == 12
# [C] Section 4, the knots of Proposition 3.2(a), (b): the band searches of plumbing_dinv.py and surgery_dinv.py (the search of
# smooth_top_gap.py) on all of them, and the representative bands rechecked by scripts/gap_new.py.
gn = load("gap_new.json")
gns = gn["summary"]
assert gns["candidates"] == gns["searched"] == len(newd) and set(gn["knots"]) <= newd
assert all(r["J"]["alexander_is_one"] for r in gn["knots"].values())
assert ps["gap"]["errors"] == [] and sdv["summary"]["gap"]["contradictions"] == 0
assert set(ps["gap"]["certified"]) == {k for k, r in gn["knots"].items() if r["source"] == "plumbing"}
cg = load("surgery_cfk_gap.json")["summary"]
assert cg["contradictions"] == 0 and set(cg["knots"]) == cob and set(cg["certified"]) == cob
assert {k for k in sdv["summary"]["gap"]["certified"] if k in sob - pob} | set(cg["certified"]) == \
    {k for k, r in gn["knots"].items() if r["source"] == "surgery"}
numbers["nGapNewCandidates"] = gns["candidates"]
numbers["nGapNew"] = len(gn["knots"])
numbers["nGapNewPlumb"], numbers["nGapNewSurg"] = gns["by_source"]["plumbing"], gns["by_source"]["surgery"]
numbers["nGapAll"] = numbers["nGapCertified"] + numbers["nGapNew"]
numbers["nGapNewKI"] = len(gns["representative_from_knotinfo_diagram"])
numbers["nGapTauZero"], numbers["nGapTauZeroJ"] = len(gns["tau_zero"]), gns["tau_zero_J_distinct"]
# the four J are pairwise distinct: the total ranks of their knot Floer homology are 97, 97, 257 and 145, and the two
# of rank 97 have Seifert genera 3 and 4
assert [(gn["knots"][k]["J"]["total_rank"], gn["knots"][k]["J"]["seifert_genus"])
        for k in ("13n_1281", "13n_1307", "13n_2010", "13n_2932")] == [(97, 3), (97, 4), (257, 4), (145, 3)]
assert gns["tau_zero"] == ["13n_1281", "13n_1307", "13n_2010", "13n_2932"]          # named in [C] Section 4
assert (numbers["nGapNew"], numbers["nGapNewPlumb"], numbers["nGapNewSurg"], numbers["nGapAll"], numbers["nGapNewKI"],
        numbers["nGapTauZeroJ"]) == (29, 19, 10, 41, 13, 4)
# KnotInfo 2026.10.5 lists a range for the topological value of every knot of [C] Theorem B
gaptop = Counter(op[k]["topological_4d_crosscap_number"].replace(" ", "")
                 for k in set(gs["certified"]) | set(gn["knots"]))
numbers["nGapTopOneThree"], numbers["nGapTopOneFive"], numbers["nGapTopOneSeven"] = gaptop["[1,3]"], gaptop["[1,5]"], gaptop["[1,7]"]
assert gaptop == Counter({"[1,3]": 35, "[1,5]": 4, "[1,7]": 2})
# the knots with gamma_4 = 2 whose upper bound comes from punctured Klein bottles built from band moves
numbers["nKleinTwo"] = numbers["nKleinAll"] + numbers["nNDOneThree"] + numbers["nDefOneThree"]
assert numbers["nKleinTwo"] == 2057
# [C] Section 3.6: Jabuka's rational Witt span ws, recomputed from the residues of V + V^T with Hilbert symbols and the
# Hasse-Minkowski theorem (witt_span.py), against the closed form of [T] Theorem 4.6(c), (d) and against h_lambda.
wsd = load("witt_span.json")
wsk, wss = wsd["knots"], wsd["summary"]
assert set(wsk) == set(lg) and wss["ws_differs_from_theorem"] == wss["h_lambda_differs_from_theorem"] == []
assert all(v["h_lambda"] == lg[k]["h_witt"] for k, v in wsk.items())
assert all(v["ws"] == v["h_lambda"] or (v["ws"], v["h_lambda"]) == (2, 3) for v in wsk.values())  # [T] Theorem 4.6(b)
assert wss["pairs"] == {"0,0": 725, "1,1": 7275, "2,2": 4916, "2,3": 17, "3,3": 32}
cert = load("hthree_certificates.json")["knots"]
only = lambda k, c: cert[k]["criteria"] == [c]
mod8 = lambda k: (cert[k]["sigma"] % 8, cert[k]["det"] % 8)
ws2 = {k for k, v in wsk.items() if v["ws"] < v["h_lambda"]}
assert ws2 == {k for k in cert if only(k, "(iii)") and mod8(k) == (0, 5)}                          # [T] Theorem 4.6(d)
assert "8_18" in ws2 and all(wsk[k]["even_binary_form"] for k in ws2)
assert wsk["8_18"]["even_binary_form"] == [[6, -3], [-3, -6]]
# rank2.py finds an even binary form of odd determinant (linking signs, no odd signs) for exactly these knots
r2all = {**load("rank2_controls.json"), **load("rank2_targets.json"), **load("rank2_checks.json")}
assert ws2 == {k for k in cert if r2all[k]["linking_signs"] and not r2all[k]["odd_signs"]}
numbers["nWsTwo"] = len(ws2); assert len(ws2) == 17
ws3 = {k for k in cert if wsk[k]["ws"] == 3}
ws3i = {k for k in ws3 if "(i)" in cert[k]["criteria"]}
ws3iii = ws3 - ws3i
assert ws3 == set(cert) - ws2 and not any("(ii)" in c["criteria"] for c in cert.values())
assert ws3iii == {"12a_569", "12n_626", "13a_1974", "13a_3121", "13n_3005", "13n_4927"}
assert all(only(k, "(iii)") and mod8(k) == (4, 1) for k in ws3iii)                                  # [T] Theorem 4.6(c)
assert not any(r2all[k]["linking_signs"] for k in ws3iii)       # and rank2.py finds no binary form of either parity
numbers["nWsThree"], numbers["nWsThreeOne"], numbers["nWsThreeIII"] = len(ws3), len(ws3i), len(ws3iii)
assert (len(ws3), len(ws3i), len(ws3iii)) == (32, 26, 6)
# Jabuka's Example 1.7 lists the knots with at most 12 crossings that satisfy his Thm 1.6(ii): a single prime p = 1 mod 4
# with a nonzero residue, of rank 2.  They are the type (i) knots of [C] Table 3 up to 12 crossings; 12a_569 and 12n_626 have
# ws = 3 without that criterion.  The values of ws in Jabuka's Example 3.5 (and Example 1.8 for 9_40) are reproduced.
crossings = lambda k: int(k.split("_")[0].rstrip("an"))
jab = {"9_49", "11n_133", "12a_561", "12a_664", "12a_780", "12a_907", "12n_276"}
assert {k for k, c in cert.items() if "(i)" in c["criteria"] and crossings(k) <= 12} == jab


def jabuka_16ii(k):
    nonzero = [int(p) for p, a in wsk[k]["a"].items() if a]
    return len(nonzero) == 1 and nonzero[0] % 4 == 1 and wsk[k]["a"][str(nonzero[0])] == 2


assert {k for k in wsk if crossings(k) <= 12 and jabuka_16ii(k)} == jab
assert {k for k in ws3iii if crossings(k) <= 12} == {"12a_569", "12n_626"} and not any(map(jabuka_16ii, ws3iii))
assert wss["jabuka_examples"] == {"9_40": 2, "7_5": 2, "9_49": 3, "11n_133": 3, "12a_561": 3, "12a_664": 3,
                                  "12a_780": 3, "12a_907": 3, "12n_276": 3, "3_1": 1, "5_1": 1, "7_1": 1}
# Fourteen crossings: the 196 knots with h_lambda = 3, from a Goeritz matrix of the SnapPy diagram.
s14, k14 = wsd["summary_ls14"], wsd["knots_ls14"]
assert set(k14) == h3ls and s14["ws_differs_from_theorem"] == s14["h_lambda_differs_from_theorem"] == []
assert s14["ws_distribution"] == {"2": 50, "3": 146}
assert all(k14[k]["criteria"] == ["(iii)"] and (k14[k]["sigma8"], k14[k]["det"] % 8) == (0, 5)
           for k in s14["ws_less_than_h_lambda"])
numbers["nLSWsTwo"] = len(s14["ws_less_than_h_lambda"]); assert numbers["nLSWsTwo"] == 50
mc = load("montesinos_dinv_check.json")["summary"]
assert mc["targets_obstructed"] and mc["controls_obstructed"] == [] and mc["d_spin_disagreements"] == []
numbers["nMontCheckControls"] = mc["controls_checked"]; assert mc["controls_checked"] == 157   # seeded run
# [C] Section 3.7: 11n_57, for which Fairchild (Prop. 3.2), the table of Lee-Sabloff and the KnotInfo package 2026.8.1
# give gamma_4 = 1.  Both correction-term computations exclude a smooth Moebius band, with the same d of the spin
# structure, and the knot bounds a punctured Klein bottle (range [1,3] in KnotInfo 2026.10.5).
ml = load("montesinos_dinv_check_literature.json")["knots"]["11n_57"]
assert ml["obstructed"] and ml["obstructed_rank1_dinv"] and ml["d_spin"] == ml["d_spin_rank1_dinv"] == "3/2"
assert kk["11n_57"]["knotinfo_gamma4"] == "[1,3]" and "11n_57" in kt and lg["11n_57"]["h_witt"] == 1
# [C] Section 2.1: the knots K and K' of [T] Examples 3.5 and 3.6 built from PD codes (scripts/small_examples.py).
with open(os.path.join(ROOT, "data", "small_examples.txt")) as f:
    se = f.read().splitlines()
assert se[-1] == "all checks pass: True" and "complement of 8_10: True" in se[2] and "complement of 8_11: True" in se[4]
assert "[[0, 3, 0, 0], [3, -3, 1, 0], [0, 1, 1, -2], [0, 0, -2, 1]]" in se[4] and "complement of 10_143: True" in se[3]
# The PD codes of K and K' quoted in [C] Section 2.1, with the checks repeated on them.
PD_K = [(20, 14, 1, 13), (14, 2, 15, 1), (2, 16, 3, 15), (5, 12, 6, 13), (11, 4, 12, 5), (3, 10, 4, 11),
        (16, 9, 17, 10), (17, 7, 18, 6), (7, 19, 8, 18), (19, 9, 20, 8)]
PD_KP = [(22, 9, 1, 10), (10, 21, 11, 22), (20, 11, 21, 12), (17, 9, 18, 8), (7, 19, 8, 18), (19, 7, 20, 6),
         (12, 5, 13, 6), (13, 16, 14, 17), (1, 4, 2, 5), (2, 15, 3, 16), (14, 3, 15, 4)]
assert (f"PD code of K (KnotInfo convention): {PD_K}; Goeritz matrix that of the paper: True (minus it: False); "
        f"complement of 8_10: True") in se
assert (f"PD code of K' (KnotInfo convention): {PD_KP}; Goeritz matrix that of the paper: True (minus it: False); "
        f"complement of 8_11: True") in se


if __name__ == "__main__":
    for key, value in numbers.items():
        print(f"{key:28s} {value:>12,}")
    print(f"{len(numbers)} counts, all equal to the values quoted in [C] and the README")
