# vendor/

Unmodified copies of eight modules from the author's earlier code: three from the code for unknotting numbers of
alternating knots, and five from the code for negative-definite plumbings (MIT License, Seong-Jin Lee; the repository
negative-definite-plumbings of arXiv:2610.25430, commit `01dca0e`).  The scripts reach them through `vendor/`.

| file | sha256 | used by |
|---|---|---|
| `code/lower_bounds/owens/owens_obstruction.py` | `66d19b2b358f2d1cf5dc6bae0a6fbba4cfcb09a7c0ebc48b3ed5f5f73bc7ede5` | `rank1_dinv`, `jk_lattice`, `rank2` |
| `code/lower_bounds/montesinos/montesinos_u1.py` | `28ba5e47a7699c3d71ab45fd901c9901d89d7cde4e5b0def1c9dc24ad07c86e0` | `rank1_dinv` (modes `mcontrols`, `mtargets`), `plumbing_dinv` (mode `controls`) |
| `code/lower_bounds/linking_pairing.py` | `b79c856105fb0ed925fb783995f1b861d4bbec5d895dda9d9bd2efbe61252b4a` | `montesinos_u1` |
| `paper_followup/negative-definite-plumbings/rank_two/bf2.py` | `945226d46ebd212b472fcbae765449ecc6f2a96ad53d95a012aee6ae94cc48fb` | `rank1_linking.brute` |
| `paper_followup/negative-definite-plumbings/plumbing/neumann/formula.py` | `f1e1d35e9b46d8b6e66f74fb75504b7be2d21a59750315a8e4fdd88d6b4313f1` | `plumbing_dinv` (Regina names of trees by the calibrated formula) |
| `paper_followup/negative-definite-plumbings/plumbing/plumbing.py` | `057c6a815accabc58f8109300f263a6f027c02e811a045d430695136b91d57d0` | `plumbing_dinv` (`plumbed_manifold`: a tree as surgery on its plumbing link) |
| `paper_followup/negative-definite-plumbings/plumbing/algebraic/algebraic_plumbing.py` | `e37fcea7be99be5658c1ffef2cb82d33c1b4fb8764f0ab61e3b3cf636f65c490` | `plumbing_dinv` (the tree of an algebraic diagram and its normal form) |
| `paper_followup/negative-definite-plumbings/plumbing/algebraic/basic_polyhedron.py` | `492ad4532dcb7189448f3f24ff0347e8b28ac566413895ce3dfcf254eb87d82f` | `plumbing_dinv`, `algebraic_plumbing` (bigon merging) |

Only the functions that the scripts call are supported.  Some other functions need files that are not included:

| function | missing files |
|---|---|
| `owens_obstruction.definite_goeritz_candidates` | `_compat`, `taitgraph` |
| the σ = 2 branch of the Owens test | `traczyk` |
| `montesinos_u1.test_knot` | `hf_red` |
| `montesinos_u1.main` | the KnotInfo reference table in `results/` |
| `plumbing.regina_name` | none, but `plumbing_dinv` uses its own `fast_regina_name` (that of `graph_u1.py`), since `simplifyExhaustive` does not terminate in practice on the larger plumbing triangulations |

The other functions of the follow-up pipeline that `plumbing_dinv.py` needs (`general.py`, `invert.py`,
`neumann/precheck.py`, `neumann/both_orientations.py`, `neumann/identify_trees.py`, `graph_u1.py`, `reach_orb.py`,
`recognise_more.py`, `recognise_retry.py`, `verify_hyp.py`, `jsj.py`, `jsj_verify.py`, `pachner_recognise.py`,
`algebraic/alg_diagrams.py`, `algebraic/crosscheck_u.py`, `algebraic/u_plumbing.py`, `algebraic/cw_general.py`) are
copied into `scripts/plumbing_dinv.py` itself, each under its original name in the docstring of the function, because
those modules import further modules of that repository at import time.
