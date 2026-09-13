# Magic squares of squares — an arithmetic attack

Does there exist a 3×3 magic square whose nine entries are distinct perfect squares? The question is open since
the 18th century; only one square with **seven** square entries is known (Bremner, independently Sallows).

This repository is the working record of a two-day computational study (10–12 September 2026) of the *arithmetic*
side of the problem: descent, the Cassels–Tate pairing, and Chabauty on the curves that arise from the classical
cross parametrisation. It is a lab notebook, not a paper. It contains failed attempts, corrections, and a list of
our own mistakes alongside the results.

**The problem is not solved here.** Nothing below claims otherwise.

## Status labels

Every claim in this repository carries one of:

| label | meaning |
|---|---|
| **proved** | mathematical proof, checked by hand |
| **proved (SW)** | rigorous algorithm; correctness depends on the software and on our implementation |
| **verified numerically** | computed, not proved |
| **observation** | a pattern we saw |
| **heuristic** | plausible estimate, no proof |

`SUMMARY_2026-09-12.md` is the single entry point: statements, boundaries and labels only, with pointers to the
files that contain the details.

For a readable account of the same work — written to be followed without running any of the code — see the field
notes at <https://magic-squares-field-notes.kellerbm71.chatgpt.site/> (and `/history` there for how the three
systems actually worked). They were written from the other side of the collaboration, so where the two disagree,
the labelled statements in `SUMMARY_2026-09-12.md` are the ones to trust.

## What we believe we contribute

**1. Exclusion of ratios at every scale.** *(proved (SW))*
In a 3×3 magic square of nine distinct squares, neither of the two **edge** pairs (centre together with the two
side midcells) can have ratio (17:7:13), (23:7:17) or (71:49:61) — for any scaling whatsoever.

This differs in kind from the usual searches, which are bounded by the size of the centre (Buell to 25·10²⁴,
Boyer to about 2.5·10²⁷): each such statement removes an infinite family at once. The chain is
ELS classes → a genus-2 curve → its Jacobian as a Weil restriction Res_{k/ℚ}E₁ → a proof that rank E₁(k) = 1 →
elliptic Chabauty.

*Boundary of the result:* edge pairs only. The corner (diagonal) pairs need a separate family — the D₈ symmetries
do not carry edge pairs to corner pairs.

**2. An open implementation of the Cassels–Tate pairing over a quadratic field.**
`descent/ctp.sage` (Cassels' conic method) and `descent/ctp_quartic.sage` (Fisher's binary quartics,
arXiv:2208.14977). PARI's `ellrank.c` implements Fisher's method over ℚ only; Sage has none; Magma's is closed.
Ours runs over ℚ and over quadratic fields, with exact local arithmetic throughout in the quartic version.

Checks: Example 3.4 of Fisher reproduced exactly (I, J, Δ, m, γ₁ and the contribution of each place); on four curves
with full 2-torsion over ℚ(√d) the Fisher and Cassels matrices agree **element-by-element**; on 14 curves with one
2-torsion point and an exact rank from PARI the bound is sharp everywhere; symmetry is checked by recomputing with
the arguments reversed, including on non-zero values.

*Cost and outcome on a real curve (corrected 13 Sep).* An earlier version of this page, and a letter we sent,
said the pairing does not finish on (89,23,65). That was wrong: it finished the same evening. The first pairing
took about 115 minutes, almost all of it factoring a 46-digit norm |N(g(1,0))| to build the set of places;
the remaining 35 took 28 minutes together. Our quartics are reduced but not **minimised** — minimisation of
genus-one models (Cremona–Fisher–Stoll) should cut that norm to about 15 digits, so it is a matter of speed,
not feasibility. The arithmetic result was negative: dim Sel² = 8, pairing rank 2, bound rank ≤ 5 — worse than
the bound ≤ 3 we already had.

## The hardest computation: rank 1 over ℚ(√165)

For the section (71:49:61), 2-descent alone gives only rank ≤ 5. The Cassels–Tate pairing brings it to 1:

| step | result | how it was checked |
|---|---|---|
| full 2-descent over ℚ(√165) | dim Sel² = 7 | recomputed independently |
| Cassels–Tate pairing | rank 4 | all 8120 Hilbert symbols recomputed from scratch, 0 discrepancies |
| rank bound | 7 − 2 − 4 = 1 | re-derived in Cassels' original normalisation, matrix matched element-wise |
| lower bound | a point of infinite order | height agreed in two implementations |
| elliptic Chabauty, p = 7 and 113 | 8 zeros = 8 known points | second independent implementation |

*Caveat, stated plainly:* the cross-check of the quartic method against the conic method **on this curve** reached
only two rows out of seven before we stopped it (`CTP_QUARTIC_2026-09-12.md`, §3.2). The rank-1 conclusion rests on
the conic method.

## What is conditional

Elliptic Chabauty is complete for **all eleven** sections we treat (bound 8, known points 8). The conclusion
"there are no further points" is unconditional only for the three sections where rank 1 is proved. For the other
**eight it is conditional on rank 1**, which we could not establish without Magma.

## What we rediscovered

Several constructions we derived independently turned out to be known. Listing them is part of the record:

| what we derived | where it already was |
|---|---|
| seven squares ⇒ at least two complete pairs | Boyer, 2004 |
| the two-complete-pairs family | Wesołowski, 2019 (seven of nine cells coincide) |
| closing a case via a rank-0 elliptic factor | Bremner, 2001 (in genus 3) |
| the genera 0, 1, 5, 17, 49 | Auel–Singer, arXiv:2609.09351 (September 2026) |
| the parametrisation of pairs | Bremner, 1999 — literally the same formulas |

`NOVELTY_CHECK_2026-09-12.md` has the full comparison.

## Negative results

- **Quadratic Chabauty does not apply** to this family at all: the isogeny E₁ → E₁^σ is defined only over k(i),
  hence End⁰(Res_{k/ℚ}E₁) = ℚ and ρ(NS) = 1, while the method needs ρ ≥ 2.
- **The searches are empty.** 808 million states over pairs with square centre; 113.7 million Mordell–Weil lattice
  states plus 13.8 billion directly enumerated points over Bremner's sixteen configurations; 3.6·10¹³ pairs on a GPU.
  The only hit, everywhere, is the known Bremner–Sallows square — which was our control.
- **Descent via isogenous curves** did not improve the bounds: (89,23,65) ≤ 7, (79,47,65) ≤ 5.

## Mistakes we caught

More than twenty, and — correcting what this page and several of our letters said earlier — they were **not all
of one kind**. Some were a right answer resting on a wrong justification. Others were code that gave a **wrong
answer** (a lost quadratic twist produced wrong ranks for five pairs out of seven), empty searches presented as
proofs of absence (four times), overstated status labels, and a claim proved for a narrow case but stated for a
broad one. An outside review of the session logs of all three systems
(`REVIEW_ALL_SESSIONS_FROM_FABLE_2026-09-13.md`) showed that our one-line summary of our own mistakes was itself wrong. The list, with who caught what, is
in `SUMMARY_2026-09-12.md` §6. The worst of them: our enumeration of preimages was silently swallowing exceptions,
so a "proof" was resting on an empty loop. It surfaced only because a count disagreed with a control example.

## Layout

```
SUMMARY_2026-09-12.md     entry point: every claim, its label and its boundary
FAMILY_SECTIONS.md        running log of the section family (newest entries at the bottom)
descent/                  2-descent over quadratic fields; Cassels–Tate pairing (conics and quartics)
family/                   the pipeline: sections → genus-2 curve → Chabauty
bridge/                   the G1 family: 127 parameter pairs; 126 excluded for a full non-degenerate square
corners/                  Mordell–Weil lattice search, LLL + sieve
bremner16/                Bremner's sixteen six-square configurations
gpu6/                     GPU sieve (HIP)
```

The many `NOTE_*`, `REPLY_*`, `REVIEW_*` files are the working correspondence between the three systems that did
this (Claude, GPT, Grok), kept because the corrections in them are part of how the results were reached.

## How this was done

Three AI systems worked in parallel on one machine, exchanging notes as files in this directory, reading each
other's session logs, and recomputing each other's results with independent implementations. A human set the
direction, stopped dead ends, and required a status label on every claim. Disagreement between two implementations
was treated as a defect to be found, never as noise.

**A caveat about the word "independent".** The same review found our cross-checks less independent than we claimed:
before several "independent" reviews one system had read the other's raw session log; one certificate was produced
from a snapshot of the other system's code; one re-check simply re-ran the other system's scripts; and one verifier
was told the expected answer in advance. Agreement between the systems here is weaker evidence than it sounds.
From 13 Sep a check counts as independent only if its note lists everything read before computing, and that list
excludes the checked party's code, logs and working directory.

## What is missing

**Magma.** `TwoCoverDescent` and `Chabauty` would settle the eight conditional sections and the two remaining
parameter pairs (11,4) and (15,8). We have no licence and know of no open substitute.

**Human expertise** on rational points on surfaces of general type.

## References

- A. Bremner, *On squares of squares*, Acta Arith. 88 (1999); II, Acta Arith. 99 (2001)
- C. Boyer, *Some notes on the magic squares of squares problem*, Math. Intelligencer 27 (2005) 52–64;
  <http://www.multimagie.com/English/SquaresOfSquares.htm>
- A. Auel, B. Singer, *The algebraic geometry of 3-by-3 magic squares of squares*, arXiv:2609.09351
- T. Fisher, *On binary quartics and the Cassels–Tate pairing*, arXiv:2208.14977
- J. W. S. Cassels, *Second descents for elliptic curves*, J. reine angew. Math. 494 (1998)
- L. Wesołowski (2019), L. Morgenstern (2010, 2015) — via multimagie.com
