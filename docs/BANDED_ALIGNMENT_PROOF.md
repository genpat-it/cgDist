# Certified banded alignment: why it is bit-identical to parasail

cgDist counts SNPs, InDel events and InDel bases from a global alignment of
two alleles. Up to 0.1.4 every pair was aligned by parasail over the full
dynamic-programming (DP) matrix (`nw_trace_striped_sat`). Since then
(`src/core/banded.rs`) cgDist first tries a **certified banded alignment**:
it fills only a diagonal band of the matrix and returns a result only if it
can **prove** that the result equals the full-matrix result. Otherwise the
pair goes to parasail as before.

This document states what is proved, what is verified, and how to rerun
every check.

## 1. The specification: parasail's alignment

Let `q` (length `n`, rows `i`) be the query and `r` (length `m`, columns
`j`). Scores come from `parasail_matrix_create("ACGT", match, mismatch)`:
case-insensitive `A/C/G/T` score `match` against themselves and `mismatch`
against each other, and every other byte (e.g. `N`) scores `0` against
anything. Gaps are affine, and a gap of length `L` costs `open + (L-1)*extend`
(`open, extend >= 0`).

Recurrences (parasail `nw_trace.c`), for `i, j >= 1`:

```
F(i,j) = max(H(i-1,j) - open, F(i-1,j) - extend)      vertical   (deletion)
E(i,j) = max(H(i,j-1) - open, E(i,j-1) - extend)      horizontal (insertion)
H(i,j) = max(H(i-1,j-1) + s(q_i, r_j), E(i,j), F(i,j))
```

The boundary is `H(0,0) = 0`, `H(0,j) = -open-(j-1)*extend`,
`H(i,0) = -open-(i-1)*extend`, and `E(i,0) = F(0,j) = -inf`.

Trace flags at each cell, which also fix the **tie-breaking**:

* `H`: `DIAG` if `H == diagonal`, else `DEL` if `H == F`, else `INS`. The
  diagonal wins ties, then deletion, then insertion.
* `E`: `DIAG_E` if `E_open > E_extend` (strict), else `INS_E`. Extension wins
  ties.
* `F`: `DIAG_F` if `F_open > F_extend` (strict), else `DEL_F`.

The traceback (parasail `cigar_template.c` / `traceback_template.c`) starts
in state `H` at `(n,m)` and follows these flags as a deterministic state
machine. Once one sequence is exhausted, the rest of the other is gaps.
cgDist counts statistics over the resulting alignment columns
(`compute_alignment_stats`): a SNP is an aligned column whose two **bytes**
differ, and every maximal run of gap columns is one InDel event.

The alignment is therefore a **deterministic function** of the trace flags
on the path that the traceback visits.

## 2. The band and the certificate

A band is the set of cells with diagonal offset `d = j - i` in
`[dlo, dhi]`, where `dlo <= min(0, m-n)` and `dhi >= max(0, m-n)`. The banded
DP uses the same recurrences, with `H, E, F = -inf` outside the band.

**Lemma 1 (banded values never exceed full values).** For every cell,
`X_band <= X_full` for `X` in {H, E, F}. Each banded value is a maximum over
a subset of the full value's candidates, the others being `-inf`. The proof
is by induction over the cells in DP order.

**Lemma 2 (banded values are in-band optima).** For a cell in the band,
`X_band` is the best score among the paths from `(0,0)` to that cell (ending
in the move type of `X`) that stay inside the band. If some optimal path to
that cell lies in the band, then `X_band = X_full`.

**Lemma 3 (paths that leave the band).** Take any global alignment path that
visits a cell with offset `d > dhi >= max(0, m-n)`.

* To reach offset `d` it needs at least `d` insertion columns. It must also
  end at offset `m-n`, so it needs at least `d-(m-n) > 0` deletion columns.
* Therefore it has at least `G = 2d-(m-n)` gap columns, including at least
  one insertion run and one deletion run.
* It has at most `m-d` aligned columns, each scoring at most
  `smax = max(match, mismatch, 0)`.
* With `k >= 2` gap runs its gap penalty is `k*open + (G-k)*extend`. This is
  at least `2*(open-extend) + G*extend` when `open >= extend`, and at least
  `G*open` otherwise.

Its score is therefore at most

```
U(d) = smax*(m-d) - mingap(G)
```

and `U` is non-increasing in `d`, so `U(dhi+1)` bounds every path that
exceeds `dhi`. The case `d < dlo` is symmetric. This is
`outside_upper_bound` in the code.

**Theorem (certificate).** Let `S = H_band(n,m)`. Suppose `S > U(dhi+1)`
(or `dhi` already reaches the last column) and `S > U(dlo-1)` (or `dlo`
already reaches the last row). Then:

1. `S = H_full(n,m)`, and **every** optimal global path lies inside the band;
2. the banded traceback visits exactly the same cells and states, and reads
   exactly the same flags, as the full-matrix traceback. It therefore
   produces the same alignment, the same SNP/InDel counts and the same score.

*Proof.* (1) In-band paths reach `S` (Lemma 2). By Lemma 3, a path leaving
the band scores at most `U < S`. So the optimum is `S`, and no optimal path
leaves the band.

(2) By induction over traceback steps we keep the following invariant: the
current state `(c, X)` lies on an optimal global path, and
`X_band(c) = X_full(c)`. It holds at `(n,m)`.

*State H at cell `c`.* Consider the three candidates of `H(c)`: diagonal,
`F(c)` and `E(c)`.

* If a candidate equals `H_full(c)` in the full matrix, it extends to an
  optimal global path (its best prefix followed by the current optimal
  suffix). That path lies in the band by (1), so by Lemma 2 the candidate
  has the same banded value, and it also equals `H_band(c)`.
* If a candidate is strictly below `H_full(c)`, then by Lemma 1 its banded
  value is also strictly below `H_full(c) = H_band(c)`.

So the pattern of equalities that determines the flag, with its fixed
priority `DIAG > DEL > INS`, is the same in both matrices, and the chosen
predecessor satisfies the invariant.

*State E (and symmetrically F).*

* If `E_open > E_extend` in the full matrix, the opening move is on an
  optimal path, so `H(i,j-1)` has the same value in the band. Its extension
  alternative can only be smaller in the band (Lemma 1), so the band also
  picks `DIAG_E`.
* If `E_open <= E_extend`, the extension is on an optimal path and is
  unchanged in the band. The opening alternative can only be smaller, so
  the band again picks `INS_E`.

The boundary moves (one sequence exhausted) depend only on positions. ∎

The theorem holds for **every** pair of sequences, every length and every
non-negative gap penalty, including inputs no test has ever seen.

**Choosing the band.** A band whose certificate can succeed is computed from
a lower bound `LB <= S_full`: the score of the best "single gap" alignment,
which is a valid alignment (`simple_lower_bound`). The smallest `w` with
`U < LB` on both sides is used. Since `S_band >= LB` once the band contains
that alignment, the certificate then holds. If the band would exceed half
of the matrix, the pair goes to parasail.

## 3. Implementation invariants

`align_in_band_dp` computes the banded recurrences on anti-diagonals
`k = i + j`, storing cell `(i,j)` at slot `(d - dlo + 2) >> 1` of the arrays
of parity `k & 1`. It relies on these invariants:

* The diagonal predecessor `(k-2, d)` is the same slot of the same array,
  updated in place. `F` reads `(k-1, d+1)` and `E` reads `(k-1, d-1)` from
  the other parity.
* The band edges `d = dlo-1` and `d = dhi+1` are fixed slots that always
  hold `-inf`. Block padding can overwrite the right edge, and it is
  restored after every anti-diagonal.
* Every predecessor of an interior cell is a matrix cell, so slots of cells
  outside the matrix, including block-padding garbage, are never read.
* Trace flags are stored at `k*stride + slot`. The traceback only reads
  interior in-band cells, all of which were written for the current pair.
* Arithmetic is on integers only: no floating point and no randomness.
  Every pair is independent and uses per-thread buffers that are fully
  (re)initialised, except the trace buffer, which is only read where it was
  written. Results therefore do not depend on the number of threads, the
  order of evaluation or the CPU.
* The AVX2 kernel performs the same operations as the portable kernel, 8
  cells at a time. The portable kernel is used on CPUs without AVX2 and on
  ARM. Both are verified to be identical.
* Pairs with a NUL byte are left to parasail, as are pairs with
  `--save-alignments`, which needs the aligned strings.

## 4. What is proved, and what is verified

* **Proved (Section 2):** a certified banded result equals the full-matrix
  result of the specification in Section 1, for all inputs.
* **Verified, not proved:**
  * that the code implements Sections 1 and 3 without bugs;
  * that parasail's SIMD kernels (striped, used up to 0.1.4, and scan, now
    the fallback) implement the scalar specification with identical
    tie-breaking. This is parasail's design contract; cgDist relied on it
    before as well.

The evidence:

| Check | Scope | Differences |
|---|---|---|
| Exhaustive: **every** pair over {A,C,G,T,N} up to length 5, all 3 presets, **every** band width; SIMD, scalar and row-major implementations, parasail scan-16, all against parasail's original `nw_trace_striped_sat` | 15.2 M pairs, 272.6 M band checks | **0** |
| Exhaustive: every pair over {A,C,G,T} up to length 6 | 29.8 M pairs, 620.1 M band checks | **0** |
| Exhaustive: every pair over {A,C,G,T,a} up to length 4 (lower-case scoring vs byte-wise SNP counting) | 0.6 M pairs, 9.1 M band checks | **0** |
| Real allele pairs, L. monocytogenes and S. enterica | 739,554 pairs | **0** |
| Adversarial fuzzing (tandem repeats, homopolymers, N, lower case, large length differences, unrelated sequences), 3 presets | 200,000 cases, 109,030 certified | **0** |
| Full runs with `--verify-alignments 1` (every pair re-checked in production) | both datasets | **0** |
| Distance matrices and cache contents vs cgdist 0.1.4 | both datasets, all modes | identical |

## 5. Runtime safeguard

`--verify-alignments <fraction>` re-checks a deterministic fraction of new
alignments (selected by allele hashes, so the selection is reproducible)
against parasail's original production kernel. The run stops with an error
naming the locus and alleles at the first difference. Use
`--verify-alignments 1` to validate a new schema or machine end to end.

## 6. Reproducing the checks

```bash
cargo test --release                                           # unit tests incl. banded
cargo run --release --example banded_exhaustive -- ACGTN 5     # exhaustive small-scope
cargo run --release --example banded_exhaustive -- ACGT 6
cargo run --release --example banded_exhaustive -- ACGTa 4
cargo run --release --example banded_fuzz -- 200000 1          # adversarial fuzzing
cargo run --release --example align_bench -- pairs.tsv         # real pairs from --save-alignments
cgdist ... --verify-alignments 1                               # production cross-check
```
