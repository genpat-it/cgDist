// banded.rs - Certified banded global alignment, bit-identical to parasail
//
// Alleles of one locus are near-identical, so the optimal global alignment
// stays close to the main diagonal. This module fills only a diagonal band of
// the dynamic-programming matrix and *proves* for every pair that the result
// equals the full-matrix result, or reports that it could not (the caller then
// uses the full parasail alignment). It never returns an uncertified answer.
//
// Exactness has two parts.
//
// 1. Same recurrences and tie-breaking as parasail's NW trace kernels
//    (nw_trace.c): affine gaps where a gap of length L costs
//    open + (L-1)*extend; H prefers diagonal, then deletion (F), then
//    insertion (E) on ties; E and F prefer extension over opening on ties;
//    the traceback state machine of cigar_template.c / traceback_template.c.
//    The scoring matrix is parasail_matrix_create("ACGT", match, mismatch):
//    case-insensitive ACGT, every other symbol scores 0 against anything.
//
// 2. A band certificate. Any alignment path that visits a cell with diagonal
//    offset d = j - i outside the band [dlo, dhi] needs enough insertions and
//    deletions to reach that offset and come back, which bounds its score
//    from above (`outside_upper_bound`). If the best in-band score S is
//    strictly greater than that bound on both sides, every optimal global
//    path lies inside the band. Then every value the traceback compares
//    (H, E, F along the optimal path and the alternatives that tie with it)
//    is identical in the banded and the full matrix, and every alternative
//    that is strictly worse in the full matrix is still strictly worse in the
//    band (banded values are never larger than full values). So the traceback
//    takes exactly the same decisions, including on ties.

use std::cell::RefCell;

const NEG_INF: i32 = i32::MIN / 2;

// parasail trace flags (parasail.h)
const INS: u8 = 1;
const DEL: u8 = 2;
const DIAG: u8 = 4;
const DIAG_E: u8 = 8;
const INS_E: u8 = 16;
const DIAG_F: u8 = 32;
const DEL_F: u8 = 64;

/// Result of a certified alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BandedStats {
    pub snps: usize,
    pub indel_events: usize,
    pub indel_bases: usize,
    pub score: i32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Scoring {
    pub match_score: i32,
    pub mismatch: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
}

/// Symbol code as in parasail_matrix_create("ACGT", ...): case-insensitive
/// A/C/G/T -> 0..=3, every other byte -> 4. A table avoids branches.
const CODE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

#[inline(always)]
fn code(c: u8) -> u8 {
    CODE[c as usize]
}

/// Lower bound on the total gap penalty of an alignment with `g` gap
/// columns that contains at least one insertion and one deletion run.
fn min_gap_cost(g: i64, s: &Scoring) -> i64 {
    let (o, e) = (s.gap_open as i64, s.gap_extend as i64);
    // k runs cost k*open + (g-k)*extend = k*(open-extend) + g*extend, k >= 2
    if o >= e {
        2 * (o - e) + g * e
    } else {
        // cheaper to split into as many runs as possible (k = g)
        g * o
    }
}

/// Upper bound on the score of any global alignment of lengths n (query,
/// rows) and m (reference, columns) that visits diagonal offset `d`, where
/// `d` lies outside [min(0, m-n), max(0, m-n)].
fn outside_upper_bound(n: i64, m: i64, d: i64, s: &Scoring) -> i64 {
    let smax = s.match_score.max(s.mismatch).max(0) as i64;
    let (aligned_max, gaps_min) = if d > 0 {
        // reach j - i = d: insertions I >= d, deletions D = I - (m - n)
        let i_min = d;
        let d_min = d - (m - n);
        (m - i_min, i_min + d_min)
    } else {
        let e = -d;
        let d_min = e;
        let i_min = e + (m - n);
        (n - d_min, i_min + d_min)
    };
    smax * aligned_max.max(0) - min_gap_cost(gaps_min, s)
}

/// Row-major reference implementation (used by tests and `verify`): same
/// recurrences and traceback as `align_in_band_dp`, written as directly as
/// possible.
#[derive(Default)]
struct Buffers {
    h: Vec<i32>,
    f: Vec<i32>,
    trace: Vec<u8>,
    qc: Vec<u8>,
    rc: Vec<u8>,
}

thread_local! {
    static BUFFERS: RefCell<Buffers> = RefCell::new(Buffers::default());
}

/// Why a band attempt produced no result.
enum Uncertified {
    /// The band certificate does not hold.
    Certificate,
    /// Structural failure (should not happen); the caller falls back.
    Invalid,
}

/// Fill the band [dlo, dhi] and trace back, or explain why not.
fn align_in_band(
    q: &[u8],
    r: &[u8],
    s: &Scoring,
    dlo: i64,
    dhi: i64,
    b: &mut Buffers,
) -> Result<BandedStats, Uncertified> {
    let n = q.len();
    let m = r.len();
    let (open, gap) = (s.gap_open, s.gap_extend);
    let width = (dhi - dlo + 1) as usize;

    b.qc.clear();
    b.qc.extend(q.iter().map(|&c| code(c)));
    b.rc.clear();
    b.rc.extend(r.iter().map(|&c| code(c)));
    let score_of = |a: u8, c: u8| -> i32 {
        if a == 4 || c == 4 {
            0
        } else if a == c {
            s.match_score
        } else {
            s.mismatch
        }
    };

    b.h.clear();
    b.h.resize(m + 1, NEG_INF);
    b.f.clear();
    b.f.resize(m + 1, NEG_INF);
    b.trace.clear();
    b.trace.resize(n * width, 0);

    // first row (i = 0): offsets 0..=m, inside the band up to dhi
    b.h[0] = 0;
    let row0_hi = (dhi.max(0) as usize).min(m);
    for j in 1..=row0_hi {
        b.h[j] = -open - (j as i32 - 1) * gap;
    }

    let h = &mut b.h;
    let f = &mut b.f;
    let trace = &mut b.trace;
    for i in 1..=n {
        let ii = i as i64;
        let jlo = (ii + dlo).max(1);
        let jhi = (ii + dhi).min(m as i64);
        if jlo > jhi {
            return Err(Uncertified::Invalid); // band misses this row (invalid band)
        }
        let (jlo, jhi) = (jlo as usize, jhi as usize);

        // column 0 of row i has offset -i
        let col0 = if -ii >= dlo {
            -open - (i as i32 - 1) * gap
        } else {
            NEG_INF
        };
        let mut nh = h[jlo - 1]; // row i-1, column jlo-1 (becomes NWH)
        let mut wh = if jlo == 1 { col0 } else { NEG_INF };
        let mut e = NEG_INF;
        h[0] = col0;

        let qa = b.qc[i - 1];
        let row = &mut trace[(i - 1) * width..i * width];
        for j in jlo..=jhi {
            let nwh = nh;
            nh = h[j];
            let f_opn = nh - open;
            let f_ext = f[j] - gap;
            let fj = f_opn.max(f_ext);
            f[j] = fj;
            let e_opn = wh - open;
            let e_ext = e - gap;
            e = e_opn.max(e_ext);
            let h_dag = nwh + score_of(qa, b.rc[j - 1]);
            let mut v = h_dag.max(e);
            v = v.max(fj);
            h[j] = v;
            wh = v;

            let mut t = if f_opn > f_ext { DIAG_F } else { DEL_F };
            t |= if e_opn > e_ext { DIAG_E } else { INS_E };
            t |= if v == h_dag {
                DIAG
            } else if v == fj {
                DEL
            } else {
                INS
            };
            row[(j as i64 - ii - dlo) as usize] = t;
        }
        // the next row reads h[jhi + 1] as its top neighbour of the new
        // rightmost cell: it must be outside the band (never written)
    }

    let score = h[m];

    // certificate
    let (nn, mm) = (n as i64, m as i64);
    if dhi < mm {
        let ub = outside_upper_bound(nn, mm, dhi + 1, s);
        if (score as i64) <= ub {
            return Err(Uncertified::Certificate);
        }
    }
    if dlo > -nn {
        let ub = outside_upper_bound(nn, mm, dlo - 1, s);
        if (score as i64) <= ub {
            return Err(Uncertified::Certificate);
        }
    }

    // traceback (0-based i, j as in parasail's cigar/traceback templates)
    let mut i = n as i64 - 1;
    let mut j = m as i64 - 1;
    let mut where_ = DIAG;
    let (mut snps, mut events, mut bases) = (0usize, 0usize, 0usize);
    let mut in_gap = false;
    let gap_col = |events: &mut usize, bases: &mut usize, in_gap: &mut bool| {
        if !*in_gap {
            *events += 1;
            *in_gap = true;
        }
        *bases += 1;
    };
    loop {
        if i < 0 && j < 0 {
            break;
        }
        if i < 0 {
            while j >= 0 {
                gap_col(&mut events, &mut bases, &mut in_gap);
                j -= 1;
            }
            break;
        }
        if j < 0 {
            while i >= 0 {
                gap_col(&mut events, &mut bases, &mut in_gap);
                i -= 1;
            }
            break;
        }
        let d = j - i;
        if d < dlo || d > dhi {
            return Err(Uncertified::Invalid); // cannot happen when certified; refuse rather than guess
        }
        let t = trace[i as usize * width + (d - dlo) as usize];
        match where_ {
            DIAG => {
                if t & DIAG != 0 {
                    in_gap = false;
                    if q[i as usize] != r[j as usize] {
                        snps += 1;
                    }
                    i -= 1;
                    j -= 1;
                } else if t & INS != 0 {
                    where_ = INS;
                } else if t & DEL != 0 {
                    where_ = DEL;
                } else {
                    return Err(Uncertified::Invalid);
                }
            }
            INS => {
                gap_col(&mut events, &mut bases, &mut in_gap);
                j -= 1;
                if t & DIAG_E != 0 {
                    where_ = DIAG;
                } else if t & INS_E != 0 {
                    where_ = INS;
                } else {
                    return Err(Uncertified::Invalid);
                }
            }
            _ => {
                gap_col(&mut events, &mut bases, &mut in_gap);
                i -= 1;
                if t & DIAG_F != 0 {
                    where_ = DIAG;
                } else if t & DEL_F != 0 {
                    where_ = DEL;
                } else {
                    return Err(Uncertified::Invalid);
                }
            }
        }
    }

    Ok(BandedStats {
        snps,
        indel_events: events,
        indel_bases: bases,
        score,
    })
}

/// Cells computed per SIMD step; anti-diagonal segments are padded to a
/// multiple of this so the inner loop always runs on full blocks.
const LANES: usize = 8;

/// Per-thread buffers of the diagonal-parity DP, for cell type `T`.
struct DpBuffers<T> {
    ha: Vec<T>,
    hb: Vec<T>,
    ea: Vec<T>,
    eb: Vec<T>,
    fa: Vec<T>,
    fb: Vec<T>,
    trace: Vec<u8>,
    qrev: Vec<u8>,
    rc: Vec<u8>,
}

impl<T> Default for DpBuffers<T> {
    fn default() -> Self {
        Self {
            ha: Vec::new(),
            hb: Vec::new(),
            ea: Vec::new(),
            eb: Vec::new(),
            fa: Vec::new(),
            fb: Vec::new(),
            trace: Vec::new(),
            qrev: Vec::new(),
            rc: Vec::new(),
        }
    }
}

thread_local! {
    static DP_BUFFERS: RefCell<DpBuffers<i32>> = RefCell::new(DpBuffers::default());
    static DP_BUFFERS16: RefCell<DpBuffers<i16>> = RefCell::new(DpBuffers::default());
}

/// Padding added to sequence code arrays: the widest block.
const MAX_LANES: usize = 16;

/// A cell kernel: the per-cell recurrence over one block run.
///
/// Every kernel performs exactly the operations of `fill_blocks_io` below.
/// The i16 kernel is used only when `fits_i16` proves that no real cell value
/// or intermediate comes near the i16 range, so its saturating arithmetic
/// never saturates on real cells and equals exact integer arithmetic; the
/// "outside the band" value -32000 stays below every real value.
trait Kernel {
    type T: Copy;
    const LANES: usize;
    const NEG: Self::T;
    fn from_i32(v: i32) -> Self::T;
    fn to_i32(v: Self::T) -> i32;
    /// # Safety
    /// The CPU must support the kernel's instruction set, and every slice
    /// must hold at least `blocks * LANES` elements.
    #[allow(clippy::too_many_arguments)]
    unsafe fn fill(
        blocks: usize,
        h_io: &mut [Self::T],
        h_f: &[Self::T],
        f_f: &[Self::T],
        h_e: &[Self::T],
        e_e: &[Self::T],
        qs: &[u8],
        rs: &[u8],
        e_out: &mut [Self::T],
        f_out: &mut [Self::T],
        t_out: &mut [u8],
        s: &Scoring,
    );
}

/// Portable kernel (any CPU), i32 cells.
struct Scalar32;
impl Kernel for Scalar32 {
    type T = i32;
    const LANES: usize = LANES;
    const NEG: i32 = NEG_INF;
    fn from_i32(v: i32) -> i32 {
        v
    }
    fn to_i32(v: i32) -> i32 {
        v
    }
    #[inline(always)]
    unsafe fn fill(
        blocks: usize,
        h_io: &mut [i32],
        h_f: &[i32],
        f_f: &[i32],
        h_e: &[i32],
        e_e: &[i32],
        qs: &[u8],
        rs: &[u8],
        e_out: &mut [i32],
        f_out: &mut [i32],
        t_out: &mut [u8],
        s: &Scoring,
    ) {
        fill_blocks_io(
            blocks, h_io, h_f, f_f, h_e, e_e, qs, rs, e_out, f_out, t_out, s,
        )
    }
}

/// AVX2 kernel, 8 x i32 cells.
#[cfg(target_arch = "x86_64")]
struct Avx2x32;
#[cfg(target_arch = "x86_64")]
impl Kernel for Avx2x32 {
    type T = i32;
    const LANES: usize = LANES;
    const NEG: i32 = NEG_INF;
    fn from_i32(v: i32) -> i32 {
        v
    }
    fn to_i32(v: i32) -> i32 {
        v
    }
    #[inline(always)]
    unsafe fn fill(
        blocks: usize,
        h_io: &mut [i32],
        h_f: &[i32],
        f_f: &[i32],
        h_e: &[i32],
        e_e: &[i32],
        qs: &[u8],
        rs: &[u8],
        e_out: &mut [i32],
        f_out: &mut [i32],
        t_out: &mut [u8],
        s: &Scoring,
    ) {
        fill_blocks_io_avx2(
            blocks, h_io, h_f, f_f, h_e, e_e, qs, rs, e_out, f_out, t_out, s,
        )
    }
}

/// "Outside the band" for i16 cells (see `fits_i16`).
const NEG_INF16: i16 = -32000;

/// True if every real cell value of any in-band DP of these lengths, and
/// every intermediate derived from it, lies strictly inside (-31000, 32000),
/// so the i16 kernel is exact. Upper bound: at most min(n, m) aligned
/// columns of score <= smax. Lower bound: some in-band path reaches every
/// in-band cell (diagonal, then one gap), scoring at least
/// min(mismatch, 0) * min(n, m) - open - extend * (n + m); E/F and their
/// intermediates subtract at most one more open and one more extend.
fn fits_i16(n: usize, m: usize, s: &Scoring) -> bool {
    let (n, m) = (n as i64, m as i64);
    let smax = s.match_score.max(s.mismatch).max(0) as i64;
    let mmin = s.mismatch.min(0) as i64;
    let hi = smax * (n.min(m) + 1);
    // real H >= mmin*min(n,m) - open - extend*(n+m); E/F and their
    // intermediates (E - extend) subtract at most one more open and extend
    let lo = mmin * n.min(m) - 2 * s.gap_open as i64 - s.gap_extend as i64 * (n + m + 1) - smax;
    hi < 30_000 && lo > -30_000
}

/// AVX2 kernel, 16 x i16 cells, saturating arithmetic (exact under `fits_i16`).
#[cfg(target_arch = "x86_64")]
struct Avx2x16;
#[cfg(target_arch = "x86_64")]
impl Kernel for Avx2x16 {
    type T = i16;
    const LANES: usize = 16;
    const NEG: i16 = NEG_INF16;
    fn from_i32(v: i32) -> i16 {
        v as i16
    }
    fn to_i32(v: i16) -> i32 {
        v as i32
    }
    #[inline(always)]
    unsafe fn fill(
        blocks: usize,
        h_io: &mut [i16],
        h_f: &[i16],
        f_f: &[i16],
        h_e: &[i16],
        e_e: &[i16],
        qs: &[u8],
        rs: &[u8],
        e_out: &mut [i16],
        f_out: &mut [i16],
        t_out: &mut [u8],
        s: &Scoring,
    ) {
        use std::arch::x86_64::*;
        let len = blocks * 16;
        debug_assert!(
            h_io.len() >= len
                && h_f.len() >= len
                && f_f.len() >= len
                && h_e.len() >= len
                && e_e.len() >= len
                && qs.len() >= len
                && rs.len() >= len
                && e_out.len() >= len
                && f_out.len() >= len
                && t_out.len() >= len
        );
        let open = _mm256_set1_epi16(s.gap_open as i16);
        let gap = _mm256_set1_epi16(s.gap_extend as i16);
        let ms = _mm256_set1_epi16(s.match_score as i16);
        let mm = _mm256_set1_epi16(s.mismatch as i16);
        let four = _mm256_set1_epi16(4);
        let c_diag_f = _mm256_set1_epi16(DIAG_F as i16);
        let c_del_f = _mm256_set1_epi16(DEL_F as i16);
        let c_diag_e = _mm256_set1_epi16(DIAG_E as i16);
        let c_ins_e = _mm256_set1_epi16(INS_E as i16);
        let c_diag = _mm256_set1_epi16(DIAG as i16);
        let c_del = _mm256_set1_epi16(DEL as i16);
        let c_ins = _mm256_set1_epi16(INS as i16);
        for c in 0..blocks {
            let o = c * 16;
            let ld = |v: &[i16]| _mm256_loadu_si256(v.as_ptr().add(o) as *const __m256i);
            let hd = _mm256_loadu_si256(h_io.as_ptr().add(o) as *const __m256i);
            let hu = ld(h_f);
            let fu = ld(f_f);
            let hl = ld(h_e);
            let el = ld(e_e);
            let qa = _mm256_cvtepu8_epi16(_mm_loadu_si128(qs.as_ptr().add(o) as *const __m128i));
            let rb = _mm256_cvtepu8_epi16(_mm_loadu_si128(rs.as_ptr().add(o) as *const __m128i));

            let unk = _mm256_or_si256(_mm256_cmpeq_epi16(qa, four), _mm256_cmpeq_epi16(rb, four));
            let eq = _mm256_cmpeq_epi16(qa, rb);
            let sc = _mm256_andnot_si256(unk, _mm256_blendv_epi8(mm, ms, eq));

            let f_opn = _mm256_subs_epi16(hu, open);
            let f_ext = _mm256_subs_epi16(fu, gap);
            let f = _mm256_max_epi16(f_opn, f_ext);
            let e_opn = _mm256_subs_epi16(hl, open);
            let e_ext = _mm256_subs_epi16(el, gap);
            let e = _mm256_max_epi16(e_opn, e_ext);
            let h_dag = _mm256_adds_epi16(hd, sc);
            let v = _mm256_max_epi16(_mm256_max_epi16(h_dag, e), f);

            _mm256_storeu_si256(h_io.as_mut_ptr().add(o) as *mut __m256i, v);
            _mm256_storeu_si256(e_out.as_mut_ptr().add(o) as *mut __m256i, e);
            _mm256_storeu_si256(f_out.as_mut_ptr().add(o) as *mut __m256i, f);

            let tf = _mm256_blendv_epi8(c_del_f, c_diag_f, _mm256_cmpgt_epi16(f_opn, f_ext));
            let te = _mm256_blendv_epi8(c_ins_e, c_diag_e, _mm256_cmpgt_epi16(e_opn, e_ext));
            let th = _mm256_blendv_epi8(
                _mm256_blendv_epi8(c_ins, c_del, _mm256_cmpeq_epi16(v, f)),
                c_diag,
                _mm256_cmpeq_epi16(v, h_dag),
            );
            let t = _mm256_or_si256(_mm256_or_si256(tf, te), th);
            // 16 x i16 (values < 128) -> 16 bytes
            let p8 = _mm256_packus_epi16(t, t);
            let lo = _mm256_extract_epi64::<0>(p8) as u64;
            let hi = _mm256_extract_epi64::<2>(p8) as u64;
            let dst = t_out.as_mut_ptr().add(o);
            std::ptr::copy_nonoverlapping(lo.to_le_bytes().as_ptr(), dst, 8);
            std::ptr::copy_nonoverlapping(hi.to_le_bytes().as_ptr(), dst.add(8), 8);
        }
    }
}

/// In-place kernel for the diagonal-parity layout: `h_io` holds H of
/// anti-diagonal k-2 on entry (the diagonal predecessor, same slot) and H of
/// k on exit.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn fill_blocks_io(
    blocks: usize,
    h_io: &mut [i32],
    h_f: &[i32],
    f_f: &[i32],
    h_e: &[i32],
    e_e: &[i32],
    qs: &[u8],
    rs: &[u8],
    e_out: &mut [i32],
    f_out: &mut [i32],
    t_out: &mut [u8],
    s: &Scoring,
) {
    // Wrapping arithmetic matches the AVX2 i32 kernel lane for lane. It only
    // matters for block-padding lanes, whose values are never read; real
    // cells stay far from the i32 limits.
    let (open, gap, ms, mm) = (s.gap_open, s.gap_extend, s.match_score, s.mismatch);
    let len = blocks * LANES;
    let (h_io, h_f, f_f, h_e, e_e) = (
        &mut h_io[..len],
        &h_f[..len],
        &f_f[..len],
        &h_e[..len],
        &e_e[..len],
    );
    let (qs, rs, e_out, f_out, t_out) = (
        &qs[..len],
        &rs[..len],
        &mut e_out[..len],
        &mut f_out[..len],
        &mut t_out[..len],
    );
    for x in 0..len {
        let sc = if qs[x] == 4 || rs[x] == 4 {
            0
        } else if qs[x] == rs[x] {
            ms
        } else {
            mm
        };
        let f_opn = h_f[x].wrapping_sub(open);
        let f_ext = f_f[x].wrapping_sub(gap);
        let f = f_opn.max(f_ext);
        let e_opn = h_e[x].wrapping_sub(open);
        let e_ext = e_e[x].wrapping_sub(gap);
        let e = e_opn.max(e_ext);
        let h_dag = h_io[x].wrapping_add(sc);
        let v = h_dag.max(e).max(f);
        h_io[x] = v;
        e_out[x] = e;
        f_out[x] = f;
        let tf = if f_opn > f_ext { DIAG_F } else { DEL_F };
        let te = if e_opn > e_ext { DIAG_E } else { INS_E };
        let th = if v == h_dag {
            DIAG
        } else if v == f {
            DEL
        } else {
            INS
        };
        t_out[x] = tf | te | th;
    }
}

/// AVX2 version of `fill_blocks_io` (8 x i32 lanes), identical operations.
#[cfg(target_arch = "x86_64")]
#[inline(always)]
#[allow(clippy::too_many_arguments)]
unsafe fn fill_blocks_io_avx2(
    blocks: usize,
    h_io: &mut [i32],
    h_f: &[i32],
    f_f: &[i32],
    h_e: &[i32],
    e_e: &[i32],
    qs: &[u8],
    rs: &[u8],
    e_out: &mut [i32],
    f_out: &mut [i32],
    t_out: &mut [u8],
    s: &Scoring,
) {
    use std::arch::x86_64::*;
    let len = blocks * LANES;
    debug_assert!(
        h_io.len() >= len
            && h_f.len() >= len
            && f_f.len() >= len
            && h_e.len() >= len
            && e_e.len() >= len
            && qs.len() >= len
            && rs.len() >= len
            && e_out.len() >= len
            && f_out.len() >= len
            && t_out.len() >= len
    );
    let open = _mm256_set1_epi32(s.gap_open);
    let gap = _mm256_set1_epi32(s.gap_extend);
    let ms = _mm256_set1_epi32(s.match_score);
    let mm = _mm256_set1_epi32(s.mismatch);
    let four = _mm256_set1_epi32(4);
    let c_diag_f = _mm256_set1_epi32(DIAG_F as i32);
    let c_del_f = _mm256_set1_epi32(DEL_F as i32);
    let c_diag_e = _mm256_set1_epi32(DIAG_E as i32);
    let c_ins_e = _mm256_set1_epi32(INS_E as i32);
    let c_diag = _mm256_set1_epi32(DIAG as i32);
    let c_del = _mm256_set1_epi32(DEL as i32);
    let c_ins = _mm256_set1_epi32(INS as i32);
    for c in 0..blocks {
        let o = c * LANES;
        let ld = |v: &[i32]| _mm256_loadu_si256(v.as_ptr().add(o) as *const __m256i);
        let hd = _mm256_loadu_si256(h_io.as_ptr().add(o) as *const __m256i);
        let hu = ld(h_f);
        let fu = ld(f_f);
        let hl = ld(h_e);
        let el = ld(e_e);
        let qa = _mm256_cvtepu8_epi32(_mm_loadl_epi64(qs.as_ptr().add(o) as *const __m128i));
        let rb = _mm256_cvtepu8_epi32(_mm_loadl_epi64(rs.as_ptr().add(o) as *const __m128i));

        let unk = _mm256_or_si256(_mm256_cmpeq_epi32(qa, four), _mm256_cmpeq_epi32(rb, four));
        let eq = _mm256_cmpeq_epi32(qa, rb);
        let sc = _mm256_andnot_si256(unk, _mm256_blendv_epi8(mm, ms, eq));

        let f_opn = _mm256_sub_epi32(hu, open);
        let f_ext = _mm256_sub_epi32(fu, gap);
        let f = _mm256_max_epi32(f_opn, f_ext);
        let e_opn = _mm256_sub_epi32(hl, open);
        let e_ext = _mm256_sub_epi32(el, gap);
        let e = _mm256_max_epi32(e_opn, e_ext);
        let h_dag = _mm256_add_epi32(hd, sc);
        let v = _mm256_max_epi32(_mm256_max_epi32(h_dag, e), f);

        _mm256_storeu_si256(h_io.as_mut_ptr().add(o) as *mut __m256i, v);
        _mm256_storeu_si256(e_out.as_mut_ptr().add(o) as *mut __m256i, e);
        _mm256_storeu_si256(f_out.as_mut_ptr().add(o) as *mut __m256i, f);

        let tf = _mm256_blendv_epi8(c_del_f, c_diag_f, _mm256_cmpgt_epi32(f_opn, f_ext));
        let te = _mm256_blendv_epi8(c_ins_e, c_diag_e, _mm256_cmpgt_epi32(e_opn, e_ext));
        let th = _mm256_blendv_epi8(
            _mm256_blendv_epi8(c_ins, c_del, _mm256_cmpeq_epi32(v, f)),
            c_diag,
            _mm256_cmpeq_epi32(v, h_dag),
        );
        let t = _mm256_or_si256(_mm256_or_si256(tf, te), th);
        let p16 = _mm256_packs_epi32(t, t);
        let p8 = _mm256_packus_epi16(p16, p16);
        let lo = _mm256_extract_epi32::<0>(p8) as u32;
        let hi = _mm256_extract_epi32::<4>(p8) as u32;
        let dst = t_out.as_mut_ptr().add(o);
        std::ptr::copy_nonoverlapping(lo.to_le_bytes().as_ptr(), dst, 4);
        std::ptr::copy_nonoverlapping(hi.to_le_bytes().as_ptr(), dst.add(4), 4);
    }
}

/// Which kernel the dispatcher may use.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum KernelChoice {
    /// Best available: AVX2 i16 when `fits_i16`, else AVX2 i32, else scalar.
    Auto,
    /// Portable i32 kernel (as on CPUs without AVX2 and on ARM).
    Scalar32,
    /// AVX2 i32 kernel (falls back to scalar without AVX2).
    Avx2x32,
}

fn align_in_band_dp(
    q: &[u8],
    r: &[u8],
    s: &Scoring,
    dlo: i64,
    dhi: i64,
    choice: KernelChoice,
) -> Result<BandedStats, Uncertified> {
    #[cfg(target_arch = "x86_64")]
    {
        if choice != KernelChoice::Scalar32 && std::arch::is_x86_feature_detected!("avx2") {
            if choice == KernelChoice::Auto && fits_i16(q.len(), r.len(), s) {
                // SAFETY: AVX2 checked at runtime.
                return DP_BUFFERS16
                    .with(|c| unsafe { dp_avx2_16(q, r, s, dlo, dhi, &mut c.borrow_mut()) });
            }
            // SAFETY: AVX2 checked at runtime.
            return DP_BUFFERS
                .with(|c| unsafe { dp_avx2_32(q, r, s, dlo, dhi, &mut c.borrow_mut()) });
        }
    }
    let _ = choice;
    DP_BUFFERS.with(|c| align_in_band_dp_impl::<Scalar32>(q, r, s, dlo, dhi, &mut c.borrow_mut()))
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn dp_avx2_32(
    q: &[u8],
    r: &[u8],
    s: &Scoring,
    dlo: i64,
    dhi: i64,
    b: &mut DpBuffers<i32>,
) -> Result<BandedStats, Uncertified> {
    align_in_band_dp_impl::<Avx2x32>(q, r, s, dlo, dhi, b)
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn dp_avx2_16(
    q: &[u8],
    r: &[u8],
    s: &Scoring,
    dlo: i64,
    dhi: i64,
    b: &mut DpBuffers<i16>,
) -> Result<BandedStats, Uncertified> {
    align_in_band_dp_impl::<Avx2x16>(q, r, s, dlo, dhi, b)
}

/// Geometry of one band: slot mapping and the cells of an anti-diagonal.
#[derive(Clone, Copy)]
struct Band {
    n: i64,
    m: i64,
    dlo: i64,
    dhi: i64,
    stride: usize,
}

impl Band {
    #[inline(always)]
    fn slot(&self, d: i64) -> usize {
        ((d - self.dlo + 2) >> 1) as usize
    }
    /// (dmin, dmax) of matrix cells on anti-diagonal k inside the band, with
    /// d of the parity of k (may be empty: dmin > dmax).
    #[inline(always)]
    fn cells(&self, k: i64) -> (i64, i64) {
        let mut dmin = self.dlo.max(-k).max(k - 2 * self.n);
        let mut dmax = self.dhi.min(k).min(2 * self.m - k);
        dmin += (dmin - k) & 1;
        dmax -= (dmax - k) & 1;
        (dmin, dmax)
    }
}

/// One anti-diagonal: interior cells in SIMD blocks, then boundary cells.
/// `own` arrays have the parity of k, `oth` arrays the other parity.
#[allow(clippy::too_many_arguments)]
#[inline(always)]
fn wave_step<K: Kernel>(
    k: i64,
    bd: &Band,
    h_own: &mut [K::T],
    e_own: &mut [K::T],
    f_own: &mut [K::T],
    h_oth: &[K::T],
    e_oth: &[K::T],
    f_oth: &[K::T],
    own_has_right_sentinel: bool,
    trace: &mut [u8],
    qrev: &[u8],
    rc: &[u8],
    s: &Scoring,
) {
    let (dmin, dmax) = bd.cells(k);
    if dmin > dmax {
        return;
    }
    let di_lo = dmin.max(2 - k);
    let di_hi = dmax.min(k - 2);
    if di_lo <= di_hi {
        let ta = bd.slot(di_lo);
        let tb = bd.slot(di_hi);
        let blocks = (tb - ta + 1).div_ceil(K::LANES);
        let len = blocks * K::LANES;
        // u = d - dlo + 2 has the parity of k - dlo on this anti-diagonal
        let (tf, te) = if (k - bd.dlo) & 1 == 0 {
            (ta, ta - 1)
        } else {
            (ta + 1, ta)
        };
        let qoff = (bd.n - (k - di_lo) / 2) as usize; // qrev index of q[i - 1]
        let roff = ((k + di_lo) / 2 - 1) as usize; // rc index of r[j - 1]
        let row = k as usize * bd.stride;
        let (h_io, e_out, f_out) = (
            &mut h_own[ta..ta + len],
            &mut e_own[ta..ta + len],
            &mut f_own[ta..ta + len],
        );
        let (h_f, f_f) = (&h_oth[tf..tf + len], &f_oth[tf..tf + len]);
        let (h_e, e_e) = (&h_oth[te..te + len], &e_oth[te..te + len]);
        let (qs, rs) = (&qrev[qoff..qoff + len], &rc[roff..roff + len]);
        let t_out = &mut trace[row + ta..row + ta + len];
        // SAFETY: K's instruction set is available (the AVX2 kernels are
        // only instantiated inside #[target_feature(enable = "avx2")]
        // functions reached after a runtime check) and every slice has
        // length len = blocks * K::LANES.
        unsafe {
            K::fill(
                blocks, h_io, h_f, f_f, h_e, e_e, qs, rs, e_out, f_out, t_out, s,
            )
        };
        // block padding may have overwritten the right band-edge sentinel
        if own_has_right_sentinel {
            let t = bd.slot(bd.dhi + 1);
            h_own[t] = K::NEG;
            e_own[t] = K::NEG;
            f_own[t] = K::NEG;
        }
    }
    // boundary cells: row 0 (i = 0, j = k, d = k), column 0 (i = k, j = 0, d = -k)
    let bval = if k == 0 {
        0
    } else {
        -s.gap_open - (k as i32 - 1) * s.gap_extend
    };
    if dmax == k {
        let t = bd.slot(k);
        h_own[t] = K::from_i32(bval);
        e_own[t] = K::NEG;
        f_own[t] = K::NEG;
    }
    if dmin == -k {
        let t = bd.slot(-k);
        h_own[t] = K::from_i32(bval);
        e_own[t] = K::NEG;
        f_own[t] = K::NEG;
    }
}

/// Band fill in the diagonal-parity layout, then traceback.
///
/// Cell (i, j) lies on anti-diagonal k = i + j at diagonal offset d = j - i
/// and is stored in the arrays of parity k & 1 at slot t = (d - dlo + 2) >> 1.
/// Its predecessors are: diagonal (k-2, d) = same slot of the same array
/// (updated in place); vertical, for F, (k-1, d+1) and horizontal, for E,
/// (k-1, d-1) in the arrays of the other parity. Band edges are the fixed
/// slots of d = dlo-1 and d = dhi+1, which hold the "outside" value
/// throughout, so no per-anti-diagonal resets are needed. Slots of cells
/// outside the matrix are never read: every predecessor of an interior cell
/// is a matrix cell. Trace flags live at trace[k * stride + slot]; only
/// interior cells of the band are ever read back.
#[inline(always)]
fn align_in_band_dp_impl<K: Kernel>(
    q: &[u8],
    r: &[u8],
    s: &Scoring,
    dlo: i64,
    dhi: i64,
    b: &mut DpBuffers<K::T>,
) -> Result<BandedStats, Uncertified> {
    let n = q.len();
    let m = r.len();
    let (nn, mm) = (n as i64, m as i64);

    let mut bd = Band {
        n: nn,
        m: mm,
        dlo,
        dhi,
        stride: 0,
    };
    let slots = bd.slot(dhi + 1) + 1 + K::LANES + 1;
    bd.stride = slots;
    for v in [
        &mut b.ha, &mut b.hb, &mut b.ea, &mut b.eb, &mut b.fa, &mut b.fb,
    ] {
        v.clear();
        v.resize(slots, K::NEG);
    }
    b.qrev.clear();
    b.qrev.extend(q.iter().rev().map(|&c| code(c)));
    b.qrev.resize(n + MAX_LANES, 4);
    b.rc.clear();
    b.rc.extend(r.iter().map(|&c| code(c)));
    b.rc.resize(m + MAX_LANES, 4);
    let kmax = nn + mm;
    let need = (kmax as usize + 1) * slots;
    if b.trace.len() < need {
        b.trace.resize(need, 0); // never cleared: only written cells are read
    }
    let right_parity_even = (dhi + 1) & 1 == 0;

    let DpBuffers {
        ha,
        hb,
        ea,
        eb,
        fa,
        fb,
        trace,
        qrev,
        rc,
    } = b;
    // a* arrays hold even anti-diagonals, b* arrays odd ones
    let mut k = 0i64;
    while k <= kmax {
        wave_step::<K>(
            k,
            &bd,
            ha,
            ea,
            fa,
            hb,
            eb,
            fb,
            right_parity_even,
            trace,
            qrev,
            rc,
            s,
        );
        if k < kmax {
            wave_step::<K>(
                k + 1,
                &bd,
                hb,
                eb,
                fb,
                ha,
                ea,
                fa,
                !right_parity_even,
                trace,
                qrev,
                rc,
                s,
            );
        }
        k += 2;
    }
    let score = K::to_i32(if kmax & 1 == 0 {
        ha[bd.slot(mm - nn)]
    } else {
        hb[bd.slot(mm - nn)]
    });

    if dhi < mm && (score as i64) <= outside_upper_bound(nn, mm, dhi + 1, s) {
        return Err(Uncertified::Certificate);
    }
    if dlo > -nn && (score as i64) <= outside_upper_bound(nn, mm, dlo - 1, s) {
        return Err(Uncertified::Certificate);
    }

    // 0 <= i < n and 0 <= j < m here, so (i + 1, j + 1) is an interior
    // matrix cell; it was computed iff it lies in the band.
    let flag_at = |i: i64, j: i64| -> Option<u8> {
        let d = j - i;
        if d < dlo || d > dhi {
            return None;
        }
        Some(trace[(i + j + 2) as usize * bd.stride + bd.slot(d)])
    };
    let mut i = nn - 1;
    let mut j = mm - 1;
    let mut where_ = DIAG;
    let (mut snps, mut events, mut bases) = (0usize, 0usize, 0usize);
    let mut in_gap = false;
    loop {
        if i < 0 && j < 0 {
            break;
        }
        if i < 0 || j < 0 {
            let rest = if i < 0 { j + 1 } else { i + 1 } as usize;
            if !in_gap {
                events += 1;
            }
            bases += rest;
            break;
        }
        let Some(t) = flag_at(i, j) else {
            return Err(Uncertified::Invalid);
        };
        match where_ {
            DIAG => {
                if t & DIAG != 0 {
                    in_gap = false;
                    if q[i as usize] != r[j as usize] {
                        snps += 1;
                    }
                    i -= 1;
                    j -= 1;
                } else if t & INS != 0 {
                    where_ = INS;
                } else if t & DEL != 0 {
                    where_ = DEL;
                } else {
                    return Err(Uncertified::Invalid);
                }
            }
            INS => {
                if !in_gap {
                    events += 1;
                    in_gap = true;
                }
                bases += 1;
                j -= 1;
                where_ = if t & DIAG_E != 0 {
                    DIAG
                } else if t & INS_E != 0 {
                    INS
                } else {
                    return Err(Uncertified::Invalid);
                };
            }
            _ => {
                if !in_gap {
                    events += 1;
                    in_gap = true;
                }
                bases += 1;
                i -= 1;
                where_ = if t & DIAG_F != 0 {
                    DIAG
                } else if t & DEL_F != 0 {
                    DEL
                } else {
                    return Err(Uncertified::Invalid);
                };
            }
        }
    }
    Ok(BandedStats {
        snps,
        indel_events: events,
        indel_bases: bases,
        score,
    })
}

/// Score of the best "single gap" global alignment: no internal gaps
/// except one gap of |n - m| columns, placed at the best position (at either
/// end or anywhere inside). It is the score of a valid alignment, hence a
/// lower bound on the optimal score. O(n) with running sums.
fn simple_lower_bound(q: &[u8], r: &[u8], s: &Scoring) -> i64 {
    let (n, m) = (q.len(), r.len());
    // score table over (code(a), code(b)); 4 = any non-ACGT symbol
    let mut table = [[0i32; 5]; 5];
    for (a, row) in table.iter_mut().enumerate().take(4) {
        for (b, v) in row.iter_mut().enumerate().take(4) {
            *v = if a == b { s.match_score } else { s.mismatch };
        }
    }
    let sc = |a: u8, b: u8| -> i64 { table[code(a) as usize][code(b) as usize] as i64 };
    if n == m {
        return q.iter().zip(r).map(|(&a, &b)| sc(a, b)).sum();
    }
    // the shorter sequence is aligned without gaps; the longer one skips
    // `delta` symbols after position p of the shorter one
    let (short, long) = if n < m { (q, r) } else { (r, q) };
    let (k, delta) = (short.len(), long.len() - short.len());
    let gap = s.gap_open as i64 + (delta as i64 - 1) * s.gap_extend as i64;
    // suffix = sum_{i >= p} sc(short[i], long[i + delta]); prefix = sum_{i < p} sc(short[i], long[i])
    let mut suffix: i64 = short
        .iter()
        .zip(&long[delta..])
        .map(|(&a, &b)| sc(a, b))
        .sum();
    let mut prefix: i64 = 0;
    let mut best = suffix;
    for p in 0..k {
        prefix += sc(short[p], long[p]);
        suffix -= sc(short[p], long[p + delta]);
        best = best.max(prefix + suffix);
    }
    best - gap
}

fn band_for(n: i64, m: i64, w: i64) -> (i64, i64) {
    ((0i64.min(m - n) - w).max(-n), (0i64.max(m - n) + w).min(m))
}

/// True if a best in-band score of at least `score` certifies half-width `w`.
fn certifies(n: i64, m: i64, w: i64, score: i64, s: &Scoring) -> bool {
    let (dlo, dhi) = band_for(n, m, w);
    (dhi >= m || outside_upper_bound(n, m, dhi + 1, s) < score)
        && (dlo <= -n || outside_upper_bound(n, m, dlo - 1, s) < score)
}

/// Smallest half-width certified by an in-band score of at least `score`.
/// The bound decreases as the band widens, so a binary search applies.
fn needed_width(n: i64, m: i64, score: i64, s: &Scoring) -> i64 {
    let (mut lo, mut hi) = (0i64, n.max(m));
    while lo < hi {
        let mid = (lo + hi) / 2;
        if certifies(n, m, mid, score, s) {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    lo
}

/// Certified banded alignment. Tries bands of half-width `w` around the
/// diagonals 0..=(m-n), widening while the band stays below `max_fraction` of
/// the full matrix. Returns None when no tried band could be certified; the
/// caller must then compute the full alignment.
pub fn align_certified(q: &[u8], r: &[u8], s: &Scoring, max_fraction: f64) -> Option<BandedStats> {
    if q.is_empty() || r.is_empty() || s.gap_open < 0 || s.gap_extend < 0 {
        return None;
    }
    let (n, m) = (q.len() as i64, r.len() as i64);
    let full = (n * m) as f64;
    {
        // The best single-gap alignment is a valid alignment, so its score
        // bounds the optimal score from below; the band it certifies is
        // therefore guaranteed to certify. On real allele pairs the bound is
        // usually exact, so one attempt with the narrowest provable band.
        let lb = simple_lower_bound(q, r, s);
        let w = needed_width(n, m, lb, s);
        let (dlo, dhi) = band_for(n, m, w);
        if (dhi - dlo + 1) as f64 * n as f64 > full * max_fraction {
            return None;
        }
        align_in_band_dp(q, r, s, dlo, dhi, KernelChoice::Auto).ok()
    }
}

/// Hooks for exhaustive verification (examples/banded_exhaustive.rs); not a
/// stable API.
#[doc(hidden)]
pub mod verify {
    use super::*;

    /// Result of each implementation for band half-width `w`: the default
    /// dispatch (AVX2 i16 when it fits), AVX2 i32, the portable scalar i32
    /// kernel, and the row-major reference. None = not certified.
    pub fn band_results(q: &[u8], r: &[u8], s: &Scoring, w: i64) -> [Option<BandedStats>; 4] {
        let (dlo, dhi) = band_for(q.len() as i64, r.len() as i64, w);
        let auto = align_in_band_dp(q, r, s, dlo, dhi, KernelChoice::Auto).ok();
        let avx32 = align_in_band_dp(q, r, s, dlo, dhi, KernelChoice::Avx2x32).ok();
        let scalar = align_in_band_dp(q, r, s, dlo, dhi, KernelChoice::Scalar32).ok();
        let row = BUFFERS.with(|c| align_in_band(q, r, s, dlo, dhi, &mut c.borrow_mut()).ok());
        [auto, avx32, scalar, row]
    }

    pub fn lower_bound(q: &[u8], r: &[u8], s: &Scoring) -> i64 {
        simple_lower_bound(q, r, s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::alignment::compute_alignment_stats;
    use parasail_rs::{Aligner, Matrix};

    const DNA: Scoring = Scoring {
        match_score: 2,
        mismatch: -1,
        gap_open: 5,
        gap_extend: 2,
    };
    const PRESETS: [Scoring; 3] = [
        DNA,
        Scoring {
            match_score: 3,
            mismatch: -2,
            gap_open: 8,
            gap_extend: 3,
        },
        Scoring {
            match_score: 1,
            mismatch: 0,
            gap_open: 3,
            gap_extend: 1,
        },
    ];

    fn parasail(q: &[u8], r: &[u8], s: &Scoring) -> BandedStats {
        let m = Matrix::create(b"ACGT", s.match_score, s.mismatch).unwrap();
        let a = Aligner::new()
            .matrix(m)
            .gap_open(s.gap_open)
            .gap_extend(s.gap_extend)
            .global()
            .use_trace()
            .build();
        let res = a.align(Some(q), r).unwrap();
        let tb = res.get_traceback_strings(q, r).unwrap();
        let (snps, indel_events, indel_bases) = compute_alignment_stats(&tb.query, &tb.reference);
        BandedStats {
            snps,
            indel_events,
            indel_bases,
            score: res.get_score(),
        }
    }

    struct Rng(u64);
    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n.max(1)
        }
    }

    fn random_pair(rng: &mut Rng) -> (Vec<u8>, Vec<u8>) {
        let alphabet: &[u8] = if rng.below(8) == 0 {
            b"ACGTNacgt"
        } else {
            b"ACGT"
        };
        let max_len = if rng.below(4) == 0 { 12 } else { 400 };
        let len = 1 + rng.below(max_len) as usize;
        let mut a = Vec::with_capacity(len);
        while a.len() < len {
            if rng.below(3) == 0 {
                let unit_len = 1 + rng.below(4);
                let unit: Vec<u8> = (0..unit_len)
                    .map(|_| alphabet[rng.below(alphabet.len() as u64) as usize])
                    .collect();
                let reps = 2 + rng.below(6);
                for _ in 0..reps {
                    a.extend_from_slice(&unit);
                }
            } else {
                a.push(alphabet[rng.below(alphabet.len() as u64) as usize]);
            }
        }
        a.truncate(len);
        let rate = [5u64, 20, 60, 200][rng.below(4) as usize];
        let mut b = Vec::new();
        let mut i = 0;
        while i < a.len() {
            if rng.below(1000) < rate {
                match rng.below(3) {
                    0 => b.push(alphabet[rng.below(alphabet.len() as u64) as usize]),
                    1 => {
                        b.push(a[i]);
                        let ins = 1 + rng.below(5);
                        for _ in 0..ins {
                            b.push(alphabet[rng.below(4) as usize]);
                        }
                    }
                    _ => i += rng.below(5) as usize,
                }
            } else {
                b.push(a[i]);
            }
            i += 1;
        }
        if b.is_empty() {
            b.push(b'A');
        }
        (a, b)
    }

    #[test]
    fn simd_layout_matches_row_reference_on_every_band() {
        let mut rng = Rng(11);
        let mut row = Buffers::default();
        for _ in 0..3000 {
            let (q, r) = random_pair(&mut rng);
            let s = PRESETS[rng.below(3) as usize];
            let (n, m) = (q.len() as i64, r.len() as i64);
            for w in [0, 1, 2, 5, 17] {
                let (dlo, dhi) = band_for(n, m, w);
                let a = align_in_band_dp(&q, &r, &s, dlo, dhi, KernelChoice::Auto);
                let b = align_in_band(&q, &r, &s, dlo, dhi, &mut row);
                match (a, b) {
                    (Ok(x), Ok(y)) => assert_eq!(x, y),
                    (Err(Uncertified::Certificate), Err(Uncertified::Certificate)) => {}
                    _ => panic!("layouts disagree on certification (w={w}, lens {n}/{m})"),
                }
            }
        }
    }

    #[test]
    fn scalar_kernel_equals_simd_kernel() {
        let mut rng = Rng(7);
        for _ in 0..1500 {
            let (q, r) = random_pair(&mut rng);
            let s = PRESETS[rng.below(3) as usize];
            for w in [0, 3, 9, 40] {
                let res = verify::band_results(&q, &r, &s, w);
                assert_eq!(res[0], res[1], "auto (i16) vs AVX2 i32, w={w}");
                assert_eq!(res[0], res[2], "auto vs scalar, w={w}");
                assert_eq!(res[0], res[3], "auto vs row reference, w={w}");
            }
        }
    }

    #[test]
    fn certified_results_equal_parasail() {
        let mut rng = Rng(12345);
        let mut certified = 0;
        for _ in 0..4000 {
            let (q, r) = random_pair(&mut rng);
            let s = PRESETS[rng.below(3) as usize];
            if let Some(got) = align_certified(&q, &r, &s, 1.0) {
                certified += 1;
                assert_eq!(
                    got,
                    parasail(&q, &r, &s),
                    "q={} r={}",
                    String::from_utf8_lossy(&q),
                    String::from_utf8_lossy(&r)
                );
            }
        }
        assert!(certified > 2000, "too few certified cases: {certified}");
    }

    #[test]
    fn edge_cases() {
        let cases: &[(&[u8], &[u8])] = &[
            (b"A", b"A"),
            (b"A", b"C"),
            (b"A", b"ACGTACGT"),
            (b"ACGTACGT", b"T"),
            (b"NNNN", b"NNNN"),
            (b"ACGTNNNNACGT", b"ACGTACGT"),
            (b"acgtacgtac", b"ACGTACGTAC"),
            (b"AAAAAAAAAA", b"AAAAAAA"),
            (b"ACACACACACACAC", b"ACACACACAC"),
            (b"ATATATATGCGCGCGC", b"ATATATGCGCGC"),
        ];
        for (q, r) in cases {
            for s in PRESETS {
                if let Some(got) = align_certified(q, r, &s, 1.0) {
                    assert_eq!(
                        got,
                        parasail(q, r, &s),
                        "q={} r={}",
                        String::from_utf8_lossy(q),
                        String::from_utf8_lossy(r)
                    );
                }
            }
        }
        // identical sequences always certify with a zero-width band
        let q = b"ACGTTGCAACGTTGCAAC";
        assert_eq!(align_certified(q, q, &DNA, 0.5), Some(parasail(q, q, &DNA)));
        assert!(align_certified(b"", b"A", &DNA, 1.0).is_none());
    }

    #[test]
    fn i16_guard_bounds() {
        assert!(fits_i16(1200, 1250, &DNA));
        assert!(fits_i16(5000, 5000, &DNA));
        assert!(!fits_i16(12000, 12000, &DNA));
        // every preset fits typical gene lengths
        for s in PRESETS {
            assert!(fits_i16(3000, 3000, &s));
        }
    }

    #[test]
    fn lower_bound_is_a_valid_alignment_score() {
        let mut rng = Rng(99);
        for _ in 0..2000 {
            let (q, r) = random_pair(&mut rng);
            let s = PRESETS[rng.below(3) as usize];
            assert!(simple_lower_bound(&q, &r, &s) <= parasail(&q, &r, &s).score as i64);
        }
    }
}
