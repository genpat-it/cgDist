// parasail_trace.rs - Global alignment with traceback through parasail's C API
//
// parasail-rs's AlignResult::get_traceback_strings (0.7.x to 0.9.1) takes the
// three strings returned by parasail_result_get_traceback with
// CString::from_raw: they are malloc'd by the C library and then freed by
// Rust's allocator (undefined behaviour unless both are the same malloc), and
// the parasail_traceback_t itself is never freed (~30 bytes leaked per call).
//
// NwTracer calls the very same C kernel parasail-rs would (looked up by the
// same name, e.g. "nw_trace_scan_16", with the same parasail-rs scoring
// matrix and gap penalties), copies the strings and releases everything with
// parasail's own functions. Scores, saturation flags and gapped strings are
// therefore identical to parasail-rs's.

use libparasail_sys::{
    parasail_function_t, parasail_lookup_function, parasail_result_free, parasail_result_get_score,
    parasail_result_get_traceback, parasail_result_is_saturated, parasail_traceback_free,
};
use parasail_rs::Matrix;
use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int};

/// Result of one traced alignment. `query`/`reference` are empty when the
/// kernel saturated (retry at a wider solution width).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Traced {
    pub score: i32,
    pub saturated: bool,
    pub query: String,
    pub reference: String,
}

/// A parasail global-alignment kernel with traceback, its scoring matrix and
/// gap penalties.
pub struct NwTracer {
    func: parasail_function_t,
    matrix: Matrix,
    gap_open: i32,
    gap_extend: i32,
}

impl NwTracer {
    /// `name` is a parasail function name such as "nw_trace_scan_16" or
    /// "nw_trace_striped_sat" (what parasail-rs builds from
    /// `global().use_trace().scan().solution_width(16)`, resp.
    /// `global().use_trace()`). None if parasail has no such function.
    pub fn new(name: &str, matrix: Matrix, gap_open: i32, gap_extend: i32) -> Option<Self> {
        let cname = CString::new(name).ok()?;
        let func = unsafe { parasail_lookup_function(cname.as_ptr()) };
        func?;
        Some(Self {
            func,
            matrix,
            gap_open,
            gap_extend,
        })
    }

    /// Align `query` (first sequence) against `reference`. Err when the input
    /// contains a NUL byte, parasail returns no result or traceback, or the
    /// traceback is not UTF-8 (the cases in which parasail-rs errs too).
    pub fn align(&self, query: &[u8], reference: &[u8]) -> Result<Traced, String> {
        let f = self.func.ok_or("parasail function missing")?;
        let q = CString::new(query).map_err(|_| "NUL byte in query".to_string())?;
        let r = CString::new(reference).map_err(|_| "NUL byte in reference".to_string())?;
        let (ql, rl) = (query.len() as c_int, reference.len() as c_int);
        let m = *self.matrix;
        unsafe {
            let res = f(
                q.as_ptr(),
                ql,
                r.as_ptr(),
                rl,
                self.gap_open,
                self.gap_extend,
                m,
            );
            if res.is_null() {
                return Err("parasail returned no result".into());
            }
            let score = parasail_result_get_score(res);
            if parasail_result_is_saturated(res) != 0 {
                parasail_result_free(res);
                return Ok(Traced {
                    score,
                    saturated: true,
                    query: String::new(),
                    reference: String::new(),
                });
            }
            let tb = parasail_result_get_traceback(
                res,
                q.as_ptr(),
                ql,
                r.as_ptr(),
                rl,
                m,
                b'|' as c_char,
                b' ' as c_char,
                b' ' as c_char,
            );
            if tb.is_null() {
                parasail_result_free(res);
                return Err("parasail returned no traceback".into());
            }
            // the strings belong to parasail: copy them, then free everything
            let copy = |p: *const c_char| {
                CStr::from_ptr(p)
                    .to_str()
                    .map(str::to_owned)
                    .map_err(|_| "traceback is not UTF-8".to_string())
            };
            let out = copy((*tb).query).and_then(|a| copy((*tb).ref_).map(|b| (a, b)));
            parasail_traceback_free(tb);
            parasail_result_free(res);
            let (query, reference) = out?;
            Ok(Traced {
                score,
                saturated: false,
                query,
                reference,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use parasail_rs::Aligner;

    /// Same kernel, matrix and gaps through parasail-rs: identical output.
    #[test]
    fn identical_to_parasail_rs() {
        let pairs: [(&[u8], &[u8]); 5] = [
            (b"ACGTACGTACGT", b"ACGTTCGTACGT"),
            (b"ACGTACGTACGTAAAA", b"ACGTACGTACGT"),
            (
                b"ATGAAACGCATTAGCACCACCATTACCACC",
                b"ATGAAACGCATTACCACCACCATTACC",
            ),
            (b"A", b"ACGT"),
            (b"GGGGGGGGGG", b"CCCCCCCCCC"),
        ];
        for (name, scan) in [("nw_trace_scan_16", true), ("nw_trace_striped_sat", false)] {
            let t = NwTracer::new(name, Matrix::create(b"ACGT", 2, -1).unwrap(), 5, 2).unwrap();
            let mut b = Aligner::new();
            b.matrix(Matrix::create(b"ACGT", 2, -1).unwrap())
                .gap_open(5)
                .gap_extend(2)
                .global()
                .use_trace();
            if scan {
                b.scan().solution_width(16);
            }
            let a = b.build();
            for (q, r) in pairs {
                let mine = t.align(q, r).unwrap();
                let res = a.align(Some(q), r).unwrap();
                assert_eq!(mine.score, res.get_score());
                assert_eq!(mine.saturated, res.is_saturated());
                let tb = res.get_traceback_strings(q, r).unwrap();
                assert_eq!(
                    (mine.query.as_str(), mine.reference.as_str()),
                    (tb.query.as_str(), tb.reference.as_str())
                );
            }
        }
        assert!(NwTracer::new(
            "no_such_kernel",
            Matrix::create(b"ACGT", 2, -1).unwrap(),
            5,
            2
        )
        .is_none());
    }
}
