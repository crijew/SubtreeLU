#pragma once
#include <memory>

namespace subtree_lu {
// parameter enums
enum {
    // input
    I_PRUNE_DENSE,          // threshold for pruning dense rows
    I_PIVTOL,               // pivoting tolerance in millionths
    I_MAX_SUPERNODE_ROWS,   // maximum number of rows in a supernode
    I_MIN_SUPERNODE_COLS,   // minimum number of columns in a supernode
    I_INIT_SUPERNODE_ROWS,  // initial number of rows in a supernode
    I_MEM_FACTOR,           // memory growth factor for supernode in percent

    // output
    O_ANALYZE_TIME,     // runtime for analysis in microseconds
    O_FACTORIZE_TIME,   // runtime for factorization or refactorization in microseconds
    O_SOLVE_TIME,       // runtime for solving in microseconds
    O_NPIVOTS,          // number of off-diagonal pivots
    O_LNNZ,             // number of non-zeros in L (including diagonal)
    O_UNNZ,             // number of non-zeros in U (excluding diagonal)
    O_NSUPERNODES,      // number of supernodes
    O_NREALLOC,         // number of reallocations
    O_FACTORIZE_FLOPS,  // FLOPs for factorization
    O_SOLVE_FLOPS,      // FLOPs for solving
    O_FACTORIZE_MEM,    // memory for factorization in bytes
    O_SINGULAR_ROW,     // row index of singularity when E_SINGULAR
};

// error enums
enum {
    E_OK = 0,
    E_BAD_ARG = -1,
    E_OOM = -2,
    E_SINGULAR = -3,
    E_NOT_ANALYZED = -4,
    E_NOT_FACTORIZED = -5,
    E_UNKNOWN = -100,
};

template <typename idx_t, typename val_t>
class SubtreeLU final {
public:
    SubtreeLU(const SubtreeLU &) = delete;

    SubtreeLU &operator=(const SubtreeLU &) = delete;

    SubtreeLU(SubtreeLU &&) = delete;

    SubtreeLU &operator=(SubtreeLU &&) = delete;

    SubtreeLU();

    ~SubtreeLU();

    int analyze(idx_t n, const idx_t *ap, const idx_t *ai, const val_t *ax, int nthreads, idx_t *rp, idx_t *cp);

    int factorize(const val_t *ax);

    int refactorize(const val_t *ax);

    int solve(const val_t *b, val_t *x);

    int extract_factors(idx_t *lp, idx_t *li, val_t *lx, idx_t *up, idx_t *ui, val_t *ux, idx_t *rp, idx_t *cp) const;

    long long *parm;

private:
    class Impl;

    std::unique_ptr<Impl> impl;
};
}  // namespace subtree_lu
