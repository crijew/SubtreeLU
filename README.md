# SubtreeLU -- Parallel Sparse LU Factorization

SubtreeLU is a high-performance parallel sparse LU factorization algorithm for SPICE-like circuit simulation.
It is designed to be used in circuit simulation software,
particularly for solving large sparse linear systems that arise in circuit analysis.

SubtreeLU adopts a novel task scheduling scheme based
on collapsing or partitioning the separator tree derived from the nested dissection ordering of matrix. It leverages the inherent subtree structure of the separator tree.

SubtreeLU ensures the compatibility with partial pivoting by constraining the pivot selection within subdomains, preserving numerical stability without compromising concurrency.

The supernodal method is integrated into all scheduling phases of SubtreeLU, so as to further accelerate the computation.

# Pre-requisites

- C++17 compatible compiler
- GLIBC 2.17 or higher
- GLIBCXX 3.4.19 or higher

# Usage of the demo program

Use `make` to build the demo program. The demo program is a simple example of how to use the SubtreeLU library. Its source code is located at `demo/demo.cpp`. After building the project, the demo program is located at `bin/subtree_lu`.

Here is the help message of the demo program:

```text
demo program for SubtreeLU
Usage:
  ./bin/subtree_lu [OPTION...]

  -m, --mat arg      Coefficient matrix file (Matrix Market)
  -p, --pivtol arg   Pivoting tolerance (default: 1e-3)
  -r, --run arg      Run test (0: fact, 1: refact, 2: solve)
  -t, --threads arg  Number of threads (default: 1)
  -n, --repeats arg  Number of repeats (default: 1)
  -s, --save         Whether to save results
  -h, --help         Print usage
```

A typical command to run the demo program is as follows:

```text
./bin/subtree_lu -m demo/rajat14.mtx -p 1e-3 -r 0 -t 16 -n 10
```

This command will read the coefficient matrix in MatrixMarket format from `demo/rajat14.mtx`, use a pivoting tolerance of 0.001, run the test for factorization with pivoting with 16 threads, and repeat the test 10 times.

# Usage of the library

All matrices passed to the interfaces use CSR (Compressed Sparse Row) format with 0-based indices.

Type aliases:
- `idx_t = int`
- `val_t = double`

Typical workflow:
(1) `analyze(...)` → (2) `factorize(...)` → (3) `solve(...)`

For a new numeric matrix with the same sparsity pattern, call `refactorize(...)` then `solve(...)`.

To export factors, call `extract_factors(...)`.

Note that caller is responsible for allocating sufficient storage for all output arrays.

## Analyze
Signature:
- `int analyze(int n, const idx_t* ap, const idx_t* ai, const val_t* ax, int nthreads, idx_t* rp, idx_t* cp)`

Purpose:
- Performs symbolic analysis, sparsity inspection, ordering, and internal preparation for subsequent factorizations.

Parameters:
- `n`: Matrix dimension (A is `n`-by-`n`).
- `ap`: CSR row pointer array of length `n`+1.
- `ai`: CSR column indices array of length `ap[n]`.
- `ax`: CSR nonzero values array of length `ap[n]`.
- `nthreads`: Number of threads used for analysis and subsequent factorization.
- `rp, cp`: Row/column permutation vectors (length `n`). The $i$-th row/column in the permuted matrix is the `rp[i]`/`cp[i]`-th row/column in the original matrix. They act as both input and output:
  - If both pointers are null, or both `rp[0] == -1` and `cp[0] == -1`, an internal ordering is computed. For any non-null pointer among `rp`/`cp`, the computed ordering is written back to that array.
  - If user-supplied permutations are provided, they are used for reordering the matrix. In this case, the factorization is serial.

Returns:
- 0 on success; non-zero on failure.

## Factorize
Signature:
- `int factorize(const val_t* ax)`

Purpose:
- LU factorization with pivoting.

Parameters:
- `ax`: CSR values array (same sparsity structure and ordering as used in analyze).

Returns:
- 0 on success; non-zero on failure.

## Refactorize
Signature:
- `int refactorize(const val_t* ax)`

Purpose:
- Recompute numeric LU factors for a new set of values with the same CSR structure (`ap, ai` unchanged) and pivoting choices.

Parameters:
- `ax`: CSR values array matching the structure provided to analyze.

Returns:
- 0 on success; non-zero on failure.

## Solve
Signature:
- `int solve(const val_t* b, val_t* x)`

Purpose:
- Solve $Ax=b$ using the most recent LU factors (from factorize/refactorize).

Parameters:
- `b`: Right-hand side vector of length `n`.
- `x`: Output solution vector of length `n`.

Returns:
- 0 on success; non-zero on failure.

## Extract the LU factors
Signature:
- `int extract_factors(idx_t* lp, idx_t* li, val_t* lx, idx_t* up, idx_t* ui, val_t* ux, idx_t* rp, idx_t* cp)`

Purpose:
- Export the computed L and U factors in CSR, along with the row/column permutation vectors.

Parameters:
- `lp`: CSR row pointer for L (length `n+1`).
- `li`: CSR column indices for L (length ≥ `lp[n]`).
- `lx`: CSR values for L (length ≥ `lp[n]`).
- `up`: CSR row pointer for U (length `n+1`).
- `ui`: CSR column indices for U (length ≥ `up[n]`).
- `ux`: CSR values for U (length ≥ `up[n]`).
- `rp`: Row permutation vector (length `n`). Factorization does not change row order, so this is identical to the `rp` determined in `analyze(...)`.
- `cp`: Final column permutation vector after numeric pivoting (length `n`). It may differ from the `cp` produced by `analyze(...)` because pivoting can change the order of columns.


Returns:
- 0 on success; non-zero on failure.

# Publications

[1] Jiawen Cheng, Yibin Zhang, and Wenjian Yu, "SubtreeLU: High-performance parallel sparse LU factorization for circuit simulation," in Proc. International Conference on Computer-Aided Design (ICCAD), Munich, Germany, Oct. 2025.
