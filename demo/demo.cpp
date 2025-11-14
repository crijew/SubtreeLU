#include <algorithm>
#include <fstream>
#include <cxxopts.h>
#include <fast_matrix_market/fast_matrix_market.hpp>
#include "subtree_lu.h"

using namespace std;
using namespace subtree_lu;
using idx_t = int;
using val_t = double;
using VecI = vector<idx_t>;
using VecV = vector<val_t>;

enum RunType { FACT, REFACT, SOLVE };

VecV spmv(const VecI &ap, const VecI &ai, const VecV &ax, const VecV &x) {
    auto n = ap.size() - 1;
    VecV y(n);
    for (idx_t i = 0; i < n; i++) {
        val_t temp{};
        for (idx_t m = ap[i]; m < ap[i + 1]; m++) {
            temp += ax[m] * x[ai[m]];
        }
        y[i] = temp;
    }
    return y;
}

double inf_norm(const VecV &x, const VecV &y) {
    double norm = 0;
    for (idx_t i = 0; i < x.size(); i++) {
        norm = max(norm, static_cast<double>(abs(x[i] - y[i])));
    }
    return norm;
}

int main(int argc, char **argv) {
    // parse arguments
    cxxopts::Options options(argv[0], "demo program for SubtreeLU");
    // clang-format off
    options.add_options()
        ("m,mat", "Coefficient matrix file (Matrix Market)", cxxopts::value<string>())
        ("p,pivtol", "Pivoting tolerance", cxxopts::value<double>()->default_value("1e-3"))
        ("r,run", "Run test (0: fact, 1: refact, 2: solve)", cxxopts::value<int>())
        ("t,threads", "Number of threads", cxxopts::value<int>()->default_value("1"))
        ("n,repeats", "Number of repeats", cxxopts::value<int>()->default_value("1"))
        ("s,save", "Whether to save results", cxxopts::value<bool>()->default_value("false"))
        ("h,help", "Print usage");
    // clang-format on
    auto args = options.parse(argc, argv);
    if (args.count("help")) {
        cout << options.help() << endl;
        return 0;
    }
    auto mat_path = args["mat"].as<string>();
    double pivtol = clamp(args["pivtol"].as<double>(), 0.0, 1.0);
    int run_type = args["run"].as<int>();
    if (run_type < FACT || run_type > SOLVE) {
        cerr << "Invalid run type" << endl;
        return 1;
    }
    int nthreads = max(args["threads"].as<int>(), 1);
    int repeats = max(args["repeats"].as<int>(), 1);
    bool save = args["save"].as<bool>();

    // read matrix
    idx_t n, ncols;
    VecI rows, cols;
    VecV values;
    ifstream fin(mat_path);
    fast_matrix_market::read_options read_options;
    read_options.num_threads = 8;
    read_matrix_market_triplet(fin, n, ncols, rows, cols, values, read_options);
    fin.close();
    auto nnz = rows.size();
    cout << "matrix loaded, m = " << n << ", n = " << ncols << ", nnz = " << nnz << endl;
    if (n != ncols) {
        cerr << "Matrix must be square" << endl;
        return 1;
    }

    // convert COO to CSR
    VecI count(n, 0);
    for (idx_t m = 0; m < nnz; m++) {
        count[rows[m]]++;
    }
    VecI ap(n + 1), ai;
    VecV ax;
    ai.reserve(nnz);
    ax.reserve(nnz);
    for (idx_t i = 0; i < n; i++) {
        ap[i + 1] = ap[i] + count[i];
    }
    fill_n(count.begin(), n, 0);
    for (idx_t m = 0; m < nnz; m++) {
        idx_t i = rows[m];
        ai[ap[i] + count[i]] = cols[m];
        ax[ap[i] + count[i]++] = values[m];
    }
    VecI{}.swap(rows);
    VecI{}.swap(cols);
    VecV{}.swap(values);

    // generate rhs
    VecV x(n, val_t{1}), b = spmv(ap, ai, ax, x), sol(n);

    // stats
    double factorize_time = 0, refactorize_time = 0, solve_time = 0;
    int factorize_runs = 0, refactorize_runs = 0, solve_runs = 0;

    // create solver
    SubtreeLU<idx_t, val_t> solver;
    solver.parm[I_PIVTOL] = static_cast<long long>(pivtol * 1e6);

    // analyze
    VecI rp(n), cp(n);
    rp[0] = -1;
    cp[0] = -1;
    int err = solver.analyze(n, &ap[0], &ai[0], &ax[0], nthreads, &rp[0], &cp[0]);
    if (err != E_OK) {
        goto FINAL;
    }
    if (save) {
        fast_matrix_market::matrix_market_header header(n, 1);
        ofstream fout("rp.txt");
        write_matrix_market_array(fout, header, rp);
        fout.close();
        fout.open("cp.txt");
        write_matrix_market_array(fout, header, cp);
        fout.close();
    }

    // run task
    if (run_type > FACT) {
        err = solver.factorize(&ax[0]);
        factorize_time += static_cast<double>(solver.parm[O_FACTORIZE_TIME]) * 1e-6;
        factorize_runs++;
    }
    for (int i = 0; i < repeats; i++) {
        if (err != E_OK) {
            break;
        }
        switch (run_type) {
            case FACT:
                err = solver.factorize(&ax[0]);
                factorize_time += static_cast<double>(solver.parm[O_FACTORIZE_TIME]) * 1e-6;
                factorize_runs++;
                break;
            case REFACT:
                err = solver.refactorize(&ax[0]);
                refactorize_time += static_cast<double>(solver.parm[O_FACTORIZE_TIME]) * 1e-6;
                refactorize_runs++;
                break;
            default:
                err = solver.solve(&b[0], &sol[0]);
                solve_time += static_cast<double>(solver.parm[O_SOLVE_TIME]) * 1e-6;
                solve_runs++;
                break;
        }
    }
    if (err != E_OK) {
        goto FINAL;
    }
    solver.solve(&b[0], &sol[0]);

    if (save) {
        VecI lp(n + 1), up(n + 1), li(solver.parm[O_LNNZ]), ui(solver.parm[O_UNNZ]);
        VecV lx(solver.parm[O_LNNZ]), ux(solver.parm[O_UNNZ]);
        solver.extract_factors(&lp[0], &li[0], &lx[0], &up[0], &ui[0], &ux[0], &rp[0], &cp[0]);
        fast_matrix_market::matrix_market_header header(n, 1);
        ofstream fout("cp1.txt");
        write_matrix_market_array(fout, header, cp);
        fout.close();
        header.ncols = n;
        fout.open("L.mtx");
        fast_matrix_market::write_options write_options;
        write_options.num_threads = 8;
        write_options.precision = 10;
        write_matrix_market_csc(fout, header, lp, li, lx, true, write_options);
        fout.close();
        fout.open("U.mtx");
        write_matrix_market_csc(fout, header, up, ui, ux, true, write_options);
        fout.close();
    }

    // print summary
FINAL:
    printf("\n=== Input summary ===\n");
    printf("run type:                   %s\n", array{"fact", "refact", "solve"}[run_type]);
    printf("matrix:                     %s\n", mat_path.c_str());
    printf("n:                          %d\n", n);
    printf("nnz:                        %lu\n", nnz);
    printf("threads:                    %d\n", nthreads);
    printf("pivtol:                     %e\n", pivtol);
    printf("repeats:                    %d\n", repeats);
    printf("\n=== Output summary ===\n");
    printf("error code:                 %d\n", err);
    if (err == E_SINGULAR) {
        printf("singular row:               %lld\n", solver.parm[O_SINGULAR_ROW]);
    }
    printf("number of off-diag pivots:  %lld\n", solver.parm[O_NPIVOTS]);
    printf("nnz (LU):                   %lld\n", solver.parm[O_LNNZ] + solver.parm[O_UNNZ]);
    printf("number of supernodes:       %lld\n", solver.parm[O_NSUPERNODES]);
    printf("factorize FLOPs:            %lld\n", solver.parm[O_FACTORIZE_FLOPS]);
    printf("solve FLOPs:                %lld\n", solver.parm[O_SOLVE_FLOPS]);
    printf("factorize memory:           %lld\n", solver.parm[O_FACTORIZE_MEM]);
    printf("inf-norm of error:          %e\n", inf_norm(sol, x));
    printf("time\n");
    printf("  analyze:                  %f\n", static_cast<double>(solver.parm[O_ANALYZE_TIME]) * 1e-6);
    printf("  factorize:                %f\n", factorize_time / factorize_runs);
    if (run_type == REFACT) {
        printf("  refactorize:              %f\n", refactorize_time / refactorize_runs);
    } else if (run_type == SOLVE) {
        printf("  solve:                    %f\n", solve_time / solve_runs);
    }
    return err;
}