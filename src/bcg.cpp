#include "abcd.h"
#include <Eigen/src/misc/blas.h>

void abcd::bcg()
{
    // s is the block size of the current run
    int s = std::max<int>(block_size, nrhs);
    if(s < 1) throw - 51;

    Eigen::MatrixXd P(n, s);
    Eigen::MatrixXd QP(n, s);
}

void abcd::gmgs(Eigen::MatrixXd &p, Eigen::MatrixXd &ap, Eigen::MatrixXd &r, int s, bool use_a)
{
    Eigen::SparseMatrix<double> g(ap.rows(), ap.rows());
    g.reserve(ap.rows());
    for(int i = 0; i < ap.rows(); i++)
        g.insert(i, i) = 1;
    g.makeCompressed();

    abcd::gmgs(p, ap, r, g, s, use_a);
}

void abcd::gqr(Eigen::MatrixXd &p, Eigen::MatrixXd &ap, Eigen::MatrixXd &r, int s, bool use_a)
{
    Eigen::SparseMatrix<double> g(ap.rows(), ap.rows());
    g.reserve(ap.rows());
    for(int i = 0; i < ap.rows(); i++)
        g.insert(i, i) = 1;
    g.makeCompressed();

    abcd::gqr(p, ap, r, g, s, use_a);
}

/** Generalized modifed Gram-Schmidt
 *
 */
void abcd::gmgs(Eigen::MatrixXd &p, Eigen::MatrixXd &ap, Eigen::MatrixXd &r,
                Eigen::SparseMatrix<double> g, int s, bool use_a)
{
    r.setZero();

    // OK!! we have our R here, lets have some fun :)
    for(int k = 0; k < s; k++) {
        r.block(k, k, 1, 1) = abcd::ddot(p.col(k), ap.col(k));

        // if it's negative, make it positive
        if(r(k, k) < 0) {
            r(k, k) = abs(r(k, k));
            cout << "PROBLEM IN GMGS : FOUND A NEGATIVE ELMENT " << r(k, k) << " postion " << k << endl;
        }
        if(r(k, k) == 0) {
            r(k, k) = 1;
            cout << "PROBLEM IN GMGS : FOUND A ZERO ELMENT " <<  r(k, k) << " postion " << k << endl;
        }
        // if all is fine, do an sqrt:
        r(k, k) = sqrt(r(k, k));
        p.col(k)  = p.col(k) / r(k, k);
        ap.col(k) = ap.col(k) / r(k, k);

        if(k < s - 1) {
            r.block(k, k + 1, 1, s - (k + 1)) = abcd::ddot(p.col(k), ap.block(0, k + 1, n, s - (k + 1)));

            Eigen::MatrixXd tcol = p.col(k);
            p.block(0, k + 1, n, s - (k + 1))   -=  tcol * r.block(k, k + 1, 1, s - (k + 1));
            tcol = ap.col(k);
            ap.block(0, k + 1, n, s - (k + 1))  -=  tcol * r.block(k, k + 1, 1, s - (k + 1));
        }
    }
}

void abcd::gqr(Eigen::MatrixXd &p, Eigen::MatrixXd &ap, Eigen::MatrixXd &r,
               Eigen::SparseMatrix<double> g, int s, bool use_a)
{
    Eigen::MatrixXd loc_p(n, s);
    Eigen::MatrixXd loc_ap(n, s);
    Eigen::MatrixXd loc_r(s, s);
    loc_ap.setZero();
    loc_p.setZero();
    loc_r.setZero();

    int pos = 0;
    if(use_a) {
        for(int i = 0; i < n; i++) {
            if(comm_map(i) == 1) {
                for(int j = 0; j < nrhs; j++) {
                    loc_p(pos, j) = p(i, j);
                    loc_ap(pos, j) = ap(i, j);
                }
                pos++;
            }
        }
    } else {
        for(int i = 0; i < n; i++) {
            if(comm_map(i) == 1) {
                for(int j = 0; j < nrhs; j++) {
                    loc_p(pos, j) = p(i, j);
                }
                pos++;
            }
        }
    }

    // R = P'AP
    if(use_a)
        loc_r = loc_p.transpose() * loc_ap;
    else
        loc_r = loc_p.transpose() * loc_p;

    const double *l_r_ptr = loc_r.data();
    double *r_ptr = r.data();
    mpi::all_reduce(inter_comm, l_r_ptr, s * s,  r_ptr, std::plus<double>());

    // P = PR^-1    <=> P^T = R^-T P^T
    int ierr;
    char up = 'U';
    char right = 'R';
    char no = 'N';
    double alpha = 1;
    double *p_ptr = p.data();
    double *ap_ptr = ap.data();
    cout.precision(15);
    cout << scientific << r << endl << endl;

    dpotrf_(&up, &s, r_ptr, &s, &ierr);
    cout << r << endl << endl;

//     For the moment if there is an error, just crash!
    if(ierr != 0) throw - 9010 + ierr;

    dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, p_ptr, &n);
    if(use_a) dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, ap_ptr, &n);

    cout << p.transpose() * p << endl;

}

