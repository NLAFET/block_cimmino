#include <abcd.h>
#include <Eigen/src/misc/blas.h>

void abcd::bcg()
{
    // s is the block size of the current run
    int s = std::max<int>(block_size, nrhs);
    if(s < 1) throw - 51;

    cout << n << endl;
    //exit(0);
    if(!use_xk) {
        xk = Eigen::MatrixXd(n, nrhs);
        xk.setZero();
    }

    Eigen::MatrixXd u(b.rows(), s);
    Eigen::MatrixXd p(n, s);
    Eigen::MatrixXd qp(n, s);
    Eigen::MatrixXd r(n, s);
    Eigen::MatrixXd gammak(s, s);
    Eigen::MatrixXd betak(s, s);
    Eigen::MatrixXd lambdak(s, nrhs);
    Eigen::MatrixXd prod_gamma(nrhs, nrhs);
    Eigen::MatrixXd e1(s, nrhs);
    double *qp_ptr = qp.data();
    double *betak_ptr = betak.data();
    double *l_ptr;

    double lnrmBs = std::pow(b.norm(), 2);
    double thresh = 1e-12;

    // Sync B norm :
    mpi::all_reduce(inter_comm, &lnrmBs, 1,  &nrmB, std::plus<double>());
    nrmB = sqrt(nrmB);

    u.setZero();
    p.setZero();
    qp.setZero();
    r.setZero();
    gammak.setZero();
    betak.setZero();
    lambdak.setZero();
    e1.setZero();
    e1.row(0).setOnes();

    char up = 'U';
    char left = 'L';
    char right = 'R';
    char tr = 'T';
    char notr = 'N';
    double alpha = 1;

    /***************************************************
     * ITERATION k = 0                                 *
     **************************************************/

    u.col(0) = b;

    if(use_xk) {
        r.block(0, 0, n, 1) = sumProject(1e0, b, -1e0, xk);
    } else {
        r.block(0, 0, n, 1) = sumProject(1e0, b, 0, xk);
    }
    
    if(s > nrhs) {
        Eigen::MatrixXd RR =  MatrixXd::Random(n, s - nrhs);
        //RR.setOnes();
        r.block(0, nrhs, n, s - nrhs) = RR;
    }

    // orthogonalize
    // r = r*gamma^-1
    if(gqr(r, r, gammak, s, false) != 0)
        gmgs(r, r, gammak, s, false);
    p = r;
    prod_gamma = gammak.block(0, 0, nrhs, nrhs);

    int it = 0;
    double rho = 1;
    std::vector<double> grho(inter_comm.size());
    double mrho;

    double ti = MPI_Wtime();

    if(inter_comm.rank() == inter_comm.size() - 1) {
        cout << "ITERATION " << 0 << endl;
    }
    rho = compute_rho(xk, b, thresh);

    while((rho > thresh) && (it < itmax)) {
        if(inter_comm.rank() == inter_comm.size() - 1) {
            // a simple hack to clear the screen on unix systems
            cout << "ITERATION " << it + 1 << endl;
        }
        it++;

        qp = sumProject(0e0, u, 1e0, p);
        if(gqr(p, qp, betak, s, true) != 0)
            gmgs(p, qp, betak, s, true);

        lambdak = e1;
        l_ptr = lambdak.data();
        dtrsm_(&left, &up, &tr, &notr, &s, &nrhs, &alpha, betak_ptr,
               &s, l_ptr, &n);

        lambdak = lambdak * prod_gamma;
        xk.block(0, 0, n, nrhs) += (p * lambdak).block(0, 0, n, nrhs);

        // R = R - QP * B^-T
        dtrsm_(&right, &up, &tr, &notr, &n, &s, &alpha, betak_ptr,
               &s, qp_ptr, &n);
        r = r - qp;

        if(gqr(r, r, gammak, s, false) != 0)
            gmgs(r, r, gammak, s, false);

        prod_gamma = gammak.block(0, 0, nrhs, nrhs).triangularView<Eigen::Upper>() * prod_gamma;
        gammak = gammak.triangularView<Eigen::Upper>();
        betak = betak.triangularView<Eigen::Upper>() * gammak.transpose() ;

        p = r + p * betak;

        rho = abcd::compute_rho(xk, b, thresh);
        normres.push_back(rho);
        //mpi::all_gather(inter_comm, rho, grho);
        //mrho = *std::max_element(grho.begin(), grho.end());
        //cout << "ITERATION " << it
                //<< " rho = " << rho << endl;
    }
    if(inter_comm.rank() == inter_comm.size() - 1) {
        //cout << xk << endl;
        cout << "TIME : " << MPI_Wtime() - ti << endl;
    }
}

double abcd::compute_rho(Eigen::MatrixXd x, Eigen::MatrixXd u, double thresh)
{
    //double nrmX = x.norm();
    int s = x.cols();
    Eigen::MatrixXd R(m, s);
    R.setZero();
    int pos = 0;
    int ci;
    double gnrmx, nrmXfmX;
    std::vector<double> vnrmx(inter_comm.size());

    for(int p = 0; p < parts.size(); p++) {
        Eigen::VectorXd compressed_x(parts[p].cols());
        for(int j = 0; j < x.cols(); j++) {
            compressed_x.setZero();

            int x_pos = 0;
            for(int i = 0; i < local_column_index[p].size(); i++) {
                int ci = local_column_index[p][i];
                compressed_x(x_pos) = x(ci, j);
                x_pos++;
            }
            R.col(j).segment(pos, parts[p].rows()) = parts[p] * compressed_x;
        }

        pos += parts[p].rows();
    }
    R = u - R;
    //
    double nrmR, nrmX, rho;
    abcd::get_nrmres(x, nrmR, nrmX, nrmXfmX);
    rho = nrmR / (nrmMtx*nrmX + nrmB);
    //cout << "X -> " << x.col(0).norm() << endl;
    //cout << "X -> " << nrmX << " " << x.col(0).norm() << endl;
    //cout << "R -> " << nrmR << endl;
    //cout << "B -> " << nrmB << endl;
    //cout << "A -> " << nrmA << endl;
    //cout << "M -> " << nrmMtx << endl;
    //return R.col(0).norm() / (nrmA * x.col(0).norm() + u.col(0).norm());
    //cout<< R.col(0).norm() / (nrmA * x.col(0).norm() + nrmB) << endl;
    if(inter_comm.rank() == inter_comm.size() - 1) {
        cout << "Rho = " << rho << endl;
        if(use_xf) cout << "Forward = " << nrmXfmX/nrmXf << endl << endl;
    }
    return rho;
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

int abcd::gqr(Eigen::MatrixXd &p, Eigen::MatrixXd &ap, Eigen::MatrixXd &r, int s, bool use_a)
{
    Eigen::SparseMatrix<double> g(ap.rows(), ap.rows());
    g.reserve(ap.rows());
    for(int i = 0; i < ap.rows(); i++)
        g.insert(i, i) = 1;
    g.makeCompressed();

    return abcd::gqr(p, ap, r, g, s, use_a);
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
        if(r(k, k) < r(0,0)*1e-16) {
            r(k, k) = 1;
            cout << "PROBLEM IN GMGS : FOUND AN EXT. SMALL ELMENT " <<  r(k, k) << " postion " << k << endl;
        }
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
        if(use_a)
            ap.col(k) = ap.col(k) / r(k, k);

        if(k < s - 1) {
            r.block(k, k + 1, 1, s - (k + 1)) = abcd::ddot(p.col(k), ap.block(0, k + 1, n, s - (k + 1)));

            Eigen::MatrixXd tcol = p.col(k);
            p.block(0, k + 1, n, s - (k + 1))   -=  tcol * r.block(k, k + 1, 1, s - (k + 1));

            if(use_a) {
                tcol = ap.col(k);
                ap.block(0, k + 1, n, s - (k + 1))  -=  tcol * r.block(k, k + 1, 1, s - (k + 1));
            }
        }
    }
}

int abcd::gqr(Eigen::MatrixXd &p, Eigen::MatrixXd &ap, Eigen::MatrixXd &r,
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
                for(int j = 0; j < s; j++) {
                    loc_p(pos, j) = p(i, j);
                    loc_ap(pos, j) = ap(i, j);
                }
                pos++;
            }
        }
    } else {
        for(int i = 0; i < n; i++) {
            if(comm_map(i) == 1) {
                for(int j = 0; j < s; j++) {
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

    dpotrf_(&up, &s, r_ptr, &s, &ierr);

//     For the moment if there is an error, just crash!
    if(ierr != 0){
        cout << "PROBLEM IN GQR"  << endl;
        return ierr;
    }

    dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, p_ptr, &n);
    if(use_a) dtrsm_(&right, &up, &no, &no, &n, &s, &alpha, r_ptr, &s, ap_ptr, &n);

    return 0;
}




