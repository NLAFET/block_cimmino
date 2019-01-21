/*!
 * \file spldlt/sumProject.cpp
 * \brief Implementation of the sum of projections with the augmented systems using SpLDLT
 * \author S. Cayrols
 * \version 1.0
 */

#include<abcd.h>
#include <spldlt.h>

/*!
 *  \brief Compute sum of projections Ai^+(alpha*Rhs - beta*Ai*X) with the augmented systems using SpLDLT
 *
 *  Compute sum of projections Ai^+(alpha*Rhs - beta*Ai*X) with the augmented systems using SpLDLT
 *
 *  \param alpha: scaling of the Rhs, if 0 then the Rhs is not taken into account
 *  \param Rhs: vector first part of the RHS
 *  \param beta: scaling of the Ai*X, if 0 then the X is not taken into account
 *  \param X: vector Ai*X second part of the RHS
 *
 */
MV_ColMat_double abcd::sumProjectSpLDLT(double            alpha,
                                        MV_ColMat_double  &Rhs,
                                        double            beta,
                                        MV_ColMat_double  &X)
{
  
  // Check size of X equals size of Rhs, if both Rhs and X provided
  if (alpha!=0 && beta!=0) assert(X.dim(1) == Rhs.dim(1));
  // s is number of columns of RhS or X
  int s = alpha != 0 ? Rhs.dim(1) : X.dim(1);

  // Initialize SpLDLT rhs
//LINFO << "Allocate inner_solver.rhs of size " << inner_solver.n * s <<
//  " = " << inner_solver.n <<" x "<< s;
  inner_solver.rhs = new double[inner_solver.n * s];
  for(int i = 0; i < inner_solver.n * s; i++)
    inner_solver.rhs[i] = 0;
  MV_ColMat_double spldlt_rhs(inner_solver.rhs, inner_solver.n, s, MV_Matrix_::ref);

  // Initialize data for the sum of the projections
  MV_ColMat_double Delta(n, s, 0); // Delta for sum of projections
  int pos     = 0; // current position in the augmented system
  int b_pos   = 0; // current position in the local rhs
  double *xpt = X.ptr(); // pointer to X
  int xlda    = X.lda();
  int dlda    = Delta.lda();

  if(beta != 0 || alpha != 0){
    // Build SpLDLT local rhs for each partition
    for(int k = 0; k < nb_local_parts; k++) {
      CompRow_Mat_double *part = &partitions[k];

      /* Build local */
      // local alpha*Rhs
      MV_ColMat_double r(part->dim(0), s, 0); // rhs
      double *rpt = r.ptr(); // pointer to rhs
      int rlda    = r.lda();
      // local beta*AiX
      MV_ColMat_double compressed_x(part->dim(1), s, 0); // compressed X
      double *cxpt  = compressed_x.ptr(); // pointer to compressed X
      int cxlda     = compressed_x.lda();

      // if X provided, rhs=Ai*X (beware compressed format)
      if(beta != 0){
        // Compress X
        int x_pos = 0;
        for(size_t i = 0; i < local_column_index[k].size(); i++) {
          int ci = local_column_index[k][i];
          for(int j = 0; j < s; j++) {
            cxpt[x_pos + j * cxlda] = xpt[ci + j * xlda];
          }
          x_pos++;
        }
        // rhs=beta*Ai*X
        r = smv(*part, compressed_x) * beta;
      }

      // If X and Rhs provided, rhs=alpha*Rhs+beta*Ai*X
      if(alpha != 0 && beta !=0){
        // Build beta*Rhs
        MV_ColMat_double rr(part->dim(0), s);
        for (int i = b_pos; i < b_pos + part->dim(0); ++i)
          for (int j = 0; j < s; ++j)
            rr(i - b_pos, j) = Rhs(i, j);
        // rhs=alpha*Rhs+beta*Ai*X
        r = r + rr * alpha;
      }

      // If only Rhs provided, rhs=alpha Rhs
      if(alpha != 0 && beta ==0){
        for (int i = b_pos; i < b_pos + part->dim(0); ++i)
          for (int j = 0; j < s; ++j)
            r(i - b_pos, j) = Rhs(i, j) * alpha;
      }

      // rhs to SpLDLT rhs
      int j = 0;
      for(int i = pos + part->dim(1); i < pos + part->dim(1) + part->dim(0); i++) {
        for(int r_p = 0; r_p < s; r_p++){
          inner_solver.rhs[i + r_p * inner_solver.n] = rpt[j + r_p * rlda];
        //LINFO << "inner_solver.rhs["<< i + r_p * inner_solver.n <<
        //  "] = "<< rpt[j + r_p * rlda];
        }
        j++;
      }

      // shift for next partition
      b_pos += part->dim(0);
      pos   += part->dim(1) + part->dim(0);

    }

    int job = 1;
    mpi::broadcast(intra_comm, job, 0);

    inner_solver.nrhs = s;
    inner_solver.lrhs = inner_solver.n;

//  double *given_rhs = NULL;
//  given_rhs = (double*)malloc(s * inner_solver.n * sizeof(double));
//  memcpy(given_rhs, inner_solver.rhs, s * inner_solver.n * sizeof(double));
  //for(int i = 0 ; i < s * inner_solver.n; i++)
  //  given_rhs[i] = inner_solver.rhs[i]

    double t = MPI_Wtime();

    // Run SpLDLT solve: Ai^+(alpha*Rhs + beta*Ai*X)
    //TODO put the solve of SpLDLT here
    spldlt_solve(0, inner_solver.nrhs, inner_solver.rhs, inner_solver.n, 
      inner_solver.akeep, inner_solver.fkeep, &inner_solver.info);

//  spldlt_chkerr(inner_solver.n, inner_solver.ptr, inner_solver.row, 
//      inner_solver.val, inner_solver.nrhs, inner_solver.rhs, given_rhs);
//  free(given_rhs);

    t = MPI_Wtime() - t;

    int x_pos   = 0;
    double *dpt = Delta.ptr(); // pointer to delta
    //TODO remove this test ? why having several local partitions ?
    //if(nb_local_parts > 1)
    //{
    //  //// Compute compressed local sum of projections (delta)
    //  //for(int k = 0; k < nb_local_parts; k++) {
    //  //    for(size_t i = 0; i < local_column_index[k].size(); i++) {
    //  //        //int ci = local_column_index[k][i];
    //  //        int ci = local_column_index[k][i];
    //  //        for(int j = 0; j < s; j++) {
    //  //            //Delta(ci, j) = Delta(ci, j) + mumps_rhs(x_pos, j) ;
    //  //            dpt[ci + j * dlda] += mumps.rhs[x_pos + j * mumps.n];
    //  //        }
    //  //        x_pos++;
    //  //    }
    //  //    x_pos += partitions[k].dim(0);
    //  //    //cout << Delta << endl;
    //  //}

    //} else
    //{
    // If only 1 partition, delta is just the compressed solution
    for(size_t i = 0; i < local_column_index[0].size(); i++) {
      int ci = local_column_index[0][i];
      for(int j = 0; j < s; j++) {
        dpt[ci + j * dlda] = inner_solver.rhs[x_pos + j * inner_solver.n];
//      if(comm.rank() == 0)
//        LINFO << i << " -> " << dpt[ci + j * dlda];
      }
      x_pos++;
    }
//  if(comm.rank() == 0)
//    LINFO << "---------";
    //}
  }

//delete[] inner_solver.rhs;
//return Delta;

// free memory and return in case of no slaves
  if(inter_comm.size() == 1) {
//  LINFO << "Only 1 worker in sumProject";
    return Delta;
  }

//LINFO << "Parallel Sum of Delta";
  /* Parallel sum of deltas */
  std::vector<double *>itcp(icntl[Controls::nbparts]); // delta rows to send
  std::vector<double *>otcp(icntl[Controls::nbparts]); // delta rows to receive

  std::vector<mpi::request> reqs; // send/recv delta requests
  std::vector<mpi::status> sts; // status of the requests

  double t = MPI_Wtime();
  double t1 = t;
  double t2 = 0;
  // For each interconnection, send/receive the corresponding rows of delta
  for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
      it != col_interconnections.end(); ++it) {

    if(it->second.size() == 0) continue;

    // Prepare the data to be sent
    //create it if does not exist
    itcp[it->first] = new double[it->second.size()*s];
    otcp[it->first] = new double[it->second.size()*s];

    // add rows of delta to send for interconnections
    double t3 = MPI_Wtime();
    int kp = 0;
    for(int j = 0; j < s; j++) {
      for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); ++i) {
        itcp[it->first][kp] = Delta(*i,j);
        kp++;
      }
    }
    t2 += MPI_Wtime() - t3;
//  if(comm.rank() == 0)
//    LINFO << "Send data to " << it->first;
    reqs.push_back(inter_comm.irecv(it->first, 31, otcp[it->first], kp));
    reqs.push_back(inter_comm.isend(it->first, 31, itcp[it->first], kp));
  }

  // wait for communication to finish
  mpi::wait_all(reqs.begin(), reqs.end());
  t1 = MPI_Wtime() - t1;

  // For each of the received rows of delta, compute the sum
  for(std::map<int, std::vector<int> >::iterator it = col_interconnections.begin();
      it != col_interconnections.end(); ++it) {

    if(it->second.size() == 0) continue;

    // sum of deltas
    int p = 0;
    for(int j = 0; j < s; j++) {
      for(std::vector<int>::iterator i = it->second.begin(); i != it->second.end(); ++i) {
        Delta(*i, j) += otcp[it->first][p];
        p++;
      }
    }

  }

  // free memory
  for(size_t i = 0; i < itcp.size(); i ++) {
    delete[] itcp[i];
    delete[] otcp[i];
  }

  return Delta;

}               /* -----  end of function abcd::sumProjectSpLDLT  ----- */
