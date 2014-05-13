!
!  This file is part of MUMPS 5.0beta, built on Wed Oct 24 12:18:13 UTC 2012
!
!
!  Copyright 1991-2012 CERFACS, CNRS, ENS Lyon, INPT(ENSEEIHT)-IRIT,
!  INRIA, and University of Bordeaux. Authors: Emmanuel Agullo, Patrick Amestoy,
!  Maurice Bremond, Alfredo Buttari, Philippe Combes, Aurelia Fevre, Abdou
!  Guermouche, Guillaume Joslin, Jacko Koster, Jean-Yves L'Excellent,
!  Stephane Pralet, Francois-Henry Rouet, Wissam Sid-Lakhdar, Tzvetomila Slavova,
!  Bora Ucar, Clement Weisbecker.
!
!  We are grateful to Caroline Bousquet, Indranil Chowdhury, Christophe
!  Daniel, Iain Duff, Vincent Espirat, Chiara Puglisi, Gregoire Richard,
!  Miroslav Tuma and Christophe Voemel who have been contributing to this
!  project.
!
!  This version of MUMPS is provided to you free of charge. It is
!  released under the CeCILL-C license (compatible with the GNU LGPL licence), 
!   http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html, except for
!  the directory PORD, which is public domain (see PORD/README).
!
!  You can acknowledge (using references [1] and [2]) the contribution of
!  this package in any scientific publication dependent upon the use of
!  the package. Please use reasonable endeavours to notify the authors
!  of the package of this publication.
!
!   [1] P. R. Amestoy, I. S. Duff, J. Koster and  J.-Y. L'Excellent,
!   A fully asynchronous multifrontal solver using distributed dynamic
!   scheduling, SIAM Journal of Matrix Analysis and Applications,
!   Vol 23, No 1, pp 15-41 (2001).
!
!   [2] P. R. Amestoy, A. Guermouche, J.-Y. L'Excellent and
!   S. Pralet, Hybrid scheduling for the parallel solution of linear
!   systems. Parallel Computing Vol 32 (2), pp 136-156 (2006).
!
!  As a counterpart to the access to the source code and rights to copy,
!  modify and redistribute granted by the license, users are provided only
!  with a limited warranty  and the software's author,  the holder of the
!  economic rights,  and the successive licensors  have only  limited
!  liability. 
!
!  In this respect, the user's attention is drawn to the risks associated
!  with loading,  using,  modifying and/or developing or reproducing the
!  software by the user in light of its specific status of free software,
!  that may mean  that it is complicated to manipulate,  and  that  also
!  therefore means  that it is reserved for developers  and  experienced
!  professionals having in-depth computer knowledge. Users are therefore
!  encouraged to load and test the software's suitability as regards their
!  requirements in conditions enabling the security of their systems and/or 
!  data to be ensured and,  more generally, to use and operate it in the 
!  same conditions as regards security. 
!
!  The fact that you are presently reading this means that you have had
!  knowledge of the CeCILL-C license and that you accept its terms.
!
      INCLUDE 'dmumps_root.h'
      TYPE DMUMPS_STRUC
        SEQUENCE
!
! This structure contains all parameters 
! for the interface to the user, plus internal
! information from the solver
!
! *****************
! INPUT PARAMETERS
! *****************
!    -----------------
!    MPI Communicator
!    -----------------
        INTEGER :: COMM
!    ------------------
!    Problem definition
!    ------------------
!    Solver (SYM=0 unsymmetric,SYM=1 symmetric Positive Definite, 
!        SYM=2 general symmetric)
!    Type of parallelism (PAR=1 host working, PAR=0 host not working)
        INTEGER ::  SYM, PAR
        INTEGER ::  JOB 
!    --------------------
!    Order of Input matrix 
!    --------------------
        INTEGER ::  N
!
!    ----------------------------------------
!    Assembled input matrix : User interface
!    ----------------------------------------
        INTEGER :: NZ
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A
        INTEGER, DIMENSION(:), POINTER :: IRN, JCN
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COLSCA, ROWSCA, pad0
!
!       ------------------------------------
!       Case of distributed assembled matrix
!       matrix on entry:
!       ------------------------------------
        INTEGER :: NZ_loc, pad1
        INTEGER, DIMENSION(:), POINTER :: IRN_loc, JCN_loc
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A_loc, pad2
!
!    ----------------------------------------
!    Unassembled input matrix: User interface
!    ----------------------------------------
        INTEGER :: NELT, pad3
        INTEGER, DIMENSION(:), POINTER :: ELTPTR
        INTEGER, DIMENSION(:), POINTER :: ELTVAR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: A_ELT, pad4
!
!    ---------------------------------------------
!    Symmetric permutation : 
!               PERM_IN if given by user (optional)
!    ---------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: PERM_IN
!
!
! ******************
! INPUT/OUTPUT data 
! ******************
!    --------------------------------------------------------
!    RHS / SOL_loc
!    -------------
!       right-hand side and solution
!    -------------------------------------------------------
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS, REDRHS
        DOUBLE PRECISION, DIMENSION(:), POINTER :: RHS_SPARSE
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SOL_loc
        INTEGER, DIMENSION(:), POINTER :: IRHS_SPARSE
        INTEGER, DIMENSION(:), POINTER :: IRHS_PTR
        INTEGER, DIMENSION(:), POINTER :: ISOL_loc
        INTEGER ::  LRHS, NRHS, NZ_RHS, LSOL_loc, LREDRHS
        INTEGER ::  pad5
!    ----------------------------
!    Control parameters,
!    statistics and output data
!    ---------------------------
        INTEGER ::  ICNTL(40)
        INTEGER ::  INFO(40) 
        INTEGER :: INFOG(40)
        DOUBLE PRECISION ::  COST_SUBTREES
        DOUBLE PRECISION ::  CNTL(15)
        DOUBLE PRECISION ::  RINFO(40)
        DOUBLE PRECISION ::  RINFOG(40)
!    ---------------------------------------------------------
!    Permutations computed during analysis:
!       SYM_PERM: Symmetric permutation 
!       UNS_PERM: Column permutations (optionnal)
!    ---------------------------------------------------------
        INTEGER, DIMENSION(:), POINTER :: SYM_PERM, UNS_PERM
! 
!    -----
!    Schur
!    -----
        INTEGER ::  NPROW, NPCOL, MBLOCK, NBLOCK
        INTEGER ::  SCHUR_MLOC, SCHUR_NLOC, SCHUR_LLD
        INTEGER ::  SIZE_SCHUR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: SCHUR_CINTERFACE
        INTEGER, DIMENSION(:), POINTER :: LISTVAR_SCHUR
!    -------------------------------------
!    Case of distributed matrix on entry:
!    DMUMPS potentially provides mapping
!    -------------------------------------
        INTEGER, DIMENSION(:), POINTER :: MAPPING
!    --------------
!    Version number
!    --------------
        CHARACTER(LEN=14) ::  VERSION_NUMBER
!    -----------
!    Out-of-core
!    -----------
        CHARACTER(LEN=255) :: OOC_TMPDIR
        CHARACTER(LEN=63) :: OOC_PREFIX
!    ------------------------------------------
!    To save the matrix in matrix market format
!    ------------------------------------------
        CHARACTER(LEN=255) ::  WRITE_PROBLEM
        CHARACTER(LEN=5) :: pad8
!
!
! **********************
! INTERNAL Working data
! *********************
        INTEGER(8) :: KEEP8(150), MAX_SURF_MASTER
        INTEGER ::  INST_Number
!       For MPI
        INTEGER ::  COMM_NODES, MYID_NODES, COMM_LOAD
        INTEGER ::   MYID, NPROCS, NSLAVES
        INTEGER ::  ASS_IRECV
        INTEGER ::  LBUFR
        INTEGER ::  LBUFR_BYTES
        INTEGER, DIMENSION(:), POINTER :: POIDS
        INTEGER, DIMENSION(:), POINTER ::  BUFR
!       IS is used for the factors + workspace for contrib. blocks
        INTEGER, DIMENSION(:), POINTER :: IS
!       IS1 (maxis1) contains working arrays computed 
!       and used only during analysis
        INTEGER, DIMENSION(:), POINTER :: IS1
!       For analysis/facto/solve phases
        INTEGER ::  MAXIS1, Deficiency
        INTEGER ::  KEEP(500)
!       The following data/arrays are computed during the analysis
!       phase and used during the factorization and solve phases.
        INTEGER ::  LNA
        INTEGER ::  NBSA
        INTEGER,POINTER,DIMENSION(:) :: STEP, NE_STEPS, ND_STEPS
!  Info for pruning tree 
        INTEGER,POINTER,DIMENSION(:) :: Step2node
!  ---------------------
        INTEGER,POINTER,DIMENSION(:) :: FRERE_STEPS, DAD_STEPS
        INTEGER,POINTER,DIMENSION(:) :: FILS, PTRAR, FRTPTR, FRTELT
        INTEGER,POINTER,DIMENSION(:) :: NA, PROCNODE_STEPS
!       The two pointer arrays computed in facto and used by the solve
!          (except the factors) are PTLUST_S and PTRFAC. 
        INTEGER, DIMENSION(:), POINTER :: PTLUST_S
        INTEGER(8), DIMENSION(:), POINTER :: PTRFAC
!       main real working arrays for factorization/solve phases
        DOUBLE PRECISION, DIMENSION(:), POINTER :: S
!       Information on mapping
        INTEGER, DIMENSION(:), POINTER :: PROCNODE
!       Input matrix ready for numerical assembly 
!           -arrowhead format in case of assembled matrix
!           -element format otherwise
        INTEGER, DIMENSION(:), POINTER :: INTARR
        DOUBLE PRECISION, DIMENSION(:), POINTER :: DBLARR
!       Element entry: internal data
        INTEGER :: NELT_loc, LELTVAR, NA_ELT, pad11
        INTEGER, DIMENSION(:), POINTER :: ELTPROC
!       Candidates and node partitionning
        INTEGER, DIMENSION(:,:), POINTER :: CANDIDATES
        INTEGER, DIMENSION(:),   POINTER :: ISTEP_TO_INIV2
        INTEGER, DIMENSION(:),   POINTER :: FUTURE_NIV2
        INTEGER, DIMENSION(:,:), POINTER :: TAB_POS_IN_PERE 
        LOGICAL, DIMENSION(:),   POINTER :: I_AM_CAND
!       For heterogeneous architecture
        INTEGER, DIMENSION(:), POINTER :: MEM_DIST
!       Compressed RHS
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP_ROW
        INTEGER, DIMENSION(:),   POINTER :: POSINRHSCOMP_COL
        LOGICAL  :: POSINRHSCOMP_COL_ALLOC, pad111
        DOUBLE PRECISION, DIMENSION(:),   POINTER :: RHSCOMP
!       Info on the subtrees to be used during factorization
        DOUBLE PRECISION, DIMENSION(:), POINTER :: MEM_SUBTREE
        DOUBLE PRECISION, DIMENSION(:), POINTER :: COST_TRAV
        INTEGER, DIMENSION(:),   POINTER :: MY_ROOT_SBTR
        INTEGER, DIMENSION(:),   POINTER :: MY_FIRST_LEAF
        INTEGER, DIMENSION(:),   POINTER :: MY_NB_LEAF
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST
        INTEGER, DIMENSION(:),   POINTER :: DEPTH_FIRST_SEQ
        INTEGER, DIMENSION(:),   POINTER :: SBTR_ID
        INTEGER, DIMENSION(:),   POINTER :: SCHED_DEP
        INTEGER, DIMENSION(:),   POINTER :: SCHED_GRP
        INTEGER, DIMENSION(:),   POINTER :: SCHED_SBTR
        INTEGER, DIMENSION(:),   POINTER :: CROIX_MANU
        DOUBLE PRECISION, DIMENSION(:), POINTER :: WK_USER
        INTEGER :: NBSA_LOCAL
        INTEGER :: LWK_USER
!    Internal control array
        DOUBLE PRECISION ::  DKEEP(130)
!    For simulating parallel out-of-core stack.
        DOUBLE PRECISION, DIMENSION(:),POINTER :: CB_SON_SIZE, pad12
!    Instance number used/managed by the C/F77 interface
        INTEGER ::  INSTANCE_NUMBER
!    OOC management data that must persist from factorization to solve.
        INTEGER ::  OOC_MAX_NB_NODES_FOR_ZONE
        INTEGER, DIMENSION(:,:),   POINTER :: OOC_INODE_SEQUENCE, pad13
        INTEGER(8),DIMENSION(:,:), POINTER :: OOC_SIZE_OF_BLOCK
        INTEGER(8), DIMENSION(:,:),   POINTER :: OOC_VADDR
        INTEGER,DIMENSION(:), POINTER :: OOC_TOTAL_NB_NODES
        INTEGER,DIMENSION(:), POINTER :: OOC_NB_FILES
        INTEGER :: OOC_NB_FILE_TYPE
        CHARACTER,DIMENSION(:,:), POINTER :: OOC_FILE_NAMES  
        INTEGER,DIMENSION(:), POINTER :: OOC_FILE_NAME_LENGTH
!    Indices of nul pivots
        INTEGER,DIMENSION(:), POINTER :: PIVNUL_LIST
!    Array needed to manage additionnal candidate processor 
        INTEGER, DIMENSION(:,:), POINTER :: SUP_PROC, pad14
!    Lists of nodes where processors work. Built/used in solve phase.
        INTEGER, DIMENSION(:), POINTER :: IPTR_WORKING, WORKING
!    Root structure(internal)
        TYPE (DMUMPS_ROOT_STRUC) :: root
!    Low-rank
        INTEGER, POINTER, DIMENSION(:) :: LRGROUPS
        INTEGER :: NBGRP
!    Multicore
        INTEGER :: LPOOL_AFTER_L0_OMP, LPOOL_BEFORE_L0_OMP
        INTEGER :: L_PHYS_L0_OMP
        INTEGER :: L_VIRT_L0_OMP                                    
        INTEGER(8) :: THREAD_LA
        INTEGER, DIMENSION(:), POINTER :: IPOOL_BEFORE_L0_OMP ! Pool before L0_OMP
        INTEGER, DIMENSION(:), POINTER :: IPOOL_AFTER_L0_OMP ! Pool after L0_OMP
        INTEGER, DIMENSION(:), POINTER :: PHYS_L0_OMP   ! Subtrees
        INTEGER, DIMENSION(:), POINTER :: VIRT_L0_OMP   ! Amalgamated subtrees
        INTEGER, DIMENSION(:), POINTER :: PERM_L0_OMP   ! From heaviest to lowest subtree
        INTEGER, DIMENSION(:), POINTER :: PTR_LEAFS_L0_OMP ! To get leafs in global pool
        INTEGER, DIMENSION(:), POINTER :: L0_OMP_MAPPING   ! Mapping of the subtrees
      END TYPE DMUMPS_STRUC
