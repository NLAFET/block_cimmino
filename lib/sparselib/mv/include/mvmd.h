
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                   MV++ Numerical Matrix/Vector C++ Library                */
/*                             MV++ Version 1.5                              */
/*                                                                           */
/*                                  R. Pozo                                  */
/*               National Institute of Standards and Technology              */
/*                                                                           */
/*                                  NOTICE                                   */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that this permission notice appear in all copies and             */
/* supporting documentation.                                                 */
/*                                                                           */
/* Neither the Institution (National Institute of Standards and Technology)  */
/* nor the author makes any representations about the suitability of this    */
/* software for any purpose.  This software is provided ``as is''without     */
/* expressed or implied warranty.                                            */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

//
//      mvmtp.h  : basic templated numerical matrix class, storage
//                  by columns (Fortran oriented.)
//
//
//


#ifndef _MV_MATRIX_double_H_
#define _MV_MATRIX_double_H_    

#include "mvvd.h"
#include "mvmrf.h"


#include <sstream>

#ifdef MV_MATRIX_BOUNDS_CHECK
#   include <assert.h>
#endif


class MV_ColMat_double
{                                                                      
    private:                                                           
           MV_Vector_double v_;
           int dim0_;   // perferred to using dim_[2]. some compilers
           int dim1_;   // refuse to initalize these in the constructor.
           int lda_;
           int ref_;   // true if this is declared as a reference vector,
                        // i.e. it does not own the memory space, but 
                        // rather it is a view to another vector or array.
    public:                                                            
                                                                       
        /*::::::::::::::::::::::::::*/                                 
        /* Constructors/Destructors */                                 
        /*::::::::::::::::::::::::::*/                                 
                                                                       
            MV_ColMat_double();                             
            MV_ColMat_double( int,  int); 

    // some compilers have difficulty with inlined 'for' statements.
    MV_ColMat_double( int,  int, const double&);   

    // usual copy by value
    // (can't use default parameter lda=m, because m is not a constant...)
    //
    MV_ColMat_double(double*,  int m,  int n);
    MV_ColMat_double(double*,  int m,  int n,  int lda);

    // the "reference" versions
    //
    //
    MV_ColMat_double(MV_ColMat_double &A, MV_Matrix_::ref_type i);
    MV_ColMat_double(double*,  int m,  int n, MV_Matrix_::ref_type i);
    MV_ColMat_double(double*,  int m,  int n,  int lda,
                MV_Matrix_::ref_type i);

    MV_ColMat_double(const MV_ColMat_double&); 
    ~MV_ColMat_double();                              
                                                                       
        /*::::::::::::::::::::::::::::::::*/                           
        /*  Indices and access operations */                           
        /*::::::::::::::::::::::::::::::::*/                           
                                                                       
    inline double&       operator()( int,  int); 
    inline const double& operator()( int,  int) const; 


    MV_ColMat_double operator()(const MV_VecIndex &I, const MV_VecIndex &J) ;
    const MV_ColMat_double operator()(const MV_VecIndex &I, const MV_VecIndex &J) const;
     int            dim(int i) const; 
     int            lda(void) const{ return lda_; }
     int            size(int i) const { return dim(i);} 
    MV_ColMat_double&        newsize( int,  int);
    const MV_Vector_double & data() const { return v_; }
    double *ptr() {double *p = v_.ptr(); return p; }
    void setData(const MV_Vector_double &V);
    int ref() const { return ref_;}
                                                                       
        /*::::::::::::::*/                                             
        /*  Assignment  */                                             
        /*::::::::::::::*/                                             
                                                                       
    MV_ColMat_double & operator=(const MV_ColMat_double&);
    MV_ColMat_double & operator=(const double&);

    MV_ColMat_double operator+(const MV_ColMat_double&);
    MV_ColMat_double &operator+=(const MV_ColMat_double &M);
    MV_ColMat_double operator-(const MV_ColMat_double&);
    MV_ColMat_double operator*(const double);
    MV_ColMat_double operator*(const MV_ColMat_double&);
    MV_ColMat_double transpose();
    double infNorm();

    MV_Vector_double operator()(int);
    void setCol(MV_Vector_double &, int);
    void setCols(MV_ColMat_double &V, int start, int nbcols);
    double squaredSum();

    friend std::ostream& operator<<(std::ostream &s, const MV_ColMat_double &A);

};                                                                     

inline void MV_ColMat_double::setData(const MV_Vector_double &V){
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(V.size() == _v.size());
#endif
    v_ = V;
}

inline MV_ColMat_double MV_ColMat_double::operator*(const double alpha)
{
    double *v = new double[dim0_*dim1_];
    for(int k = 0; k < dim0_*dim1_; k++){
        v[k] = v_[k] * alpha;
    }
    MV_ColMat_double R(v, dim0_, dim1_);
    delete[] v;
    return R;
}

inline MV_ColMat_double &MV_ColMat_double::operator+=(const MV_ColMat_double &M)
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(dim(0) ==  M.dim(0) && dim(1) == M.dim(1));
#endif
    MV_Vector_double u = M.data();
    for(int k = 0; k < v_.size(); k++){
        v_[k] += u[k];
    }
    return *this;
}


inline MV_ColMat_double MV_ColMat_double::operator+(const MV_ColMat_double &M)
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(dim(0) ==  M.dim(0) && dim(1) == M.dim(1));
#endif
    MV_Vector_double v = v_;
    MV_Vector_double u = M.data();
    for(int k = 0; k < v_.size(); k++){
        v[k] += u[k];
    }
    MV_ColMat_double R(dim(0), dim(1));
    R.setData(v);
    return R;
}

inline MV_ColMat_double MV_ColMat_double::operator-(const MV_ColMat_double &M)
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(dim(0) ==  M.dim(0) && dim(1) == M.dim(1));
#endif
    MV_Vector_double v = v_;
    MV_Vector_double u = M.data();
    for(int k = 0; k < v_.size(); k++){
        v[k] -= u[k];
    }
    MV_ColMat_double R(dim(0), dim(1));
    R.setData(v);
    return R;
}

inline double MV_ColMat_double::infNorm(){
    double s = 0;
    for(int i=0; i<v_.size(); i++) 
        if(abs(v_[i]) > s) s = abs(v_[i]);
    return s;
}

inline double MV_ColMat_double::squaredSum(){
    double s = 0;
    for(int i=0; i<v_.size(); i++) s+= v_[i]*v_[i];
    return s;
}

inline void MV_ColMat_double::setCols(MV_ColMat_double &V, int start, int nbcols){
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=start && start<dim(1));
    assert(0<=start+nbcols && start+nbcols<dim(1));
#endif

    v_(MV_VecIndex(start*lda_, (start+nbcols)*lda_ -1)) = V.data();
}

inline void MV_ColMat_double::setCol(MV_Vector_double &V, int j){
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=j && j<dim(1));
#endif

    v_(MV_VecIndex(j*lda_, (j+1)*lda_ -1)) = V(MV_VecIndex());
}

inline MV_Vector_double MV_ColMat_double::operator()(int j){
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=j && j<dim1_);
#endif

    MV_Vector_double v = v_(MV_VecIndex(j*lda_, (j+1)*lda_ -1));
    return v;
}

inline double& MV_ColMat_double::operator()( int i,  int j)
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=i && i<dim(0));
    assert(0<=j && j<dim(1));
#endif
    return v_(j*lda_ + i);      // could use indirect addressing
                                // instead...
}

inline const double& MV_ColMat_double::operator()
                    ( int i,  int j) const
{
#ifdef MV_MATRIX_BOUNDS_CHECK
    assert(0<=i && i<dim(0));
    assert(0<=j && j<dim(1));
#endif
    return v_(j*lda_ + i);
}

inline MV_ColMat_double::MV_ColMat_double(double* d,  int m, 
         int n, MV_Matrix_::ref_type i ):
            v_(d,m*n, MV_Vector_::ref), dim0_(m), dim1_(n), lda_(m), ref_(i) {}

inline MV_ColMat_double::MV_ColMat_double( MV_ColMat_double &A, 
                MV_Matrix_::ref_type i ):
                v_(&A(0,0), A.dim(0)*A.dim(1), MV_Vector_::ref), 
                dim0_(A.dim(0)), dim1_(A.dim(1)), lda_(A.lda()), ref_(i) {}

inline MV_ColMat_double::MV_ColMat_double(double* d,  int m,  int n,
             int lda, MV_Matrix_::ref_type i) :
            v_(d, lda*n, MV_Vector_::ref), dim0_(m), dim1_(n), lda_(lda),
            ref_(i) {}


#endif

