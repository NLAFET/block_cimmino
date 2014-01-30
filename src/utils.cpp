#include <abcd.h>
#include <Eigen/src/misc/blas.h>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

using namespace std;
using namespace boost::lambda;

void abcd::partitioning(std::vector<vector<int> > &parts, std::vector<int> weights, int nb_parts){
    int total_size = std::accumulate(weights.begin(), weights.end(), 0);
    int mean = floor((double)total_size / nb_parts);
    vector<int> sets(nb_parts);
    map<int, vector<int> > pts;

    for(int i = 0; i < nb_parts; i++){
        sets[i] = 0;
        pts[i];
    }
    
    int cur = 0;
    int ls = 0;
    while(cur != weights.size()){
        int sm = ls;
        for(int i = 0; i < nb_parts; i++){
            if(sets[i] < sets[sm]) sm = i;
        }
        sets[sm] = sets[sm] + weights[cur];
        pts[sm].push_back(cur);
        cur++;
        ls = sm;
    }

    for(int i = 0; i < nb_parts; i++){
        parts.push_back(pts[i]);
    }
}

/// Partition weigts
void abcd::partitionWeights(std::vector<int> &parts, std::vector<int> weights, int nb_parts)
{
    int total_size = std::accumulate(weights.begin(), weights.end(), 0);
    int mean = floor((double)total_size / nb_parts);
    int cum = 0;
    int precum = 0;

    if(nb_parts == weights.size()){
        for(int i = 0; i < weights.size(); i++){
            parts.push_back(i); 
        }
        return;
    }

    for(int c = 0; c < weights.size(); c++) {
        precum = cum;
        cum += weights[c];

        if(cum > mean) {
            if((mean - precum) > 1.5*(cum - mean)) {
                parts.push_back(c);
                cum = 0;
            } else {
                parts.push_back(c - 1);
                cum = weights[c];
            }
        }
    }
    // the last partitions were not taken as they are less than the mean
    if(parts.size() == nb_parts) {
        if(parts[nb_parts - 1 ] < weights.size())
            parts[nb_parts - 1 ] = weights.size() - 1;
    } else {
        parts.push_back(weights.size() - 1);
    }
}
///DDOT
double abcd::ddot(VECTOR_double &p, VECTOR_double &ap)
{
    int lm = p.size();
    int rm = ap.size();
    if(lm != rm) throw - 800;

    VECTOR_double loc_p(lm, 0);
    VECTOR_double loc_ap(rm, 0);
    double loc_r, r;

    int pos = 0;
    for(int i = 0; i < lm; i++) {
        if(comm_map[i] == 1) {
            loc_p(pos) = p(i);
            loc_ap(pos) = ap(i);
            pos++;
        }
    }

    // R = P'AP
    for(int i = 0; i < lm; i++){
        loc_r += loc_p(i) * loc_ap(i);
    }

    mpi::all_reduce(inter_comm, loc_r, r, std::plus<double>());

    return r;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  abcd::get_nrmres
 *  Description:  Computes ||X_k|| and ||X_f - X_k||/||X_f||
 * =====================================================================================
 */
void abcd::get_nrmres(MV_ColMat_double &x, MV_ColMat_double &b, double &nrmR, double &nrmX, double &nrmXfmX)
{
    mpi::communicator world;
    int rn = 1;
    int rm = x.dim(0);

    nrmX = 0;
    nrmR = 0;

    VECTOR_double nrmXV(rn, 0);
    VECTOR_double nrmRV(rn, 0);

    MV_ColMat_double loc_x(rm, rn, 0);
    MV_ColMat_double loc_r(m, rn, 0);
    MV_ColMat_double loc_xfmx(rm, rn, 0);

    int pos = 0;
    for(int i = 0; i < rm; i++) {
        if(comm_map[i] == 1) {
            for(int j = 0; j < rn; j++) {
                nrmXV(j) += abs(x(i, j));
            }
            pos++;
        }
    }

    pos = 0;

    double loc_nrmxfmx = 0;

    for(int p = 0; p < nb_local_parts; p++) {
        for(int j = 0; j < block_size; j++) {
            VECTOR_double compressed_x = VECTOR_double((partitions[p].dim(1)), 0);

            int x_pos = 0;
            for(int i = 0; i < local_column_index[p].size(); i++) {
                int ci = local_column_index[p][i];
                compressed_x(x_pos) = x(ci, j);

                x_pos++;
            }
            VECTOR_double vj = loc_r(j);
            vj(MV_VecIndex(pos, pos+partitions[p].dim(0) - 1)) = partitions[p] * compressed_x;
            loc_r.setCol(vj, j);
        }

        pos += partitions[p].dim(0);
    }

    loc_r  = b - loc_r;

    for(int j = 0; j < rn ; j++){
        VECTOR_double loc_r_j = loc_r(j);
        double loc_nrm = infNorm(loc_r_j);

        mpi::all_reduce(inter_comm, &loc_nrm, 1, &nrmR, mpi::maximum<double>());
        mpi::all_reduce(inter_comm, &nrmXV[j], 1, &nrmX, std::plus<double>());

    }

}

/// Compair pairs
bool ip_comp(const dipair &l, const dipair &r)
{
    return l.first > r.first;
}

template <class K, class V>
std::vector<K> get_keys(std::map<K,V> my_map){
    std::vector<K> keys;
    for(typename std::map<K,V>::iterator it = my_map.begin(); it != my_map.end(); it++){
        keys.push_back(it->first);
    }
    return keys;
}

double or_bin(double &a, double &b){
    if(a!=0) return a;
    else if(b!=0) return b;
    else return 0;
}

void setVal(int *lst, int sz, int ival) {
        int i;
        for (i=0;i < sz; i++) lst[i] = ival;
}

vector<int> sort_indexes(const int *v, const int nb_el) {

    typedef std::pair<int,int> pair_type;
    std::vector< std::pair<int,int> > vp;
    vp.reserve(nb_el);

    for(int i = 0; i < nb_el; i++)
        vp.push_back( std::make_pair<int,int>(v[i], i) );


    // sort indexes based on comparing values in v
    sort(vp.begin(), vp.end(), bind(&pair_type::first, _1) < bind(&pair_type::first, _2));

    std::vector<int> idx(vp.size());
    transform(vp.begin(), vp.end(), idx.begin(), bind(&pair_type::second, _1));
    return idx;
}

template <typename T>
vector<int> sort_indexes(const vector<T> &v) {

    typedef std::pair<T,int> pair_type;
    std::vector< std::pair<T,int> > vp;
    vp.reserve(v.size());
    for(int i = 0; i < v.size(); i++)
        vp.push_back( std::make_pair<T,int>(v[i], i) );

    // sort indexes based on comparing values in v
    //sort(vp.begin(), vp.end(), bind(&pair_type::first, _1) < bind(&pair_type::first, _2));

    std::vector<int> idx;
    //transform(vp.begin(), vp.end(), idx.begin(), bind(&pair_type::second, _1));
    return idx;
}
