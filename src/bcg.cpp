#include "abcd.h"

void abcd::bcg()
{
    // s is the block size of the current run
    int s = std::max<int>(block_size, nrhs);
    if(s < 1) throw - 51;
    
    Eigen::MatrixXd P(n, s);
    Eigen::MatrixXd QP(n, s);
}
