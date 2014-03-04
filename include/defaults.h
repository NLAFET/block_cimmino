#ifndef _DEFAULTS_H_
#define _DEFAULTS_H_
///\xrefitem exp "Experimental" "Experimental List

/// Defines the control parameters indices in a safe way
namespace Controls{
    enum icontrols {
        nbparts         = 1, ///< The number of partitions
        part_type       = 2, ///< The partitioning type
        part_guess      = 4, ///< Guess the number of partitions
        scaling         = 5, ///< The scaling type
        itmax           = 6, ///< The max number of iterations
        block_size      = 7, ///< The Block-CG block-size
        verbose_level   = 8, ///< The verbose level

        aug_type        = 10, ///< The augmentation type
        aug_blocking    = 11, ///< The blocking factor when building S
        aug_analysis    = 12, ///< Analyse the augmentation process

        exploit_sparcity= 13, ///< Exploit the sparcity in MUMPS

        aug_iterative   = 14, ///< \exp Enable or disable iterative solving of Sz=f

        aug_project     = 15, ///< \deprecated Compute the projection only
        aug_dense       = 16, ///< \deprecated Use dense RHS when doing the computation
    };
    enum dcontrols {
        part_imbalance  = 1, ///< The imbalance factor in PaToH case
        threshold       = 2, ///< The stoping threshold

        aug_filter      = 10, ///< \deprecated The filtering value
        aug_precond     = 11, ///< \exp The preconditioner criteria of selection
    };
    enum info {
        status          = 1, ///< Exit status
    };

    enum dinfo {
        residual        = 1, ///< The resulting residual
        forward_error   = 2, ///< The resulting forward error
    };

}

/// Get some goodies from C++11
#define nullptr 0

#endif // _DEFAULTS_H_
