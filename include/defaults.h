// Copyright Institut National Polytechnique de Toulouse (2014) 
// Contributor(s) :
// M. Zenadi <mzenadi@enseeiht.fr>
// D. Ruiz <ruiz@enseeiht.fr>
// R. Guivarch <guivarch@enseeiht.fr>

// This software is governed by the CeCILL-C license under French law and
// abiding by the rules of distribution of free software.  You can  use, 
// modify and/ or redistribute the software under the terms of the CeCILL-C
// license as circulated by CEA, CNRS and INRIA at the following URL
// "http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html"

// As a counterpart to the access to the source code and  rights to copy,
// modify and redistribute granted by the license, users are provided only
// with a limited warranty  and the software's author,  the holder of the
// economic rights,  and the successive licensors  have only  limited
// liability. 

// In this respect, the user's attention is drawn to the risks associated
// with loading,  using,  modifying and/or developing or reproducing the
// software by the user in light of its specific status of free software,
// that may mean  that it is complicated to manipulate,  and  that  also
// therefore means  that it is reserved for developers  and  experienced
// professionals having in-depth computer knowledge. Users are therefore
// encouraged to load and test the software's suitability as regards their
// requirements in conditions enabling the security of their systems and/or 
// data to be ensured and,  more generally, to use and operate it in the 
// same conditions as regards security. 

// The fact that you are presently reading this means that you have had
// knowledge of the CeCILL-C license and that you accept its terms.

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
        use_gmgs2       = 9, ///< Force GMGS^2 in B-CG

        aug_type        = 10, ///< The augmentation type
        aug_blocking    = 11, ///< The blocking factor when building S
        aug_analysis    = 12, ///< Analyse the augmentation process

        exploit_sparcity= 13, ///< Exploit the sparcity in MUMPS

#ifdef WIP
        aug_iterative   = 14, ///< \exp Enable or disable iterative solving of Sz=f

        aug_project     = 15, ///< \deprecated Compute the projection only
        aug_dense       = 16, ///< \deprecated Use dense RHS when doing the computation
#endif //WIP
    };
    enum dcontrols {
        part_imbalance  = 1, ///< The imbalance factor in PaToH case
        threshold       = 2, ///< The stoping threshold

#ifdef WIP
        aug_filter      = 10, ///< \deprecated The filtering value
        aug_precond     = 11, ///< \exp The preconditioner criteria of selection
#endif //WIP
    };
    enum info {
        status          = 1, ///< Exit status
        nb_iter         = 2, ///< Number of iterations after CG
    };

    enum dinfo {
        residual        = 1, ///< The resulting residual
        forward_error   = 2, ///< The resulting forward error
        backward        = 3, ///< The resulting residual
        scaled_residual = 4, ///< The resulting residual
    };

}

/// Get some goodies from C++11
#define nullptr 0

#endif // _DEFAULTS_H_
