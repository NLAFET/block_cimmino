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

/*!
 * \file logger.cpp
 * \brief Implementation of the configuration of the easyloggingpp logger
 * \author R. Guivarch, P. Leleux, D. Ruiz, S. Torun, M. Zenadi
 * \version 1.0
 */

#include "abcd.h"

using namespace easyloggingpp;

/*!
 *  \brief Configure logger
 *
 *  Configure the logger by setting format of the lines, and parameters
 *
 *  \param log_output: log filename
 *
 */
void configure_logger(std::string log_output)
{
    Configurations abcd_conf;

    abcd_conf.setToDefault();
    abcd_conf.set(Level::Info,
          ConfigurationType::Format,
          "[%time] %log");
    abcd_conf.set(Level::Debug,
          ConfigurationType::Format,
          "[%level - %time] %log");

    abcd_conf.set(Level::Error,
          ConfigurationType::Format,
          "[%level - %time] %log");

    // reset each time
    abcd_conf.setAll(ConfigurationType::RollOutSize, "1");

    Loggers::setDefaultConfigurations(abcd_conf);
    Loggers::reconfigureAllLoggers(abcd_conf);
}               /* -----  end of function configure_logger  ----- */

/*!
 *  \brief logger_set_filename
 *
 *  Configure the logger to output logs to the file
 *
 *  \param log_output: log filename
 *
 */
void logger_set_filename(std::string log_output)
{
    Configurations c;
    if (log_output != "") {
        c.setAll(ConfigurationType::Filename, log_output);
        c.setAll(ConfigurationType::ToFile, "true");
    } else {
        c.setAll(ConfigurationType::ToFile, "false");
    }
    Loggers::reconfigureAllLoggers(c);
}               /* -----  end of function logger_set_filename  ----- */
