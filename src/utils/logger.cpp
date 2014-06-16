#include "abcd.h"

using namespace easyloggingpp;

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
}

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
}

