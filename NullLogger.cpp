#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include "NullLogger.hpp"

NullLogger::~NullLogger() {}

void NullLogger::level(enum NullLogger::level)
{
}

enum NullLogger::level NullLogger::level() const
{
    return L_OFF;
}

void NullLogger::logv(enum level lv, char const* format, va_list ap)
{
}

void NullLogger::flush()
{
}

NullLoggerFactory::~NullLoggerFactory() {}

void NullLoggerFactory::level(enum Logger::level)
{
}

Logger* NullLoggerFactory::operator()(char const* logger_name) const
{
    return new NullLogger(logger_name);
}

char const* NullLoggerFactory::get_name() const
{
    return "NullLogger";
}
