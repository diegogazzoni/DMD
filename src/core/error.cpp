#include "error.h"

ForceException::ForceException(ForceError e, const char* comp, const char* msg)
    : std::runtime_error(msg), err(e), component(comp) {}

ForceException::ForceException(ForceError e, const char* comp, const std::string& msg)
    : std::runtime_error(msg), err(e), component(comp) {}
