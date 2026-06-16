#pragma once

#include <stdexcept>
#include <string>

enum class ForceError {
    NaN,
    BondExplosion,
    Mismatch,
    ParmError,
    GPUError
};

class ForceException : public std::runtime_error {
public:
    ForceError err;
    std::string component;

    ForceException(ForceError e, const char* comp, const char* msg);
    ForceException(ForceError e, const char* comp, const std::string& msg);
};
