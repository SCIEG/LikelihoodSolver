//==================================================================================================
// Name        : StringUtil.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include "StringUtil.h"

#include <cerrno>
#include <cstdlib>

namespace LabRetriever {

bool ToDouble(const std::string& input, double* output) {
    char* endPtr;
    errno = 0;
    *output = strtod(input.c_str(), &endPtr);
    return *endPtr == '\0' && errno == 0;
}

bool ToInt(const std::string& input, int* output) {
    char* endPtr;
    errno = 0;
    *output = static_cast<int>(strtol(input.c_str(), &endPtr, 10 /* base */));
    return *endPtr == '\0' && errno == 0;
}

}  // namespace LabRetriever
