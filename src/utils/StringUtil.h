//==================================================================================================
// Name        : StringUtil.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include <string>

namespace LabRetriever {

// Returns true on success.
bool ToDouble(const std::string& input, double* output);
bool ToInt(const std::string& input, int* output);

}  // namespace LabRetriever
