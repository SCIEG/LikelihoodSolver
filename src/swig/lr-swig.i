//==================================================================================================
// Name        : lr-swig.i
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

%module LabRetriever
%{
#include "lr-swig.h"
%}

%include "std_string.i"
%include "std_vector.i"
namespace std {
    %template(StringVector) vector<string>;
    %template(StringVectorVector) vector<vector<string> >;
}

%include "lr-swig.h"
