//==================================================================================================
// Name        : lr-swig.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : This file contains the interface for possible SWIG'ing to other languages.
//==================================================================================================

#include <string>
#include <vector>

double CalculateLikelihood(bool calculateSuspect, unsigned int numUnknowns,
			   const std::string& alleleFrequencyTableDirectory,
                           const std::string& locusName,
			   const std::string& race,
                           const std::vector<std::string>& suspectAlleles,
                           const std::vector<std::vector<std::string> >& assumedAlleles,
                           const std::vector<std::vector<std::string> >& unattributedAlleles,
                           double zeroAllelesInCommonProb, double oneAlleleInCommonProb,
                           double bothAllelesInCommonProb, double dropoutRate, double dropinRate,
                           double alpha, double fst);
