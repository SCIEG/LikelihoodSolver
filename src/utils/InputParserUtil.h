//==================================================================================================
// Name        : InputParserUtil.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#ifndef INPUTPARSERUTIL_H_
#define INPUTPARSERUTIL_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "../Configuration.h"
#include "FileReaderUtil.h"

namespace LabRetriever {

inline std::string GetAlleleFrequencyTableFileName(const std::string& alleleFrequencyTablePath, const std::string locus) {
  return alleleFrequencyTablePath + locus + "_B.count.csv";
}

std::vector<Race> GetRaces(const Race& race, const std::string& alleleFrequencyTablePath,
                           const std::set<std::string>& lociToRun);

std::map<Race, std::map<string, double> > GetRaceToAlleleProportions(
        const std::string& alleleFrequencyTableFileName,
        const std::vector<std::string>& suspectAlleles,
        const std::vector<std::set<std::string> >& assumedAlleles,
        const std::vector<std::set<std::string> >& unattributedAlleles, double fst);

void RetrieveDataFromCSV(
        const std::string& inputFileName, double* alpha, double* dropinRate, double* dropoutRate,
        double* fst, Race* race, IdenticalByDescentProbability* identicalByDescentProbability,
        std::map<std::string, std::vector<std::string> >* locusToSuspectAlleles,
        std::map<std::string, std::vector<std::set<std::string> > >* locusToAssumedAlleles,
        std::map<std::string, std::vector<std::set<std::string> > >* locusToUnattributedAlleles,
        std::map<std::string, double>* locusSpecificDropout,
        std::map<std::string, double>* locusSpecificDropin, std::set<std::string>* lociToRun);

}  // namespace LabRetriever

#endif /* INPUTPARSERUTIL_H_ */
