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

#include <vector>
#include <map>

#include "../Configuration.h"
#include "FileReaderUtil.h"

using namespace std;

namespace LabRetriever {

vector<Race> GetRaces(const Race& race, const string& alleleFrequencyTablePath,
		      const set<string>& lociToRun);

void RetrieveDataFromCSV(const string& inputFileName, double* alpha, double* dropinRate,
                         double* dropoutRate, double* fst, Race* race,
                         IdenticalByDescentProbability* identicalByDescentProbability,
                         map<string, vector<string> >* locusToSuspectAlleles,
                         map<string, vector<set<string> > >* locusToAssumedAlleles,
                         map<string, vector<set<string> > >* locusToUnattributedAlleles,
                         map<string, double>* locusSpecificDropout,
                         map<string, double>* locusSpecificDropin,
                         set<string>* lociToRun);

Configuration CreateConfiguration(
  const string& alleleFrequencyTablePath,
    const string& locus, const Race& race, double alpha, double dropinRate, double dropoutRate,
    double fst, const IdenticalByDescentProbability& identicalByDescentProbability,
    const map<string, vector<string> >& locusToSuspectAlleles,
    const map<string, vector<set<string> > >& locusToAssumedAlleles,
    const map<string, vector<set<string> > >& locusToUnattributedAlleles,
    const map<string, double>& locusSpecificDropout,
    const map<string, double>& locusSpecificDropin);
}  // namespace LabRetriever

#endif /* INPUTPARSERUTIL_H_ */
