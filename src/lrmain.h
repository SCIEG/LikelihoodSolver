//==================================================================================================
// Name        : lrmain.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

#ifndef LRMAIN_H
#define LRMAIN_H

#include "Configuration.h"
#include "utils/FileReaderUtil.h"
#include "LikelihoodSolver/LikelihoodSolver.h"
#include "LikelihoodSolver/UnknownsSolver.h"
#include "utils/ProbabilityUtil.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <float.h>
#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;
using namespace LabRetriever;

void outputData(const set<string>& lociToCheck, const vector<LikelihoodSolver*>& likelihoodSolvers,
        const map<Race, vector<map<string, double> > >& raceToSolverIndexToLocusLogProb,
        const map<Race, vector<double> >& raceToSolverIndexToLogProb,
        const vector<Race> races,
        const string& outputFileName);

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

map<Race, vector<double> > run(const string& inputFileName, const string& outputFileName,
        vector<LikelihoodSolver*> likelihoodSolvers);


#endif // LRMAIN_H
