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

map<Race, vector<double> > run(const string& inputFileName, const string& outputFileName,
        vector<LikelihoodSolver*> likelihoodSolvers);


#endif // LRMAIN_H
