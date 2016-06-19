//==================================================================================================
// Name        : lr-swig.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

#include "lr-swig.h"

#include <algorithm>
#include <cmath>
#include <set>

#include "../Configuration.h"
#include "../LikelihoodSolver/UnknownsSolver.h"
#include "../utils/FileReaderUtil.h"
#include "../utils/InputParserUtil.h"

using namespace LabRetriever;

double CalculateLikelihoodImpl(LikelihoodSolver* solver, const Configuration& config) {
    return exp(solver->getLogLikelihood(config));
}

double CalculateLikelihood(bool calculateSuspect, unsigned int numUnknowns,
                           const std::string& alleleFrequencyTableDirectory,
                           const std::string& locusName,
			   const std::string& race,
                           const std::vector<std::string>& suspectAlleles,
                           const std::vector<std::vector<std::string> >& assumedAlleles,
                           const std::vector<std::vector<std::string> >& unattributedAlleles,
                           double zeroAllelesInCommonProb, double oneAlleleInCommonProb,
                           double bothAllelesInCommonProb, double dropoutRate, double dropinRate,
                           double alpha, double fst) {
    std::vector<std::set<std::string> > assumedAllelesSet;
    assumedAllelesSet.reserve(assumedAlleles.size());
    for (unsigned int i = 0; i < assumedAlleles.size(); ++i) {
        assumedAllelesSet.push_back(
                std::set<std::string>(assumedAlleles[i].begin(), assumedAlleles[i].end()));
    }
    std::vector<std::set<std::string> > unattributedAllelesSet;
    unattributedAllelesSet.reserve(unattributedAlleles.size());
    for (unsigned int i = 0; i < unattributedAlleles.size(); ++i) {
        unattributedAllelesSet.push_back(std::set<std::string>(unattributedAlleles[i].begin(),
                                                               unattributedAlleles[i].end()));
    }
    
    const std::string alleleFrequencyTableFileName =
	GetAlleleFrequencyTableFileName(alleleFrequencyTableDirectory, locusName);

    Race r(race);
    std::map<std::string, double> alleleProportions =
            GetRaceToAlleleProportions(alleleFrequencyTableFileName, suspectAlleles,
                                       assumedAllelesSet, unattributedAllelesSet, fst)[r];

    // They should be equal, but take the min, just in case.
    const unsigned int size = min(assumedAllelesSet.size(), unattributedAllelesSet.size());
    std::vector<ReplicateData> replicateData;
    replicateData.reserve(size);
    for (unsigned int i = 0; i < size; ++i) {
        replicateData.push_back(ReplicateData::fromUnattributedAndMaskedAlleles(
                unattributedAllelesSet[i], assumedAllelesSet[i]));
    }

    IdenticalByDescentProbability identicalByDescentProbability(
            zeroAllelesInCommonProb, oneAlleleInCommonProb, bothAllelesInCommonProb);

    Configuration config(suspectAlleles, assumedAllelesSet, unattributedAllelesSet,
                         alleleProportions, identicalByDescentProbability, dropoutRate, dropinRate,
                         alpha);
    // Note that the solvers expect number of unknown alleles; not the number of unknown people.
    // TODO: Only one of these needs to be constructed.
    SuspectUnknownsSolver suspectSolver(numUnknowns * 2);
    UnknownsSolver unknownsSolver(numUnknowns * 2);
    LikelihoodSolver* solver = calculateSuspect ? static_cast<LikelihoodSolver*>(&suspectSolver)
                                                : static_cast<LikelihoodSolver*>(&unknownsSolver);
    return CalculateLikelihoodImpl(solver, config);
}
