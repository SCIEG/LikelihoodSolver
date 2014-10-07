//==================================================================================================
// Name        : OneSuspectOneUnkownLikelihoodSolver.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================


#include "CachingSolver.h"
#include "../utils/LikelihoodUtil.h"
#include "../utils/ProbabilityUtil.h"
#include <iostream>
#include <algorithm>
#include <cassert>


using namespace std;
namespace LabRetriever {
    class OneSuspectOneUnknownLikelihoodSolver: public CachingSolver {
        public:
            double calculateLogLikelihood(const Configuration& config);
        private:
            static bool cmpFn(const Configuration& left, const Configuration& right);
            OneSuspectOneUnknownLikelihoodSolver();


            const static OneSuspectOneUnknownLikelihoodSolver EXEMPLAR;
    };

    START_COMPARE_FUNCTION(OneSuspectOneUnknownLikelihoodSolver, cmpFn) {
        COMPARE_CONFIGURATION_ELEMENT(suspectProfile);
        COMPARE_CONFIGURATION_ELEMENT(data);
        COMPARE_CONFIGURATION_ELEMENT(alleleProportions);
        COMPARE_CONFIGURATION_ELEMENT(dropoutRate);
        COMPARE_CONFIGURATION_ELEMENT(dropinRate);
        COMPARE_CONFIGURATION_ELEMENT(alpha);
    } END_COMPARE_FUNCTION;

    OneSuspectOneUnknownLikelihoodSolver::OneSuspectOneUnknownLikelihoodSolver() :
            CachingSolver(ONE_SUSPECT_ONE_UNKNOWN, "One-Suspect_One-Unknown", cmpFn) {}

    double OneSuspectOneUnknownLikelihoodSolver::calculateLogLikelihood(
            const Configuration& config) {
        double dropinRate = config.dropinRate;
        double noDropinProb = calculateNoDropinProbability(config);
        const map<string, double>& alleleProportions = config.alleleProportions;
        const vector<ReplicateData>& data = config.data;
        const AlleleProfile& suspectProfile = config.suspectProfile;
        double alpha = config.alpha;
        double dropoutRate = config.dropoutRate;
        numComplete = 0; totalToComplete = pow((double) alleleProportions.size(), 2);

        double logLikelihood = LOG_ZERO;

        BEGIN_CHOOSE_RANDOM_ALLELES(alleleProportions, allele1, logProbRandomAllele1, allele2,
                logProbRandomAllele2) {
            double logProbHavingTheseAlleles = logProbRandomAllele1 + logProbRandomAllele2;
            if (allele1 != allele2) {
                logProbHavingTheseAlleles += LOG_TWO;
            }

            // Create a suspect profile + unknown.
            AlleleProfile newSuspectProfile = suspectProfile;
            newSuspectProfile.addAllele(allele1).addAllele(allele2);

            double logProbGivenAlleles = calculateLogProbability(newSuspectProfile, data,
                    alleleProportions, alpha, dropoutRate, dropinRate, noDropinProb);
            double partialLogLikelihood = logProbGivenAlleles + logProbHavingTheseAlleles;

            logLikelihood = addLogProbability(logLikelihood, partialLogLikelihood);
            numComplete += 1;
        } END_CHOOSE_RANDOM_ALLELES;

        return logLikelihood;
    }

    const OneSuspectOneUnknownLikelihoodSolver OneSuspectOneUnknownLikelihoodSolver::EXEMPLAR;
} /* namespace LabRetriever */

