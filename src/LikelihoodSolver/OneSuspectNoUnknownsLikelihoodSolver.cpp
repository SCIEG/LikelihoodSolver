//==================================================================================================
// Name        : OneSuspectLikelihoodSolver.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

//#include "DebugUtil.h"

#include "CachingSolver.h"
#include "../utils/LikelihoodUtil.h"
#include "../utils/ProbabilityUtil.h"
#include <iostream>
#include <algorithm>
#include <cassert>


using namespace std;
namespace LabRetriever {

    class OneSuspectLikelihoodSolver: public LabRetriever::CachingSolver {
        public:
            double calculateLogLikelihood(const Configuration& config);
            double getLikelihood(const Configuration& config);
        private:
            static bool cmpFn(const Configuration& left, const Configuration& right);
            OneSuspectLikelihoodSolver();


            const static OneSuspectLikelihoodSolver EXEMPLAR;
    };

    START_COMPARE_FUNCTION(OneSuspectLikelihoodSolver, cmpFn) {
        COMPARE_CONFIGURATION_ELEMENT(suspectProfile);
        COMPARE_CONFIGURATION_ELEMENT(data);
        COMPARE_CONFIGURATION_ELEMENT(alleleProportions);
        COMPARE_CONFIGURATION_ELEMENT(dropoutRate);
        COMPARE_CONFIGURATION_ELEMENT(dropinRate);
        COMPARE_CONFIGURATION_ELEMENT(alpha);
    } END_COMPARE_FUNCTION;

    OneSuspectLikelihoodSolver::OneSuspectLikelihoodSolver() : CachingSolver(
            ONE_SUSPECT_NO_UNKNOWNS, "One-Suspect_No-Unknowns", cmpFn) {}

    double OneSuspectLikelihoodSolver::calculateLogLikelihood(const Configuration& config) {
        numComplete = 0; totalToComplete = 1;
        double dropinRate = config.dropinRate;
        double noDropinProb = calculateNoDropinProbability(config);
        const map<string, double>& alleleProportions = config.alleleProportions;
        const vector<ReplicateData>& data = config.data;
        const AlleleProfile& suspectProfile = config.suspectProfile;
        double alpha = config.alpha;
        double dropoutRate = config.dropoutRate;

        double logLikelihood = calculateLogProbability(suspectProfile, data, alleleProportions,
                alpha, dropoutRate, dropinRate, noDropinProb);
        numComplete = 1;
        return logLikelihood;
    }

//    double OneSuspectLikelihoodSolver::getLikelihood(const Configuration& config) {
//        numComplete = 0; totalToComplete = 1;
//        double dropinRate = config.dropinRate;
//        double noDropinProb = calculateNoDropinProbability(config);
//        const map<string, double>& alleleProportions = config.alleleProportions;
//        const vector<ReplicateData>& data = config.data;
//        const AlleleProfile& suspectProfile = config.suspectProfile;
//        double alpha = config.alpha;
//        double dropoutRate = config.dropoutRate;
//
//        double likelihood = calculateProbability(suspectProfile, data, alleleProportions,
//                alpha, dropoutRate, dropinRate, noDropinProb);
//        numComplete = 1;
//        return likelihood;
//    }

    const OneSuspectLikelihoodSolver OneSuspectLikelihoodSolver::EXEMPLAR;
} /* namespace LabRetriever */
