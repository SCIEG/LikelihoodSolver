//==================================================================================================
// Name        : ProbabilityUtil.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include "ProbabilityUtil.h"
#include <cmath>

namespace LabRetriever {
    // TODO: log(0) isn't always defined in all implementations; replace with -inf.
    const double LOG_ZERO = log(0);
    const double LOG_TWO = log(2);

    double addLogProbability(double logProb1, double logProb2) {
        if (logProb1 == LOG_ZERO) return logProb2;
        else if (logProb2 == LOG_ZERO) return logProb1;

        return logProb1 + log(1 + exp(logProb2 - logProb1));
    }

    double complementLogProbability(double logProb) {
//        static const double INV_E = exp(-1);
//        return -log(INV_E + exp(logProb));
        return log(1 - exp(logProb));
    }

    double calculateKDropoutsLogProbability(double alpha, double dropoutRate, unsigned int k) {
        if (dropoutRate == 0) return LOG_ZERO;
        if (k == 0) return 0;
        double logAlpha = log(alpha);
        double logDropoutProb = log(dropoutRate);
        double logFinalProb = (k - 1) * (logAlpha + logDropoutProb) + logDropoutProb;
        return logFinalProb;
    }

    double calculateKDropoutsProbability(double alpha, double dropoutRate, unsigned int k) {
        if (k == 1 || dropoutRate == 0) return dropoutRate;
        if (k == 0) return 1;
        double logFinalProb = calculateKDropoutsLogProbability(alpha, dropoutRate, k);
        return exp(logFinalProb);
    }

    double calculateNoDropinProbability(const Configuration& config) {
       double dropinRate = config.dropinRate;
       int numOfAlleles = config.alleleProportions.size();
       return (1 - 2 * dropinRate + pow(dropinRate, numOfAlleles + 1)) / (1 - dropinRate);
    }

    double calculateNoDropinLogProbability(const Configuration& config) {
        // Since the drop-in rate is relatively small, the regular function on normal inputs is
        // numerically stable. There is no need to do anything special to calculate the log
        // probability.
        return log(calculateNoDropinProbability(config));
    }
}


