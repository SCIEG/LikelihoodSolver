//==================================================================================================
// Name        : ProbabilityUtil.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : Defines some useful functions to help calculate probabilities.
//==================================================================================================

#ifndef PROBABILITYUTIL_H_
#define PROBABILITYUTIL_H_

#include "../Configuration.h"
#include <cmath>
#include <map>
#include <string>

using namespace std;

namespace LabRetriever {
    extern const double LOG_ZERO;
    extern const double LOG_TWO;

    double addLogProbability(double logProb1, double logProb2);

    /*
     * Returns the log of (1 - e^logProb).
     */
    double complementLogProbability(double logProb);

    /*
     * Returns a mapping from alleles to the proportion that the allele appears in people, after
     * sampling adjustment and the Balding-Nichols FST correction.
     *
     * alleleCounts - a mapping from alleles to a count of how many times that allele appeared in
     *     the sample.
     * suspectProfile - the current suspect's profile.
     * samplingAdjustment -
     * fst - See Balding-Nichols FST correction. Balding recommends that fst should be at least
     *     0.01. It can be as high as 0.05, depending on how representative the sample population
     *     is of the whole target population.
     */
    map<string, double> getAlleleProportionsFromCounts(
            const map<string, unsigned int>& alleleCounts, const AlleleProfile& suspectProfile,
            double fst = 0.01, unsigned int samplingAdjustment = 2);

    /*
     * Calculates the probability of k alleles dropping out, given some extra data. The current
     * implementation computes it like this:
     *
     *     alpha ^ (k - 1) * dropoutRate ^ k
     *
     * which is suggested by Balding & Buckleton (2009). The Configuration object is passed in,
     * just in case a different model is needed with different parameters.
     */
    double calculateKDropoutsProbability(double alpha, double dropoutRate, unsigned int k);
    double calculateKDropoutsLogProbability(double alpha, double dropoutRate, unsigned int k);

    /*
     * Calculates the probability of none of the alleles dropping out. Since drop-in can only be
     * estimated by sampling, this function attempts to correct some error using this calculation:
     *
     *     (1 - 2 * dropoutRate + dropoutRate ^ (1 + numOfAlleles)) / (1 - dropoutRate)
     */
    double calculateNoDropinProbability(const Configuration& config);
    double calculateNoDropinLogProbability(const Configuration& config);
}

#endif /* PROBABILITYUTIL_H_ */
