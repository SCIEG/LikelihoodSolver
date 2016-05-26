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
#include <algorithm>
#include <cmath>
#include <cstdlib>
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

template <typename Type>
class DistributionSampler {
   public:
    // Assumes distribution is normalized and lives longer than this object.
    DistributionSampler(const std::map<Type, double>& distribution) : distribution(distribution) {
        double total_prob_so_far = 0;
        for (typename std::map<Type, double>::const_iterator iter = distribution.begin();
             iter != distribution.end(); ++iter) {
            if (iter->second == 0) continue;
            total_prob_so_far += iter->second;
            cdf.push_back(std::pair<int, const Type*>(
                    static_cast<int>(total_prob_so_far * RAND_MAX), &iter->first));
        }
    }

    const Type& sample() const {
        int value = rand();
	typename std::vector<std::pair<int, const Type*> >::const_iterator iter =
                std::lower_bound(cdf.begin(), cdf.end(), value, comp_less);
        return iter == cdf.end() ? *cdf.back().second : *iter->second;
    }

   private:
    static bool comp_less(const std::pair<double, const Type*>& left, int right) {
        return left.first < right;
    }
    const std::map<Type, double>& distribution;
    std::vector<std::pair<int, const Type*> > cdf;
};
}

#endif /* PROBABILITYUTIL_H_ */
