//==================================================================================================
// Name        : UnknownsSolver.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include "UnknownsSolver.h"

#include "../Configuration.h"
#include "../utils/ProbabilityUtil.h"

// #include <iostream>
#include <sstream>

using namespace std;

namespace LabRetriever {

class UnknownsSolverImpl {
   public:
    static double getUnknownProbability(const AlleleProfile& alleleProfile,
                                        const std::vector<ReplicateData>& data,
                                        const std::map<std::string, double>& alleleProportions,
                                        double alpha, double dropoutRate, double dropinRate,
                                        double noDropinProb, unsigned int numUnknownAlleles);

   private:
    // TODO: Fold renormalizationFactors with data.
    UnknownsSolverImpl(const std::vector<std::string>& alleles, const AlleleProfile& alleleProfile,
                       const std::vector<ReplicateData>& data,
                       const std::map<std::string, double>& alleleProportions,
                       const std::vector<double>& renormalizationFactors, double alpha, double dropoutRate,
                       double dropinRate, double noDropinProb, unsigned int numUnknownAlleles);
    virtual ~UnknownsSolverImpl();

    double solve() const;

    double f(const std::string& allele, unsigned int count) const;
    double X(int numAlleles, int alleleStartIndex, int alleleEndIndex) const;
    double Y(int numAlleles, int alleleStartIndex, int alleleEndIndex,
             const std::vector<unsigned int>& preallocatedCounts) const;

    const std::vector<std::string>& alleles;
    const AlleleProfile& alleleProfile;
    const std::vector<ReplicateData>& data;
    const std::map<std::string, double>& alleleProportions;
    const std::vector<double>& renormalizationFactors;
    unsigned int numReplicates;
    double alpha;
    double dropoutRate;
    double dropinRate;
    double noDropinProb;
    unsigned int numUnknownAlleles;
};

UnknownsSolverImpl::UnknownsSolverImpl(const std::vector<std::string>& alleles,
                                       const AlleleProfile& alleleProfile,
                                       const std::vector<ReplicateData>& data,
                                       const std::map<std::string, double>& alleleProportions,
                                       const std::vector<double>& renormalizationFactors, double alpha,
                                       double dropoutRate, double dropinRate, double noDropinProb,
                                       unsigned int numUnknownAlleles)
    : alleles(alleles),
      alleleProfile(alleleProfile),
      data(data),
      alleleProportions(alleleProportions),
      renormalizationFactors(renormalizationFactors),
      numReplicates(data.size()),
      alpha(alpha),
      dropoutRate(dropoutRate),
      dropinRate(dropinRate),
      noDropinProb(noDropinProb),
      numUnknownAlleles(numUnknownAlleles) {
    // TODO: Check that data.size() == renormalizationFactors.size().
}

UnknownsSolverImpl::~UnknownsSolverImpl() {
    // TODO Auto-generated destructor stub
}

namespace {

unsigned int choose(unsigned int n, unsigned int r) {
    static map<pair<unsigned int, unsigned int>, unsigned int> memoMap;
    if (r == 0 || r == n) return 1;
    if (r == 1 || r == n - 1) return n;
    pair<unsigned int, unsigned int> index(n, min(r, n-r));
    map<pair<unsigned int, unsigned int>, unsigned int>::iterator it = memoMap.find(index);
    if (it != memoMap.end()) {
        return it->second;
    }
    // Memoize and return.
    return (memoMap[index] = choose(n-1, r-1) + choose(n-1, r));
}

// TODO find stabler version of permutation function
unsigned int perm(const map<string, unsigned int>& b) {
    unsigned int totalNum = 0;
    for (map<string, unsigned int>::const_iterator iter = b.begin(); iter != b.end(); iter++) {
        totalNum += iter->second;
    }

    unsigned int retVal = 1;
    for (unsigned int i = 2; i <= totalNum; i++) {
        retVal *= i;
    }

    for (map<string, unsigned int>::const_iterator iter = b.begin(); iter != b.end(); iter++) {
        for (unsigned int i = 2; i <= iter->second; i++) {
            retVal /= i;
        }
    }

    return retVal;
}

}  // namespace

double UnknownsSolverImpl::f(const string& allele, unsigned int count) const {
    const double oneAlleleFrequency = alleleProportions.find(allele)->second;
    const double countAllelesFrequency = pow(oneAlleleFrequency, (double)count);
    const unsigned int adjCount = count + alleleProfile.getAlleleCounts(allele);

    double prob = 1;
    for (unsigned int index = 0; index < numReplicates; ++index) {
        const ReplicateData& replicate_data = data[index];
        bool alleleIsMasked = replicate_data.maskedAlleles.count(allele) != 0;

        if (alleleIsMasked) {
	    continue;
        }

        // Note if count == 0, countAllelesFrequency is 1.
        bool alleleIsPresent = replicate_data.unattributedAlleles.count(allele) != 0;
        if (alleleIsPresent) {
            if (adjCount == 0) {
                prob *= dropinRate * oneAlleleFrequency / renormalizationFactors[index];
            } else {
                prob *= 1 - calculateKDropoutsProbability(alpha, dropoutRate, adjCount);
            }
        } else {
            if (adjCount == 0) {
                continue;
            } else {
                prob *= calculateKDropoutsProbability(alpha, dropoutRate, adjCount);
            }
        }
    }
    return countAllelesFrequency * prob;
}

double UnknownsSolverImpl::getUnknownProbability(const AlleleProfile& alleleProfile,
        const vector<ReplicateData>& data,
        const map<string, double>& alleleProportions, double alpha,
        double dropoutRate, double dropinRate, double noDropinProb,
        unsigned int numUnknownAlleles) {
//    const map<string, unsigned int>& b = alleleProfile.getAlleleCounts();
    vector<string> alleles;
    for (map<string, double>::const_iterator it = alleleProportions.begin();
            it != alleleProportions.end(); it++) {
        alleles.push_back(it->first);
    }

    // Parallel vector to data.
    vector<double> renormalizationFactors(data.size());
    for (unsigned int i = 0; i < data.size(); i++) {
        double inverseNormalizationFactor = 0;
        const set<string>& maskedAlleles = data[i].maskedAlleles;
        for (set<string>::iterator it = maskedAlleles.begin(); it != maskedAlleles.end(); it++) {
            const string& maskedAllele = *it;
            double maskedAlleleFrequency = alleleProportions.find(maskedAllele)->second;
            inverseNormalizationFactor += maskedAlleleFrequency;
        }
        renormalizationFactors[i] = (1 - inverseNormalizationFactor);
    }

    UnknownsSolverImpl solver(alleles, alleleProfile, data, alleleProportions,
                              renormalizationFactors, alpha, dropoutRate, dropinRate, noDropinProb,
                              numUnknownAlleles);
    return solver.solve();
}

namespace {

class CombinationIterator {
   public:
    CombinationIterator(unsigned int n, unsigned int k) : n(n) {
        currentCombination.reserve(k);
        for (unsigned int i = 0; i < k; ++i) {
            currentCombination.push_back(i);
        }
        // Decrement the last value so that next() can be called first.
        --currentCombination.back();
    }

    const std::vector<unsigned int>& getCurrentCombination() const { return currentCombination; }

    // Changes the current combination, and returns if combination is valid.
    bool next() {
        ++currentCombination.back();

        // Check for ripple.
        unsigned int limit = n;
        unsigned int index = currentCombination.size() - 1;
        while (currentCombination[index] >= limit) {
            if (index == 0) { return false; }
            // Ripple previous digit, then check for ripple on digit before.
            ++currentCombination[index - 1];
            --index;
            --limit;
        }
        // Fix all digits past the "okay" digit.
        for (++index; index < currentCombination.size(); ++index) {
            currentCombination[index] = currentCombination[index - 1] + 1;
        }
        return true;
    }

   private:
    const unsigned int n;
    std::vector<unsigned int> currentCombination;
};
}

double UnknownsSolverImpl::solve() const {
    const double x = this->X(numUnknownAlleles, 0, alleles.size());

    // Use inclusion-exclusion principle over replicates to correct drop-in. That is, correct
    // drop-in for all cases where at least one replicate has drop-in, then correct the correction
    // for all cases where at least two replicate has drop-in, etc.
    double prob = 0;
    double previousFactor = 0;
    for (unsigned int numPossibleDropInCases = 1; numPossibleDropInCases <= numReplicates;
         ++numPossibleDropInCases) {
        CombinationIterator combinationIterator(numReplicates, numPossibleDropInCases);
        double y = 0;
        while (combinationIterator.next()) {
            const std::vector<unsigned int> currentCombination =
                    combinationIterator.getCurrentCombination();
            // Find which alleles in the CSP are not explained by the suspect. If all of these
            // alleles appear in some suspect-unknowns combination, we have a case of drop-in.
            set<string> leftoverCSPAlleles;
            for (std::vector<unsigned int>::const_iterator iter = currentCombination.begin();
                 iter != currentCombination.end(); ++iter) {
                unsigned int replicateIndex = *iter;
                const ReplicateData& replicateData = data[replicateIndex];
                for (set<string>::const_iterator iter = replicateData.unattributedAlleles.begin();
                     iter != replicateData.unattributedAlleles.end(); iter++) {
                    const string& allele = *iter;
                    if (!alleleProfile.contains(allele)) {
                        leftoverCSPAlleles.insert(allele);
                    }
                }
            }

            if (numUnknownAlleles >= leftoverCSPAlleles.size()) {
                vector<unsigned int> preallocatedCounts(alleles.size() + 1);
                preallocatedCounts[0] = 0;
                for (unsigned int i = 0; i < alleles.size(); i++) {
                    const string& allele = alleles[i];
                    preallocatedCounts[i + 1] =
                            preallocatedCounts[i] + (leftoverCSPAlleles.count(allele) ? 1 : 0);
                }
                y += this->Y(numUnknownAlleles - leftoverCSPAlleles.size(), 0, alleles.size(),
                             preallocatedCounts);
            }
        }
        int n = numPossibleDropInCases;  // For arithmetic purposes
        double newFactor = 1 - pow(noDropinProb, n);
        double tmp = y * (n * previousFactor - newFactor);
        previousFactor = newFactor;
        prob += tmp;
    }
    // std::cout << "x: " << x << std::endl;
    // std::cout << "y: " << y << std::endl;
    // std::cout << noDropinProb << std::endl;
    return x + prob;
}

// alleleStartIndex is inclusive, alleleEndIndex is exclusive.
double UnknownsSolverImpl::X(int numAlleles, int alleleStartIndex, int alleleEndIndex) const {
    const vector<unsigned int> empty(alleleEndIndex + 1, 0);
    return Y(numAlleles, alleleStartIndex, alleleEndIndex, empty);
}

// alleleStartIndex is inclusive, alleleEndIndex is exclusive.
double UnknownsSolverImpl::Y(int numAlleles, int alleleStartIndex, int alleleEndIndex,
        const vector<unsigned int>& preallocatedCounts) const {
//    std::cout << "calling y with numAlleles " << numAlleles << " from " << alleleStartIndex << " to " << alleleEndIndex << std::endl;
    if (numAlleles < 0) return 0;

    unsigned int numPreallocatedAllelesInCurrentSet =
            preallocatedCounts[alleleEndIndex] - preallocatedCounts[alleleStartIndex];
    if (alleleStartIndex == alleleEndIndex - 1) {
        const string& allele = alleles[alleleStartIndex];
        // BASE CASE
        return f(allele, numAlleles + numPreallocatedAllelesInCurrentSet);
    }

    int alleleMidpointIndex = (alleleStartIndex + alleleEndIndex) / 2;
    unsigned int numPreallocatedAllelesInLowerHalf =
            preallocatedCounts[alleleMidpointIndex] - preallocatedCounts[alleleStartIndex];
    double likelihood = 0;
    for (int i = 0; i <= numAlleles; i++) {
        likelihood += choose(numAlleles + numPreallocatedAllelesInCurrentSet,
                             i + numPreallocatedAllelesInLowerHalf)
                * Y(i, alleleStartIndex, alleleMidpointIndex, preallocatedCounts)
                * Y(numAlleles - i, alleleMidpointIndex, alleleEndIndex, preallocatedCounts);
    }
    return likelihood;
}

////////////////////////////////////////////////////////////////////////////////
// UnknownsSolver

UnknownsSolver::UnknownsSolver(unsigned int numUnknownAlleles) :
    LikelihoodSolver(ARBITRARY_UNKNOWNS, static_cast<ostringstream*>(
            &(ostringstream() << (numUnknownAlleles/2) << "-UNKNOWNS"))->str()),
    numUnknownAlleles(numUnknownAlleles) {}

double UnknownsSolver::getLogLikelihood(const Configuration& config) {
    AlleleProfile emptyProfile;
    // TODO: If IBD probabilities are set, check that there are only two alleles in the profile.
    if (config.identicalByDescentProbability.zeroAllelesInCommonProb == 1) {
        return log(UnknownsSolverImpl::getUnknownProbability(emptyProfile, config.data,
                config.alleleProportions, config.alpha, config.dropoutRate, config.dropinRate,
                calculateNoDropinProbability(config), numUnknownAlleles));
    } else {
        double prob = 0;
        if (config.identicalByDescentProbability.zeroAllelesInCommonProb != 0) {
            prob += config.identicalByDescentProbability.zeroAllelesInCommonProb *
                    UnknownsSolverImpl::getUnknownProbability(emptyProfile,
                            config.data, config.alleleProportions, config.alpha,
                            config.dropoutRate, config.dropinRate,
                            calculateNoDropinProbability(config),
                            numUnknownAlleles);
        }
        if (config.identicalByDescentProbability.oneAlleleInCommonProb != 0) {
            const set<string>& suspectAlleles = config.suspectProfile.getAlleles();
            for (set<string>::const_iterator iter = suspectAlleles.begin();
                    iter != suspectAlleles.end(); iter++) {
                AlleleProfile partialAlleleProfile;
                partialAlleleProfile.addAllele(*iter);
                
                double probabilityFactor =
                    config.identicalByDescentProbability.oneAlleleInCommonProb;
                if (config.suspectProfile.getAlleleCounts(*iter) == 1) {
                    // Divide the probability by 2, since you can choose one of two alleles.
                    probabilityFactor /= 2;
                }
                prob += probabilityFactor *
                        UnknownsSolverImpl::getUnknownProbability(
                                        partialAlleleProfile, config.data,
                                        config.alleleProportions, config.alpha,
                                        config.dropoutRate, config.dropinRate,
                                        calculateNoDropinProbability(config),
                                        numUnknownAlleles - 1);
            }
        }
        if (config.identicalByDescentProbability.bothAllelesInCommonProb != 0) {
            prob += config.identicalByDescentProbability.bothAllelesInCommonProb *
                    UnknownsSolverImpl::getUnknownProbability(
                            config.suspectProfile, config.data,
                            config.alleleProportions, config.alpha,
                            config.dropoutRate, config.dropinRate,
                            calculateNoDropinProbability(config),
                            numUnknownAlleles - 2);
        }
        return log(prob);
    }
}

SuspectUnknownsSolver::SuspectUnknownsSolver(unsigned int numUnknownAlleles) :
    LikelihoodSolver(SUSPECT_PLUS_ARBITRARY_UNKNOWNS,  "SUSPECT-PLUS-" +
            static_cast<ostringstream*>(
                    &(ostringstream() << (numUnknownAlleles/2) << "-UNKNOWNS"))->str()),
    numUnknownAlleles(numUnknownAlleles) {}

double SuspectUnknownsSolver::getLogLikelihood(const Configuration& config) {
    return log(UnknownsSolverImpl::getUnknownProbability(config.suspectProfile, config.data,
            config.alleleProportions, config.alpha, config.dropoutRate, config.dropinRate,
            calculateNoDropinProbability(config), numUnknownAlleles));
}

}  // namespace LabRetriever
