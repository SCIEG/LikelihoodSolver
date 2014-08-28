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

//#include <iostream>
#include <sstream>

using namespace std;

namespace LabRetriever {

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

UnknownsSolverImpl::UnknownsSolverImpl(const std::vector<std::string>& alleles,
        const AlleleProfile& alleleProfile, const ReplicateData& data,
        const std::map<std::string, double>& alleleProportions,
        double renormalizationFactor,
        double alpha, double dropoutRate, double dropinRate, double noDropinProb,
        unsigned int numUnknownAlleles)
    : alleles(alleles), alleleProfile(alleleProfile), data(data),
      alleleProportions(alleleProportions), renormalizationFactor(renormalizationFactor),
      alpha(alpha), dropoutRate(dropoutRate), dropinRate(dropinRate), noDropinProb(noDropinProb),
      numUnknownAlleles(numUnknownAlleles) {
}

UnknownsSolverImpl::~UnknownsSolverImpl() {
    // TODO Auto-generated destructor stub
}

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

double UnknownsSolverImpl::f(const string& allele, unsigned int count) const {
    const double oneAlleleFrequency = alleleProportions.find(allele)->second;
    const double countAllelesFrequency = pow(oneAlleleFrequency, (double) count);
    const unsigned int adjCount = count + alleleProfile.getAlleleCounts(allele);

    bool alleleIsMasked = data.maskedAlleles.count(allele) != 0;

    if (alleleIsMasked) {
        return countAllelesFrequency;
    }

    // Note if count == 0, countAllelesFrequency is 1.
    bool alleleIsPresent = data.unattributedAlleles.count(allele) != 0;
    if (alleleIsPresent) {
        if (adjCount == 0) {
            return countAllelesFrequency * dropinRate * oneAlleleFrequency / renormalizationFactor;
        } else {
            return countAllelesFrequency *
                    (1 - calculateKDropoutsProbability(alpha, dropoutRate, adjCount));
        }
    } else {
        if (adjCount == 0) {
            return countAllelesFrequency;  // Note, this is equal to 1.
        } else {
            return countAllelesFrequency *
                    calculateKDropoutsProbability(alpha, dropoutRate, adjCount);
        }
    }
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
    for (int i = 0; i < data.size(); i++) {
        double inverseNormalizationFactor = 0;
        const set<string>& maskedAlleles = data[i].maskedAlleles;
        for (set<string>::iterator it = maskedAlleles.begin(); it != maskedAlleles.end(); it++) {
            const string& maskedAllele = *it;
            double maskedAlleleFrequency = alleleProportions.find(maskedAllele)->second;
            inverseNormalizationFactor += maskedAlleleFrequency;
        }
        renormalizationFactors[i] = (1 - inverseNormalizationFactor);
    }

    double probability = 1;
    for (int i = 0; i < data.size(); i++) {
        UnknownsSolverImpl solver(alleles, alleleProfile, data[i], alleleProportions,
                renormalizationFactors[i], alpha, dropoutRate, dropinRate, noDropinProb,
                numUnknownAlleles);
        probability *= solver.solve();
    }

    return probability;
}

double UnknownsSolverImpl::solve() const {
    if (alleleProfile.getAlleles().empty()) {
        unsigned int numCSPAlleles = data.unattributedAlleles.size();
        double x = this->X(numUnknownAlleles, 0, alleles.size());
        double y = 0;
        if (numUnknownAlleles >= numCSPAlleles) {
            // TODO: consider putting this as a member.
            vector<unsigned int> CSPCounts(alleles.size() + 1);
            CSPCounts[0] = 0;
            for (unsigned int i = 0; i < alleles.size(); i++) {
                CSPCounts[i + 1] = CSPCounts[i] +
                        (data.unattributedAlleles.count(alleles[i]) > 0 ? 1 : 0);
            }
            y = this->Y(numUnknownAlleles - numCSPAlleles, 0, alleles.size(), CSPCounts);
        }
        return x - (1 - noDropinProb) * y;
    } else {
//        vector<unsigned int> suspectCounts(alleles.size() + 1);
//        suspectCounts[0] = 0;
//        for (unsigned int i = 0; i < alleles.size(); i++) {
//            const string& allele = alleles[i];
//            suspectCounts[i + 1] = suspectCounts[i] + alleleProfile.getAlleleCounts(allele);
//            std::cout << suspectCounts[i+1] << " ";
//        }
        double x = this->X(numUnknownAlleles, 0, alleles.size());

        set<string> leftoverCSPAlleles;
        for (set<string>::const_iterator iter = data.unattributedAlleles.begin();
                iter != data.unattributedAlleles.end(); iter++) {
            const string& allele = *iter;
            if (!alleleProfile.contains(allele)) {
                leftoverCSPAlleles.insert(allele);
            }
        }

        double y = 0;
        if (numUnknownAlleles >= leftoverCSPAlleles.size()) {
            vector<unsigned int> preallocatedCounts(alleles.size() + 1);
            preallocatedCounts[0] = 0;
            for (unsigned int i = 0; i < alleles.size(); i++) {
                const string& allele = alleles[i];
                preallocatedCounts[i + 1] = preallocatedCounts[i] +
                        (leftoverCSPAlleles.count(allele) ? 1 : 0);
            }
            y = this->Y(numUnknownAlleles - leftoverCSPAlleles.size(), 0, alleles.size(),
                    preallocatedCounts);
        }
//        std::cout << "x: " << x << std::endl;
//        std::cout << "y: " << y << std::endl;
//        std::cout << noDropinProb << std::endl;
        return x - (1 - noDropinProb) * y;
    }
}

// alleleStartIndex is inclusive, alleleEndIndex is exclusive.
double UnknownsSolverImpl::X(int numAlleles, int alleleStartIndex, int alleleEndIndex) const {
    if (alleleStartIndex == alleleEndIndex - 1) {
        // BASE CASE
        return f(alleles[alleleStartIndex], numAlleles);
    }
    int alleleMidpointIndex = (alleleStartIndex + alleleEndIndex) / 2;
    double likelihood = 0;
    for (int i = 0; i <= numAlleles; i++) {
        likelihood += choose(numAlleles, i)
                * X(i, alleleStartIndex, alleleMidpointIndex)
                * X(numAlleles - i, alleleMidpointIndex, alleleEndIndex);
    }
    return likelihood;
}

// alleleStartIndex is inclusive, alleleEndIndex is exclusive.
double UnknownsSolverImpl::Y(int numAlleles, int alleleStartIndex, int alleleEndIndex,
        const vector<unsigned int>& preallocatedCounts) const {
//    std::cout << "calling y with numAlleles " << numAlleles << " from " << alleleStartIndex << " to " << alleleEndIndex << std::endl;
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
}
