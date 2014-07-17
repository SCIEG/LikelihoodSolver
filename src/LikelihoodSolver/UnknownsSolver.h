//==================================================================================================
// Name        : UnknownsSolver.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#ifndef UNKNOWNSSOLVER_H_
#define UNKNOWNSSOLVER_H_

#include <map>
#include <string>
#include <vector>

#include "LikelihoodSolver.h"

namespace LabRetriever {

class AlleleProfile;
class ReplicateData;

class UnknownsSolver : public LikelihoodSolver {
public:
    UnknownsSolver(unsigned int numUnknownAlleles);
    virtual ~UnknownsSolver() {}

    virtual double getLogLikelihood(const Configuration& config);

private:
    unsigned int numUnknownAlleles;
};

class SuspectUnknownsSolver :public LikelihoodSolver {
public:
    SuspectUnknownsSolver(unsigned int numUnknownAlleles);
    virtual ~SuspectUnknownsSolver() {}

    virtual double getLogLikelihood(const Configuration& config);

private:
    unsigned int numUnknownAlleles;
};

class UnknownsSolverImpl {
public:
    static double getUnknownProbability(const AlleleProfile& alleleProfile,
            const std::vector<ReplicateData>& data,
            const std::map<std::string, double>& alleleProportions,
            double alpha, double dropoutRate, double dropinRate, double noDropinProb,
            unsigned int numUnknownAlleles);
private:
    UnknownsSolverImpl(const std::vector<std::string>& alleles, const AlleleProfile& alleleProfile,
            const ReplicateData& data,
            const std::map<std::string, double>& alleleProportions,
            double renormalizationFactors,
            double alpha, double dropoutRate, double dropinRate, double noDropinProb,
            unsigned int numUnknownAlleles);
    virtual ~UnknownsSolverImpl();

    double solve() const;

    double f(const std::string& allele, unsigned int count) const;
    double X(int numAlleles, int alleleStartIndex, int alleleEndIndex) const;
    double Y(int numAlleles, int alleleStartIndex, int alleleEndIndex,
            const std::vector<unsigned int>& preallocatedCounts) const;

    const std::vector<std::string>& alleles;
    const AlleleProfile& alleleProfile;
    const ReplicateData& data;
    const std::map<std::string, double>& alleleProportions;
    double renormalizationFactor;
    double alpha;
    double dropoutRate;
    double dropinRate;
    double noDropinProb;
    unsigned int numUnknownAlleles;
};

}

#endif /* UNKNOWNSSOLVER_H_ */
