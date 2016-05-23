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

class SuspectUnknownsSolver : public LikelihoodSolver {
public:
    SuspectUnknownsSolver(unsigned int numUnknownAlleles);
    virtual ~SuspectUnknownsSolver() {}

    virtual double getLogLikelihood(const Configuration& config);

private:
    unsigned int numUnknownAlleles;
};

}

#endif /* UNKNOWNSSOLVER_H_ */
