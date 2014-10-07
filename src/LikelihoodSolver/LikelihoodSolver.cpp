//==================================================================================================
// Name        : LikelihoodSolver.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include "LikelihoodSolver.h"

namespace LabRetriever {
    std::map<LikelihoodSolver::ScenarioType, LikelihoodSolver*>&
            LikelihoodSolver::getExemplarMap() {
        static std::map<LikelihoodSolver::ScenarioType, LikelihoodSolver*> EXEMPLAR_MAP;
        return EXEMPLAR_MAP;
    }

    LikelihoodSolver* LikelihoodSolver::getSolver(ScenarioType type) {
        return getExemplarMap()[type];
    }

    LikelihoodSolver::LikelihoodSolver(ScenarioType type, const string& name) : name(name) {
        getExemplarMap()[type] = this;
        numComplete = 0;
        totalToComplete = 0;
    }
}


