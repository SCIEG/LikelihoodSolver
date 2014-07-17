//==================================================================================================
// Name        : LikelihoodSolver.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : Functions that calculate likelihood are inherently stateless, so it generally does
//   not make sense to make a class for each situation. However, to avoid function pointers all
//   over the place, all derived LikelihoodSolvers will follow the singleton pattern, and they will
//   be accessed via the exemplar pattern, which takes advantage of static code. To create a new
//   LikelihoodSolver, first you need a new ScenarioType, which is just an ID. When you create the
//   class, you statically create the singleton, a.k.a. the exemplar, and add it to the EXEMPLAR_MAP
//   defined below. For example, if I'm creating a FooLikelihoodSolver, then in the ScenarioType
//   enum, I would add a new enum value like:
//
//       enum ScenarioType {
//           ...,
//           FOO_TYPE
//       };
//
//   and the class would look something like:
//
//       class FooLikelihoodSolver : LikelihoodSolver {
//           FooLikelihoodSolver() : LikelihoodSolver(FOO_TYPE, "Foo") {}
//           ...
//           static FooLikelihoodSolver EXEMPLAR;
//       };
//
//       /* Create the singleton exemplar and run the initialization code. */
//       FooLikelihoodSolver FooLikelihoodSolver::EXEMPLAR;
//
//   This avoids the tediousness of creating a pointer to each object and instead allows enums to
//   represent LikelihoodSolvers. Thus it allows code like:
//
//       /* Read in input */
//       ScenarioType type = getScenarioType();
//       Configuration config = getConfig();
//
//       double likelihood = LikelihoodSolver::getSolver(type)->getLikelihood(config);
//==================================================================================================

#ifndef LIKELIHOODSOLVER_H_
#define LIKELIHOODSOLVER_H_

#include "../Configuration.h"
#include <map>


namespace LabRetriever {

    class LikelihoodSolver {
        public:
            enum ScenarioType {
                ONE_SUSPECT_NO_UNKNOWNS,
                ONE_SUSPECT_ONE_UNKNOWN,
                ONE_SUSPECT_TWO_UNKNOWNS,
                ONE_SUSPECT_THREE_UNKNOWNS,
                NO_SUSPECT_ONE_UNKNOWN,
                NO_SUSPECT_TWO_UNKNOWNS,
                NO_SUSPECT_THREE_UNKNOWNS,
                NO_SUSPECT_FOUR_UNKNOWNS,
                ARBITRARY_UNKNOWNS,
                SUSPECT_PLUS_ARBITRARY_UNKNOWNS
            };

            virtual ~LikelihoodSolver() {}
            virtual double getLogLikelihood(const Configuration& config) = 0;
            static LikelihoodSolver* getSolver(ScenarioType type);
            double numComplete;
            double totalToComplete;
            const string name;
        protected:
            LikelihoodSolver(ScenarioType type, const string& name);
        private:
            /*
             * This function exists so that EXEMPLAR_MAP will always be initialized before it is
             * used. See <http://stackoverflow.com/questions/1005685/c-static-initialization-order>.
             */
            static std::map<ScenarioType, LikelihoodSolver*>& getExemplarMap();
            static std::map<ScenarioType, LikelihoodSolver*> EXEMPLAR_MAP;
    };
} /* namespace LabRetriever */
#endif /* LIKELIHOODSOLVER_H_ */
