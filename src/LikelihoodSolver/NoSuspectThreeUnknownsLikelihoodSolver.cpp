//==================================================================================================
// Name        : NoSuspectThreeUnknownsLikelihoodSolver.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================


#include "CachingSolver.h"
#include "../utils/LikelihoodUtil.h"
#include "../utils/ProbabilityUtil.h"
#include <iostream>
#include <algorithm>
#include <cassert>


using namespace std;
namespace LabRetriever {
    class NoSuspectThreeUnknownsLikelihoodSolver: public CachingSolver {
        public:
            double calculateLogLikelihood(const Configuration& config);
        private:
            static bool cmpFn(const Configuration& left, const Configuration& right);
            NoSuspectThreeUnknownsLikelihoodSolver();


            const static NoSuspectThreeUnknownsLikelihoodSolver EXEMPLAR;
    };

    START_COMPARE_FUNCTION(NoSuspectThreeUnknownsLikelihoodSolver, cmpFn) {
        COMPARE_CONFIGURATION_ELEMENT(suspectProfile);
        COMPARE_CONFIGURATION_ELEMENT(data);
        COMPARE_CONFIGURATION_ELEMENT(alleleProportions);
        COMPARE_CONFIGURATION_ELEMENT(identicalByDescentProbability);
        COMPARE_CONFIGURATION_ELEMENT(dropoutRate);
        COMPARE_CONFIGURATION_ELEMENT(dropinRate);
        COMPARE_CONFIGURATION_ELEMENT(alpha);
    } END_COMPARE_FUNCTION;

    NoSuspectThreeUnknownsLikelihoodSolver::NoSuspectThreeUnknownsLikelihoodSolver() :
            CachingSolver(NO_SUSPECT_THREE_UNKNOWNS, "No-Suspect_Three-Unknowns", cmpFn) {}

    double NoSuspectThreeUnknownsLikelihoodSolver::calculateLogLikelihood(
            const Configuration& config) {
        double dropinRate = config.dropinRate;
        double noDropinProb = calculateNoDropinProbability(config);
        const map<string, double>& alleleProportions = config.alleleProportions;
        const IdenticalByDescentProbability& ibdProbability =
                config.identicalByDescentProbability;
        const vector<ReplicateData>& data = config.data;
        const AlleleProfile& suspectProfile = config.suspectProfile;
        double alpha = config.alpha;
        double dropoutRate = config.dropoutRate;
        double numAlleles = alleleProportions.size();
        numComplete = 0; totalToComplete =  pow(numAlleles, 4) *
            (2 * numAlleles + numAlleles * numAlleles + 1);

        double zeroIBDLogLikelihood = LOG_ZERO;
        if (ibdProbability.zeroAllelesInCommonProb != 0) {
            BEGIN_CHOOSE_RANDOM_ALLELES(alleleProportions, allele1, logProbRandomAllele1, allele2,
                    logProbRandomAllele2) {
                BEGIN_CHOOSE_RANDOM_ALLELES(alleleProportions, allele3, logProbRandomAllele3,
                        allele4, logProbRandomAllele4) {
                    BEGIN_CHOOSE_RANDOM_ALLELES(alleleProportions, allele5, logProbRandomAllele5,
                            allele6, logProbRandomAllele6) {
                        double logProbHavingTheseAlleles =
                                logProbRandomAllele1 + logProbRandomAllele2 +
                                logProbRandomAllele3 + logProbRandomAllele4 +
                                logProbRandomAllele5 + logProbRandomAllele6;
                        if (allele1 != allele2) {
                            logProbHavingTheseAlleles += LOG_TWO;
                        }
                        if (allele3 != allele4) {
                            logProbHavingTheseAlleles += LOG_TWO;
                        }
                        if (allele5 != allele6) {
                            logProbHavingTheseAlleles += LOG_TWO;
                        }

                        // Create a new suspect profile.
                        AlleleProfile newSuspectProfile;
                        newSuspectProfile.addAllele(allele1).addAllele(allele2)
                                .addAllele(allele3).addAllele(allele4)
                                .addAllele(allele5).addAllele(allele6);

                        double logProbGivenAlleles = calculateLogProbability(newSuspectProfile,
                                data, alleleProportions, alpha, dropoutRate, dropinRate,
                                noDropinProb);
                        double partialZeroIBDLogLikelihood = logProbGivenAlleles +
                                logProbHavingTheseAlleles;

                        zeroIBDLogLikelihood = addLogProbability(zeroIBDLogLikelihood,
                                partialZeroIBDLogLikelihood);
                        numComplete += 1;
                    } END_CHOOSE_RANDOM_ALLELES;
                } END_CHOOSE_RANDOM_ALLELES;
            } END_CHOOSE_RANDOM_ALLELES;
        }
        numComplete = pow(numAlleles, 6);

        double oneIBDLogLikelihood = LOG_ZERO;
        if (ibdProbability.oneAlleleInCommonProb != 0) {
            // Choose one allele from the suspect:
            const set<string>& suspectAlleles = suspectProfile.getAlleles();
            BEGIN_CHOOSE_ONE_RANDOM_ALLELE(alleleProportions, suspectAlleles, suspectAllele,
                    otherAllele, logProbHavingRandomAllele) {
                BEGIN_CHOOSE_RANDOM_ALLELES(alleleProportions, allele1, logProbRandomAllele1,
                        allele2, logProbRandomAllele2) {
                    BEGIN_CHOOSE_RANDOM_ALLELES(alleleProportions, allele3, logProbRandomAllele3,
                            allele4, logProbRandomAllele4) {
                        double logProbHavingTheseAlleles =
                                logProbHavingRandomAllele + logProbRandomAllele1 +
                                logProbRandomAllele2 + logProbRandomAllele3 + logProbRandomAllele4;
                        if (allele1 == allele2) {
                            logProbHavingTheseAlleles += LOG_TWO;
                        }
                        if (allele3 == allele4) {
                            logProbHavingTheseAlleles += LOG_TWO;
                        }

                        AlleleProfile relatedSuspectProfile;
                        relatedSuspectProfile.addAllele(suspectAllele).addAllele(otherAllele)
                                .addAllele(allele1).addAllele(allele2)
                                .addAllele(allele3).addAllele(allele4);

                        double logProbGivenAllele = calculateLogProbability(relatedSuspectProfile,
                                data, alleleProportions, alpha, dropoutRate, dropinRate,
                                noDropinProb);
                        double partialOneIBDLogLikelihood = logProbGivenAllele +
                                logProbHavingRandomAllele;

                        oneIBDLogLikelihood = addLogProbability(oneIBDLogLikelihood,
                                partialOneIBDLogLikelihood);
                        numComplete += 1;
                    } END_CHOOSE_RANDOM_ALLELES;
                } END_CHOOSE_RANDOM_ALLELES;
            } END_CHOOSE_ONE_RANDOM_ALLELE;
        }
        // Divide by two for the two alleles you can choose from the suspect.
        oneIBDLogLikelihood -= LOG_TWO;
        numComplete = totalToComplete - pow(numAlleles, 4);

        double bothIBDLogLikelihood = (ibdProbability.bothAllelesInCommonProb == 0) ?
                LOG_ZERO :
                LikelihoodSolver::getSolver(ONE_SUSPECT_TWO_UNKNOWNS)->getLogLikelihood(config);

        double logLikelihood = addLogProbability(
                log(ibdProbability.zeroAllelesInCommonProb) + zeroIBDLogLikelihood,
                addLogProbability(
                        log(ibdProbability.oneAlleleInCommonProb) + oneIBDLogLikelihood,
                        log(ibdProbability.bothAllelesInCommonProb) + bothIBDLogLikelihood));
        numComplete = totalToComplete;
        return logLikelihood;
    }

    const NoSuspectThreeUnknownsLikelihoodSolver NoSuspectThreeUnknownsLikelihoodSolver::EXEMPLAR;
} /* namespace LabRetriever */

