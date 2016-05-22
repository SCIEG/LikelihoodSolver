//==================================================================================================
// Name        : lrmain.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

#include "lrmain.h"

#include "LikelihoodSolver/UnknownsSolver.h"
#include "utils/InputParserUtil.h"
#include "utils/StlUtil.h"

using namespace LabRetriever;

void outputData(const set<string>& lociToRun, const vector<LikelihoodSolver*>& likelihoodSolvers,
        const map<Race, vector<map<string, double> > >& raceToSolverIndexToLocusLogProb,
        const map<Race, vector<double> >& raceToSolverIndexToLogProb,
        const vector<Race> races,
        const string& outputFileName) {
    // TODO: Output in proper format to file:
    stringstream outputStringStream;
    outputStringStream.setf(ios::scientific);
    outputStringStream.precision(4);

    for (unsigned int raceIndex = 0; raceIndex < races.size(); raceIndex++) {
        Race curRace = races[raceIndex];

        const vector<map<string, double> >& solverIndexToLocusLogProb =
                raceToSolverIndexToLocusLogProb.find(curRace)->second;
        const vector<double>& solverIndexToLogProb =
                raceToSolverIndexToLogProb.find(curRace)->second;

        // Output race
        outputStringStream << curRace << endl;

        // Output probability header
        outputStringStream << "Probabilities, total";
        for (set<string>::const_iterator iter = lociToRun.begin(); iter != lociToRun.end();
                iter++) {
            const string& locus = *iter;
            outputStringStream << ", " << locus;
        }
        outputStringStream << endl;

        // Output probability data
        for (unsigned int solverIndex = 0; solverIndex < likelihoodSolvers.size(); solverIndex++) {
            map<string, double> locusToLogProb = solverIndexToLocusLogProb[solverIndex];
            outputStringStream << likelihoodSolvers[solverIndex]->name << ", "
                    << exp(solverIndexToLogProb[solverIndex]);
            for (set<string>::const_iterator iter = lociToRun.begin(); iter != lociToRun.end();
                    iter++) {
                const string& locus = *iter;
                double logProb = locusToLogProb[locus];
                outputStringStream << ", " << exp(logProb);
            }
            outputStringStream << endl;
        }

        outputStringStream << endl;

        // Output probability ratio header
        outputStringStream << "Probabilities Ratios, total";
        for (set<string>::const_iterator iter = lociToRun.begin(); iter != lociToRun.end();
                iter++) {
            const string& locus = *iter;
            outputStringStream << ", " << locus;
        }

        outputStringStream << endl;

        // Output probability ratios
        for (unsigned int i = 0; i < likelihoodSolvers.size(); i++) {
            map<string, double> locusToLogProb_i = solverIndexToLocusLogProb[i];
            for (unsigned int j = i + 1; j < likelihoodSolvers.size(); j++) {
                map<string, double> locusToLogProb_j = solverIndexToLocusLogProb[j];
                string ratioName = likelihoodSolvers[i]->name + " to " + likelihoodSolvers[j]->name;
                double diff = solverIndexToLogProb[i] - solverIndexToLogProb[j];
                outputStringStream << ratioName << ", " << exp(diff);
                for (set<string>::const_iterator iter = lociToRun.begin(); iter != lociToRun.end();
                        iter++) {
                    const string& locus = *iter;
                    double logProbDiff = locusToLogProb_i[locus] - locusToLogProb_j[locus];
                    outputStringStream << ", " << exp(logProbDiff);
                }
                outputStringStream << endl;
            }
        }
        outputStringStream << endl;
    }


    string dataToOutput = outputStringStream.str();

    cout << dataToOutput;

    ofstream myFileStream;
    myFileStream.open(outputFileName.c_str());
    myFileStream << dataToOutput;
    myFileStream.close();
}


map<Race, vector<double> > run(const string& executablePath, const string& inputFileName,
                               const string& outputFileName,
                               vector<LikelihoodSolver*> likelihoodSolvers) {
    // These are defaults; if specified, will be in the input file.
    double alpha = 0.5;
    double dropinRate = 0.01;
    double dropoutRate = 0.01;
    double fst = 0.01;
    Race race = ALL_RACE;
    IdenticalByDescentProbability identicalByDescentProbability(1, 0, 0);
    map<string, vector<string> > locusToSuspectAlleles;
    map<string, vector<set<string> > > locusToAssumedAlleles;
    map<string, vector<set<string> > > locusToUnattributedAlleles;
    map<string, double> locusSpecificDropout;
    map<string, double> locusSpecificDropin;

    set<string> lociToRun;

    RetrieveDataFromCSV(inputFileName, &alpha, &dropinRate, &dropoutRate, &fst, &race,
                        &identicalByDescentProbability, &locusToSuspectAlleles,
                        &locusToAssumedAlleles, &locusToUnattributedAlleles,
                        &locusSpecificDropout, &locusSpecificDropin, &lociToRun);

    const string alleleFrequencyTablePath = executablePath + "Allele Frequency Tables/";
    vector<Race> races = GetRaces(race, alleleFrequencyTablePath, lociToRun);

    map<Race, vector<double> > raceToSolverIndexToLogProb;
    map<Race, vector<map<string, double> > > raceToSolverIndexToLocusLogProb;

    for (unsigned int raceIndex = 0; raceIndex < races.size(); raceIndex++) {
        Race r = races[raceIndex];
        raceToSolverIndexToLogProb.insert(
                pair<Race, vector<double> >(r, vector<double>(likelihoodSolvers.size(), 0)));
        raceToSolverIndexToLocusLogProb.insert(
                pair<Race, vector<map<string, double> > >(r,
                        vector<map<string, double> >(likelihoodSolvers.size())));
    }

    // Create configurations and run on likelihood solvers.
    for (set<string>::const_iterator iter = lociToRun.begin();
            iter != lociToRun.end(); iter++) {
        string locus = *iter;
        // const string alleleFrequencyTableFileName =
        //         alleleFrequencyTablePath + locus + "_B.count.csv";
        const string alleleFrequencyTableFileName =
                GetAlleleFrequencyTableFileName(alleleFrequencyTablePath, locus);
        const vector<string>& suspectAlleles = locusToSuspectAlleles[locus];
	const vector<set<string> >& assumedAlleles = locusToAssumedAlleles[locus];
	const vector<set<string> >& unattributedAlleles = locusToUnattributedAlleles[locus];
        map<Race, map<string, double> > raceToAlleleProportions =
                GetRaceToAlleleProportions(alleleFrequencyTableFileName, suspectAlleles,
                                           assumedAlleles, unattributedAlleles, fst);
        for (unsigned int raceIndex = 0; raceIndex < races.size(); raceIndex++) {
            Race curRace = races[raceIndex];
	    const map<string, double>& alleleProportions = raceToAlleleProportions[curRace];
            Configuration config(suspectAlleles, assumedAlleles, unattributedAlleles,
                                 alleleProportions, identicalByDescentProbability,
                                 GetValueOrDefault(locusSpecificDropout, locus, dropoutRate),
                                 GetValueOrDefault(locusSpecificDropin, locus, dropinRate), alpha);
            for (unsigned int solverIndex = 0; solverIndex < likelihoodSolvers.size();
                 solverIndex++) {
                LikelihoodSolver* solver = likelihoodSolvers[solverIndex];
                double logLikelihood = solver->getLogLikelihood(config);
                raceToSolverIndexToLocusLogProb[curRace][solverIndex][locus] = logLikelihood;
                raceToSolverIndexToLogProb[curRace][solverIndex] += logLikelihood;
            }
        }
    }

    // TODO: Output in proper format to file:
    outputData(lociToRun, likelihoodSolvers, raceToSolverIndexToLocusLogProb,
            raceToSolverIndexToLogProb, races, outputFileName);
    return raceToSolverIndexToLogProb;
}

int main(int argc, char *argv[]) {
    vector<LikelihoodSolver*> solversToUse;

    if (argc < 2) {
        std::cout << "Usage is <inputfile> <outputfile> [<[su][0123]>, ...]\n";
        std::cout << "first character is with [s]uspect or [u]nknowns only, second is number of unknowns.\n";
        std::cout << "\nexample: for no suspect, one unknown and one suspect, one unknown\n";
        std::cout << argv[0] << " input.csv output.csv 01 11\n\n";
        return -1;
    }
    if (argc > 2) {
        for (int i = 3; i < argc; i++) {
            if (argv[i][0] == 'u') {
                if (argv[i][1] <= '9' && argv[i][1] >= '0') {
                    unsigned int numUnknowns = argv[i][1] - '0';
                    // TODO fix memory leak
                    LikelihoodSolver* solver = new UnknownsSolver(numUnknowns * 2);
                    solversToUse.push_back(solver);
                }
            } else if (argv[i][0] == 's') {
                if (argv[i][1] <= '9' && argv[i][1] >= '0') {
                    unsigned int numUnknowns = argv[i][1] - '0';
                    // TODO fix memory leak
                    LikelihoodSolver* solver = new SuspectUnknownsSolver(numUnknowns * 2);
                    solversToUse.push_back(solver);
                }
            } else if (argv[i][0] == '0') {
                switch(argv[i][1]) {
                case '4':
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::NO_SUSPECT_FOUR_UNKNOWNS));
                    break;
                case '3':
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::NO_SUSPECT_THREE_UNKNOWNS));
                    break;
                case '2':
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::NO_SUSPECT_TWO_UNKNOWNS));
                    break;
                default:
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::NO_SUSPECT_ONE_UNKNOWN));
                }
            } else { // one suspect
                switch(argv[i][1]) {
                case '3':
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::ONE_SUSPECT_THREE_UNKNOWNS));
                    break;
                case '2':
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::ONE_SUSPECT_TWO_UNKNOWNS));
                    break;
                case '1':
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::ONE_SUSPECT_ONE_UNKNOWN));
                    break;
                default:
                    solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::ONE_SUSPECT_NO_UNKNOWNS));
                }
            }
        }
    } else {
        solversToUse.push_back(LikelihoodSolver::getSolver(LikelihoodSolver::NO_SUSPECT_ONE_UNKNOWN));
    }

    string executablePath = string(argv[0]);

    // The executable file name is 'labr' but on windows is 'labr.exe'
    int fileNameLength = 4;
    if (executablePath.find(".exe") != string::npos) {
        fileNameLength += 4;
    }

    run(executablePath.substr(0, executablePath.length() - fileNameLength), argv[1], argv[2], solversToUse);
}
