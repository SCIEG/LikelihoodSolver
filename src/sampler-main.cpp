//==================================================================================================
// Name        : sampler-main.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Configuration.h"
#include "LikelihoodSolver/UnknownsSolver.h"
#include "utils/FileReaderUtil.h"
#include "utils/InputParserUtil.h"
#include "utils/StringUtil.h"

double GetDemonimatorLogLikelihood(int numUnknowns, const LabRetriever::Configuration& config) {
  LabRetriever::UnknownsSolver solver(2 * numUnknowns /* numUnknownAlleles */);
  return solver.getLogLikelihood(config);
}

// TODO: Add a --full flag.

int main(int argc, char* argv[]) {
    const std::string prog = argv[0];
    if (argc != 4) {
        std::cerr << "Usage:" << std::endl
                  << "    " << prog << " input.csv N samples" << std::endl
                  << std::endl
                  << " - N is the number of unknowns in the numerator of the LR. Thus, N+1"
                  << std::endl
                  << "   is the number of unknowns in the denominator." << std::endl
                  << " - samples is the number of suspects to sample." << std::endl
		  << std::endl
		  << "This program calculates the likelihood ratio of the suspect (as denoted in"
		  << std::endl
		  << "the input file) and of a number of randomly sampled allele profiles, and"
		  << std::endl
		  << "compares the suspect's LR to the sampled distribution of LRs. The LR"
		  << std::endl
		  << "calculations are made under the same distribution (corrected for fst and"
		  << std::endl
		  << "rare alleles), and thus is not the same as manually calling labr multiple"
		  << std::endl
		  << "times.";
        return -1;
    }

    // The executable file name is 'sampler' but on windows is 'sampler.exe'
    int fileNameLength = 7;
    if (prog.find(".exe") != std::string::npos) {
        fileNameLength += 4;
    }
    const std::string executablePath = prog.substr(0, prog.length() - fileNameLength);

    const std::string inputFile = argv[1];

    int n;
    if (!LabRetriever::ToInt(argv[2], &n) || n < 0) {
      std::cerr << "N must be a non-negative integer!" << std::endl;
      return -1;
    }

    int samples;
    if (!LabRetriever::ToInt(argv[3], &samples) || samples <= 0) {
      std::cerr << "samples must be a positive integer!" << std::endl;
      return -1;
    }

    double alpha;
    double dropinRate;
    double dropoutRate;
    double fst;
    LabRetriever::Race race;
    LabRetriever::IdenticalByDescentProbability identicalByDescentProbability(1, 0, 0);
    std::map<std::string, std::vector<std::string> > locusToSuspectAlleles;
    std::map<std::string, std::vector<std::set<std::string> > > locusToAssumedAlleles;
    std::map<std::string, std::vector<std::set<std::string> > > locusToUnattributedAlleles;
    std::map<std::string, double> locusSpecificDropout;
    std::map<std::string, double> locusSpecificDropin;
    std::set<std::string> lociToRun;
    LabRetriever::RetrieveDataFromCSV(inputFile, &alpha, &dropinRate, &dropoutRate, &fst, &race,
                                      &identicalByDescentProbability, &locusToSuspectAlleles,
                                      &locusToAssumedAlleles, &locusToUnattributedAlleles,
                                      &locusSpecificDropout, &locusSpecificDropin, &lociToRun);

    const std::vector<LabRetriever::Race> races =
            LabRetriever::GetRaces(race, executablePath, lociToRun);

    const std::string alleleFrequencyTableFileName;
    
    return 0;
}
