//==================================================================================================
// Name        : sampler-main.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description :
//==================================================================================================

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "Configuration.h"
#include "LikelihoodSolver/UnknownsSolver.h"
#include "utils/FileReaderUtil.h"
#include "utils/InputParserUtil.h"
#include "utils/ProbabilityUtil.h"
#include "utils/StlUtil.h"
#include "utils/StringUtil.h"

class Person {
   public:
    void setLocusAlleles(const std::string& locus, const std::string& allele1,
                         const std::string& allele2) {
        std::pair<std::string, std::string>& p = locusToAlleles[locus];
	bool oneSmallerThanTwo = allele1 < allele2;
        p.first = oneSmallerThanTwo ? allele1 : allele2;
        p.second = oneSmallerThanTwo ? allele2 : allele1;
    }

    std::string ToString() const {
        if (locusToAlleles.empty()) return "";
        std::stringstream ss;
        std::map<std::string, std::pair<std::string, std::string> >::const_iterator iter =
                locusToAlleles.begin();
        appendString(iter->first, iter->second, &ss);
        ++iter;
        for (; iter != locusToAlleles.end(); ++iter) {
            ss << ";";
            appendString(iter->first, iter->second, &ss);
        }
        return ss.str();
    }

   private:
    static void appendString(const std::string& locus,
                             const std::pair<std::string, std::string>& alleles,
                             std::stringstream* ss) {
        *ss << locus << ": " << alleles.first << ", " << alleles.second;
    }
    std::map<std::string, std::pair<std::string, std::string> > locusToAlleles;
};

bool cmp(const std::pair<double, Person>& left, const std::pair<double, Person>& right) {
  return left.first < right.first;
}

bool cmp2(double left, const std::pair<double, Person>& right) {
  return left < right.first;
}

int calculatePercentileIndex(double percentile, int total) {
    return static_cast<int>((total - 1) * percentile / 100);
}

// TODO: Add a --full flag.

int main(int argc, char* argv[]) {
    const std::string prog = argv[0];
    if (argc < 4 || argc > 5) {
        std::cerr << "Usage:" << std::endl
                  << "    " << prog << " input.csv N samples [SEED]" << std::endl
                  << std::endl
                  << " - N is the number of unknowns in the numerator of the LR. Thus, N+1"
                  << std::endl
                  << "   is the number of unknowns in the denominator." << std::endl
                  << " - samples is the number of suspects to sample." << std::endl
                  << " - SEED is optional, which controls the random sampling. If SEED is missing"
                  << std::endl
		  << "   then the seed is set using the current time. Otherwise if it is a number,"
		  << std::endl
		  << "   the seed will be that number."
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
		  << "times." << std::endl;
        return -1;
    }

    // The executable file name is 'sampler' but on windows is 'sampler.exe'
    int fileNameLength = 7;
    if (prog.find(".exe") != std::string::npos) {
        fileNameLength += 4;
    }
    const std::string executablePath = prog.substr(0, prog.length() - fileNameLength);
    const std::string alleleFrequencyTablePath = executablePath + "Allele Frequency Tables/";

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

    int seed;
    if (argc == 5) {
	if (!LabRetriever::ToInt(argv[4], &seed)) {
	    std::cerr<< "SEED must be a positive integer!" << std::endl;
	}
    } else {
	seed = time(NULL);
    }
    srand(seed);

    LabRetriever::SuspectUnknownsSolver numeratorSolver(2 * n /* numUnknownAlleles */);
    LabRetriever::UnknownsSolver denominatorSolver(2 * (n + 1) /* numUnknownAlleles */);

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
            LabRetriever::GetRaces(race, alleleFrequencyTablePath, lociToRun);

    std::map<LabRetriever::Race, std::map<std::string, std::map<std::string, double> > >
            raceAndLocusToAlleleProportions;

    for (std::set<std::string>::const_iterator locusIter = lociToRun.begin();
         locusIter != lociToRun.end(); ++locusIter) {
        const std::string locus = *locusIter;
        const std::map<LabRetriever::Race, std::map<std::string, double> > raceToAlleleProportions =
                LabRetriever::GetRaceToAlleleProportions(
                        LabRetriever::GetAlleleFrequencyTableFileName(alleleFrequencyTablePath, locus),
                        GetValueOrDie(locusToSuspectAlleles, locus),
                        GetValueOrDie(locusToAssumedAlleles, locus),
                        GetValueOrDie(locusToUnattributedAlleles, locus), fst);
        for (std::map<LabRetriever::Race, std::map<std::string, double> >::const_iterator iter =
                     raceToAlleleProportions.begin();
             iter != raceToAlleleProportions.end(); ++iter) {
            const LabRetriever::Race race = iter->first;
            raceAndLocusToAlleleProportions[race][locus] = iter->second;
        }
    }

    // For each race
    for (std::vector<LabRetriever::Race>::const_iterator raceIter = races.begin();
         raceIter != races.end(); ++raceIter) {
        const LabRetriever::Race race = *raceIter;
        const std::map<std::string, std::map<std::string, double> >& locusToAlleleProportions =
                GetValueOrDie(raceAndLocusToAlleleProportions, race);
        double logSuspectNumer = 0;
        double logDenom = 0;
	//  Calculate the suspect's likelihood
	//  Calculate the denominator's likelihood
        for (std::set<std::string>::const_iterator locusIter = lociToRun.begin();
             locusIter != lociToRun.end(); ++locusIter) {
            const std::string locus = *locusIter;
            LabRetriever::Configuration config(
                    GetValueOrDie(locusToSuspectAlleles, locus),
                    GetValueOrDie(locusToAssumedAlleles, locus),
                    GetValueOrDie(locusToUnattributedAlleles, locus),
                    GetValueOrDie(locusToAlleleProportions, locus),
                    identicalByDescentProbability,
                    GetValueOrDefault(locusSpecificDropout, locus, dropoutRate),
                    GetValueOrDefault(locusSpecificDropin, locus, dropinRate), alpha);
            logSuspectNumer += numeratorSolver.getLogLikelihood(config);
            logDenom += denominatorSolver.getLogLikelihood(config);
        }
	const double logSuspectLR = logSuspectNumer - logDenom;
	const double suspectLR = exp(logSuspectLR);

	// Note: This can be put inside the above loop for efficiency, but is separate for clarity.
	std::map<std::string, LabRetriever::DistributionSampler<std::string> > locusToAlleleSampler;
        for (std::set<std::string>::const_iterator locusIter = lociToRun.begin();
             locusIter != lociToRun.end(); ++locusIter) {
	  const string& locus = *locusIter;
          locusToAlleleSampler.insert(
                  std::pair<std::string, LabRetriever::DistributionSampler<std::string> >(
                          locus, LabRetriever::DistributionSampler<std::string>(
                                         GetValueOrDie(locusToAlleleProportions, locus))));
        }

        //  for 1 to number of samples
        std::vector<std::pair<double, Person> > sampledPeople;
        sampledPeople.reserve(samples);
        for (int i = 0; i < samples; ++i) {
            //   generate and store a person's full profile
            Person p;
            //   calculate likelihood
            double logSampleSuspectNumer = 0;
            for (std::set<std::string>::const_iterator locusIter = lociToRun.begin();
                 locusIter != lociToRun.end(); ++locusIter) {
                const std::string& locus = *locusIter;
                const LabRetriever::DistributionSampler<std::string>& sampler =
                        GetValueOrDie(locusToAlleleSampler, locus);
                const std::string& allele1 = sampler.sample();
                const std::string& allele2 = sampler.sample();
                p.setLocusAlleles(locus, allele1, allele2);

                std::vector<std::string> sampledAlleles;
                sampledAlleles.reserve(2);
                sampledAlleles.push_back(allele1);
                sampledAlleles.push_back(allele2);
                LabRetriever::Configuration config(
                        sampledAlleles, GetValueOrDie(locusToAssumedAlleles, locus),
                        GetValueOrDie(locusToUnattributedAlleles, locus),
                        GetValueOrDie(locusToAlleleProportions, locus),
                        identicalByDescentProbability,
                        GetValueOrDefault(locusSpecificDropout, locus, dropoutRate),
                        GetValueOrDefault(locusSpecificDropin, locus, dropinRate), alpha);
                logSampleSuspectNumer += numeratorSolver.getLogLikelihood(config);
            }
            sampledPeople.push_back(
                    std::pair<double, Person>(exp(logSampleSuspectNumer - logDenom), p));
        }
        //  sort samples by likelihood
	std::sort(sampledPeople.begin(), sampledPeople.end(), cmp);
	//  compute statistics and print samples with those statistics
        // suspect percentile
        std::vector<std::pair<double, Person> >::iterator percentile_iter =
                std::upper_bound(sampledPeople.begin(), sampledPeople.end(), suspectLR, cmp2);
        const double percentile = (100 * std::distance(sampledPeople.begin(), percentile_iter)) /
                                  sampledPeople.size();
        // min
        const double min = sampledPeople.front().first;
        // 25%
        const double pc25 = sampledPeople[calculatePercentileIndex(25, sampledPeople.size())].first;
        // 50%
        const double pc50 = sampledPeople[calculatePercentileIndex(50, sampledPeople.size())].first;
        // 75%
        const double pc75 = sampledPeople[calculatePercentileIndex(75, sampledPeople.size())].first;
        // max
        const double max = sampledPeople.back().first;

        std::cout << "Race: " << race << std::endl
                  << std::endl
                  << "Suspect LR: " << suspectLR << std::endl
                  << "Suspect Percentile: " << percentile << "%" << std::endl
                  << "Min: " << min << std::endl
                  << "25%: " << pc25 << std::endl
                  << "50%: " << pc50 << std::endl
                  << "75%: " << pc75 << std::endl
                  << "max: " << max << std::endl
                  << std::endl;
        // Top 5 profiles
        const int num_profiles_to_show = (5 < sampledPeople.size()) ? 5 : sampledPeople.size();
        std::cout << "Top " << num_profiles_to_show << " Sampled Profile(s):" << std::endl;
        for (int i = 0; i < num_profiles_to_show; ++i) {
            const int index = sampledPeople.size() - i - 1;
            const std::pair<double, Person>& val = sampledPeople[index];
            std::cout << "LR: " << val.first << "; Alleles: " << val.second.ToString() << std::endl;
        }
	std::cout << endl;
    }
    return 0;
}
