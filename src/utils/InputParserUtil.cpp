//==================================================================================================
// Name        : InputParserUtil.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include "InputParserUtil.h"

#include <cstdlib>

#include "ProbabilityUtil.h"
#include "StringUtil.h"

namespace LabRetriever {

void RetrieveDataFromCSV(const string& inputFileName, double* alpha, double* dropinRate,
                         double* dropoutRate, double* fst, Race* race,
                         IdenticalByDescentProbability* identicalByDescentProbability,
                         map<string, vector<string> >* locusToSuspectAlleles,
                         map<string, vector<set<string> > >* locusToAssumedAlleles,
                         map<string, vector<set<string> > >* locusToUnattributedAlleles,
                         map<string, double>* locusSpecificDropout,
                         map<string, double>* locusSpecificDropin, set<string>* lociToRun) {
    vector<vector<string> > inputData = readRawCsv(inputFileName);
    unsigned int csvIndex = 0;
    for (; csvIndex < inputData.size(); csvIndex++) {
        const vector<string>& row = inputData[csvIndex];
        if (row.size() == 0) continue;

        const string& header = row[0];
        if (header == "alpha") {
            if (row.size() <= 1) continue;
	    double value;
	    if (!ToDouble(row[1], &value)) continue;
            *alpha = value;
        } else if (header == "Drop-in rate") {
            if (row.size() <= 1) continue;
	    double value;
            if (!ToDouble(row[1], &value)) continue;
            *dropinRate = value;
        } else if (header == "Drop-out rate") {
            if (row.size() <= 1) continue;
            double value;
            if (!ToDouble(row[1], &value)) continue;
            *dropoutRate = value;
        } else if (header == "fst") {
            if (row.size() <= 1) continue;
            double value;
            if (!ToDouble(row[1], &value)) continue;
            *fst = value;
        } else if (header == "Race") {
            if (row.size() <= 1) continue;
            // Currently only expect one race or ALL. Change this to vector if multiple races
            // are needed but not all.
            *race = row[1];
        } else if (header == "IBD Probs") {
            if (row.size() <= 3 || row[1].size() == 0 || row[2].size() == 0 || row[3].size() == 0)
                continue;

            double zeroAllelesInCommonProb, oneAlleleInCommonProb, twoAllelesInCommonProb;
            if (!ToDouble(row[1], &zeroAllelesInCommonProb) ||
                !ToDouble(row[2], &oneAlleleInCommonProb) ||
                !ToDouble(row[3], &twoAllelesInCommonProb)) {
                continue;
            }

            identicalByDescentProbability->zeroAllelesInCommonProb = zeroAllelesInCommonProb;
            identicalByDescentProbability->oneAlleleInCommonProb = oneAlleleInCommonProb;
            identicalByDescentProbability->bothAllelesInCommonProb = twoAllelesInCommonProb;
        } else {
            size_t index = header.find("-");
            if (index == string::npos) continue;
            string locus = header.substr(0, index);
            string locusType = header.substr(index + 1, header.size());
            if (locusType == "Drop-in" || locusType == "Drop-out") {
                double value;
                if (!ToDouble(row[1], &value)) continue;
                if (locusType == "Drop-in") {
                    (*locusSpecificDropin)[locus] = value;
                } else if (locusType == "Drop-out") {
                    (*locusSpecificDropout)[locus] = value;
                }
                continue;
            }
            set<string> alleleSet;
            vector<string> alleles;
            for (unsigned int i = 1; i < row.size(); i++) {
                string data = row[i];
                if (data.length() != 0) {
                    alleles.push_back(data);
                    alleleSet.insert(data);
                }
            }
            if (locusType == "Assumed") {
                (*locusToAssumedAlleles)[locus].push_back(alleleSet);
            } else if (locusType == "Unattributed") {
                (*locusToUnattributedAlleles)[locus].push_back(alleleSet);
            } else if (locusType == "Suspected") {
                // If there are no suspected alleles, then there's no point of calculating this.
                if (alleles.size() == 0) continue;
                (*locusToSuspectAlleles)[locus] = alleles;
            }
        }
    }

    // Check for loci where there is a suspect allele, and the assumed and unattributed
    // rows are present (but possibly empty).
    for (map<string, vector<string> >::const_iterator iter = locusToSuspectAlleles->begin();
         iter != locusToSuspectAlleles->end(); iter++) {
        const string& locus = iter->first;
        if (locusToAssumedAlleles->find(locus) != locusToAssumedAlleles->end() &&
            locusToUnattributedAlleles->find(locus) != locusToUnattributedAlleles->end()) {
            lociToRun->insert(locus);
        }
    }

    // If the number of Assumed and Unattributed cases don't match, extend one of them to match the
    // other.
    for (set<string>::const_iterator iter = lociToRun->begin(); iter != lociToRun->end(); iter++) {
        string allele = *iter;
        vector<set<string> >& assumedAlleles = (*locusToAssumedAlleles)[allele];
        vector<set<string> >& unattributedAlleles = (*locusToUnattributedAlleles)[allele];
        int len = max(assumedAlleles.size(), unattributedAlleles.size());
        assumedAlleles.resize(len);
        unattributedAlleles.resize(len);
    }
}

vector<Race> GetRaces(const Race& race, const string& alleleFrequencyTablePath,
                      const set<string>& lociToRun) {
    vector<Race> races;
    // Determine the races to run.
    if (race == ALL_RACE) {
        set<Race> allRaces;
        for (set<string>::const_iterator iter = lociToRun.begin(); iter != lociToRun.end();
             iter++) {
            string locus = *iter;
            map<Race, map<string, unsigned int> > raceToAlleleCounts =
                getAlleleCountsFromFile(alleleFrequencyTablePath + locus + "_B.count.csv");
            for (map<Race, map<string, unsigned int> >::const_iterator race_iter =
                     raceToAlleleCounts.begin();
                 race_iter != raceToAlleleCounts.end(); race_iter++) {
                allRaces.insert(race_iter->first);
            }
        }
        for (set<Race>::const_iterator iter = allRaces.begin(); iter != allRaces.end(); iter++) {
            races.push_back(*iter);
        }
    } else {
        races.push_back(race);
    }
    return races;
}

Configuration CreateConfiguration(
        const string& alleleFrequencyTablePath, const string& locus, const Race& race, double alpha,
        double dropinRate, double dropoutRate, double fst,
        const IdenticalByDescentProbability& identicalByDescentProbability,
        const map<string, vector<string> >& locusToSuspectAlleles,
        const map<string, vector<set<string> > >& locusToAssumedAlleles,
        const map<string, vector<set<string> > >& locusToUnattributedAlleles,
        const map<string, double>& locusSpecificDropout,
        const map<string, double>& locusSpecificDropin) {
    const vector<set<string> >& unattributedAlleles =
            locusToUnattributedAlleles.find(locus)->second;
    const vector<set<string> >& assumedAlleles = locusToAssumedAlleles.find(locus)->second;
    const vector<string>& suspectAlleles = locusToSuspectAlleles.find(locus)->second;

    // TODO: Consider building these data structures outside of this function and pass them into
    // this function.
    set<string> allAlleles;
    allAlleles.insert(suspectAlleles.begin(), suspectAlleles.end());
    for (unsigned int i = 0; i < unattributedAlleles.size(); i++) {
        allAlleles.insert(unattributedAlleles[i].begin(), unattributedAlleles[i].end());
    }
    for (unsigned int i = 0; i < assumedAlleles.size(); i++) {
        allAlleles.insert(assumedAlleles[i].begin(), assumedAlleles[i].end());
    }

    AlleleProfile suspectProfile;
    vector<ReplicateData> replicateDatas;

    // If only one suspect allele is present, then assume it means that there are two of them.
    // TODO: fix quick hack for this:
    if (suspectAlleles.size() == 1) {
        suspectProfile.addAllele(suspectAlleles[0]);
    }

    for (unsigned int i = 0; i < suspectAlleles.size(); i++) {
        suspectProfile.addAllele(suspectAlleles[i]);
    }
    for (unsigned int index = 0; index < unattributedAlleles.size(); index++) {
        replicateDatas.push_back(ReplicateData::fromUnattributedAndMaskedAlleles(
                unattributedAlleles[index], assumedAlleles[index]));
    }

    map<Race, map<string, unsigned int> > raceToAlleleCounts =
            getAlleleCountsFromFile(alleleFrequencyTablePath + locus + "_B.count.csv");

    map<string, unsigned int> alleleCounts = raceToAlleleCounts[race];

    // Edit the allele counts so that every allele is given at least 5 counts, even if the
    // allele does not appear in the table.
    for (set<string>::const_iterator iter = allAlleles.begin(); iter != allAlleles.end(); iter++) {
        const string& allele = *iter;
        if (alleleCounts.count(allele) == 0 || alleleCounts[allele] < 5) {
            alleleCounts[allele] = 5;
        }
    }

    map<string, double> alleleProp =
            getAlleleProportionsFromCounts(alleleCounts, suspectProfile, fst);

    Configuration config(
            suspectProfile, replicateDatas, alleleProp, identicalByDescentProbability,
            locusSpecificDropout.count(locus) == 0 ? dropoutRate
                                                   : locusSpecificDropout.find(locus)->second,
            locusSpecificDropin.count(locus) == 0 ? dropinRate
                                                  : locusSpecificDropin.find(locus)->second,
            alpha);
    return config;
}

}  // namespace LabRetriever
