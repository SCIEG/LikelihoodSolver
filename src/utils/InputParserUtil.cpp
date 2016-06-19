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

#include "StringUtil.h"

using namespace std;

namespace LabRetriever {

string GetAlleleFrequencyTableFileName(const string& alleleFrequencyTablePath, const string locus) {
    const string path = (alleleFrequencyTablePath.size() > 0 &&
                         (alleleFrequencyTablePath[alleleFrequencyTablePath.size() - 1] != '/' ||
                          alleleFrequencyTablePath[alleleFrequencyTablePath.size() - 1] != '\\'))
                                ? alleleFrequencyTablePath + "/"
                                : alleleFrequencyTablePath;
    return path + locus + "_B.count.csv";
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
            map<Race, map<string, unsigned int> > raceToAlleleCounts = getAlleleCountsFromFile(
                    GetAlleleFrequencyTableFileName(alleleFrequencyTablePath, locus));
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

namespace {

set<string> AllAlleles(const vector<string>& suspectAlleles,
                       const vector<set<string> >& assumedAlleles,
                       const vector<set<string> >& unattributedAlleles) {
    set<string> allAlleles;
    allAlleles.insert(suspectAlleles.begin(), suspectAlleles.end());
    for (unsigned int i = 0; i < unattributedAlleles.size(); i++) {
        allAlleles.insert(unattributedAlleles[i].begin(), unattributedAlleles[i].end());
    }
    for (unsigned int i = 0; i < assumedAlleles.size(); i++) {
        allAlleles.insert(assumedAlleles[i].begin(), assumedAlleles[i].end());
    }
    return allAlleles;
}

/*
 * Returns a mapping from alleles to the proportion that the allele appears in people, after
 * sampling adjustment and the Balding-Nichols FST correction.
 *
 * alleleCounts - a mapping from alleles to a count of how many times that allele appeared in
 *     the sample.
 * suspectProfile - the current suspect's profile.
 * samplingAdjustment -
 * fst - See Balding-Nichols FST correction. Balding recommends that fst should be at least
 *     0.01. It can be as high as 0.05, depending on how representative the sample population
 *     is of the whole target population.
 */
map<string, double> getAlleleProportionsFromCounts(const map<string, unsigned int>& alleleCounts,
                                                   const AlleleProfile& suspectProfile, double fst,
                                                   unsigned int samplingAdjustment = 2) {
    const map<string, unsigned int>& suspectAlleleCounts = suspectProfile.getAlleleCounts();
    double totalCounts = 0, fstCorrection = (1 - fst) / (1 + fst);
    for (map<string, unsigned int>::const_iterator iter = alleleCounts.begin();
         iter != alleleCounts.end(); iter++) {
        totalCounts += iter->second;
    }

    double totalSuspectAlleleCounts = 0;
    for (map<string, unsigned int>::const_iterator iter = suspectAlleleCounts.begin();
         iter != suspectAlleleCounts.end(); iter++) {
        totalSuspectAlleleCounts += iter->second;
    }

    // // Find the total counts in the alleleCounts, plus extra per suspect allele.
    // totalCounts += totalSuspectAlleleCounts * samplingAdjustment;
    double multiplierCorrection = fstCorrection / totalCounts;
    map<string, double> alleleProportions;
    for (map<string, unsigned int>::const_iterator iter = alleleCounts.begin();
         iter != alleleCounts.end(); iter++) {
        const string& allele = iter->first;
        alleleProportions[allele] = iter->second * multiplierCorrection;
    }

    // For each suspect allele, add in the proportions that the adjustments would have added.
    for (map<string, unsigned int>::const_iterator iter = suspectAlleleCounts.begin();
         iter != suspectAlleleCounts.end(); iter++) {
        const string& allele = iter->first;
        unsigned int count = iter->second;
        alleleProportions[allele] += count * (fst / (1 + fst));
    }
    return alleleProportions;
}

}  // namespace

map<Race, map<string, double> > GetRaceToAlleleProportions(
        const string& alleleFrequencyTableFileName, const vector<string>& suspectAlleles,
        const vector<set<string> >& assumedAlleles, const vector<set<string> >& unattributedAlleles,
        double fst) {
  const set<string> allAlleles = AllAlleles(suspectAlleles, assumedAlleles, unattributedAlleles);
    map<Race, map<string, unsigned int> > raceToAlleleCounts =
            getAlleleCountsFromFile(alleleFrequencyTableFileName);

    map<Race, map<string, double> > raceToAlleleProp;
    for (map<Race, map<string, unsigned int> >::iterator iter = raceToAlleleCounts.begin();
         iter != raceToAlleleCounts.end(); ++iter) {
        const Race& race = iter->first;
        map<string, unsigned int>& alleleCounts = iter->second;
        // Edit the allele counts so that every allele is given at least 5 counts, even if the
        // allele does not appear in the table.
        for (set<string>::const_iterator allele_iter = allAlleles.begin();
             allele_iter != allAlleles.end(); allele_iter++) {
            const string& allele = *allele_iter;
            if (alleleCounts.count(allele) == 0 || alleleCounts[allele] < 5) {
                alleleCounts[allele] = 5;
            }
        }
        raceToAlleleProp[race] =
                getAlleleProportionsFromCounts(alleleCounts, AlleleProfile(suspectAlleles), fst);
    }
    return raceToAlleleProp;
}

}  // namespace LabRetriever
