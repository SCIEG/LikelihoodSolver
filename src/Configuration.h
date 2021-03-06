//==================================================================================================
// Name        : Configuration.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : Contains classes that hold input data.
//==================================================================================================

#ifndef CONFIGURATION_H_
#define CONFIGURATION_H_

#include <map>
#include <set>
#include <string>
#include <vector>

namespace LabRetriever {
using namespace std;

/*
 * Represents the alleles belonging to the suspect.
 */
class AlleleProfile {
   public:
    AlleleProfile() {};
    AlleleProfile(const vector<string>& alleles);

    /* Returns a reference to this AlleleProfile, for chaining. */
    AlleleProfile& addAllele(const string& allele);

    bool contains(const string& allele) const;
    const map<string, unsigned int>& getAlleleCounts() const;
    unsigned int getAlleleCounts(const string& allele) const;
    const set<string>& getAlleles() const;

   private:
    map<string, unsigned int> alleleCounts;
    set<string> alleles;
};

/*
 * Represents data gathered from a replicate of LTDNA.
 */
class ReplicateData {
   public:
    set<string> unattributedAlleles;
    set<string> maskedAlleles;

    /*
     * Static factory methods.
     */
    static ReplicateData fromUnattributedAndMaskedAlleles(const set<string>& unattributedAlleles,
                                                          const set<string>& maskedAlleles);

   private:
    ReplicateData(const set<string>& unattributedAlleles, const set<string>& maskedAlleles);
};

struct IdenticalByDescentProbability {
   public:
    IdenticalByDescentProbability(double oneAlleleInCommonProb, double bothAllelesInCommonProb);
    IdenticalByDescentProbability(double zeroAllelesInCommonProb, double oneAlleleInCommonProb,
                                  double bothAllelesInCommonProb);
    double zeroAllelesInCommonProb;
    double oneAlleleInCommonProb;
    double bothAllelesInCommonProb;
};

/*
 * Contains data with which the solvers will use to calculate likelihoods.
 */
struct Configuration {
   public:
    AlleleProfile suspectProfile;
    vector<ReplicateData> data;
    map<string, double> alleleProportions;
    IdenticalByDescentProbability identicalByDescentProbability;
    double dropoutRate;
    double dropinRate;
    double alpha;

    Configuration(const AlleleProfile& suspectProfile, const vector<ReplicateData>& data,
                  const map<string, double>& alleleProportions,
                  const IdenticalByDescentProbability& identicalByDescentProbability,
                  double dropoutRate, double dropinRate, double alpha)
        : suspectProfile(suspectProfile),
          data(data),
          alleleProportions(alleleProportions),
          identicalByDescentProbability(identicalByDescentProbability),
          dropoutRate(dropoutRate),
          dropinRate(dropinRate),
          alpha(alpha) {};

    Configuration(const vector<string>& suspectAlleles, const vector<set<string> >& assumedAlleles,
                  const vector<set<string> >& unattributedAlleles,
                  const map<string, double>& alleleProportions,
                  const IdenticalByDescentProbability& identicalByDescentProbability,
                  double dropoutRate, double dropinRate, double alpha);

    Configuration& setSuspectProfile(const AlleleProfile& suspectProfile);
    Configuration& setData(const vector<ReplicateData>& data);
    Configuration& setAlleleProportions(const map<string, double>& alleleProportions);
    Configuration& setIdenticalByDescentProbability(const IdenticalByDescentProbability&);
    Configuration& setDropoutRate(double dropoutRate);
    Configuration& setDropinRate(double dropinRate);
    Configuration& setAlpha(double alpha);
};

} /* namespace LabRetriever */
#endif /* CONFIGURATION_H_ */
