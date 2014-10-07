//==================================================================================================
// Name        : Configuration.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include "Configuration.h"
#include <algorithm>
#include <cassert>


namespace LabRetriever {
    /**********************************************************************************************
     * SuspectProfile
     **********************************************************************************************/
    AlleleProfile& AlleleProfile::addAllele(const string& allele) {
        if (alleleCounts.count(allele) == 0) {
            alleleCounts[allele] = 1;
        } else {
            alleleCounts[allele]++;
        }
        alleles.insert(allele);
        return *this;
    }

    bool AlleleProfile::contains(const string& allele) const {
        return alleleCounts.count(allele) != 0;
    }

    const map<string, unsigned int>& AlleleProfile::getAlleleCounts() const {
        return alleleCounts;
    }

    unsigned int AlleleProfile::getAlleleCounts(const string& allele) const {
        map<string, unsigned int>::const_iterator iter = alleleCounts.find(allele);
        if (iter == alleleCounts.end()) {
            return 0;
        }
        return iter->second;
    }

    const set<string>& AlleleProfile::getAlleles() const {
        return alleles;
    }

    /**********************************************************************************************
     * ReplicateData
     **********************************************************************************************/
    ReplicateData ReplicateData::fromUnattributedAndMaskedAlleles(
            const set<string>& unattributedAlleles, const set<string>& maskedAlleles) {
        // Check that the intersection of the two vectors are empty.
        // Choose the smaller set for efficiency (though in practice, it shouldn't matter too much)
        const set<string> *smallerSet, *largerSet;
        if (unattributedAlleles.size() < maskedAlleles.size()) {
            smallerSet = &unattributedAlleles;
            largerSet = &maskedAlleles;
        } else {
            smallerSet = &maskedAlleles;
            largerSet = &unattributedAlleles;
        }

        for (set<string>::const_iterator iter = smallerSet->begin(); iter != smallerSet->begin();
                iter++) {
            const string& allele = *iter;
            assert(largerSet->count(allele) == 0);
        }
        return ReplicateData(unattributedAlleles, maskedAlleles);
    }

    ReplicateData fromDetectedAndMaskedAlleles(const set<string>& detectedAlleles,
            const set<string>& maskedAlleles) {
        // TODO: implement this.
        throw "fromDetectedAndMaskedAlleles not yet implemented!";
    }

    ReplicateData::ReplicateData(const set<string>& unattributedAlleles,
            const set<string>& maskedAlleles) : unattributedAlleles(unattributedAlleles),
                    maskedAlleles(maskedAlleles) {};

    /**********************************************************************************************
     * UnknownSharedAlleleProbability
     **********************************************************************************************/
    IdenticalByDescentProbability::IdenticalByDescentProbability(double oneAlleleInCommonProb,
        double bothAllelesInCommonProb) :
            zeroAllelesInCommonProb(1 - oneAlleleInCommonProb - bothAllelesInCommonProb),
            oneAlleleInCommonProb(oneAlleleInCommonProb),
            bothAllelesInCommonProb(bothAllelesInCommonProb) {
        assert(zeroAllelesInCommonProb >= 0 &&
                oneAlleleInCommonProb >= 0 &&
                bothAllelesInCommonProb >= 0);
    }

    IdenticalByDescentProbability::IdenticalByDescentProbability(double zeroAllelesInCommonProb,
        double oneAlleleInCommonProb, double bothAllelesInCommonProb) :
            zeroAllelesInCommonProb(zeroAllelesInCommonProb),
            oneAlleleInCommonProb(oneAlleleInCommonProb),
            bothAllelesInCommonProb(bothAllelesInCommonProb) {
        assert(zeroAllelesInCommonProb >= 0 &&
                oneAlleleInCommonProb >= 0 &&
                bothAllelesInCommonProb >= 0);
        double shouldBeCloseToZero = zeroAllelesInCommonProb + oneAlleleInCommonProb +
                bothAllelesInCommonProb - 1;
        if (shouldBeCloseToZero < 0) {
            shouldBeCloseToZero = -shouldBeCloseToZero;
        }
        assert(shouldBeCloseToZero < 1.0e-5);
    }

    /**********************************************************************************************
     * Configuration
     **********************************************************************************************/
    Configuration& Configuration::setSuspectProfile(const AlleleProfile& suspectProfile) {
        this->suspectProfile = suspectProfile;
        return *this;
    }

    Configuration& Configuration::setData(const vector<ReplicateData>& data) {
        this->data = data;
        return *this;
    }

    Configuration& Configuration::setAlleleProportions(
            const map<string, double>& alleleProportions) {
        this->alleleProportions = alleleProportions;
        return *this;
    }

    Configuration& Configuration::setDropoutRate(double dropoutRate) {
        this->dropoutRate = dropoutRate;
        return *this;
    }

    Configuration& Configuration::setDropinRate(double dropinRate) {
        this->dropinRate = dropinRate;
        return *this;
    }

    Configuration& Configuration::setAlpha(double alpha) {
        this->alpha = alpha;
        return *this;
    }

} /* namespace LabRetriever */
