//==================================================================================================
// Name        : FileReaderUtil.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#ifndef FILEREADERUTIL_H_
#define FILEREADERUTIL_H_

#include <map>
#include <string>
#include <vector>

using namespace std;

namespace LabRetriever {
    typedef string Race;
    extern const Race ALL_RACE;

    // TODO: Change name?
    map<Race, map<string, unsigned int> > getAlleleCountsFromFile(const string& fileName);

    vector< vector<string> > readRawCsv(const string& fileName);

    vector<string> makeTokenList(const string& line);

    string getToken(string& line);
}

#endif /* FILEREADERUTIL_H_ */
