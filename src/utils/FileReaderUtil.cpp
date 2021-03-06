//==================================================================================================
// Name        : FileReaderUtil.cpp
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#include "FileReaderUtil.h"
#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>

#include <boost/tokenizer.hpp>

using namespace std;

namespace LabRetriever {
    const Race ALL_RACE = "ALL";

    // trim from start
    static inline string &ltrim(string &s) {
            s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
            return s;
    }

    // trim from end
    static inline std::string &rtrim(std::string &s) {
            s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
            return s;
    }

    // trim from both ends
    static inline std::string &trim(std::string &s) {
            return ltrim(rtrim(s));
    }

    // TODO: optimize for 'ALL' case
    map<Race, map<string, unsigned int> > getAlleleCountsFromFile(const string& fileName) {
        map<Race, map<string, unsigned int> > retVal;
        vector< vector<string> > rawCsv = readRawCsv(fileName);

        // The header, starting at column 1, lists the races which the frequency tables are for.
        vector<string> header = rawCsv[0];

        for (unsigned int rowNum = 1; rowNum < rawCsv.size(); rowNum++) {
            vector<string> row = rawCsv[rowNum];
            string allele = row[0];

            for (unsigned int colNum = 1; colNum < row.size(); colNum++) {
                Race race = header[colNum];
                unsigned int val;
                istringstream(row[colNum]) >> val;
                retVal[race][allele] = val;
            }
        }
        return retVal;
    }

    vector< vector<string> > readRawCsv(const string& fileName) {
        ifstream file(fileName.c_str());

        if (!file.is_open()) {
            // TODO: Error here!
        }

        string line;
        vector< vector<string> > retVal;

        while (getline(file, line)) {
            vector<string> tokenList = makeTokenList(line);
            retVal.push_back(tokenList);
        }

        file.close();

        return retVal;
    }

    /*
     * Quick hack. Does not handle all csv formats.
     */
    vector<string> makeTokenList(const string& line) {
        vector<string> retVal;
        typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

        Tokenizer tok(line);
        for (Tokenizer::iterator iter = tok.begin(); iter != tok.end(); iter++) {
            string str = *iter;
            retVal.push_back(trim(str));
        }
//        while (true) {
//            // Trimming is not exactly correct, but there should not be whitespace.
//            string nextToken = trim(getToken(readInVal));
//            retVal.push_back(nextToken);

//            if (readInVal.size() == 0)
//                break;
//            else {
//                while (readInVal[0] == ' ') {
//                    readInVal.erase(0, 1);
//                }
//                readInVal.erase(0, 1);
//            }
//        }
//        return retVal;
//    }

//    /*
//     * Gets the next token (an element) from a csv string and advances the string to the next comma.
//     */
//    string getToken(string& line) {
//        while (line.size() != 0 && line[0] == ' ') {
//            line.erase(0, 1);
//        }

//        if (line.size() == 0) return "";

//        string retVal;
//        if (line[0] == '"') {
//            // Continue until you find an unescaped double quote.
//            unsigned int index = 1;
//            bool isDone = false;
//            for (; index < line.size() - 1; index++) {
//                char curChar = line[index];
//                char nextChar = line[index + 1];
//                switch (curChar) {
//                    case '\\':
//                        switch (nextChar) {
//                            case '"':
//                            case '\\':
//                                // If it's an escaped character, use the escaped character and
//                                // advance the index.
//                                retVal += nextChar;
//                                index++;
//                                break;
//                            default:
//                                retVal += '\\';
//                                break;
//                        }
//                        break;
//                    case '"':
//                        switch (nextChar) {
//                            case '"':
//                                retVal += '"';
//                                index++;
//                                break;
//                            case ',':
//                                // If it's not a double quote, then it must be the end.
//                                isDone = true;
//                                break;
//                            default:
//                                // Should not be here in a valid csv!
//                                assert(false);
//                                break;
//                        }
//                        break;
//                    default:
//                        retVal += curChar;
//                        break;
//                }

//                if (isDone) {
//                    break;
//                }
//            }

//            if (!isDone) {
//                // This must be the last token.
//                if (index != line.size() - 1 || line[index] != '"') {
//                    // If the index points past the end of the string or if the last character is
//                    // not a double quote, then it was an invalid parse.
//                    assert(false);
//                }
//            }

//            // At this point, index should point to the ending double quote.
//            line.erase(0, index + 1);
//        } else {
//            // If unquoted, then there shouldn't be a comma as an element. Just find the next comma
//            // and split there.
//            int index = line.find(',');
//            if (index == string::npos) {
//                index = line.size();
//            }
//            retVal = line.substr(0, index);
//            line.erase(0, index);
//        }

        return retVal;
    }
}
