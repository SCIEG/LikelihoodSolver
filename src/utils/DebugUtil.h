//==================================================================================================
// Name        : DebugUtil.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : 
//==================================================================================================

#ifndef DEBUGUTIL_H_
#define DEBUGUTIL_H_

#include "../Configuration.h"
#include <map>
#include <vector>
#include <set>

using namespace std;

namespace LabRetriever {
    template <class A>
    void debug(const A& a);

    template <class A, class B>
    void debugToken(const map<A, B>& m);
    template <class A>
    void debugToken(const vector<A>& v);
    template <class A>
    void debugToken(const set<A>& s);
    void debugToken(const AlleleProfile& a);
    void debugToken(const ReplicateData& r);
    void debugToken(const Configuration& c);

    // Default.
    template <class A>
    void debugToken(const A& a);
}

#endif /* DEBUGUTIL_H_ */
