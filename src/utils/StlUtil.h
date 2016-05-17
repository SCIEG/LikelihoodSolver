//==================================================================================================
// Name        : StlUtil.h
// Author      : Ken Cheng
// Copyright   : This work is licensed under the Creative Commons
//     Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this
//     license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
// Description : Contains utility functions for STL libraries.
//==================================================================================================

#include <map>

template <class K, class V>
const V& GetValueOrDie(const std::map<K, V>& m, const K& key) {
  return m.find(key)->second;
}

template <class K, class V>
const V& GetValueOrDefault(const std::map<K, V>& m, const K& key, const V& defaultValue) {
    typename std::map<K, V>::const_iterator iter = m.find(key);
    return iter == m.end() ? defaultValue : iter->second;
}
