#pragma once

#include "types.h"
#include <mutex>
#include <string.h>
#include <vector>
#include <set>

namespace hnswlib
{
class SearchCondition
{
    condition_t _conditions;

public: 
    SearchCondition(condition_t &conditions){
        _conditions = conditions;
    }

    /**
     * Returns `true` if element contains tag
     * 
     */
    inline bool checkTag(tagtype tag, const tagcontainer &element_tags) const {
        auto search = element_tags.find(tag);
        return search != element_tags.end();
    }


    /**
     * Checks if element passes condition.
     * Returns `true` if passed.
     * 
     */
    inline bool checkCondition(const tagcontainer &element_tags) const {
        if(_conditions.empty())
            return true;
        bool result = true;
        for(const or_condition_t &or_condition: _conditions) {
            bool or_result = false;
            for(const negable_tag_t& negable_tag : or_condition) {
                bool is_contain = checkTag(negable_tag.second, element_tags);
                bool check = negable_tag.first ? !is_contain : is_contain;
                or_result |= check;
            }
            result &= or_result;
        }
        return result;
    }

    /**
     * Should return list of possible categories where entry points could be found.
     */
    std::vector<tagtype> entryCandidates() {
        std::vector<tagtype> candidates;
        for(or_condition_t &or_condition: _conditions) {
            for(const negable_tag_t& negable_tag : or_condition) {
                if (!negable_tag.first) {
                    candidates.push_back(negable_tag.second);
                }
            }
        }
        return candidates;
    }
};
} // namespace hnswlib