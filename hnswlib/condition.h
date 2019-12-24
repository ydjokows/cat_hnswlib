#pragma once

#include "types.h"
#include <mutex>
#include <string.h>
#include <vector>
#include <set>

namespace hnswlib
{

class ConditionOrClause {
public:
    tagcontainer required_tags;
    tagcontainer excluded_tags;

    ConditionOrClause(const or_condition_t &or_condition) {
        for(const negable_tag_t& negable_tag : or_condition) {
            if(negable_tag.first) {
                excluded_tags.insert(negable_tag.second);
            }else{
                required_tags.insert(negable_tag.second);
            }
        }
    }

    bool checkClause(const tagcontainer &element_tags) const {
        // If there is at least one intersection with required tags set, condition id true
        if (required_tags.size() > element_tags.size()) {
            for(tagtype tag : element_tags) {
                if(required_tags.find(tag) != required_tags.end()) return true;
            }
        } else {
            for(tagtype tag : required_tags) {
                if(element_tags.find(tag) != element_tags.end()) return true;
            }
        }

        // If element does not contain any excluded element, condition is true
        if (excluded_tags.size() > element_tags.size()) {
            for(tagtype tag : element_tags) {
                if(excluded_tags.find(tag) == excluded_tags.end()) return true;
            }
        } else {
            for(tagtype tag : excluded_tags) {
                if(element_tags.find(tag) == element_tags.end()) return true;
            }
        }
        
        return false;
    }
};


class SearchCondition
{
    condition_t _conditions;
    std::vector<ConditionOrClause> clauses;

public:
    SearchCondition(condition_t &conditions){
        _conditions = conditions;
        for(const or_condition_t &or_condition: _conditions) {
            clauses.push_back(ConditionOrClause(or_condition));
        }
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
        for(const ConditionOrClause &or_condition: clauses) {
            result &= or_condition.checkClause(element_tags);
            if(!result) break;
        }
        return result;
    }

    /**
     * Should return list of possible categories where entry points could be found.
     */
    std::vector<tagtype> entryCandidates() const {
        std::vector<tagtype> candidates;
        for(const or_condition_t &or_condition: _conditions) {
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