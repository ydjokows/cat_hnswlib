#pragma once

#include "types.h"
#include "condition.h"
#include "element_container.h"
#include <mutex>
#include <string.h>
#include <vector>
#include <set>

namespace hnswlib {

    class TagsStore : public ElementContainer<tagcontainer> {
    public:

        std::unordered_map<tagtype, std::vector<tableint> > tag_mapping;

        tagcontainer *getTags(tableint idx) {
            return &data[idx];
        }

        void reset()
        {
            _max_elements = 0;
            data.clear();
            tag_mapping.clear();
        }

        virtual void addElement(tableint idx, tagtype value) {
            ElementContainer::addElement(idx, value);
            tag_mapping[value].push_back(idx);
        }

        void addTag(tableint idx, tagtype tag) {
            addElement(idx, tag);
        }

        bool checkCondition(tableint id, const SearchCondition &condition) const {
            return condition.checkCondition(data[id]);
        }
    };
}
