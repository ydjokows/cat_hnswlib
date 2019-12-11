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

        tagcontainer *getTags(tableint idx) {
            return &data[idx];
        }

        void addTag(tableint id, tagtype tag) {
            tagcontainer *element_tags = &data[id];
            element_tags->insert(element_tags->end(), tag);
        }

        void setTags(tableint id, tagcontainer &new_element_tags) {
            tagcontainer *element_tags = &data[id];
            element_tags->clear();
            for(auto &tag: new_element_tags) {
                element_tags->insert(element_tags->end(), tag);
            }
        }

        bool checkCondition(tableint id, const SearchCondition &condition) const {
            return condition.checkCondition(data[id]);
        }
    };
}
