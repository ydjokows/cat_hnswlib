#pragma once

#include "types.h"
#include "condition.h"
#include <mutex>
#include <string.h>
#include <vector>
#include <set>

namespace hnswlib {

    class TagsStore {
    public:
        std::vector<tagcontainer> tags = std::vector<tagcontainer>();
        size_t _max_elements = 0;

        TagsStore(size_t max_elements = 0) {
            tags.resize(max_elements);
            _max_elements = max_elements;
        }

        ~TagsStore() {
        }

        void resize(size_t new_size) {
            _max_elements = new_size;
            tags.resize(new_size);
        }

        void reset() {
            _max_elements = 0;
            tags.clear();
        }

        void serialize(std::ostream &out) {
            for(size_t i = 0; i < _max_elements; i++) {
                auto *element_tags = &tags[i];
                size_t length = element_tags->size();
                writeBinaryPOD(out, length);
                for(const tagtype &tag : *element_tags) {
                    writeBinaryPOD(out, tag);
                }
            }
        }

        void deserialize(std::istream &in) {
            for(size_t i = 0; i < _max_elements; i++) {
                auto *element_tags = &tags[i];
                size_t length = 0;
                readBinaryPOD(in, length);
                for(size_t j = 0; j < length; j++) {
                    tagtype tag = 0;
                    readBinaryPOD(in, tag);
                    element_tags->insert(element_tags->end(), tag);
                }
            }
        }

        tagcontainer getTags(tableint idx) {
            return tags[idx];
        }

        void addTag(tableint id, tagtype tag) {
            tagcontainer *element_tags = &tags[id];
            element_tags->insert(element_tags->end(), tag);
        }

        void setTags(tableint id, tagcontainer &new_element_tags) {
            tagcontainer *element_tags = &tags[id];
            element_tags->clear();
            for(auto &tag: new_element_tags) {
                element_tags->insert(element_tags->end(), tag);
            }
        }

        bool checkCondition(tableint id, const SearchCondition &condition) const {
            return condition.checkCondition(tags[id]);
        }
    };
}
