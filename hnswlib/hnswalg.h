#pragma once

#include "visited_list_pool.h"
#include "hnswlib.h"
#include "node_tags.h"
#include "graph_layer.h"
#include <random>
#include <stdlib.h>
#include <unordered_set>
#include <list>
#include <iostream>

#include <fstream>


namespace hnswlib {
    typedef unsigned int linklistsizeint;

    template<typename dist_t>
    class HierarchicalNSW : public AlgorithmInterface<dist_t> {
    public:

        HierarchicalNSW(SpaceInterface<dist_t> *s) {

        }

        HierarchicalNSW(SpaceInterface<dist_t> *s, const std::string &location, bool nmslib = false, size_t max_elements=0) {
            loadIndex(location, s, max_elements);
        }

        HierarchicalNSW(SpaceInterface<dist_t> *s, size_t max_elements, size_t M = 16, size_t ef_construction = 200, size_t random_seed = 100, char *prefilled_data = nullptr) :
                link_list_locks_(max_elements) {
            max_elements_ = max_elements;

            has_deletions_=false;
            s_ = s;
            data_size_ = s->get_data_size();
            fstdistfunc_ = s->get_dist_func();
            dist_func_param_ = s->get_dist_func_param();
            M_ = M;
            maxM_ = M_;
            maxM0_ = M_ * 2;
            
            maxSearchM_ = maxM_;
            maxSearchM0_ = maxM0_ * 2;

            ef_construction_ = std::max(ef_construction,M_);
            ef_ = 10;
            random_seed_ = random_seed;
            level_generator_.seed(random_seed);

            size_data_per_element_ = data_size_ + sizeof(labeltype);
            label_offset_ = data_size_;

            if (prefilled_data == nullptr) 
            {   
                do_free = true;
                data_level0_memory_ = (char *) malloc(max_elements_ * size_data_per_element_);
                if (data_level0_memory_ == nullptr)
                    throw std::runtime_error("Not enough memory");
            } else {
                do_free = false;
                data_level0_memory_ = prefilled_data;
            }
            cur_element_count = 0;

            visited_list_pool_ = new VisitedListPool(1, max_elements);

            tags.resize(max_elements_);
            layer0.resize(max_elements_, maxM0_);
            layers.resize(max_elements_);

            //initializations for special treatment of the first node
            enterpoint_node_ = -1;
            maxlevel_ = -1;

            mult_ = 1 / log(1.0 * M_);
            revSize_ = 1.0 / mult_;
        }

        struct CompareByFirst {
            constexpr bool operator()(std::pair<dist_t, tableint> const &a,
                                      std::pair<dist_t, tableint> const &b) const noexcept {
                return a.first < b.first;
            }
        };

        ~HierarchicalNSW() {
            if (do_free)
                free(data_level0_memory_);
            
            delete visited_list_pool_;
        }

        SpaceInterface<dist_t> *s_;

        GraphLayer0 layer0;
        GraphLayers layers;
        TagsStore tags;

        bool do_free;
        size_t random_seed_;

        size_t max_elements_;
        size_t cur_element_count;
        size_t size_data_per_element_;

        size_t M_;
        size_t maxM_;
        size_t maxM0_;

        size_t maxSearchM_;
        size_t maxSearchM0_;

        size_t ef_construction_;

        double mult_, revSize_;
        int maxlevel_;


        VisitedListPool *visited_list_pool_;
        std::mutex cur_element_count_guard_;

        std::vector<std::mutex> link_list_locks_;
        tableint enterpoint_node_;

        char *data_level0_memory_;

        size_t data_size_;

        bool has_deletions_;


        size_t label_offset_;
        DISTFUNC<dist_t> fstdistfunc_;
        void *dist_func_param_;
        std::unordered_map<labeltype, tableint> label_lookup_;

        std::unordered_map<tagtype, tableint> tag_to_entrypoint_;

        std::default_random_engine level_generator_;

        inline labeltype getExternalLabel(tableint internal_id) const {
            labeltype return_label;
            memcpy(&return_label,(data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), sizeof(labeltype));
            return return_label;
        }

        inline void setExternalLabel(tableint internal_id, labeltype label) const {
            memcpy((data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), &label, sizeof(labeltype));
        }

        inline labeltype *getExternalLabeLp(tableint internal_id) const {
            return (labeltype *) (data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_);
        }

        inline char *getDataByInternalId(tableint internal_id) const {
            return (data_level0_memory_ + internal_id * size_data_per_element_);
        }

        size_t getRandomLevel(double reverse_size) {
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            double r = -log(distribution(level_generator_)) * reverse_size;
            return (size_t) r;
        }

        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
        searchBaseLayer(tableint ep_id, const void *data_point, size_t layer) {
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidateSet;

            dist_t lowerBound;
            if (!isMarkedDeleted(ep_id)) {
                dist_t dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
                top_candidates.emplace(dist, ep_id);
                lowerBound = dist;
                candidateSet.emplace(-dist, ep_id);
            } else {
                lowerBound = std::numeric_limits<dist_t>::max();
                candidateSet.emplace(-lowerBound, ep_id);
            }
            visited_array[ep_id] = visited_array_tag;

            while (!candidateSet.empty()) {
                std::pair<dist_t, tableint> curr_el_pair = candidateSet.top();
                if ((-curr_el_pair.first) > lowerBound) {
                    break;
                }
                candidateSet.pop();

                tableint curNodeNum = curr_el_pair.second;

                std::unique_lock <std::mutex> lock(link_list_locks_[curNodeNum]);

                linkcontainer *data;
                if (layer == 0) {
                    data = layer0.getData(curNodeNum);
                } else {
                    data = layers.getLinks(curNodeNum, layer);
                }
                size_t size = data->size();
#ifdef USE_SSE
                if (size > 1){
                    _mm_prefetch((char *) (visited_array + data->at(0)), _MM_HINT_T0);
                    _mm_prefetch((char *) (visited_array + data->at(0) + 64), _MM_HINT_T0);
                    _mm_prefetch(getDataByInternalId(data->at(0)), _MM_HINT_T0);
                    _mm_prefetch(getDataByInternalId(data->at(1)), _MM_HINT_T0);
                }
#endif

                for (size_t j = 0; j < size; j++) {
                    tableint candidate_id = data->at(j);
#ifdef USE_SSE
                    if (size > j + 1) {
                        _mm_prefetch((char *) (visited_array + data->at(j + 1)), _MM_HINT_T0);
                        _mm_prefetch(getDataByInternalId(data->at(j + 1)), _MM_HINT_T0);
                    }
#endif
                    if (visited_array[candidate_id] == visited_array_tag) continue;
                    visited_array[candidate_id] = visited_array_tag;
                    char *currObj1 = (getDataByInternalId(candidate_id));

                    dist_t dist1 = fstdistfunc_(data_point, currObj1, dist_func_param_);
                    if (top_candidates.size() < ef_construction_ || lowerBound > dist1) {
                        candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                        _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                        if (!isMarkedDeleted(candidate_id))
                            top_candidates.emplace(dist1, candidate_id);

                        if (top_candidates.size() > ef_construction_)
                            top_candidates.pop();

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }
            visited_list_pool_->releaseVisitedList(vl);

            return top_candidates;
        }

        template <bool has_deletions>
        std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst>
        searchBaseLayerST(tableint ep_id, const void *data_point, size_t ef, const SearchCondition& condition) const {
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidate_set;

            dist_t lowerBound;
            if (!has_deletions || !isMarkedDeleted(ep_id)) {
                dist_t dist = fstdistfunc_(data_point, getDataByInternalId(ep_id), dist_func_param_);
                lowerBound = dist;
                top_candidates.emplace(dist, ep_id);
                candidate_set.emplace(-dist, ep_id);
            } else {
                lowerBound = std::numeric_limits<dist_t>::max();
                candidate_set.emplace(-lowerBound, ep_id);
            }

            visited_array[ep_id] = visited_array_tag;

            while (!candidate_set.empty()) {

                std::pair<dist_t, tableint> current_node_pair = candidate_set.top();

                if ((-current_node_pair.first) > lowerBound) {
                    break;
                }
                candidate_set.pop();

                size_t num_checked = 0;

                tableint current_node_id = current_node_pair.second;
                const linkcontainer *data = layer0.getData(current_node_id);
                size_t size = data->size();
//                bool cur_node_deleted = isMarkedDeleted(current_node_id);

#ifdef USE_SSE
                if (size > 1) {
                    _mm_prefetch((char *) (visited_array + data->at(0)), _MM_HINT_T0);
                    _mm_prefetch((char *) (visited_array + data->at(0) + 64), _MM_HINT_T0);
                    _mm_prefetch(data_level0_memory_ + data->at(0) * size_data_per_element_, _MM_HINT_T0);
                    _mm_prefetch((char *) &(data->at(1)), _MM_HINT_T0);
                }
#endif

                for (size_t j = 0; j < size; j++) {
                    if (num_checked > maxSearchM0_) {
                        break;
                    }
                    int candidate_id = data->at(j);
                    if(!tags.checkCondition(candidate_id, condition)) {
                        continue;
                    }
                    num_checked += 1;
//                    if (candidate_id == 0) continue;
#ifdef USE_SSE
                    if (size > j + 1) {
                        _mm_prefetch((char *) (visited_array + data->at(j + 1)), _MM_HINT_T0);
                        _mm_prefetch(data_level0_memory_ + data->at(j + 1) * size_data_per_element_,
                                    _MM_HINT_T0);////////////
                    }
#endif
                    if (!(visited_array[candidate_id] == visited_array_tag)) {

                        visited_array[candidate_id] = visited_array_tag;

                        char *currObj1 = (getDataByInternalId(candidate_id));
                        dist_t dist = fstdistfunc_(data_point, currObj1, dist_func_param_);

                        if (top_candidates.size() < ef || lowerBound > dist) {
                            candidate_set.emplace(-dist, candidate_id);
#ifdef USE_SSE
                            _mm_prefetch(data_level0_memory_ + candidate_set.top().second * size_data_per_element_,///////////
                                         _MM_HINT_T0);////////////////////////
#endif

                            if (!has_deletions || !isMarkedDeleted(candidate_id))
                                top_candidates.emplace(dist, candidate_id);

                            if (top_candidates.size() > ef)
                                top_candidates.pop();

                            if (!top_candidates.empty())
                                lowerBound = top_candidates.top().first;
                        }
                    }
                }
            }

            visited_list_pool_->releaseVisitedList(vl);
            return top_candidates;
        }

        void getNeighborsByHeuristic2(
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> &top_candidates,
                const size_t M) {
            if (top_candidates.size() < M) {
                return;
            }
            std::priority_queue<std::pair<dist_t, tableint>> queue_closest;
            std::vector<std::pair<dist_t, tableint>> return_list;
            while (top_candidates.size() > 0) {
                queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
                top_candidates.pop();
            }

            while (queue_closest.size()) {
                if (return_list.size() >= M)
                    break;
                std::pair<dist_t, tableint> curent_pair = queue_closest.top();
                dist_t dist_to_query = -curent_pair.first;
                queue_closest.pop();
                bool good = true;
                for (std::pair<dist_t, tableint> second_pair : return_list) {
                    dist_t curdist =
                            fstdistfunc_(getDataByInternalId(second_pair.second),
                                         getDataByInternalId(curent_pair.second),
                                         dist_func_param_);;
                    if (curdist < dist_to_query) {
                        good = false;
                        break;
                    }
                }
                if (good) {
                    return_list.push_back(curent_pair);
                }
            }

            for (std::pair<dist_t, tableint> curent_pair : return_list) {

                top_candidates.emplace(-curent_pair.first, curent_pair.second);
            }
        }

        size_t get_node_level(tableint internal_id) const {
            return layers.getData(internal_id)->size();
        }

        void mutuallyConnectNewElement(const void *data_point, tableint cur_c,
                                       std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates,
                                       size_t level) {
            size_t Mcurmax = level ? maxM_ : maxM0_;
            getNeighborsByHeuristic2(top_candidates, M_);
            if (top_candidates.size() > M_)
                throw std::runtime_error("Should be not be more than M_ candidates returned by the heuristic");

            std::vector<tableint> selectedNeighbors;
            selectedNeighbors.reserve(M_);
            while (top_candidates.size() > 0) {
                selectedNeighbors.push_back(top_candidates.top().second);
                top_candidates.pop();
            }

            {
                linkcontainer *ll_cur;
                if (level == 0)
                    ll_cur = layer0.getData(cur_c);
                else
                    ll_cur = layers.getLinks(cur_c, level);

                if (!ll_cur->empty()) {
                    throw std::runtime_error("The newly inserted element should have blank link list");
                }
                ll_cur->resize(selectedNeighbors.size());

                for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
                    if (level > get_node_level(selectedNeighbors[idx]))
                        throw std::runtime_error("Trying to make a link on a non-existent level");
                    (*ll_cur)[idx] = selectedNeighbors[idx];
                }
            }

            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {

                std::unique_lock <std::mutex> lock(link_list_locks_[selectedNeighbors[idx]]);

                linkcontainer *ll_other;
                if (level == 0)
                    ll_other = layer0.getData(selectedNeighbors[idx]);
                else
                    ll_other = layers.getLinks(selectedNeighbors[idx], level);

                size_t sz_link_list_other = ll_other->size();

                // if (sz_link_list_other > Mcurmax)
                //     throw std::runtime_error("Bad value of sz_link_list_other");
                if (selectedNeighbors[idx] == cur_c)
                    throw std::runtime_error("Trying to connect an element to itself");
                if (level > get_node_level(selectedNeighbors[idx]))
                    throw std::runtime_error("Trying to make a link on a non-existent level");

                if (sz_link_list_other < Mcurmax) {
                    ll_other->push_back(cur_c);
                } else {
                    // finding the "weakest" element to replace it with the new one
                    dist_t d_max = fstdistfunc_(getDataByInternalId(cur_c), getDataByInternalId(selectedNeighbors[idx]),
                                                dist_func_param_);
                    // Heuristic:
                    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> candidates;
                    candidates.emplace(d_max, cur_c);

                    // Allow to change only first Mcurmax links, others should be considered special
                    for (size_t j = 0; j < std::min(sz_link_list_other, Mcurmax); j++) {
                        candidates.emplace(
                                fstdistfunc_(getDataByInternalId(ll_other->at(j)), getDataByInternalId(selectedNeighbors[idx]),
                                             dist_func_param_), ll_other->at(j));
                    }

                    getNeighborsByHeuristic2(candidates, Mcurmax);

                    int indx = 0;
                    while (candidates.size() > 0) {
                        (*ll_other)[indx] = candidates.top().second;
                        candidates.pop();
                        indx++;
                    }
                }
            }
        }

        std::mutex global;
        size_t ef_;

        void setEf(size_t ef) {
            ef_ = ef;
        }

        void resizeIndex(size_t new_max_elements){
            if (new_max_elements<cur_element_count)
                throw std::runtime_error("Cannot resize, max element is less than the current number of elements");


            delete visited_list_pool_;
            visited_list_pool_ = new VisitedListPool(1, new_max_elements);

            tags.resize(new_max_elements);
            layer0.resize(new_max_elements);
            layers.resize(new_max_elements);

            std::vector<std::mutex>(new_max_elements).swap(link_list_locks_);


            // Reallocate base layer
            char * data_level0_memory_new = (char *) malloc(new_max_elements * size_data_per_element_);
            if (data_level0_memory_new == nullptr)
                throw std::runtime_error("Not enough memory: resizeIndex failed to allocate base layer");
            memcpy(data_level0_memory_new, data_level0_memory_,cur_element_count * size_data_per_element_);
            free(data_level0_memory_);
            data_level0_memory_=data_level0_memory_new;

            max_elements_=new_max_elements;

        }

        void saveIndex(const std::string &location) {
            std::ofstream output(location, std::ios::binary);
            std::streampos position;

            writeBinaryPOD(output, max_elements_);
            writeBinaryPOD(output, cur_element_count);
            writeBinaryPOD(output, size_data_per_element_);
            writeBinaryPOD(output, label_offset_);
            writeBinaryPOD(output, maxlevel_);
            writeBinaryPOD(output, enterpoint_node_);
            writeBinaryPOD(output, maxM_);

            writeBinaryPOD(output, maxM0_);
            writeBinaryPOD(output, M_);
            writeBinaryPOD(output, random_seed_);
            // writeBinaryPOD(output, maxSearchM_);
            // writeBinaryPOD(output, maxSearchM0_);
            writeBinaryPOD(output, mult_);
            writeBinaryPOD(output, ef_construction_);

            writeMap(output, tag_to_entrypoint_);
            tags.serialize(output);
            layer0.serialize(output);
            layers.serialize(output);
            

            output.write(data_level0_memory_, cur_element_count * size_data_per_element_);

            output.close();
        }

        void loadIndex(const std::string &location, SpaceInterface<dist_t> *s, size_t max_elements_i=0) {


            std::ifstream input(location, std::ios::binary);

            if (!input.is_open())
                throw std::runtime_error("Cannot open file");


            // get file size:
            input.seekg(0,input.end);
            std::streampos total_filesize=input.tellg();
            input.seekg(0,input.beg);

            readBinaryPOD(input, max_elements_);
            readBinaryPOD(input, cur_element_count);

            size_t max_elements=max_elements_i;
            if(max_elements < cur_element_count)
                max_elements = max_elements_;
            max_elements_ = max_elements;
            readBinaryPOD(input, size_data_per_element_);
            readBinaryPOD(input, label_offset_);
            readBinaryPOD(input, maxlevel_);
            readBinaryPOD(input, enterpoint_node_);

            readBinaryPOD(input, maxM_);
            readBinaryPOD(input, maxM0_);
            readBinaryPOD(input, M_);
            readBinaryPOD(input, random_seed_);
            // readBinaryPOD(input, maxSearchM_);
            // readBinaryPOD(input, maxSearchM0_);
            readBinaryPOD(input, mult_);
            readBinaryPOD(input, ef_construction_);

            readMap(input, tag_to_entrypoint_);
            tags.reset();
            tags.resize(max_elements_);
            tags.deserialize(input);

            layer0.reset();
            layer0.resize(max_elements_);
            layer0.deserialize(input);
            
            layers.reset();
            layers.resize(max_elements_);
            layers.deserialize(input);

            level_generator_.seed(random_seed_);

            s_ = s;

            data_size_ = s->get_data_size();
            fstdistfunc_ = s->get_dist_func();
            dist_func_param_ = s->get_dist_func_param();

            do_free = true;
            data_level0_memory_ = (char *) malloc(max_elements * size_data_per_element_);
            if (data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");
            input.read(data_level0_memory_, cur_element_count * size_data_per_element_);

            std::vector<std::mutex>(max_elements).swap(link_list_locks_);

            visited_list_pool_ = new VisitedListPool(1, max_elements);

            revSize_ = 1.0 / mult_;
            ef_ = 10;

            for (size_t i = 0; i < cur_element_count; i++) {
                label_lookup_[getExternalLabel(i)]=i;
            }

            has_deletions_=false;

            for (size_t i = 0; i < cur_element_count; i++) {
                if(isMarkedDeleted(i))
                    has_deletions_=true;
            }

            if(input.tellg()!=total_filesize)
                throw std::runtime_error("Index seems to be corrupted or unsupported");
            
            input.close();

            return;
        }

        template<typename data_t>
        std::vector<data_t> getDataByLabel(labeltype label)
        {
            tableint label_c = getInterenalIdByLabel(label);
           
            char* data_ptrv = getDataByInternalId(label_c);
            size_t dim = *((size_t *) dist_func_param_);
            std::vector<data_t> data;
            data_t* data_ptr = (data_t*) data_ptrv;
            for (size_t i = 0; i < dim; i++) {
                data.push_back(*data_ptr);
                data_ptr += 1;
            }
            return data;
        }

        inline tableint getInterenalIdByLabel(labeltype label)
        {
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end()) {
                throw std::runtime_error("Label not found");
            }
            return search->second;
        }

        /**
         * Assigns tag to the datapoint
         * @param label
         * @param tag
         */
        void addTagByLabel(labeltype label, tagtype tag)
        {   
            addTagByIdx(getInterenalIdByLabel(label), tag);
        }

        void addTagByIdx(tableint idx, tagtype tag)
        {
            if (tag_to_entrypoint_.find(tag) == tag_to_entrypoint_.end()){
                tag_to_entrypoint_[tag] = idx;
            }else{
                size_t currect_level = get_node_level(tag_to_entrypoint_[tag]);
                size_t new_level = get_node_level(idx);
                if (new_level > currect_level) 
                    tag_to_entrypoint_[tag] = idx;
            }
            tags.addTag(idx, tag);
        }

        /**
         * Returns assigned category
         * @param label
         */
        const tagcontainer *getTagsByLabel(labeltype label)
        {
            return tags.getTags(getInterenalIdByLabel(label));
        }

        /**
         * Bulk assign categories to datapoints. 
         * @param lables
         * @param tag
         * @param index do additional indexing inside category?
         */
        void addTags(std::vector<labeltype> &lables, const tagtype tag)
        {
            for (labeltype label : lables) {
                tableint point_id = getInterenalIdByLabel(label);
                addTagByIdx(point_id, tag);
            }
        }

        void indexTag(tagtype tag, size_t m = 0) {
            auto tag_points = tags.tag_mapping.find(tag);
            if (tag_points == tags.tag_mapping.end())
                throw std::runtime_error("No points with this tag found");

            tableint tag_entrypoint = additinalIndex(tag_points->second, m);
            tag_to_entrypoint_[tag] = tag_entrypoint;
        }


        /**
         * Creates additional links between given tags.
         * Should ensure connectivity.
         */
        void indexTags(std::vector<tagtype> &indexing_tags, size_t m = 0) {
            if (m == 0) m = M_;

            HierarchicalNSW *temp_hnsw = new HierarchicalNSW(s_, max_elements_, m, ef_construction_, random_seed_, data_level0_memory_);

            std::unordered_set<tableint> indexed_points;
            for(tagtype tag : indexing_tags) {
                auto tag_points = tags.tag_mapping.find(tag);
                if (tag_points == tags.tag_mapping.end())
                    throw std::runtime_error("No points with tag " + std::to_string(tag) + " found");

                for(tableint idx : tag_points->second){
                    if( indexed_points.find(idx) != indexed_points.end()){
                        indexed_points.insert(idx);
                        temp_hnsw->linkNewPoint(idx, get_node_level(idx));
                    }
                }
            }

            layer0.mergeOther(temp_hnsw->layer0);
            layers.mergeOther(temp_hnsw->layers);

            delete temp_hnsw;
        }

        /**
         * Remove all tag information as well as additional link, created for tags support.
         * 
         */
        void resetTags() {
            tags.reset();
            tags.resize(max_elements_);
            layer0.shrink(maxM0_);
            layers.shrink(maxM_);
        }
        /**
         * Perform additional indexing over subset of points
         * 
         */
        tableint additinalIndex(std::set<tableint> &ids, size_t m = 0)
        {
            // Create a new Index with same data, but empty links
            if (m == 0) m = M_;

            HierarchicalNSW *temp_hnsw = new HierarchicalNSW(s_, max_elements_, m, ef_construction_, random_seed_, data_level0_memory_);

            for(tableint idx : ids)
                temp_hnsw->linkNewPoint(idx, get_node_level(idx));
            
            layer0.mergeOther(temp_hnsw->layer0);
            layers.mergeOther(temp_hnsw->layers);

            tableint entrypoint = temp_hnsw->enterpoint_node_;

            delete temp_hnsw;

            return entrypoint;
        }

        /**
         * Perform additional indexing over subset of points
         * 
         */
        void additinalIndexByLables(std::vector<labeltype> &labels)
        {
            std::set<tableint> ids = std::set<tableint>();
            std::transform(labels.begin(), labels.end(), std::inserter(ids, ids.begin()), getInterenalIdByLabel);
            additinalIndex(ids);
        }

        static const unsigned char DELETE_MARK = 0x01;
//        static const unsigned char REUSE_MARK = 0x10;
        /**
         * Marks an element with the given label deleted, does NOT really change the current graph.
         * @param label
         */
        void markDelete(labeltype label)
        {
            has_deletions_=true;
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end()) {
                throw std::runtime_error("Label not found");
            }
            markDeletedInternal(search->second);
        }

        /**
         * Uses the first 8 bits of the memory for the linked list to store the mark,
         * whereas maxM0_ has to be limited to the lower 24 bits, however, still large enough in almost all cases.
         * @param internalId
         */
        void markDeletedInternal(tableint internalId) {
            // ToDo: Make this via external lookup table
        }

        /**
         * Remove the deleted mark of the node.
         * @param internalId
         */
        void unmarkDeletedInternal(tableint internalId) {
            // ToDo: Make this via external lookup table
        }

        /**
         * Checks the first 8 bits of the memory to see if the element is marked deleted.
         * @param internalId
         * @return
         */
        bool isMarkedDeleted(tableint internalId) const {
            return false; // ToDo: Make this via external lookup table
        }

        void addPoint(const void *data_point, labeltype label) {
            size_t level = getRandomLevel(mult_);
            addPoint(data_point, label, level);
        }

        /**
         * Add new data point into vectors memory and label lookup
         * 
         */
        tableint addPointData(const void *data_point, labeltype label) {
            tableint cur_c = 0;
            {
                std::unique_lock <std::mutex> lock(cur_element_count_guard_);
                if (cur_element_count >= max_elements_) {
                    throw std::runtime_error("The number of elements exceeds the specified limit");
                };

                cur_c = cur_element_count;
                cur_element_count++;

                auto search = label_lookup_.find(label);
                if (search != label_lookup_.end()) {
                    std::unique_lock <std::mutex> lock_el(link_list_locks_[search->second]);
                    has_deletions_ = true;
                    markDeletedInternal(search->second);
                }
                label_lookup_[label] = cur_c;
            }

            memset(data_level0_memory_ + cur_c * size_data_per_element_, 0, size_data_per_element_);

            // Initialisation of the data and label
            memcpy(getExternalLabeLp(cur_c), &label, sizeof(labeltype));
            memcpy(getDataByInternalId(cur_c), data_point, data_size_);

            return cur_c;
        }

        /**
         * Link new point with other points in graph.
         * 
         */
        void linkNewPoint(tableint cur_c, size_t curlevel){
            const void *data_point = getDataByInternalId(cur_c);

            std::unique_lock <std::mutex> lock_el(link_list_locks_[cur_c]);

            std::unique_lock <std::mutex> templock(global);
            size_t maxlevelcopy = maxlevel_;
            if (curlevel <= maxlevelcopy) // We need this lock only if new levels will be added
                templock.unlock();
            tableint currObj = enterpoint_node_;
            tableint enterpoint_copy = enterpoint_node_;

            if (curlevel) {
                layers.setNumLayers(cur_c, curlevel, maxM_);
            }

            if ((signed)currObj != -1) {

                if (curlevel < maxlevelcopy) {

                    dist_t curdist = fstdistfunc_(data_point, getDataByInternalId(currObj), dist_func_param_);
                    for (size_t level = maxlevelcopy; level > curlevel; level--) {


                        bool changed = true;
                        while (changed) {
                            changed = false;
                            std::unique_lock <std::mutex> lock(link_list_locks_[currObj]);
                            const linkcontainer *links = layers.getLinks(currObj, level);

                            for (tableint cand: *links) {
                                if (cand < 0 || cand > max_elements_)
                                    throw std::runtime_error("cand error");
                                dist_t d = fstdistfunc_(data_point, getDataByInternalId(cand), dist_func_param_);
                                if (d < curdist) {
                                    curdist = d;
                                    currObj = cand;
                                    changed = true;
                                }
                            }
                        }
                    }
                }

                bool epDeleted = isMarkedDeleted(enterpoint_copy);
                for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--) {
                    std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates = searchBaseLayer(
                            currObj, data_point, (size_t)level);
                    if (epDeleted) {
                        top_candidates.emplace(fstdistfunc_(data_point, getDataByInternalId(enterpoint_copy), dist_func_param_), enterpoint_copy);
                        if (top_candidates.size() > ef_construction_)
                            top_candidates.pop();
                    }
                    mutuallyConnectNewElement(data_point, cur_c, top_candidates, (size_t)level);

                    currObj = top_candidates.top().second;
                }
            } else {
                // Do nothing for the first element
                enterpoint_node_ = cur_c;
                maxlevel_ = curlevel;

            }

            //Releasing lock for the maximum level
            if (curlevel > maxlevelcopy) {
                enterpoint_node_ = cur_c;
                maxlevel_ = curlevel;
            }
        }

        tableint addPoint(const void *data_point, labeltype label, size_t level) {
            tableint cur_c = addPointData(data_point, label);
            linkNewPoint(cur_c, level);
            return cur_c;
        };

        tableint getConditionedEntrypoint(const SearchCondition &condition) const
        {
            if (tags.checkCondition(enterpoint_node_, condition))
            {
                return enterpoint_node_;
            }

            auto candidates = condition.entryCandidates();

            for (tagtype tag : candidates)
            {
                auto tag_search = tag_to_entrypoint_.find(tag);
                if (tag_search != tag_to_entrypoint_.end()) {
                    tableint tag_entry = tag_search->second;
                    if (tags.checkCondition(tag_entry, condition))
                        return tag_entry;
                }
            }

            for (tagtype tag : candidates)
            {
                if (tag_to_entrypoint_.find(tag) == tag_to_entrypoint_.end()) continue;
                const std::set<tableint> *tag_nodes = &(tags.tag_mapping.at(tag));
                for (tableint idx : *tag_nodes)
                    if (tags.checkCondition(idx, condition))
                        return idx;
            }
            return 0;
        }

        std::priority_queue<std::pair<dist_t, labeltype >>
        searchKnn(const void *query_data, size_t k, SearchCondition &condition) const {
            std::priority_queue<std::pair<dist_t, labeltype >> result;
            if (cur_element_count == 0) return result;

            tableint currObj = getConditionedEntrypoint(condition);

            if (!tags.checkCondition(currObj, condition)) return result; //happens if nothing fits the condition.

            dist_t curdist = fstdistfunc_(query_data, getDataByInternalId(enterpoint_node_), dist_func_param_);

            size_t entry_level = get_node_level(currObj);

            for (size_t level = entry_level; level > 0; level--) {
                bool changed = true;
                while (changed) {
                    changed = false;
                    size_t num_checked = 0;
                    const linkcontainer *data = layers.getLinks(currObj, level);
                    for (tableint cand: *data) {
                        if (num_checked > maxSearchM_) {
                            break;
                        }
                        if(!tags.checkCondition(cand, condition)) {
                            continue;
                        }
                        if (cand < 0 || cand > max_elements_)
                            throw std::runtime_error("cand error");
                        dist_t d = fstdistfunc_(query_data, getDataByInternalId(cand), dist_func_param_);
                        num_checked += 1;
                        if (d < curdist) {
                            curdist = d;
                            currObj = cand;
                            changed = true;
                        }
                    }
                }
            }

            std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates;
            if (has_deletions_) {
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates1=searchBaseLayerST<true>(
                        currObj, query_data, std::max(ef_, k), condition);
                top_candidates.swap(top_candidates1);
            }
            else{
                std::priority_queue<std::pair<dist_t, tableint>, std::vector<std::pair<dist_t, tableint>>, CompareByFirst> top_candidates1=searchBaseLayerST<false>(
                        currObj, query_data, std::max(ef_, k), condition);
                top_candidates.swap(top_candidates1);
            }
            while (top_candidates.size() > k) {
                top_candidates.pop();
            }
            while (top_candidates.size() > 0) {
                std::pair<dist_t, tableint> rez = top_candidates.top();
                result.push(std::pair<dist_t, labeltype>(rez.first, getExternalLabel(rez.second)));
                top_candidates.pop();
            }
            return result;
        };
    };

}
