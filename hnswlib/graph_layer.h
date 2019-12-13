#pragma once

#include "types.h"
#include "element_container.h"
#include <mutex>
#include <string.h>
#include <vector>
#include <set>

namespace hnswlib
{

void mergeLinks(linkcontainer &links, linkcontainer &other_links)
{
    auto links_start = links.begin();
    auto links_end = links.end();
    for (tableint idx : other_links)
    {
        // links is expected to be small. So linear search is optimal here
        auto search = std::find(links_start, links_end, idx);
        if (search == links_end)
        {
            links.push_back(idx);
        }
    }
}

class GraphLayer0 : public ElementContainer<linkcontainer>
{
public:
    void resize(size_t new_size, size_t reserve_size = 0)
    {
        this->data.resize(new_size);
        if ((reserve_size > 0) && (new_size > this->_max_elements))
        {
            for (size_t i = this->_max_elements; i < new_size; i++)
            {
                this->data[i].reserve(reserve_size);
            }
        }
        this->_max_elements = new_size;
    }

    void mergeOther(GraphLayer0 &other)
    {
        for (size_t element_id = 0; element_id < data.size(); element_id++)
        {
            mergeLinks(data[element_id], other.data[element_id]);
        }
    }

    void shrink(size_t max_m)
    {
        for (linkcontainer &links : data)
        {
            links.resize(std::min(max_m, links.size()));
        }
    }
};

class GraphLayers : public ElementContainer<layercontainer>
{
public:
    void setNumLayers(tableint idx, size_t num_layers, size_t max_m = 0)
    {
        layercontainer *element_layers = getData(idx);
        element_layers->resize(num_layers);
        if (max_m > 0)
        {
            for (linkcontainer &layer : *element_layers)
            {
                layer.reserve(max_m);
            }
        }
    }

    void shrink(size_t max_m)
    {
        for (layercontainer &element_layers : data)
        {
            for (linkcontainer &links : element_layers)
            {
                links.resize(std::min(max_m, links.size()));
            }
        }
    }

    linkcontainer *getLinks(tableint idx, size_t layer)
    {
        return &data[idx][layer - 1];
    }

    const linkcontainer *getLinks(tableint idx, size_t layer) const
    {
        return &data[idx][layer - 1];
    }

    /**
     * Update this object with links from other.
     * 
     */
    void mergeOther(GraphLayers &other)
    {
        if (data.size() != other.data.size())
            throw std::runtime_error("Can only merge layers with equal number of elements");

        for (size_t element_id = 0; element_id < data.size(); element_id++)
        {
            layercontainer *cur_layers = getData(element_id);
            layercontainer *other_layers = other.getData(element_id);
            size_t other_layers_count = other_layers->size();

            // If other graph has more layers - just push them back
            for (size_t layer_idx = 0; layer_idx < other_layers_count; layer_idx++)
            {
                if (layer_idx >= cur_layers->size())
                {
                    cur_layers->push_back(other_layers->at(layer_idx));
                }
                else
                {
                    mergeLinks(cur_layers->at(layer_idx), other_layers->at(layer_idx));
                }
            }
        }
    }

    void serialize(std::ostream &out)
    {
        for (layercontainer &element_layers : data)
        {
            size_t num_layers = element_layers.size();
            writeBinaryPOD(out, num_layers);
            for (const auto &element_data : element_layers)
            {
                size_t length = element_data.size();
                writeBinaryPOD(out, length);
                for (const auto &record : element_data)
                {
                    writeBinaryPOD(out, record);
                }
            }
        }
    }

    void deserialize(std::istream &in)
    {
        for (size_t i = 0; i < this->_max_elements; i++)
        {
            size_t num_layers = 0;
            readBinaryPOD(in, num_layers);
            layercontainer layer = layercontainer();
            layer.reserve(num_layers);
            for (size_t l = 0; l < num_layers; l++)
            {
                size_t length = 0;
                readBinaryPOD(in, length);
                typename layercontainer::value_type element_data;
                element_data.reserve(length); // ???
                for (size_t j = 0; j < length; j++)
                {
                    typename linkcontainer::value_type record = 0;
                    readBinaryPOD(in, record);
                    element_data.insert(element_data.end(), record);
                }
                layer.push_back(std::move(element_data));
            }
            data[i].swap(layer);
        }
    }
};
} // namespace hnswlib