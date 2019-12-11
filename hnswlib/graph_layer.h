#pragma once

#include "types.h"
#include "element_container.h"
#include <mutex>
#include <string.h>
#include <vector>
#include <set>

namespace hnswlib
{

class GraphLayer0 : public ElementContainer<linkcontainer>
{
};

class GraphLayers : public ElementContainer<layercontainer>
{
public:
    void setNumLayers(tableint idx, size_t num_layers)
    {
        getData(idx)->resize(num_layers);
    }

    linkcontainer *getLinks(tableint idx, size_t layer)
    {
        return &data[idx][layer - 1];
    }

    const linkcontainer *getLinks(tableint idx, size_t layer) const
    {
        return &data[idx][layer - 1];
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
            for (size_t l = 0; l < num_layers; l++)
            {
                size_t length = 0;
                readBinaryPOD(in, length);
                typename layercontainer::value_type element_data;
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