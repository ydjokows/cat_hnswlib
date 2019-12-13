#pragma once

#include "types.h"
#include <mutex>
#include <string.h>
#include <vector>
#include <set>

namespace hnswlib
{

template <typename DataType>
class ParamsContainer
{
public:
    std::vector<DataType> data = std::vector<DataType>();
    size_t _max_elements = 0;


    ParamsContainer(size_t max_elements = 0)
    {
        data.resize(max_elements);
        _max_elements = max_elements;
    }

    ~ParamsContainer()
    {
    }

    void resize(size_t new_size)
    {
        _max_elements = new_size;
        data.resize(new_size);
    }

    void reset()
    {
        _max_elements = 0;
        data.clear();
    }

    const DataType *getData(tableint idx) const
    {
        return &data[idx];
    }

    DataType *getData(tableint idx)
    {
        return &data[idx];
    }

    void serialize(std::ostream &out)
    {
        for (size_t i = 0; i < _max_elements; i++)
        {
            auto *element_data = &data[i];
            writeBinaryPOD(out, *element_data);
        }
    }

    void deserialize(std::istream &in)
    {
        for (size_t i = 0; i < _max_elements; i++)
        {   
            DataType record = 0;
            readBinaryPOD(in, record);
            data[i] = record;
        }
    }
};

template <typename Container>
class ElementContainer : public ParamsContainer<Container>
{
public:

    void serialize(std::ostream &out)
    {
        for (size_t i = 0; i < this->_max_elements; i++)
        {
            auto *element_data = &this->data[i];
            size_t length = element_data->size();
            writeBinaryPOD(out, length);
            for (const auto &record : *element_data)
            {
                writeBinaryPOD(out, record);
            }
        }
    }

    void deserialize(std::istream &in)
    {
        for (size_t i = 0; i < this->_max_elements; i++)
        {
            auto *element_data = &this->data[i];
            size_t length = 0;
            readBinaryPOD(in, length);
            element_data->reserve(length);
            for (size_t j = 0; j < length; j++)
            {
                typename Container::value_type record = 0;
                readBinaryPOD(in, record);
                element_data->insert(element_data->end(), record);
            }
        }
    }


};

} // namespace hnswlib