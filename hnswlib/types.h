#pragma once

#include <fstream>

namespace hnswlib
{
typedef unsigned int tableint;
typedef size_t labeltype;
typedef size_t tagtype;

typedef std::unordered_set<tagtype> tagcontainer;

typedef std::vector<tableint> linkcontainer;
typedef std::vector<linkcontainer> layercontainer;

/**
 * Used for defining conditions.
 * 
 * (A | !B) & C -> [[(0, A), (1, B)], [(0, C)]]
 */
typedef std::pair<bool, tagtype> negable_tag_t;
typedef std::vector<negable_tag_t> or_condition_t;
typedef std::vector<or_condition_t> condition_t;

template <typename T>
static void writeBinaryPOD(std::ostream &out, const T &podRef)
{
    out.write((char *)&podRef, sizeof(T));
}

template <typename T>
static void readBinaryPOD(std::istream &in, T &podRef)
{
    in.read((char *)&podRef, sizeof(T));
}

template <typename T1, typename T2>
static void writeMap(std::ostream &out, const std::unordered_map<T1, T2> &podRef)
{   
    writeBinaryPOD(out, podRef.size());
    for(const auto &pair : podRef) {
        writeBinaryPOD(out, pair.first);
        writeBinaryPOD(out, pair.second);
    }
}

template <typename T1, typename T2>
static void readMap(std::istream &in, std::unordered_map<T1, T2> &podRef)
{   
    size_t size = 0;
    readBinaryPOD(in, size);
    for (size_t i = 0; i < size ; i++) {
        T1 key;
        T2 value;
        readBinaryPOD(in, key);
        readBinaryPOD(in, value);
        podRef[key] = value;
    }
}


} // namespace hnswlib