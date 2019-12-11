#pragma once

#include <fstream>

namespace hnswlib
{
typedef unsigned int tableint;
typedef size_t labeltype;
typedef size_t tagtype;

typedef std::set<tagtype> tagcontainer;

typedef std::vector<tableint> linkcontainer;
typedef std::vector<linkcontainer> layercontainer;

/**
 * Used for defining conditions.
 * 
 * (A | !B) & C -> [[(0, A), (1, B)], [(0, C)]]
 */
typedef std::pair<bool, tagtype> negable_tag_t;
typedef std::set<negable_tag_t> or_condition_t;
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

} // namespace hnswlib