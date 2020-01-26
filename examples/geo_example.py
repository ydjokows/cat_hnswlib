import numpy as np
import tqdm
import geohash
import hnswlib
import random
import sys
from collections import defaultdict

base_alphabet = '0123456789abcdefghijklmnopqrstuv'
geo_alphabet =  '0123456789bcdefghjkmnpqrstuvwxyz' 
trantab = str.maketrans(geo_alphabet, base_alphabet)


def cosine_similarity(vector, matrix):
    return (np.sum(vector * matrix, axis=1) / (
            np.sqrt(np.sum(matrix ** 2, axis=1)) * np.sqrt(np.sum(vector ** 2))))


# The library can only use int as a tag. So we need to convert geohash into integer first
def geohash2int(geo: str) -> int:
    """
    Converts geohash string into integer
    """
    return int(geo.translate(trantab), 32)


def get_random_vector(dim):
    return np.float32(np.random.random((1, dim)))


def get_random_point(from_lat, to_lat, from_lon, to_lon):
    lat = random.uniform(from_lat, to_lat)
    lon = random.uniform(from_lon, to_lon)
    return lat, lon


def get_random_data(num_points, dim, from_lat, to_lat, from_lon, to_lon):
    points = np.random.rand(num_points, dim)
    geo_points = [get_random_point(from_lat, to_lat, from_lon, to_lon) for _ in range(num_points)]
    return points, geo_points


if __name__ == "__main__":
    from_lat, to_lat = 52.4245, 52.6176
    from_lon, to_lon = 13.1870, 13.5997

    dim = 25
    elements = 100_000
    max_precision = 6  # Minimal searchable precision. Precision of 6 is ~ 0.61 km 
    # https://en.wikipedia.org/wiki/Geohash#Number_of_geohash_characters_and_precision_in_km
    
    hnsw = hnswlib.Index(space='cosine', dim=dim)
    hnsw.init_index(max_elements = elements, M = 16, random_seed=45)
    hnsw.set_num_threads(2)

    # Generate random vectors and geo points
    points, geo_points = get_random_data(elements, dim, from_lat, to_lat, from_lon, to_lon)

    hnsw.add_items(points)

    tags_to_index = defaultdict(int)
    tags_to_ids = defaultdict(list)

    # Collect geohashes for indexing
    for idx, geo_point in enumerate(geo_points):
        lat, lon = geo_point
        ghsh = geohash.encode(lat, lon, precision=max_precision)
        # List all hashes in hierarchy: 'u337jk' -> ['u', 'u3', 'u33', 'u337', 'u337j', 'u337jk'] 
        tags = [ghsh[:i + 1] for i in range(max_precision)]

        # Save small geohash indexes with further indexing
        tags_to_index[ghsh[:max_precision]] += 1
        tags_to_index[ghsh[:max_precision - 1]] += 1
        
        # Assign geotags to points
        for tag in tags:
            tags_to_ids[tag].append(idx)
            hnsw.add_tags([idx], geohash2int(tag))

    # Additionally index points inside small regions 
    for tag in tqdm.tqdm(tags_to_index):
        # This will create additional links in a graph for each geohash region.
        # So search should work on nodes inside this region only.
        hnsw.index_tagged(geohash2int(tag))
        # With M=16 additional indexing is only required for regions containing less than ~5% of all points
        # Additional info here: https://comprehension.ml/posts/categorical-hnsw/

    for tag in tqdm.tqdm(tags_to_index):
        # This code will also create additional connections between points in neighbor regions.
        # So search in multiple neighbor regions will also work
        neighbors = [geohash2int(ntag) for ntag in geohash.neighbors(tag) if ntag in tags_to_index]
        hnsw.index_cross_tagged(neighbors)

    # Performing query

    target_query = get_random_vector(dim)
    # Hash precision defines radius of a seearch. Precision of 5 is ~ 2.4Km 
    # https://en.wikipedia.org/wiki/Geohash#Number_of_geohash_characters_and_precision_in_km
    target_precision = 5  
    target_lat, target_lon = 52.5175, 13.3937 

    # Generate integer tag from geohash
    target_ghsh = geohash.encode(target_lat, target_lon, precision=target_precision)
    target_tag = geohash2int(target_ghsh)

    # Obtain search condition from geohash
    # You can also search in multiple squares with conjunction: 
    # [[(False, hash1), (False, hash2), ..., (False, hashN)]]
    condition = [[(False, target_tag)]]

    found, dist = hnsw.knn_query(target_query, k=3, conditions=condition)
    
    print(found, dist)

    # Check search precision with brutforce approach

    true_distance = 1 - cosine_similarity(target_query, points)
    mask = np.zeros(elements, dtype=bool)
    mask[tags_to_ids[target_ghsh]] = True  # Search in given geo-region only
    np.putmask(true_distance, ~mask, 1_000_000)  
    closest = list(np.argsort(true_distance))  # Closest by mask

    print(closest[:3], true_distance[closest[:3]])
