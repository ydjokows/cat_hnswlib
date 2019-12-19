import unittest
from typing import List
import tqdm
import numpy as np
from statsmodels.stats.proportion import proportion_confint
from collections import defaultdict


def calc_precision_at(found_pos: List[int], limit):
    hits = np.array(found_pos) < limit
    return np.mean(hits), proportion_confint(sum(hits), len(found_pos))

def cosine_similarity(vector, matrix):
    return (np.sum(vector * matrix, axis=1) / (
            np.sqrt(np.sum(matrix ** 2, axis=1)) * np.sqrt(np.sum(vector ** 2))))

def get_closest_brut(target, data, mask):
    true_distance = 1 - cosine_similarity(target, data)
    np.putmask(true_distance, ~mask, 1_000_000)
    closest = list(np.argsort(true_distance))
    return closest

def get_constraints(num_elements, parts_count, part_id):
    mask = np.arange(0, num_elements) % parts_count == part_id

    tags = defaultdict(list)
    for i in range(num_elements):
        tags[i % parts_count].append(i)

    condition = [[(False, part_id)]]

    return tags, mask, condition

def get_random_vector(dim):
    return np.float32(np.random.random((1, dim)))

def get_top_fount(true_labels, found_labels):
    found_top = []
    for found_label in found_labels:
        found_top.append(true_labels.index(found_label))
    return found_top


class ConditionalSeachTestCase(unittest.TestCase):

    def test_random_subsample(self):
        import hnswlib        
        dim = 50
        elements = 10_000
        attempts = 100

        points = np.random.rand(elements, dim)

        tags, mask, condition = get_constraints(elements, 100, 66)

        hnsw = hnswlib.Index(space='cosine', dim=dim)
        hnsw.init_index(max_elements = elements, ef_construction = 10, M = 16, random_seed=45)
        hnsw.set_num_threads(1)

        hnsw.add_items(points)

        for tag, ids in tqdm.tqdm(tags.items()):
            hnsw.add_tags(ids, tag)
            hnsw.index_tagged(tag)

        top_hits = []
        for _ in range(attempts):
            target = get_random_vector(dim)
            true_closest = get_closest_brut(target, points, mask)
            found, _ = hnsw.knn_query(target, k=10, conditions=condition)
            top_found = get_top_fount(true_closest, found[0])

            top_hit = top_found[0]
            top_hits.append(top_hit)
        
        precision, conf_interval = calc_precision_at(top_hits, 1)
        print("precision:", precision)

        self.assertGreater(precision, 0.9)


if __name__ == "__main__":
    unittest.main()