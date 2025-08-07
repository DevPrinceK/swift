'''
SWIFT (Smith-Waterman with Intelligent Filtering and Tiling)
An extended version of the Smith-Waterman algorithm 
that incorporates k-mer filtering and tiling for 
improved performance in local sequence alignment.
'''

from decorators import timer
import numpy as np
from collections import defaultdict
from typing import List, Tuple, Dict

from test_data import generate_random_dna_sequence

def rolling_kmers(sequence: str, k: int) -> Dict[str, List[int]]:
    """
    Extract k-mers with their positions from a sequence.
    """
    kmers = defaultdict(list)
    for i in range(len(sequence) - k + 1):
        kmers[sequence[i:i+k]].append(i)
    return kmers

def find_tile_candidates(query_kmers: Dict[str, List[int]], target: str, k: int, tile_size: int, threshold: int) -> List[Tuple[int, int]]:
    """
    Find candidate tile regions in target sequence based on k-mer hits.
    """
    hits = defaultdict(int)
    for kmer in query_kmers:
        for i in range(len(target) - k + 1):
            if target[i:i+k] == kmer:
                tile_start = (i // tile_size) * tile_size
                hits[tile_start] += 1

    # Filter tiles that exceed the threshold number of k-mer hits
    candidate_tiles = [(start, start + tile_size) for start, count in hits.items() if count >= threshold]
    return candidate_tiles

def smith_waterman_tile(query: str, target: str, match_score: int = 2, mismatch_penalty: int = -1, gap_penalty: int = -2) -> Tuple[int, str, str]:
    """
    Perform Smith-Waterman alignment on a tile (local alignment).
    """
    m, n = len(query), len(target)
    score_matrix = np.zeros((m+1, n+1), dtype=int)
    max_score = 0
    max_pos = (0, 0)

    # Fill the score matrix
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score_matrix[i-1][j-1] + (match_score if query[i-1] == target[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(0, match, delete, insert)

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback
    aligned_query = []
    aligned_target = []
    i, j = max_pos
    while score_matrix[i][j] != 0:
        current_score = score_matrix[i][j]
        if current_score == score_matrix[i-1][j-1] + (match_score if query[i-1] == target[j-1] else mismatch_penalty):
            aligned_query.append(query[i-1])
            aligned_target.append(target[j-1])
            i -= 1
            j -= 1
        elif current_score == score_matrix[i-1][j] + gap_penalty:
            aligned_query.append(query[i-1])
            aligned_target.append('-')
            i -= 1
        else:
            aligned_query.append('-')
            aligned_target.append(target[j-1])
            j -= 1

    return max_score, ''.join(reversed(aligned_query)), ''.join(reversed(aligned_target))

@timer
def swift_alignment(query: str, target: str, k: int = 3, tile_size: int = 32, hit_threshold: int = 2) -> List[Tuple[int, str, str, Tuple[int, int]]]:
    """
    SWIFT alignment across multiple tiles based on k-mer filtering.
    """
    query_kmers = rolling_kmers(query, k)
    candidate_tiles = find_tile_candidates(query_kmers, target, k, tile_size, hit_threshold)

    results = []
    for start, end in candidate_tiles:
        tile = target[start:end]
        score, aligned_q, aligned_t = smith_waterman_tile(query, tile)
        if score > 0:
            results.append((score, aligned_q, aligned_t, (start, end)))

    results.sort(reverse=True, key=lambda x: x[0])  # Sort by score
    return results



if __name__ == "__main__":
    query_sequence = generate_random_dna_sequence(1000)
    target_sequence = generate_random_dna_sequence(5000)
    results = swift_alignment(query_sequence, target_sequence)

    for score, aligned_q, aligned_t, tile in results:
        print(f"Score: {score}, Aligned Query: {aligned_q}, Aligned Target: {aligned_t}, Tile: {tile}")