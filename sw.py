from decorators import timer
import numpy as np # type: ignore
from swift import swift_alignment
from collections import defaultdict
from typing import List, Tuple, Dict

from test_data import generate_random_dna_sequence

@timer
def smith_waterman_full(query: str, target: str, match_score: int = 2, mismatch_penalty: int = -1, gap_penalty: int = -2) -> Tuple[int, str, str]:
    """
    Traditional Smith-Waterman algorithm for full matrix alignment.
    Returns the best alignment score and the aligned sequences.
    """
    m, n = len(query), len(target)
    score_matrix = np.zeros((m+1, n+1), dtype=int)
    max_score = 0
    max_pos = (0, 0)

    # Fill score matrix
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



if __name__ == "__main__":
    query_sequence = generate_random_dna_sequence(1000)
    target_sequence = generate_random_dna_sequence(5000)
    results = swift_alignment(query_sequence, target_sequence)

    for score, aligned_q, aligned_t, tile in results:
        print(f"Score: {score}, Aligned Query: {aligned_q}, Aligned Target: {aligned_t}, Tile: {tile}")