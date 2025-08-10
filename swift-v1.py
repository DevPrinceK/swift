'''
SWIFT (Smith-Waterman with Intelligent Filtering and Tiling)
An extended version of the Smith-Waterman algorithm 
that incorporates k-mer filtering and tiling for 
improved performance in local sequence alignment.
'''

from decorators import timer
import numpy as np # type: ignore
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


def smith_waterman(query: str, target: str, match_score: int = 2, mismatch_penalty: int = -1, gap_penalty: int = -2) -> Tuple[int, int, int, int, int]:
    """
    Perform Smith-Waterman local alignment.

    Returns:
        A tuple (query_start, query_end, target_start, target_end, score)
    """
    m, n = len(query), len(target)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i-1][j-1] + (match_score if query[i-1] == target[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            score_matrix[i][j] = max(0, match, delete, insert)
            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_pos = (i, j)

    # Traceback
    i, j = max_pos
    query_end, target_end = i, j
    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        current = score_matrix[i][j]
        diag = score_matrix[i-1][j-1]
        up = score_matrix[i-1][j]
        left = score_matrix[i][j-1]

        if current == diag + (match_score if query[i-1] == target[j-1] else mismatch_penalty):
            i -= 1
            j -= 1
        elif current == up + gap_penalty:
            i -= 1
        elif current == left + gap_penalty:
            j -= 1
        else:
            break

    query_start, target_start = i, j
    return query_start, query_end, target_start, target_end, max_score


def merge_and_refine_alignments(alignments: List[Tuple[int, int, int, int, int]], query: str, target: str) -> List[Tuple[int, int, int, int, int]]:
    """
    Merge and refine overlapping or adjacent alignments.

    Parameters:
        alignments (List[Tuple[int, int, int, int, int]]): List of local alignments from tiles.
        query (str): The full query sequence.
        target (str): The full target sequence.

    Returns:
        List[Tuple[int, int, int, int, int]]: Refined list of alignments.
    """
    if not alignments:
        return []

    # Step 1: Sort alignments by query_start and target_start
    alignments.sort(key=lambda x: (x[0], x[2]))

    # Step 2: Merge overlapping or adjacent alignments
    merged = []
    current = list(alignments[0])

    for aln in alignments[1:]:
        q_start, q_end, t_start, t_end, score = aln
        if q_start <= current[1] and t_start <= current[3]:  # Overlapping or adjacent
            current[1] = max(current[1], q_end)
            current[3] = max(current[3], t_end)
            current[4] += score  # Aggregate scores (can be replaced by better logic)
        else:
            merged.append(tuple(current))
            current = list(aln)
    merged.append(tuple(current))

    # Step 3: Refine merged regions with Smith-Waterman
    refined_alignments = []
    for q_start, q_end, t_start, t_end, _ in merged:
        refined_aln = smith_waterman(query[q_start:q_end], target[t_start:t_end])
        refined_alignments.append((
            q_start + refined_aln[0],
            q_start + refined_aln[1],
            t_start + refined_aln[2],
            t_start + refined_aln[3],
            refined_aln[4]
        ))

    return refined_alignments


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
    query_sequence = ""
    target_sequence = ""

    # read sequences from file
    with open("test_data.txt", "r") as f:
        lines = f.readlines()
        for i in range(len(lines)):
            lines[i] = lines[i].strip()
        query_sequence = lines[1]
        target_sequence = lines[3]

    # Step 1: Run tile-wise SWIFT alignment
    raw_results = swift_alignment(query_sequence, target_sequence)

    # Step 2: Convert string-aligned results to coordinate results for merging
    alignment_spans = []
    for score, aligned_q, aligned_t, (tile_start, tile_end) in raw_results:
        # Rough estimation of aligned region
        q_aln_len = len(aligned_q.replace("-", ""))
        t_aln_len = len(aligned_t.replace("-", ""))
        alignment_spans.append((
            0, q_aln_len,  # assuming query spans from 0 to aligned length
            tile_start, tile_start + t_aln_len,
            score
        ))

    # Step 3: Merge and refine
    refined_alignments = merge_and_refine_alignments(alignment_spans, query_sequence, target_sequence)

    # Step 4: Display refined alignments
    for q_start, q_end, t_start, t_end, score in refined_alignments:
        print(f"Refined Alignment - Score: {score}")
        print(f"Query [{q_start}:{q_end}]: {query_sequence[q_start:q_end]}")
        print(f"Target [{t_start}:{t_end}]: {target_sequence[t_start:t_end]}\n")
