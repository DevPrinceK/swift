'''
SWIFT (Smith-Waterman with Intelligent Filtering and Tiling)
An extended version of the Smith-Waterman algorithm 
that incorporates k-mer filtering and tiling for 
improved performance in local sequence alignment.
EXTRAS:
- Parallel tile alignment using multiprocessing
- Merge, refine (with Smith-Waterman), and export results (JSON + SAM)
- ASCII + matplotlib visualization
- Performance comparison with traditional Smith-Waterman
'''


from time import perf_counter
from multiprocessing import Pool, cpu_count
from collections import defaultdict
from typing import List, Tuple, Dict, Optional
import numpy as np
import json
import matplotlib.pyplot as plt
import os

# ---------------------------
# Timer decorator (keeps your original behavior)
# ---------------------------
def timer(func):
    def wrapper(*args, **kwargs):
        start_time = perf_counter()
        result = func(*args, **kwargs)
        end_time = perf_counter()
        print(f"Execution time for ({func.__name__}): {end_time - start_time:.4f} seconds")
        return result
    return wrapper

# ---------------------------
# Basic utilities: k-mers and tile candidate finding
# ---------------------------
def rolling_kmers(sequence: str, k: int) -> Dict[str, List[int]]:
    kmers = defaultdict(list)
    for i in range(len(sequence) - k + 1):
        kmers[sequence[i:i+k]].append(i)
    return kmers

def find_tile_candidates(query_kmers: Dict[str, List[int]], target: str, k: int, tile_size: int, threshold: int) -> List[Tuple[int, int]]:
    hits = defaultdict(int)
    # iterate target once and check k-mer membership for speed
    query_kmer_set = set(query_kmers.keys())
    for i in range(len(target) - k + 1):
        t_kmer = target[i:i+k]
        if t_kmer in query_kmer_set:
            tile_start = (i // tile_size) * tile_size
            hits[tile_start] += 1
    candidate_tiles = [(start, min(start + tile_size, len(target))) for start, count in hits.items() if count >= threshold]
    candidate_tiles.sort()
    return candidate_tiles

# ---------------------------
# Smith-Waterman for a tile (returns max local alignment + aligned strings)
# ---------------------------
def smith_waterman_tile(query: str, target: str, match_score: int = 2, mismatch_penalty: int = -1, gap_penalty: int = -2) -> Tuple[int, str, str]:
    m, n = len(query), len(target)
    score_matrix = np.zeros((m+1, n+1), dtype=int)
    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score_matrix[i-1][j-1] + (match_score if query[i-1] == target[j-1] else mismatch_penalty)
            delete = score_matrix[i-1][j] + gap_penalty
            insert = score_matrix[i][j-1] + gap_penalty
            val = max(0, match, delete, insert)
            score_matrix[i][j] = val
            if val > max_score:
                max_score = val
                max_pos = (i, j)

    # Traceback
    aligned_query = []
    aligned_target = []
    i, j = max_pos
    while i > 0 and j > 0 and score_matrix[i][j] != 0:
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

# ---------------------------
# Full Smith-Waterman returning coordinates & score (used for refine / comparison)
# ---------------------------
def smith_waterman(query: str, target: str, match_score: int = 2, mismatch_penalty: int = -1, gap_penalty: int = -2) -> Tuple[int, int, int, int, int]:
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

# ---------------------------
# Helper to map aligned strings to coordinates (when the query used to align was the *full* query)
# ---------------------------
def tile_alignment_to_coords(score: int, aligned_q: str, aligned_t: str, tile_start: int) -> Tuple[int, int, int, int, int]:
    # find first and length of non-gap portions
    def first_non_gap(s):
        for idx, ch in enumerate(s):
            if ch != '-':
                return idx
        return len(s)

    def non_gap_count(s):
        return sum(1 for ch in s if ch != '-')

    q_first = first_non_gap(aligned_q)
    t_first = first_non_gap(aligned_t)
    q_non = non_gap_count(aligned_q)
    t_non = non_gap_count(aligned_t)

    query_start = q_first
    query_end = q_first + q_non
    target_start = tile_start + t_first
    target_end = tile_start + t_first + t_non
    return query_start, query_end, target_start, target_end, score

# ---------------------------
# Merge and refine alignments (your earlier logic, cleaned)
# ---------------------------
def merge_and_refine_alignments(alignments: List[Tuple[int, int, int, int, int]], query: str, target: str) -> List[Tuple[int, int, int, int, int]]:
    if not alignments:
        return []
    # Sort by query_start then target_start
    alignments.sort(key=lambda x: (x[0], x[2]))

    merged = []
    current = list(alignments[0])

    for aln in alignments[1:]:
        q_start, q_end, t_start, t_end, score = aln
        # if overlapping/adjacent in both sequences
        if q_start <= current[1] and t_start <= current[3]:
            current[1] = max(current[1], q_end)
            current[3] = max(current[3], t_end)
            current[4] += score
        else:
            merged.append(tuple(current))
            current = list(aln)
    merged.append(tuple(current))

    # Refine each merged region with Smith-Waterman (local) on the spans
    refined = []
    for q_start, q_end, t_start, t_end, _ in merged:
        # Ensure bounds
        q_sub = query[q_start:q_end]
        t_sub = target[t_start:t_end]
        if len(q_sub) == 0 or len(t_sub) == 0:
            continue
        refined_qs, refined_qe, refined_ts, refined_te, refined_score = smith_waterman(q_sub, t_sub)
        refined.append((
            q_start + refined_qs,
            q_start + refined_qe,
            t_start + refined_ts,
            t_start + refined_te,
            refined_score
        ))
    return refined

# ---------------------------
# Multiprocessing wrapper for tiles
# ---------------------------
def _process_tile(args):
    query, tile, start = args
    score, aligned_q, aligned_t = smith_waterman_tile(query, tile)
    if score > 0:
        return (score, aligned_q, aligned_t, (start, start + len(tile)))
    return None

@timer
def swift_alignment_parallel(query: str, target: str, k: int = 3, tile_size: int = 64, hit_threshold: int = 2, workers: Optional[int] = None) -> List[Tuple[int, str, str, Tuple[int, int]]]:
    if workers is None:
        workers = max(1, cpu_count() - 0)
    query_kmers = rolling_kmers(query, k)
    candidate_tiles = find_tile_candidates(query_kmers, target, k, tile_size, hit_threshold)
    if not candidate_tiles:
        return []

    args = [(query, target[start:end], start) for start, end in candidate_tiles]
    with Pool(processes=workers) as pool:
        results = pool.map(_process_tile, args)

    results = [r for r in results if r]
    results.sort(reverse=True, key=lambda x: x[0])
    return results

# ---------------------------
# ASCII alignment + matplotlib plot
# ---------------------------
def ascii_alignment(aligned_q: str, aligned_t: str) -> str:
    match_line = "".join("|" if a == b and a != '-' else (" " if a != '-' or b != '-' else " ") for a, b in zip(aligned_q, aligned_t))
    return f"{aligned_q}\n{match_line}\n{aligned_t}"

def plot_alignment_identity(aligned_q: str, aligned_t: str, title: str = "Alignment Identity"):
    identity = [1 if a == b and a != '-' else 0 for a, b in zip(aligned_q, aligned_t)]
    plt.figure(figsize=(max(6, len(identity)/20), 2))
    plt.bar(range(len(identity)), identity)
    plt.xlabel("Alignment Position")
    plt.ylabel("Match (1) / Mismatch (0)")
    plt.title(title)
    plt.tight_layout()
    plt.show()


def export_alignments_json(alignments, file_path="alignments.json"):
    def convert(obj):
        if isinstance(obj, (np.integer,)):
            return int(obj)
        elif isinstance(obj, (np.floating,)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()
        elif isinstance(obj, (tuple, list)):
            return [convert(x) for x in obj]
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        return obj

    alignments_info = []
    for aln in alignments:
        converted = {k: convert(v) for k, v in aln.items()}
        alignments_info.append(converted)

    with open(file_path, "w") as f:
        json.dump(alignments_info, f, indent=2)

    print(f"[INFO] Alignments exported to {file_path}")



def make_cigar(aligned_q: str, aligned_t: str) -> str:
    # Build simple CIGAR string: M (both non-gap), I (gap in target -> insertion in target), D (gap in query -> deletion from target)
    cigar_parts = []
    prev_op = None
    run_len = 0
    for a, b in zip(aligned_q, aligned_t):
        if a != '-' and b != '-':
            op = 'M'
        elif a != '-' and b == '-':
            op = 'I'  # insertion to reference (consumes read/query)
        elif a == '-' and b != '-':
            op = 'D'  # deletion from read (consumes ref/target)
        else:
            # both gaps shouldn't happen in valid alignments
            op = 'M'
        if op == prev_op:
            run_len += 1
        else:
            if prev_op is not None:
                cigar_parts.append(f"{run_len}{prev_op}")
            prev_op = op
            run_len = 1
    if prev_op is not None and run_len > 0:
        cigar_parts.append(f"{run_len}{prev_op}")
    return "".join(cigar_parts) if cigar_parts else f"{len(aligned_q)}M"

def export_alignments_sam(alignments_info: List[dict], target_name: str = "target", file_path: str = "alignments.sam"):
    with open(file_path, "w") as f:
        # Header: length unknown, but SAM header needs @SQ lines; leave LN unknown or compute maximum target coordinate
        max_coord = 0
        for ai in alignments_info:
            if ai['target_end'] > max_coord:
                max_coord = ai['target_end']
        if max_coord == 0:
            f.write(f"@SQ\tSN:{target_name}\tLN:0\n")
        else:
            f.write(f"@SQ\tSN:{target_name}\tLN:{max_coord}\n")
        for idx, ai in enumerate(alignments_info, start=1):
            qname = ai.get('query_name', f"query_{idx}")
            pos = ai['target_start'] + 1  # SAM is 1-based
            seq = ai['aligned_query'].replace('-', '')
            cigar = make_cigar(ai['aligned_query'], ai['aligned_target'])
            # Basic flags: 0 (mapped forward)
            f.write(f"{qname}\t0\t{target_name}\t{pos}\t255\t{cigar}\t*\t0\t0\t{seq}\t*\n")
    print(f"[INFO] Exported SAM -> {file_path}")

# ---------------------------
# Performance comparison helpers
# ---------------------------
@timer
def run_traditional_sw_tilelike(query: str, target: str):
    # Run the smith_waterman_tile on the whole target (this is "traditional" in our comparison)
    return smith_waterman_tile(query, target)



def compare_performance(query: str, target: str, swift_fn, workers: Optional[int] = None):
    print("\n--- Performance comparison ---")
    # Traditional
    sw_score, sw_q_aln, sw_t_aln = run_traditional_sw_tilelike(query, target)
    print(f"Traditional SW (single-run) best score: {sw_score}")

    # SWIFT (parallel)
    swift_results = swift_fn(query, target) if workers is None else swift_fn(query, target, workers=workers)
    swift_best = swift_results[0][0] if swift_results else 0
    print(f"SWIFT best tile score: {swift_best}")
    return {
        'traditional_score': sw_score,
        'swift_best_score': swift_best,
        'swift_num_tiles': len(swift_results)
    }

# ---------------------------
# Main executable flow
# ---------------------------
if __name__ == "__main__":
    # Read sequences from test_data.txt (same simple format you provided)
    if not os.path.exists("test_data.txt"):
        raise FileNotFoundError("Put your sequences in test_data.txt (same format as your example)")

    with open("test_data.txt", "r") as f:
        lines = [l.strip() for l in f.readlines()]
    # expected format: >query, <sequence>, >target, <sequence>
    try:
        query_sequence = lines[1]
        target_sequence = lines[3]
    except Exception as e:
        raise ValueError("test_data.txt not in expected format. Use the >query / >target layout like in your original file.") from e

    # Parameters:
    K = 3
    TILE_SIZE = 64
    HIT_THRESHOLD = 2
    WORKERS = None  # None -> cpu_count()

    # 1) Run SWIFT in parallel
    raw_results = swift_alignment_parallel(query_sequence, target_sequence, k=K, tile_size=TILE_SIZE, hit_threshold=HIT_THRESHOLD, workers=WORKERS)

    # 2) Map tile string alignments to coordinates
    alignment_spans = []
    alignments_info = []
    for score, aligned_q, aligned_t, (tile_start, tile_end) in raw_results:
        q_start, q_end, t_start, t_end, sc = tile_alignment_to_coords(score, aligned_q, aligned_t, tile_start)
        alignment_spans.append((q_start, q_end, t_start, t_end, sc))
        alignments_info.append({
            'score': sc,
            'aligned_query': aligned_q,
            'aligned_target': aligned_t,
            'query_start': q_start,
            'query_end': q_end,
            'target_start': t_start,
            'target_end': t_end
        })

    # 3) Merge & refine
    refined = merge_and_refine_alignments(alignment_spans, query_sequence, target_sequence)

    # Build final info objects (include refined aligned strings by re-running SW on refined span for display)
    final_alignments = []
    for idx, (q_start, q_end, t_start, t_end, score) in enumerate(refined, start=1):
        # Re-run tile alignment on the refined region to retrieve aligned strings for output/visualization
        score2, a_q, a_t = smith_waterman_tile(query_sequence[q_start:q_end], target_sequence[t_start:t_end])
        # map local aligned strings back to full coordinates: first non-gap in local q + q_start etc.
        # Compute absolute coordinates for accurate JSON/SAM
        q_first = next((i for i,ch in enumerate(a_q) if ch != '-'), 0)
        t_first = next((i for i,ch in enumerate(a_t) if ch != '-'), 0)
        abs_q_start = q_start + q_first
        abs_t_start = t_start + t_first
        # remove leading and trailing gaps to produce SEQ for SAM
        final_alignments.append({
            'index': idx,
            'score': score2,
            'aligned_query': a_q,
            'aligned_target': a_t,
            'query_start': abs_q_start,
            'query_end': abs_q_start + sum(1 for ch in a_q if ch != '-'),
            'target_start': abs_t_start,
            'target_end': abs_t_start + sum(1 for ch in a_t if ch != '-'),
            'query_name': f"query_{idx}"
        })

    # 4) Export results
    if final_alignments:
        export_alignments_json(final_alignments, file_path="alignments.json")
        export_alignments_sam(final_alignments, target_name="target", file_path="alignments.sam")
    else:
        print("[INFO] No final alignments to export.")

    # 5) Visualize top alignment (ASCII + matplotlib)
    if final_alignments:
        top = final_alignments[0]
        print("\n=== Top Refined Alignment (ASCII) ===")
        print(ascii_alignment(top['aligned_query'], top['aligned_target']))
        print("\nPlotting alignment identity (matplotlib)...")
        plot_alignment_identity(top['aligned_query'], top['aligned_target'], title=f"Alignment Identity (score={top['score']})")
    else:
        # If no refined alignments, fallback to best raw tile alignment (if exists)
        if raw_results:
            score, a_q, a_t, span = raw_results[0]
            print("\n=== Best Raw Tile Alignment (ASCII) ===")
            print(ascii_alignment(a_q, a_t))
            plot_alignment_identity(a_q, a_t, title=f"Raw Tile Alignment Identity (score={score})")
        else:
            print("[INFO] No alignments found to visualize.")

    # 6) Performance comparison
    # perf_summary = compare_performance(query_sequence, target_sequence)
    perf_summary = compare_performance(query_sequence, target_sequence, swift_alignment_parallel)
    print("\nPerformance summary:", perf_summary)
    # print("\nPerformance summary:", perf_summary)

    print("\nDone. Files produced (if any): alignments.json, alignments.sam")
