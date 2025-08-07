import random

def generate_random_dna_sequence(length: int) -> str:
    """
    Generate a random DNA sequence of a given length.
    """
    return ''.join(random.choices('ACGT', k=length))


if __name__ == "__main__":
    # Generate a realistic query (1,000 bp) and target (5,000 bp)
    long_query = generate_random_dna_sequence(1000)
    long_target = (
        generate_random_dna_sequence(2000) +
        long_query[200:700] +  # Insert partial match from query into target
        generate_random_dna_sequence(2300)
    )

    # save in txt file
    with open("test_data.txt", "w") as f:
        f.write(f">query\n{long_query}\n")
        f.write(f">target\n{long_target}\n")