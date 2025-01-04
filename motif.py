from replication_ori import hamming_distance, neighbors

# IMPLANTED MOTIF PROBLEM: Find all (k, d)-motifs in a collection of strings
    # Input: A collection of strings Dna, and integers k and d
    # Output: All (k, d)-motifs in Dna
# (k,d)-motif: given a collection of strings Dna and an integer d, a k-mer that appears in every string from Dna with at most d mismatches
# Brute force (i.e. exhaustive search) is an easy algorithm to implement but takes enormous amounts of time
"""
MotifEnumeration(Dna, k, d)
    Patterns ← an empty set
    for each k-mer Pattern in Dna
        for each k-mer Pattern' differing from Pattern by at most d mismatches
            if Pattern' appears in each string from Dna with at most d mismatches
                add Pattern' to Patterns
    remove duplicates from Patterns
    return Patterns
"""
def motif_enumeration(dna: list[str], k: int, d: int) -> list[str]:
    patterns = set()
    # For each string in DNA collection
    for string in dna:
        # For each k-mer in this string
        for i in range(len(string) - k + 1):
            pattern = string[i:i+k]
            # Get all possible k-mers within d mismatches
            neighborhood = neighbors(pattern, d)
            # For each potential pattern
            for pattern_prime in neighborhood:
                # Check if it appears in all strings with at most d mismatches
                appears_in_all = True
                for text in dna:
                    found_in_text = False
                    # Look for pattern_prime in this text with up to d mismatches
                    for j in range(len(text) - k + 1):
                        if hamming_distance(pattern_prime, text[j:j+k]) <= d:
                            found_in_text = True
                            break
                    if not found_in_text:
                        appears_in_all = False
                        break
                if appears_in_all:
                    patterns.add(pattern_prime)
    return list(patterns)

print(motif_enumeration(["ATTTGGC", "TGCCTTA", "CGGTATC", "GAAAATT"], 3, 1))