# HIDDEN MESSAGE PROBLEM: find a “hidden message” (frequently recurring strings) in the replication origin
    # Input: A string Text (representing the replication origin of a genome)
    # Output: A hidden message in Text
# k-mer: a substring of length k
# Count(Text, Pattern): number of times that a k-mer Pattern appears as a substring of Text
def pattern_count(text: str, pattern: str) -> int:
    count = 0
    # Slide through the DNA sequence one base at a time
    for i in range(0, len(text) - len(pattern) + 1):
        # Check if we found our DNA pattern at this position
        if text[i:i + len(pattern)] == pattern:
            count += 1
    return count


# FREQUENT WORDS PROBLEM: Find the most frequent k-mers in a string
    # Input: A string Text and an integer k
    # Output: All most frequent k-mers in Text
# most frequent k-mer: the Pattern that maximizes Count(Text, Pattern) in Text
# FrequencyTable(Text, k): a map of k-mers to their counts in the text
# FrequentWords(Text, k): a list of the most frequent k-mers in the text
def frequency_table(text: str, k: int) -> dict[str, int]:
    freq_map = {}
    n = len(text)
    # Look at every possible k-length DNA sequence in the text
    for i in range(0, n - k + 1):
        # Extract the k-mer (DNA substring of length k)
        pattern = text[i:i + k]
        # Count how many times we see each k-mer
        if pattern not in freq_map:
            freq_map[pattern] = 1
        else:
            freq_map[pattern] += 1
    return freq_map

def frequent_words(text: str, k: int) -> list[str]:
    frequent_patterns = []
    # Get counts of all k-mers in the DNA sequence
    freq_map = frequency_table(text, k)
    # Find which k-mers appear most often (could be important regulatory sequences)
    max_count = max(freq_map.values())
    for pattern in freq_map:
        if freq_map[pattern] == max_count:
            frequent_patterns.append(pattern)
    return frequent_patterns


# REVERSE COMPLEMENT PROBLEM: Find the reverse complement of a DNA string
    # Input: A DNA string Pattern
    # Output: Pattern, the reverse complement of Pattern_rc
def reverse_complement(pattern: str) -> str:
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(pattern))


# PATTERN MATCHING PROBLEM: Find all occurrences of a pattern in a string
    # Input: Strings Pattern and Genome
    # Output: All starting positions in Genome where Pattern appears as a substring
def pattern_matching(pattern: str, genome: str) -> list[int]:
    positions = []
    # Search through the genome for a specific DNA sequence
    for i in range(0, len(genome) - len(pattern) + 1):
        # Record where we find matches (could be binding sites)
        if genome[i:i + len(pattern)] == pattern:
            positions.append(i)
    return positions


# Excercise: Return a space-separated list of starting positions (in increasing order) where CTTGATCAT appears as a substring in the Vibrio cholerae genome
def exercise_1():
    pattern = "CTTGATCAT"
    genome = open("genomes/vibrio_cholerae.txt").read()
    positions = pattern_matching(pattern, genome)
    return " ".join(map(str, positions))

print(exercise_1())


# CLUMP FINDING PROBLEM: Find patterns forming clumps in a string
    # Input: A string Genome, and integers k, L, and t
    # Output: All distinct k-mers forming (L, t)-clumps in Genome
# (L, t)-clump: a k-mer Pattern that appears in at least t different locations within every interval of length L in the genome
def find_clumps(genome: str, k: int, l: int, t: int) -> list[str]:
    patterns = set()
    n = len(genome)
    if l > n:
        return list(patterns)
    
    # Count k-mers in the first window of length l
    freq_map = {}
    for i in range(0, l - k + 1):
        pattern = genome[i:i + k]
        freq_map[pattern] = freq_map.get(pattern, 0) + 1
        # If we see a k-mer enough times, it might be biologically significant
        if freq_map[pattern] >= t:
            patterns.add(pattern)
    
    # Slide the window along the genome
    for i in range(1, n - l + 1):
        # Update counts as we move the window
        first_pattern = genome[i-1:i-1 + k]
        freq_map[first_pattern] -= 1
        
        last_pattern = genome[i+l-k:i+l]
        freq_map[last_pattern] = freq_map.get(last_pattern, 0) + 1
        
        # Track k-mers that appear frequently in this window
        if freq_map[last_pattern] >= t:
            patterns.add(last_pattern)
    
    return list(patterns)


# Excercise: How many different 9-mers form (500,3)-clumps in the E. coli genome? (In other words, do not count a 9-mer more than once.)
def exercise_2():
    genome = open("genomes/e_coli.txt").read()
    k = 9
    l = 500
    t = 3
    clumps = find_clumps(genome, k, l, t)
    return len(clumps)

# print(exercise_2())


# MINIMUM SKEW PROBLEM: Find a position in a genome where the skew diagram attains a minimum
    # Input: A DNA string Genome
    # Output: All integer(s) i minimizing Skew(Genome)[i] among all values of i (from 0 to |Genome|)
def minimum_skew(genome: str) -> list[int]:
    skew = 0
    min_skew = 0
    positions = []
    # Track G-C composition along the genome
    for i in range(0, len(genome)):
        # C decreases skew, G increases it
        # This helps find the DNA replication origin
        if genome[i] == 'C':
            skew -= 1
        elif genome[i] == 'G':
            skew += 1
        # Keep track of positions where G-C skew is minimal
        if skew < min_skew:
            min_skew = skew
            positions = [i + 1]
        elif skew == min_skew:
            positions.append(i + 1)
    return positions


# HAMMING DISTANCE PROBLEM: Compute the Hamming distance between two strings
    # Input: Two strings of equal length
    # Output: The Hamming distance between these strings
# Hamming distance: number of mismatches between strings p and q; denoted HammingDistance(p, q)
def hamming_distance(p: str, q: str) -> int:
    distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance


# APPROXIMATE PATTERN MATCHING PROBLEM: Find all approximate occurrences of a pattern in a string
    # Input: Strings *Pattern* and *Text* along with an integer *d*
    # Output: All starting positions where *Pattern* appears as a substring of *Text* with at most *d* mismatches
def approximate_pattern_matching(pattern: str, text: str, d: int) -> list[int]:
    positions = []
    for i in range(0, len(text) - len(pattern) + 1):
        if hamming_distance(pattern, text[i:i + len(pattern)]) <= d:
            positions.append(i)
    return positions


# APPROXIMATE PATTERN COUNT PROBLEM: Count the occurrences of a *Pattern* in a *Text*, allowing for up to *d* mismatches
    # Input: Strings *Pattern* and *Text* as well as an integer *d*
    # Output: Count_d(Text, Pattern)
def approximate_pattern_count(text: str, pattern: str, d: int) -> int:
    count = 0
    for i in range(0, len(text) - len(pattern) + 1):
        if hamming_distance(pattern, text[i:i + len(pattern)]) <= d:
            count += 1
    return count


# Exercise: Compute Count_2(AACAAGCTGATAAACATTTAAAGAG, AAAAA)
# Count_d(Text, Pattern) is the number of times that a k-mer Pattern appears as a substring of Text with at most d mismatches
def exercise_3():
    text = "AACAAGCTGATAAACATTTAAAGAG"
    pattern = "AAAAA"
    d = 2
    return approximate_pattern_count(text, pattern, d)

print(exercise_3())


# FREQUENT WORDS WITH MISMATCHES PROBLEM: Find the most frequent k-mers with mismatches in a string
    # Input: A string *Text* as well as integers *k* and *d*
    # Output: All most frequent *k-*mers with up to *d* mismatches in *Text*
# d-neighborhood: all k-mers within Hamming distance d from Pattern, denoted Neighbors(Pattern, d)
def neighbors(pattern: str, d: int) -> set[str]:
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {'A', 'C', 'G', 'T'}
    
    # Recursively generate all possible DNA sequences
    # within d mutations of our pattern
    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)
    
    for text in suffix_neighbors:
        if hamming_distance(pattern[1:], text) < d:
            # Try all possible nucleotides if we can still make changes
            for nucleotide in ['A', 'C', 'G', 'T']:
                neighborhood.add(nucleotide + text)
        else:
            # Keep original nucleotide if we've used all our mutations
            neighborhood.add(pattern[0] + text)
            
    return neighborhood

def frequent_words_with_mismatches(text: str, k: int, d: int) -> list[str]:
    patterns = []
    freq_map = {}
    n = len(text)
    
    # Look for k-mers that appear frequently, allowing for up to d mutations
    for i in range(0, n - k + 1):
        pattern = text[i:i + k]
        # Get all similar sequences within d mutations
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            freq_map[neighbor] = freq_map.get(neighbor, 0) + 1
    
    # Find which patterns (including mutations) appear most often
    if freq_map:
        max_count = max(freq_map.values())
        for pattern in freq_map:
            if freq_map[pattern] == max_count:
                patterns.append(pattern)
                
    return patterns


# FREQUENT WORDS WITH MISMATCHES AND REVERSE COMPLEMENTS PROBLEM: Find the most frequent k-mers (with mismatches and reverse complements) in a string
    # Input: A DNA string *Text* as well as integers k and d
    # Output: All k-mers Pattern maximizing the sum Count_d(Text, Pattern) + Count_d(Text, Pattern_{rc}) over all possible k-mers
def frequent_words_mismatches_reverse_complements(text: str, k: int, d: int) -> list[str]:
    patterns = []
    freq_map = {}
    n = len(text)
    
    # Look at each k-mer in the text
    for i in range(0, n - k + 1):
        pattern = text[i:i + k]
        # Get all similar sequences within d mutations
        neighborhood = neighbors(pattern, d)
        for neighbor in neighborhood:
            # Count both the neighbor and its reverse complement
            rc = reverse_complement(neighbor)
            freq_map[neighbor] = freq_map.get(neighbor, 0) + 1
            freq_map[rc] = freq_map.get(rc, 0) + 1
    
    # Find patterns with maximum combined frequency
    if freq_map:
        max_count = max(freq_map.values())
        for pattern in freq_map:
            if freq_map[pattern] == max_count:
                patterns.append(pattern)
    
    return patterns
