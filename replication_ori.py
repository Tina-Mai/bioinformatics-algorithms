# HIDDEN MESSAGE PROBLEM: find a “hidden message” (frequently recurring strings) in the replication origin
# k-mer: a substring of length k
# Count(Text, Pattern): number of times that a k-mer Pattern appears as a substring of Text
def pattern_count(text: str, pattern: str) -> int:
    count = 0
    for i in range(0, len(text) - len(pattern) + 1):
        if text[i:i + len(pattern)] == pattern:
            count += 1
    return count


# FREQUENT WORDS PROBLEM: Find the most frequent k-mers in a string
# most frequent k-mer: the Pattern that maximizes Count(Text, Pattern) in Text
# FrequencyTable(Text, k): a map of k-mers to their counts in the text
# FrequentWords(Text, k): a list of the most frequent k-mers in the text
def frequency_table(text: str, k: int) -> dict[str, int]:
    freq_map = {}
    n = len(text)
    for i in range(0, n - k + 1):
        pattern = text[i:i + k]
        if pattern not in freq_map:
            freq_map[pattern] = 1
        else:
            freq_map[pattern] += 1
    return freq_map

def frequent_words(text: str, k: int) -> list[str]:
    frequent_patterns = []
    freq_map = frequency_table(text, k)
    max_count = max(freq_map.values())
    for pattern in freq_map:
        if freq_map[pattern] == max_count:
            frequent_patterns.append(pattern)
    return frequent_patterns