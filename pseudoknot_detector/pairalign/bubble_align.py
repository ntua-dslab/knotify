from pseudoknot_detector.consts import ACCEPTABLE_PAIRS as PAIRS

# score should be proportional to bubble distance and length somehow ....


def string_align(str1, str2):
    """Given two strings it aligns them in
    a custom way to minimize holes and
    hole lengths
    """

    # length align the two strings
    common_length()

    align1 = []
    align2 = []
    for i, j, k in enumerate(zip(str1, str2)):
        if i == j:
            align1.append(i)
            align2.append(j)
        else:
            next_candidates = [
                (str1[k:], str2[k + 1 :]),
                (str1[k + 1 :], str2[k:]),
                (str1[k + 1 :], str2[k + 1 :]),
            ]
            partial_alignments = [
                string_align(*substrings) for substrings in next_candidates
            ]

            return (
                "".join(align1) + partial_alignments[0][1],
                "".join(align2) + partial_alignments[0][1],
            )


# split the string into segments
