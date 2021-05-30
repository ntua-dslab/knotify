def consecutive_pairalign(reversed_part, identical_part):
    count = 0
    lrev = len(reversed_part)
    no_iter = min(len(reversed_part), len(identical_part))
    i = 0
    for i in range(0, no_iter):
        if (reversed_part[lrev - i - 1] + identical_part[i]) in [
            "au",
            "ua",
            "gc",
            "cg",
            "ug",
            "gu",
        ]:
            count += 1
        else:
            break
    return (lrev - count, i), (
        reversed_part[lrev - count :],
        identical_part[:count],
    )  # noqa
