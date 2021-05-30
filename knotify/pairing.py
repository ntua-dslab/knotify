# import alignment methodologies
from knotify.pairalign import consecutive_pairalign as cpairalign


def get_left_stem_aligned_indices(
    string, inner_potential_indices, outer_potential_indices
):
    """Return the indices of the aligned strings

    :param str string: the original string
    :param tuple inner_potential_indices: [start, end) of inner part
    :param tuple outer_potential_indices: [start, end] of outer part

    :return: (outer_part, inner_part), (outer_stems, inner_stems)
    :rtype: tuple(list, list), tuple(list, list)
    """
    # get substring of inner outer parts
    inner_string = string[inner_potential_indices[0] : inner_potential_indices[1]]
    outer_string = string[outer_potential_indices[0] : outer_potential_indices[1]]

    prune_indices, alignments = cpairalign(outer_string, inner_string)
    # check if there is no alignment at all
    if not alignments or not alignments[0]:
        return (
            (
                list(range(*inner_potential_indices)),
                list(range(*outer_potential_indices)),
            ),
            ([], []),
        )

    # deflate alignments to the the corresponding variables
    outer_alignment, inner_alignment = alignments

    # deflate indices to the corresponding variables
    relative_outer_first_index, relative_inner_last_index = prune_indices

    # we need to update them to the global string indexing space
    outer_first_index = outer_potential_indices[0] + relative_outer_first_index
    inner_last_index = inner_potential_indices[0] + relative_inner_last_index

    # let's find all matching indices
    (
        outer_free_indices,
        inner_free_indices,
        outer_bound_indices,
        inner_bound_indices,
    ) = get_bound_indices(
        outer_alignment, inner_alignment
    )  # noqa

    # get and update indices based on any related offset
    free_inner_indices = [
        i + inner_potential_indices[0] for i in inner_free_indices
    ]  # noqa
    free_outer_indices = [
        i + outer_potential_indices[0] for i in outer_free_indices
    ]  # noqa
    bound_inner_indices = [
        i + inner_potential_indices[0] for i in inner_bound_indices
    ]  # noqa
    bound_outer_indices = [i + outer_first_index for i in outer_bound_indices]  # noqa

    # append any defacto free index
    free_outer_indices.extend(range(outer_potential_indices[0], outer_first_index))
    free_inner_indices.extend(range(inner_last_index + 1, inner_potential_indices[1]))

    return (
        (free_inner_indices, free_outer_indices),
        (bound_inner_indices, bound_outer_indices),
    )  # noqa


def get_right_stem_aligned_indices(
    string, inner_potential_indices, outer_potential_indices
):
    """Return the indices of the aligned strings

    :param str string: the original string
    :param tuple inner_potential_indices: [start, end) of inner part
    :param tuple outer_potential_indices: [start, end] of outer part

    :return: (outer_part, inner_part), (outer_stems, inner_stems)
    :rtype: tuple(list, list), tuple(list, list)
    """

    # get substring of inner outer parts
    inner_string = string[
        inner_potential_indices[0] : inner_potential_indices[1]
    ]  # noqa
    outer_string = string[
        outer_potential_indices[0] : outer_potential_indices[1]
    ]  # noqa
    # this is the fucking catch here!
    # Caution: lower level MAGIC NUMBERS are following
    prune_indices, alignments = cpairalign(inner_string, outer_string)

    if not alignments or not alignments[0]:  # no alignmens have been found
        return (
            (
                list(range(*inner_potential_indices)),
                list(range(*outer_potential_indices)),
            ),
            ([], []),
        )

    relative_inner_first_index, relative_outer_last_index = prune_indices

    # deflate alignments to the corresponding variables
    try:
        inner_alignment, outer_alignment = alignments
    except:
        import ipdb

        ipdb.set_trace()

    # update indices to the global string indexing mapping
    inner_first_index = inner_potential_indices[0] + relative_inner_first_index
    outer_last_index = outer_potential_indices[0] + relative_outer_last_index

    (
        outer_free_indices,
        inner_free_indices,
        outer_bound_indices,
        inner_bound_indices,
    ) = get_bound_indices(
        outer_alignment, inner_alignment
    )  # noqa

    # get and update indices based on any related offset
    free_inner_indices = [
        i + inner_potential_indices[0] for i in inner_free_indices
    ]  # noqa
    free_outer_indices = [
        i + outer_potential_indices[0] for i in outer_free_indices
    ]  # noqa
    bound_inner_indices = [i + inner_first_index for i in inner_bound_indices]  # noqa
    bound_outer_indices = [
        i + outer_potential_indices[0] for i in outer_bound_indices
    ]  # noqa

    # append any defacto free index
    free_outer_indices.extend(range(outer_last_index + 1, outer_potential_indices[1]))
    free_inner_indices.extend(range(inner_potential_indices[0], inner_first_index))

    return (
        (free_inner_indices, free_outer_indices),
        (bound_inner_indices, bound_outer_indices),
    )  # noqa


def get_bound_indices(str1, str2):
    """Returns the aligned and not bound indices of str1, str2

    :param str str1: the first string to align
    :param str str2: the second string to align

    Caution: order matters, a lot! --> str1 is going to be traversed the other way around  # noqa

    :return: the free along with the bound indices
    :rtype: tuple(list, list), tuple(list, list)
    """
    free_str1_indices = []
    free_str2_indices = []
    bound_str1_indices = []
    bound_str2_indices = []

    # we need to traverse str1 reversed and use the reversed indices
    for i in range(0, len(str1)):
        if str1[i] != "-":
            bound_str1_indices.append(i)
        else:
            free_str1_indices.append(i)

    for i in range(0, len(str2)):
        if str2[i] != "-":
            bound_str2_indices.append(i)
        else:
            free_str2_indices.append(i)

    return (
        free_str1_indices,
        free_str2_indices,
        bound_str1_indices,
        bound_str2_indices,
    )  # noqa
