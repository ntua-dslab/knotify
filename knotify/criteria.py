import pandas as pd

from knotify.energy.vienna import ViennaEnergy


def apply_free_energy_and_stems_criterion(
    data: pd.DataFrame,
    sequence: str,
    max_stem_allow_smaller: int = 2,
    energy_eval: callable = ViennaEnergy().energy_eval,
):
    """
    Returns the best result for a Pandas data frame.

    Expected shape for the Pandas DataFrame rows:
    {
        "left_loop_stems": len(knot.left_loop_stems[0]),
        "right_loop_stems": len(knot.right_loop_stems[0]),
        "dot_bracket": visualize_knot(knot),
        "dd": str(len(knot.get_dd_seg()))
    }
    """
    # TODO(akolaitis): Consider supporting "strategies" where smaller number of stems
    # are allowed but are discarded based on energy. Also consider favoring pseudoknots
    # with larger loop sizes, which allows for more hairpin/bulge options.
    #
    # Note that calculating energy for all dot brackets is a very expensive
    # operation. However, we have seen dot brackets with few stems that give
    # very good energy results (and are part of the ground truth, because of
    # bulges and/or hairpins in the pseudoknot loops). For an example, see case
    # with identifier "knotify/case/2021/paper1/bulges/1".
    #
    # As a future reference, benchmarking energy calculation:
    # Sequence: guuucuccuaucgccaucuggaugggauuuaagagacuuaugccaaacucuuugaguuugaugccaauucaguuuucgucugaauugaugcccgaaaggauggc    #  noqa
    # calculating energy for 85554 rows, finished after 7.016129016876221 sec
    # calculating energy for 2949 rows, finished after 0.2533254623413086 sec

    # max stems
    data["stems"] = data["left_loop_stems"] + data["right_loop_stems"]
    data = data[
        data["stems"] >= data["stems"].max() - max_stem_allow_smaller
    ].reset_index()

    # min energy
    data["energy"] = data["dot_bracket"].apply(lambda r: energy_eval(sequence, r))
    data.sort_values(
        ["energy", "stems", "dd"], ascending=(True, False, True), inplace=True
    )

    data = data.reset_index()

    return data
