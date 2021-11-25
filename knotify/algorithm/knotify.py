import pandas as pd

from knotify import rna_analysis
from knotify.algorithm.base import BaseAlgorithm
from knotify.criteria import apply_free_energy_and_stems_criterion
from knotify.energy.base import BaseEnergy
from knotify.energy.vienna import ViennaEnergy
from knotify import hairpin
from knotify.parsers.base import BaseParser


class Knotify(BaseAlgorithm):
    def get_results(
        self,
        sequence: str,
        parser: BaseParser,
        allow_skip_final_au: bool = False,
        csv: str = None,
        max_stem_allow_smaller: int = 1,
        prune_early: bool = False,
        hairpin_grammar: str = None,
        hairpin_allow_ug: bool = False,
        min_hairpin_size: int = hairpin.MIN_HAIRPIN_SIZE,
        min_hairpin_stems: int = hairpin.MIN_HAIRPIN_STEMS,
        max_hairpins_per_loop: int = hairpin.MAX_HAIRPINS_PER_LOOP,
        max_hairpin_bulge: int = hairpin.MAX_HAIRPIN_BULGE,
        energy: BaseEnergy = ViennaEnergy(),
        *args,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Analyze RNA sequence, and predict structure. Return data frame of results
        """
        sequence = sequence.lower()
        knot_dict_list = rna_analysis.StringAnalyser(
            input_string=sequence,
            parser=parser,
        ).get_pseudoknots(
            max_stem_allow_smaller=max_stem_allow_smaller,
            prune_early=prune_early,
            allow_skip_final_au=allow_skip_final_au,
        )

        data = pd.DataFrame(knot_dict_list)
        if csv is not None:
            data.to_csv(csv)

        data = apply_free_energy_and_stems_criterion(
            data,
            sequence,
            max_stem_allow_smaller=max_stem_allow_smaller,
            energy=energy,
        )

        if hairpin_grammar is None:
            return data

        data = hairpin.find_hairpins(
            sequence,
            data,
            hairpin_grammar,
            allow_ug=hairpin_allow_ug,
            min_size=min_hairpin_size,
            min_stems=min_hairpin_stems,
            max_bulge=max_hairpin_bulge,
            max_per_loop=max_hairpins_per_loop,
        )
        data = apply_free_energy_and_stems_criterion(
            data,
            sequence,
            max_stem_allow_smaller=max_stem_allow_smaller,
            energy=energy,
        )

        return data
