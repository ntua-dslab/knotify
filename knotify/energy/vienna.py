import RNA

from knotify.energy.base import BaseEnergy


class ViennaEnergy(BaseEnergy):
    """
    Calcuate MFE using the Vienna RNA library.
    """

    def eval(self, sequence: str, dot_bracket: str) -> float:
        return RNA.energy_of_struct(sequence, dot_bracket)
