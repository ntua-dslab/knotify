import RNA


class ViennaEnergy:
    def __init__(self, **kwargs):
        pass

    def energy_eval(self, sequence: str, dot_bracket: str) -> float:
        return RNA.energy_of_struct(sequence, dot_bracket)
