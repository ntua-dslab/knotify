class BaseEnergy:
    """
    Base class for psuedoknot MFE calculation.
    """

    def eval(self, sequence: str, dot_bracket: str) -> float:
        raise NotImplementedError
