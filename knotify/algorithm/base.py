import pandas as pd


class BaseAlgorithm:
    """
    Base class for all algorithm implementations
    """

    def get_results(self, sequence: str, *args, **kwargs) -> pd.DataFrame:
        raise NotImplementedError
