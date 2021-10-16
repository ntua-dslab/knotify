class BaseParser:
    """
    Base RNA parser class. Defines the interface for using the RNA parser, and
    common configuration options for all parsers. Parser implementations are expected
    to read and respect these configurations.

    Implementations should
    """

    def __init__(
        self,
        max_dd_size: int = 2,
        allow_ug: bool = False,
        min_dd_size: int = 0,
        max_window_size: int = 100,
        min_window_size: int = 6,
    ):
        self.allow_ug = allow_ug
        self.max_dd_size = max_dd_size
        self.min_dd_size = min_dd_size
        self.max_window_size = max_window_size
        self.min_window_size = min_window_size

    def detect_pseudoknots(self, sequence: str) -> list:
        """
        Return a list of pseudoknots for the given RNA sequence (and any subsequences).
        The format should be:
        [
            "<left1>,<size1>,<left_loop_size1>,<dd_size1>",
            "<left2>,<size2>,<left_loop_size2>,<dd_size2>",
            ...
        ]
        """
        raise NotImplementedError
