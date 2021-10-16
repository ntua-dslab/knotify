class BaseParser:
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
        raise NotImplementedError
