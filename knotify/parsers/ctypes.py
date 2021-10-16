import ctypes

from knotify.parsers.base import BaseParser


class CTypesParser(BaseParser):
    """
    Base class for ctypes-based parsers. The library_path argument should be
    a shared library exposing the following functions with C linkage:

    ```
    // Will be called once to initialize the library. options is a string
    // with implementation-specific configuration.
    void initialize(
        char *options,
        int allow_ug,
        int min_dd_size,
        int max_dd_size,
        int min_window_size,
        int max_window_size
    );

    // Parse RNA sequence and return positions of possible core stems pseudoknots.
    // Returned as a char *buffer with the following format:
    //
    // ```
    // <left1>,<size1>,<leftloopsize1>,<ddsize1>
    // <left2>,<size2>,<leftloopsize2>,<ddsize2>
    // ...
    // ```
    char *detect_pseudoknots(char *sequence);
    ```
    """

    def get_options(self) -> str:
        raise NotImplementedError

    def __init__(self, library_path: str, *args, **kwargs):
        super(CTypesParser, self).__init__(*args, **kwargs)

        self.lib = ctypes.CDLL(library_path)
        self.lib.detect_pseudoknots.restype = ctypes.c_char_p
        self.lib.initialize(
            ctypes.c_char_p(self.get_options().encode()),
            ctypes.c_int(self.allow_ug),
            ctypes.c_int(self.min_dd_size),
            ctypes.c_int(self.max_dd_size),
            ctypes.c_int(self.min_window_size),
            ctypes.c_int(self.max_window_size),
        )

    def detect_pseudoknots(self, sequence: str) -> list:
        return (
            self.lib.detect_pseudoknots(ctypes.c_char_p(sequence.lower().encode()))
            .decode()
            .rstrip("\n")
            .split("\n")
        )
