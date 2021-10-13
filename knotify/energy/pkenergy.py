import ctypes


def cstr(s: str) -> str:
    return ctypes.c_char_p(s.encode())


class PKEnergy:
    """
    Load MFE calculator from a dynamic library. The library should export:

    ```
    // will be called once during instantiation of the class. use to load any
    // configuation files (e.g. parameters).
    void initialize(char *config_dir);

    // will be called for each sequence.
    float get_energy(char *sequence, char *structure);
    ```
    """

    def __init__(self, library: str, config_dir: str):
        self._lib = ctypes.CDLL(library)
        self._lib.get_energy.restype = ctypes.c_double
        self._lib.initialize(ctypes.c_char_p(config_dir.encode()))

    def energy_eval(self, sequence: str, dot_bracket: str) -> float:
        return self._lib.get_energy(
            ctypes.c_char_p(sequence.encode()),
            ctypes.c_char_p(dot_bracket.encode()),
        )
