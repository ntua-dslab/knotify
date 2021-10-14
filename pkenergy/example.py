#!/usr/bin/env python3

import ctypes

c = ctypes.CDLL("./hotknots/LE/libpkenergy.so")

c.get_energy.restype = ctypes.c_double


def cstr(s: str):
    return ctypes.c_char_p(s.encode())


c.initialize(cstr("hotknots/params"))

for x in range(1000):
    print(
        c.get_energy(
            cstr(
                "AAUGCAACUUUUAAAUAGUUUAUCUGUUAAGAUAAACCACCUAGGUUGCAUAUAUAAAAAAUAAAAGGUGCC"
            ),
            cstr(
                ".(((((((((...........................[[[[[)))))))))..............]]]]].."
            ),
        )
    )
