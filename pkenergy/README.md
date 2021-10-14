# HotKnots

This repository contains parts of the [HotKnots v2.0](http://www.cs.ubc.ca/labs/beta/Software/HotKnots/) source code.

Additionally, it has been edited to build as a shared library that can be loaded from Python code to calculate PK energy.

I do not claim any ownership on the HotKnots original source code, neither on the parameters present in this repository used to calculate energy.

Work from this repository is used in [knotify](https://github.com/ntua-dslab/knotify.git).

## Build

```bash
(cd hotknots/LE && make -j)
```

For usage, you need `./hotknots/LE/libpkenergy.so`, `./hotknots/params` and `example.py`.
