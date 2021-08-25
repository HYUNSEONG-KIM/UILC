# UILC
This repostiory contains calculation code for building linear and rectangular LEDs array for uniform irradiation pattern for given area. 

### Language Frameworks
**Language** 
* python3
* C
* Wolfram Language

**pre-request library for each implementation**

1. python<br>
* numpy
* scipy
* matplotlib

2. Wolframe language: complete itself

3. C

* GSL (Gnu Scientific Library)

**Warning**

The GNU Scientific Library is distributed under GNU General Public License. Following this License, the C implementation is distributed under GPL license v 3.0.




### /src

* uic.py : python implementation
* uic.wl : Wolfram language implementation
* uic.h, uic.c: C implementation


## Description of routines

There are sevaral differences in three types of implementations.

The python type covers most wide region of routines in papers.
### Python implementation

---
> esc
---

## Refereneces and Further Reading

If you are interesting in building flat, uniform light device from LED units, Next papers are recommended

This routines are based on "". Actually, it is written for such paper. This routines are calculates boundary propagation with fit boundary intensity to center point. This routines can be combined with multiple solutions not only analytical methods(e.g Sparrow's criterion) but also various stochastic methods.

