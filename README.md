# SCL-ABE: Adaptively Secure Lattice Attribute-Based Encryption with Sub-Linear Ciphertexts

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24214538.svg)](https://doi.org/10.5281/zenodo.24214538)

Code and data for the paper *SCL-ABE: Adaptively Secure Lattice Attribute-Based Encryption with Sub-Linear Ciphertexts* (Khaled Riad, King Faisal University).
Repository: [github.com/khaledxriad/SCL-ABE](https://github.com/khaledxriad/SCL-ABE). Archived release (v1.0.0): Zenodo, [doi:10.5281/zenodo.24214538](https://doi.org/10.5281/zenodo.24214538).

This repository holds the code, the stored data, and the independent check of every reported number.
The manuscript and its TikZ figure sources are not part of this repository.

## Layout

| Path | Content |
|---|---|
| `code/est.py` | parameter solver and primal-uSVP core-SVP estimator (0.292 beta bits, classical) |
| `code/cover.py` | complete-subtree cover and cover-size statistics |
| `code/sizes.py` | object sizes (pk, msk, keys, node ciphertext, linear baseline) for N = 2^6 ... 2^20 |
| `code/proto.py` | prototype: dual-Regev node IBE with gadget trapdoor, the ABE layer and the DEM |
| `code/bench.py` | prototype timings at toy dimensions n in {8, 16, 24, 32} |
| `data/` | `params.json`, `cover_stats.csv`, `bench.json`, `bench.log` (the stored results behind the reported tables) |
| `tests/verify_results.py` | independent re-derivation of the reported numbers (does not import `code/`) |

## Quick start

```
pip install -r requirements.txt        # numpy
make verify                            # prints PASS/FAIL per check; ends with ALL CHECKS PASSED
```

Recomputing the data from scratch (`make params`, `make cover`, `make bench`) takes longer and overwrites `data/`.
Timings depend on the machine; the paper's timings were measured on one core of a 2.1 GHz Intel Xeon.

## Scope and caveats

* The prototype runs at **toy dimensions** (n up to 32, power-of-two modulus 2^62, rounded continuous Gaussians). It is a correctness and scaling test, not a deployable implementation. Sizes for full-scale parameters are computed from formulas, not measured.
* The security numbers are a core-SVP primal-attack estimate (0.292 beta); they exclude hybrid, dual and algebraic refinements.
* The linear-size baseline is an analytical model with the same lattice parameters, not an implementation of a published scheme.

## License

MIT (see `LICENSE`).


## How to cite

If you use this code or data, please cite the archived release:

Khaled Riad. *SCL-ABE: Adaptively Secure Lattice Attribute-Based Encryption with Sub-Linear Ciphertexts* (v1.0.0). Zenodo, 2026. https://doi.org/10.5281/zenodo.24214538

A citation for the paper will be added once it is published.
