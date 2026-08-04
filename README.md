# Julia-Scientific

> Numerical analysis and applied statistics written in Julia for coursework.

<!-- badges -->
![Julia](https://img.shields.io/badge/Julia-9558B2?style=for-the-badge&logo=julia&logoColor=white)

![Express](https://img.shields.io/badge/Express-000000?style=for-the-badge&logo=express&logoColor=white)

![last commit](https://img.shields.io/github/last-commit/adzetto/Julia-Scientific?style=flat-square&color=informational) ![repo size](https://img.shields.io/github/repo-size/adzetto/Julia-Scientific?style=flat-square&color=informational) ![top language](https://img.shields.io/github/languages/top/adzetto/Julia-Scientific?style=flat-square) ![language count](https://img.shields.io/github/languages/count/adzetto/Julia-Scientific?style=flat-square) ![license](https://img.shields.io/github/license/adzetto/Julia-Scientific?style=flat-square&color=informational)

---

## What this is

Julia scripts written for a numerical analysis course, plus a set of applied
statistics algorithms. Fifty-three `.jl` files in two groups.

| Directory | Contents |
|---|---|
| `JuliaScientificSimulations` | Numerical methods and simulation scripts |
| `Applied Statistics Algorithms` | Statistical procedures implemented from scratch |

The point of writing these from scratch rather than calling a package is that
the failure modes become visible: where an iteration stops converging, where
conditioning ruins a result, where a textbook formula is numerically the wrong
way to compute a correct expression.

## Running

```bash
julia script_name.jl
```

Or from the REPL, which is faster when you are iterating because the JIT warm-up
is paid once:

```julia
include("script_name.jl")
```

If a script needs a package it will say so on the first run. Add it with:

```julia
using Pkg; Pkg.add("PackageName")
```

## Licence

See `LICENSE`.

## Repository layout

```text
Julia-Scientific/
├── Applied Statistics Algorithms/
├── JuliaScientificSimulations/
├── LICENSE
├── README.md
```

---

**59** tracked files · **0.2 MB** · **1** languages · last pushed **2024-09-25**
