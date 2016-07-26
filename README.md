# ParalogMatching.jl

| **Documentation**                       | **Build Status**                                                                                |
|:---------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |

This module implements the paralog matching technique presented in the paper
[Simultaneous identification of specifically interacting paralogs and
inter-protein contacts by Direct-Coupling Analysis](http://arxiv.org/abs/1605.03745)
by Thomas Gueudré, Carlo Baldassi, Marco Zamparo, Martin Weigt and Andrea Pagnani.

The main idea of the method is to perform a statistical analysis of two given
multiple sequence alignments, each containing one protein family. Each familiy should
comprise several species, and each species may have several sequences belonging to the
family. The algorithm tries to associate (match) interacting partners from the two families
within each species.

The underlying main assumption is that the proper matching is the one maximizing the
co-evolution signal. Such maximization is performed over the Bayesian inference of a
Gaussian model, by inverting the correlation matrix.

The code is written in [Julia](http://julialang.org), and the functions are called
from within Julia. However, a command-line interface is also provided for
those unfamiliar with the language (see the documentation).

The current code was tested on Julia versions 0.4 and 0.5.

## Installation

The package is not registered; it can be installed with `Pkg.clone`:

```
julia> Pkg.clone("https://github.com/Mirmu/ParalogMatching.jl")
```

Dependencies will be installed automatically.

However, you will also need to install at least one linear programming solver supported by
[MathProgBase](http://mathprogbasejl.readthedocs.io/en/latest/).
See the list of available solvers at the [JuliaOpt page](http://www.juliaopt.org/#packages).
Note that the solver efficiency is not particularly important for paralog matching, whose computational time
is dominated by matrix inversion operations, therefore you don't need a particularly fast solver. If unsure,
use `Pkg.add("Clp")` or `Pkg.add("GLPKMathProgInterface")`, which are free and open-source solvers.

## Documentation

- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

## Project Status

The package is tested against Julia `0.4` and *current* `0.5-dev` on Linux, OS X, and Windows.

## Contributing and Questions

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://Mirmu.github.io/ParalogMatching.jl/latest

[travis-img]: https://travis-ci.org/Mirmu/ParalogMatching.jl.svg?branch=master
[travis-url]: https://travis-ci.org/Mirmu/ParalogMatching.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/x9jkws1l4xd8q4wy/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/Mirmu/paralogmatching-jl/branch/master

[codecov-img]: https://codecov.io/gh/Mirmu/ParalogMatching.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/Mirmu/ParalogMatching.jl

[issues-url]: https://github.com/Mirmu/ParalogMatching.jl/issues
