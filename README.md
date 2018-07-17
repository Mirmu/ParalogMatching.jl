# ParalogMatching.jl

| **Documentation**                       | **Build Status**                                                                                |
|:---------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][codecov-img]][codecov-url] |

This package implements the paralog matching technique presented in the paper
*"Simultaneous identification of specifically interacting paralogs and
interprotein contacts by Direct Coupling Analysis"*
by Thomas Gueudré, Carlo Baldassi, Marco Zamparo, Martin Weigt and Andrea Pagnani,
Proc. Natl. Acad. Sci. U.S.A. 113, 12186–12191 (2016), [doi:10.1073/pnas.1607570113][paper].

The main idea of the method is to perform a statistical analysis of two given
multiple sequence alignments, each containing one protein family. Each familiy should
comprise several species, and each species may have several sequences belonging to the
family. The algorithm tries to associate (match) interacting partners from the two families
within each species. It belongs to the more general class of
[Direct Coupling Analysis methods][dca-wiki].

The underlying main assumption is that the proper matching is the one maximizing the
co-evolution signal. Such maximization is performed over the Bayesian inference of a
Gaussian model, by inverting the correlation matrix.

The code is written in [Julia][julia], and the functions are called
from within Julia. However, a command-line interface is also provided for
those unfamiliar with the language (see the documentation).

The package is tested against Julia `0.6` and `0.7` on Linux, OS X, and Windows.

## Installation

The package is not registered; it can be installed with `Pkg.clone`:

```
julia> Pkg.clone("https://github.com/Mirmu/ParalogMatching.jl")
```

Dependencies will be installed automatically.

The package requires to install at least one linear programming solver supported by
[MathProgBase][mathprogbase].
By default, it uses [GLPK][glpk], which is free and open source, but you can choose any another:
see the list of available solvers at the [JuliaOpt page][solvers].
However, note that the solver efficiency is not particularly important for paralog matching,
whose computational time is dominated by matrix inversion operations, therefore it's likely that
you won't need a particularly fast solver.

## Documentation

- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

## Contributing and Questions

Contributions are very welcome, as are feature requests and suggestions. Please open an
[issue][issues-url] if you encounter any problems or would just like to ask a question.

[paper]: http://dx.doi.org/10.1073/pnas.1607570113
[dca-wiki]: https://en.wikipedia.org/wiki/Direct_coupling_analysis
[julia]: https://julialang.org

[mathprogbase]: http://mathprogbasejl.readthedocs.io/en/latest/
[glpk]: https://github.com/JuliaOpt/GLPK.jl
[solvers]: http://www.juliaopt.org/#packages

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://Mirmu.github.io/ParalogMatching.jl/latest

[travis-img]: https://travis-ci.org/Mirmu/ParalogMatching.jl.svg?branch=master
[travis-url]: https://travis-ci.org/Mirmu/ParalogMatching.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/x9jkws1l4xd8q4wy/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/Mirmu/paralogmatching-jl/branch/master

[codecov-img]: https://codecov.io/gh/Mirmu/ParalogMatching.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/Mirmu/ParalogMatching.jl

[issues-url]: https://github.com/Mirmu/ParalogMatching.jl/issues
