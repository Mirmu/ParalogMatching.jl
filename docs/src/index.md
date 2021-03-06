# ParalogMatching.jl documentation

## Introduction

This package implements the paralog matching technique presented in the paper
*"Simultaneous identification of specifically interacting paralogs and
interprotein contacts by Direct Coupling Analysis"*
by Thomas Gueudré, Carlo Baldassi, Marco Zamparo, Martin Weigt and Andrea Pagnani,
Proc. Natl. Acad. Sci. U.S.A. 113, 12186–12191 (2016), [doi:10.1073/pnas.1607570113](http://dx.doi.org/10.1073/pnas.1607570113)

The main idea of the method is to perform a statistical analysis of two given
multiple sequence alignments, each containing one protein family. Each familiy should
comprise several species, and each species may have several sequences belonging to the
family. The algorithm tries to associate (match) interacting partners from the two families
within each species.

The underlying main assumption is that the proper matching is the one maximizing the
co-evolution signal. Such maximization is performed over the Bayesian inference of a
Gaussian model, by inverting the correlation matrix.

The code is written in [Julia](http://julialang.org), and the functions are called
from within Julia. However, a [command-line interface](@ref CLI) is also provided for
those unfamiliar with the language.

The current code was tested on Julia versions 0.6 and 0.7.

### Installation

To install the module, use this command from within Julia:

```
julia> Pkg.clone("https://github.com/Mirmu/ParalogMatching.jl")
```

Dependencies will be installed automatically.

The package requires to install at least one linear programming solver supported by
[MathProgBase](http://mathprogbasejl.readthedocs.io/en/latest/).
By default, it uses [GLPK](https://github.com/JuliaOpt/GLPK.jl), which is free and open source,
but you can choose any another: see the list of available solvers at the
[JuliaOpt page](http://www.juliaopt.org/#packages).
However, note that the solver efficiency is not particularly important for paralog matching,
whose computational time is dominated by matrix inversion operations, therefore it's likely that
you won't need a particularly fast solver.

### Usage

The module is loaded as any other Julia module:

```
julia> using ParalogMatching
```

The module provides a [high-level interface](@ref high_level), with one function
which performs all operations and writes a file in output, and the
[low-level interface](@ref low_level) functions which perform the single
steps of the algorithm.

A typical run of the algorithm could look like this:

```
julia> TrpAB, match = paralog_matching("TrpA.fasta", "TrpB.fasta", "Match_TrpAB.fasta.gz", batch=5);
```

## [High-level interface](@id high_level)

```@docs
paralog_matching
```

## [Low-level interface](@id low_level)

The low-level interface functions are called by [`paralog_matching`](@ref) in the order
in which they are documented here, but they can be called individually in order to
inspect/manipulate the intermediate results if desired.

```@docs
read_fasta_alignment
```

```@docs
prepare_alignments
```

```@docs
run_matching
```

```@docs
write_fasta_match
```

## [Command-line interface](@id CLI)

In the `utils` directory there is a script which can be called from the command line,
called `paralog_matching.jl`. Try calling it with:

```
$ julia paralog_matching.jl --help
```

to see the details. The arguments and options are basically the same as those of the
[`paralog_matching`](@ref) function.
