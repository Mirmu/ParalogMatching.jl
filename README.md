# Matching_Paralogs
The matching routine for interologs identification based on co-evolution.

Package for matching paralogs species between two physically interacting proteins, based on their co-evolution signals.
The underlying main assumption is that the proper matching maximizing the co-evolution signal.
Such maximization is performed over the Bayesian inference of a Gaussian model, by inverting the correlation matrix.

More details can be found in:
CITATION

## Installation

Use `Pkg.clone`:

```
julia> Pkg.clone("https:////github.com/Mirmu/Matching_Paralogs")
```

Dependencies will be installed automatically.

However, you will also need to install at least one linear programming solver supported by
[MathProgBase](http://mathprogbasejl.readthedocs.io/en/latest/).
See the list of available solvers at the [JuliaOpt page](http://www.juliaopt.org/#packages).
Note that the solver efficiency is not particularly important for paralog matching, whose computational time
is dominated by matrix inversion operations, therefore you don't need a particularly fast solver. If unsure,
use `Pkg.add("Clp")` or `Pkg.add("GLPKMathProgBase")`, which are free and open-source solvers.

## Usage

Given two FASTA files, run the following:

```
julia> using ParaMatch # NOT REALLY, NEEDS FIXING!!!!!!!!!!!!!!!!!

julia> X1 = read_fasta_alignment("file1.fasta");

julia> X2 = read_fasta_alignment("file2.fasta");

julia> init = initialize_matching(X1, X2);

julia> results = run_matching(init, 5);
```

## Code organization

The package is organized as follows: in the source folder, routines can be found:

* `readdata.jl` contains the reader for converting FASTA files in Julia types called Alignments

* `manip_fasta.jl` contains all the routines necessary to match the two alignements. Namely it normalizes the names of the species and orders the sequences in contiguous blocks

* `utils_fasta.jl` contains some helpers for the taks of manipulating FASTA

* `types.jl` contains the types used in the optimization/matching procedure

* `utils.jl` contains the routines performing the optimization task for a given species

* `matching_fasta.jl` contains the routines performing the whole matching of both Alignments

* `plotting.jl` contains various scoring/plotting functions to assess the quality of the solution

* `ParaMatch.jl` is a wrapper Module for Julia Package.
