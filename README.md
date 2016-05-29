# Matching_Paralogs
The matching routine for interologs identification based on co-evolution

Package for matching paralogs species between two physically interacting proteins, based on their co-evolution signals.
The underlying main assumption is that the proper matching maximizing the co-evolution signal.
Such maximization is performed over the Bayesian inference of a Gaussian model, by inverting the correlation matrix.

More details can be found in:
CITATION

##Architecture

The package is organized as follows: in the source folder, routines can be found:
1.Readdata.jl contains the reader for converting FASTA files in Julia types called Alignments
2.manip_fasta.jl contains all the routines necessary to match the two alignements. Namely it normalizes the names of the species and orders the sequences in contiguous blocks
3.utils_fasta.jl contains some helpers for the taks of manipulating FASTA
4.types.jl contains the types used in the optimization/matching procedure
5.utils.jl contains the routines performing the optimization task for a given species
6.matching_fasta.jl contains the routines performing the whole matching of both Alignments
7.Plotting.jl contains various scoring/plotting functions to assess the quality of the solution
8. ParaMatch.jl is a wrapper Module for Julia Package.

##Installation

This Package requires the following dependencies: FastaIO, MacroUtils, Iterators, Compat, Distributions, Gurobi. They can be installed either from Package Rep. or from cloning GitHub reps.
For installing the Package, just run the following command lines:





##Use

Given two FASTA files, run the following:







