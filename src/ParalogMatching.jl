# This file is a part of ParalogMatching.jl. License is GPL3+: http://github.com/Mirmu/ParalogMatching.jl/LICENCE.md

module ParalogMatching

using Compat
using Compat.LinearAlgebra, Compat.Random, Compat.SparseArrays
using FastaIO, ExtractMacro
using Distributions
using MathProgBase
using GLPKMathProgInterface
const default_lpsolver = GLPKSolverLP()

export read_fasta_alignment, prepare_alignments, run_matching,
       write_fasta_match, paralog_matching

# Including the types used for the Matching
include("types.jl")

# Including the reader for the FASTA
include("readdata.jl")

# Including some utilities for the FASTA/Alignment manip
include("utils_fasta.jl")

# Including the routines for manipulation/Preparing for merging the Alignment
include("manip_fasta.jl")

# utils.jl contains the matching/optimization routines
include("utils.jl")

# Performing the iterative/progressive matching task
include("matching_fasta.jl")

# Including the routines for evaluating/plotting the solution
include("plotting.jl")

end
