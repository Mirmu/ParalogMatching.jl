module ParaMatch

using FastaIO, MacroUtils, Iterators, Compat
using Distributions
using Gurobi

export readdata, initialize, run_matching

#Including the reader for the FASTA
include("readdata.jl")

#Including some utilities for the FASTA/Alignment manip
include("utils_fasta.jl")

#Including the routines for manipulation/Preparing for merging the Alignment
include("manip_fasta.jl")

#Including the types used for the Matching
include("types.jl")

#Utils.jl contains the matching/optimization routines
include("utils.jl")

#Performing the iterative/progressive matching task
include("matching_fasta.jl")

#Including the routines for evaluating/plotting the solution
include("plotting.jl")

end
