#!/usr/bin/env julia

# This file is a part of ParalogMatching.jl. License is GPL3+: http://github.com/Mirmu/ParalogMatching.jl/LICENCE.md

include(joinpath(dirname(@__FILE__), "../src/ParalogMatching.jl"))
using .ParalogMatching
using ArgParse

function main()
    s = ArgParseSettings("Match paralog species based on co-evolution signals.")

    valid_strats = ["covariation", "genetic", "random", "greedy"]

    @add_arg_table s begin
        "--cutoff", "-c"
            arg_type = Int
            default = 500
            range_tester = (x->x≥0)
            help = """maximum number of sequences per species
                      (discards species with more sequences than the cutoff),
                      use 0 to disable"""
        "--batch", "-b"
            arg_type = Int
            default = 1
            range_tester = (x->x≥1)
            help = """batch size, i.e. number of species to match before updating
                      the underlying alignment model (lower is more accurate but slower)"""
        "--strategy", "-s"
            default = "covariation"
            range_tester = (s->s ∈ valid_strats)
            help = "matching strategy. Allowed strategies are: $(join(valid_strats, ", ", " and "))."
        "infile1"
            help = "first input alignment"
            required = true
        "infile2"
            help = "second input alignment"
            required = true
        "outfile"
            help = "output matched alignment"
            required = true
    end

    parsed_args = parse_args(ARGS, s)

    cutoff = parsed_args["cutoff"]
    batch = parsed_args["batch"]
    infile1 = parsed_args["infile1"]
    infile2 = parsed_args["infile2"]
    outfile = parsed_args["outfile"]

    paralog_matching(infile1, infile2, outfile, cutoff=cutoff, batch=batch)
end

main()
