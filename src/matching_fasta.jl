# This file is a part of ParalogMatching.jl. License is GPL3+: http://github.com/Mirmu/ParalogMatching.jl/LICENCE.md

########################## PRE PROCESSING THE MATCHING ##############################

"""
    prepare_alignments(X1::Alignment, X2::Alignment; cutoff::Integer = 500)

Prepares two alignments (as returned by [`read_fasta_alignment`](@ref)) and returns
a `HarmonizedAlignments` object. This object contains two filtered version of the
original alignments, in which only the sequences belonging to species which exist
in both alignments are kept, and which is ready to be passed to [`run_matching`](@ref).

The `cutoff` keyword argument can be used to discard all species for which there
are more than a certain number of sequences in either alignment. Use `0` to disable this
filter entirely.
"""
function prepare_alignments(Xi1::Alignment, Xi2::Alignment; cutoff::Integer = 500)
    cutoff == 0 && (cutoff = typemax(Int))
    println("initializing the matching",
            cutoff < typemax(Int) ? ", removing families larger than $cutoff" : "",
            "...")
    return harmonize_fasta(order_and_cut(Xi1, cutoff), order_and_cut(Xi2, cutoff))
end

# Returns the initial matching between single species
# Works only for Harmonized FASTA
function start_matching(X12::HarmonizedAlignments)
    @extract X12 : X1 X2
    @extract X1  : spec_id1=spec_id
    @extract X2  : spec_id2=spec_id
    match = fill!(similar(spec_id1), 0)

    # Finds the indices of the species with one single sequence
    ind1 = index_of_unique(spec_id1)
    ind2 = index_of_unique(spec_id2)

    candi = intersect(spec_id1[ind1], spec_id2[ind2])

    for el in candi
        a1 = findfirst(isequal(el), spec_id1)::Int
        a2 = findfirst(isequal(el), spec_id2)::Int
        match[a1] = a2
    end
    return match
end

# Initializes the problem from the output of prepare_alignments,
# allocating frequency matrix and correlation matrices, inverting
# them and returning them
function initialize_matching(X12::HarmonizedAlignments, pseudo_count::Float64)
    @extract X12 : X1 X2

    # Match by uniqueness
    match = start_matching(X12)

    # Computing the prior correlation matrix and interaction matrix
    freq = FreqC(X1, X2)
    corr = FastC(freq)

    # First compute corr from single matched families
    single = X1.spec_id[findall(match .≠ 0)]

    # Computes the freq matrix for the given matched species "single"
    unitFC!(X1, X2, match, single, freq)

    # Finally compute the inverse of the corr matrix
    if isempty(findall(match .≠ 0))
	println("WARNING ! 0 sequence matched by uniqueness. No covariation strategy possible.")
	invC = zeros(size(corr.Cij))
    else
	invC = inverse_with_pseudo!(corr, freq, pseudo_count)
    end
    println("Setup computed")
    return X1, X2, match, freq, corr, invC
end

# par_corr gathers for each species in specl, the matching obtained by the "strategy" strategy
# And returns an array of those matchings
function par_corr(X1::Alignment, X2::Alignment, freq::FreqC, invC::Matrix{Float64},
                  specl::Vector{Int}, strategy::AbstractString,
                  lpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver)
    return [give_correction(X1, X2, freq, invC, i, strategy, lpsolver) for i in specl]
end

# spec_entropy computes the number of potential matchings for each species
# And returns the ordered list, from the easiest fams to the hardest
# WARNING: Works for harmonized Fasta only
function spec_entropy(X1::Alignment, X2::Alignment)
    @extract X1 : spec_id1=spec_id
    @extract X2 : spec_id2=spec_id
    bib1 = tally(spec_id1)
    bib2 = tally(spec_id2)

    entropy = Tuple{Int,Float64}[]

    for i in 1:length(bib1)
        bib1[i][1] == bib2[i][1] || error("non harmonized fasta")

        mini = min(bib1[i][2], bib2[i][2])
        maxi = max(bib1[i][2], bib2[i][2])

        ent = sum([log(i) for i in (maxi-mini+1):maxi])
        push!(entropy, (bib1[i][1], ent))
    end
    return [a[1] for a in filter(x->x[2]!=0, sort(entropy, by=x->x[2]))]
end

############################## MAIN FUNCTION ####################################

# First a helper...
# given a list of specs with their respective matching within species, it updates the global match vector
function apply_matching!(X1, X2, match, lspec, lmatch)
    @extract X1 : spec_id1=spec_id
    @extract X2 : spec_id2=spec_id
    length(lspec) == length(lmatch) || error("data non compatible")

    for (i,el) in enumerate(lspec)
        ind1 = findall(spec_id1 .== el)
        ind2 = findall(spec_id2 .== el)
        match[ind1[lmatch[i][1]]] = ind2[lmatch[i][2]]
    end
    return nothing
end

"""
    run_matching(X12::HarmonizedAlignments;
                 batch = 1,
                 strategy = "covariation",
                 lpsolver = GLPKSolverLP())

Returns the matching of the two alignments contained in `X12`, which need to be obtained by [`prepare_alignments`](@ref).
The returned matching is a Vector in which the entry `i` determines which sequence of the second alignment
is the partner to sequence `i` in the first alignment.

The keywords are:

* `batch`: how many species should be matched before updating the underlying model; smaller values give
           better results but increase the computational time

* `strategy`: the strategy to use when computing the matching. See below.

* `pseudo_count`: gives the amount of regularization used for the inversion of the correlation matrix.
                  A defaut of 0.5 gives sensible results, its value cannot be above 1.0.

* `lpsolver`: linear programming solver used when performing the matching with the `"covariation"` strategy.
              The default uses a solver provided bu the `GLPK` library.
              You can override this by passing e.g. `lpsolver = GurobiSolver(OutputFlag=false)` or similar
              (see the documentation for `MathProgBase`).

Available strategies are:

+ `"covariation"`: uses the Gaussian model and performs a matching on the scores (default).

+ `"greedy"`: same as `"covariation"` but performs a greedy matching.

+ `"random"`: produces a random matching, only useful to produce null models.

+ `"genetic"`: tries to use the Uniprot ID information to determine which sequences belong to the
               same operon (only used for testing, not a valid general strategy)

"""
function run_matching(X12::HarmonizedAlignments;
                      batch::Integer = 1,
                      pseudo_count::Float64 = 0.5,
                      strategy::AbstractString = "covariation",
                      lpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver = default_lpsolver)

    X1, X2, match, freq, corr, invC = initialize_matching(X12, pseudo_count)

    valid_strats = ["covariation", "genetic", "random", "greedy"]
    strategy ∈ valid_strats ||
        throw(ArgumentError("unknown strategy: $strategy. Must be one of: $(join(valid_strats, ", ", " or "))"))
    pseudo_count>0.0 && (pseudo_count < 1.0) || throw(ArgumentError("invalid value of the pseudo_count, must be a real between 0.0 and 1.0"))

    # Computes the entropy of the families and batch them from easiest to hardest
    spec = spec_entropy(X1, X2)
    len = length(spec)
    batchl = [spec[i*batch+1:min((i+1)*batch, len)] for i in 0:div(len, batch)]

    # For each batch...
    for el in batchl
        isempty(el) && continue

        # Performs the matching for each species of the batch
        res = par_corr(X1, X2, freq, invC, el, strategy, lpsolver)
        println("batch of species")

        for (i,spids) in enumerate(el)
            println(X1.spec_name[spids]," ")
            println(res[i])
        end
        println()

        # Applies the matching to the global matching vector
        apply_matching!(X1, X2, match, el, res)

        if strategy == "covariation" || strategy == "greedy"
            println("Recomputing the model")
            # Updates the freq and corr matrices, and its inverse
            unitFC!(X1, X2, match, el, freq)
            invC = inverse_with_pseudo!(corr, freq, pseudo_count)
        end

    end
    clear_inverse_mem()
    #return X1, X2, match, freq, corr, invC, savematch
    return match
end

"""
    paralog_matching(infile1::AbstractString, infile2::AbstractString, outfile::AbstractString;
                     cutoff = 500,
                     batch = 1,
                     strategy = "covariation",
                     lpsolver = GLPKSolverLP())

This function performs the paralog matching from two given FASTA files containing the alignments for
two different protein families. When the matching is done, it writes the result in a new FASTA file,
in which each sequence is the concatenation of two mathing sequences in the original file.

The input format is standard FASTA, but the headers are parsed as explained in [`prepare_alignments`](@ref).

The output file format is documented in [`write_fasta_match`](@ref).

The keyword arguments are documented in [`prepare_alignments`](@ref) and [`run_matching`](@ref).

Besides writing the `outfile`, this function also returns two values: the harmonized alignment (produced
by [`prepare_alignments`](@ref)) and the resulting match (as returned by [`run_matching`](@ref)).
"""
function paralog_matching(infile1::AbstractString,
                          infile2::AbstractString,
                          outfile::AbstractString;
                          cutoff::Integer = 500,
                          batch::Integer = 1,
                          strategy::AbstractString = "covariation",
                          pseudo_count::Float64= 0.5,
                          lpsolver::MathProgBase.SolverInterface.AbstractMathProgSolver = default_lpsolver)

    X1 = read_fasta_alignment(infile1)
    X2 = read_fasta_alignment(infile2)

    X12 = prepare_alignments(X1, X2, cutoff=cutoff)
    match = run_matching(X12, batch=batch, strategy=strategy, pseudo_count=pseudo_count, lpsolver=lpsolver)

    write_fasta_match(X12, match, outfile)

    println("done")

    return X12, match
end
