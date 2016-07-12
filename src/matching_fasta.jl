###########################PRE PROCESSING THE MATCHING###############################

# Returns the initial matching between single species
# Works only for Harmonized FASTA
function start_matching(X1, X2)
    @extract X1 : SpecId1=SpecId
    @extract X2 : SpecId2=SpecId
    match = zeros(SpecId1)
    #finds the indices of the species with one single sequence

    ind1 = index_of_unique(SpecId1)
    ind2 = index_of_unique(SpecId2)

    candi = intersect(SpecId1[ind1], SpecId2[ind2])

    for el in candi
	a1 = findfirst(SpecId1, el)
	a2 = findfirst(SpecId2, el)
	match[a1] = a2
    end
    return match
end

# Initializes the problem by preparing the fasta, allocating frequency matrix
# and correlation matrices, inverting them and returning them
function initialize(Xi1, Xi2, cut)
    println("initializing the matching, removing families larger than $(cut)...")
    X1, X2 = harmonize_fasta(order_and_cut(Xi1, cut), order_and_cut(Xi2, cut))
    match = start_matching(X1, X2)

    # Computing the prior correlation matrix and interaction matrix
    freq = nullF(X1, X2)
    corr = nullC(freq)

    # First compute corr from single matched families
    single = X1.SpecId[find(match)]

    # Computes the freq matrix for the given matched species "single"
    unitFC!(X1, X2, match, single, freq)

    # Finally compute the inverse of the corr matrix
    invC = inverse_with_pseudo!(corr, freq, 0.8)

    println("Setup computed")
    return X1, X2, match, freq, corr, invC
end

# par_corr gathers for each species in specl, the matching obtained by the "strat" strategy
# And returns an array of those matchings
function par_corr(X1, X2, freq::FreqC, invertC, specl::Vector{Int}, strat)
    res = [give_correction(X1, X2, freq, invertC, i, strat) for i in specl]
    return res
end

# spec_entropy computes the number of potential matchings for each species
# And returns the ordered list, from the easiest fams to the hardest
# WARNING: Works for harmonized Fasta only
function spec_entropy(X1, X2)
    @extract X1 : SpecId1=SpecId
    @extract X2 : SpecId2=SpecId
    bib1 = tally(SpecId1)
    bib2 = tally(SpecId2)

    entropy = Tuple{Int64,Float64}[]

    for i in 1:length(bib1)
	bib1[i][1] == bib2[i][1] || error("non harmonized fasta")

	mini = min(bib1[i][2], bib2[i][2])
	maxi = max(bib1[i][2], bib2[i][2])

	ent = sum([log(i) for i in (maxi-mini+1):maxi])
	push!(entropy, (bib1[i][1], ent))
    end
    return [a[1] for a in filter(x->x[2]!=0, sort(entropy, by=x->x[2]))]
end

###############################MAIN FUNCTION#####################################

# First a helper...
# given a list of specs with their respective matching within species, it updates the global match vector
function apply_matching!(X1, X2, match, lspec, lmatch)
    @extract X1 : SpecId1=SpecId
    @extract X2 : SpecId2=SpecId
    length(lspec) == length(lmatch) || error("data non compatible")

    for (i,el) in enumerate(lspec)
	ind1 = find(SpecId1 .== el)
	ind2 = find(SpecId2 .== el)
	match[ind1[lmatch[i][1]]] = ind2[lmatch[i][2]]
    end
    return nothing
end

# The main function that runs the matching
# Works for harmonized Fasta
# options "strat" for the matching are :
#   "covariation": computes the matching from co evolution signal
#   "genetic":     computes the matching from genetic proximity (if FASTA contains genetic position info)
#   "random":      computes a random matching for null hypothesis
#   "greedy":      computes a matching from a greedy strategy with the co evolution signal
# the argument "a" should be the output of the initialize function that can be found in Fasta_Manip.jl

function run_matching(a, batch, strat::AbstractString)
    # Takes the ouput of initialize
    X1, X2, match, freq, corr, invC = a

    # Computes the entropy of the families and batch them from easiest to hardest
    spec = spec_entropy(X1, X2)
    len = length(spec)
    batchl = [spec[i*batch+1:min((i+1)*batch, len)] for i in 0:div(len, batch)]

    savematch = Tuple{Vector{Int64},Vector{Int64}}[]

    # For each batch...
    for el in batchl
	isempty(el) && continue

	# Performs the matching for each species of the batch
	res = par_corr(X1, X2, freq, invC, el, strat)
	println("batch of species")
	println(el)
	println(res)

	# Applies the matching to the global matching vector
	apply_matching!(X1, X2, match, el, res)

	if strat == "covariation" || strat == "greedy"
	    println("Recomputing the model")
	    # Updates the freq and corr matrices, and its inverse
	    unitFC!(X1, X2, match, el, freq)
	    invC = inverse_with_pseudo!(corr, freq, 0.8)
	end

	# Takes a snapshot of the matching being built
	push!(savematch, (deepcopy(match), deepcopy(freq.specs)))
    end
    return X1, X2, match, freq, corr, invC, savematch
end
