######################BLOCK FOR MANIPULATION OF THE ALIGNMENTS###################

# Orders and removes families that are too large
function order_and_cut(Xtot::Alignment, cutoff::Int64)
    @extract Xtot : SpecName N q split Z Sequence UniprotId SpecId Header
    cand = tally(SpecName)
    ncand = [a[1] for a in filter(x->x[2]<=cutoff, cand)]
    kept = Int64[]

    # You may add more criteria for filtering here...

    # Orders so that families have contiguous sequences organized in blocks
    for fam in ncand
	append!(kept, find(x->x==fam, SpecName))
    end

    return Alignment(N, length(kept), q, split, Z[kept,:],
		     Sequence[kept], Header[kept], SpecName[kept], SpecId[kept], UniprotId[kept])
end

# Harmonizing consists of creating the following structure for a FASTA:
# Filter the species in common between both FASTA and remove the others, that cannot be matched
# Giving those species a common labeling
# This routine is a bit inefficient
function harmonize_fasta(X1, X2)
    @extract X1 : SpecName1=SpecName N1=N q1=q Z1=Z Sequence1=Sequence UniprotId1=UniprotId SpecId1=SpecId Header1=Header
    @extract X2 : SpecName2=SpecName N2=N q2=q Z2=Z Sequence2=Sequence UniprotId2=UniprotId SpecId2=SpecId Header2=Header
    _s1 = unique(SpecName1)
    _s2 = unique(SpecName2)
    kept = intersect(_s1, _s2)

    ind1 = Vector{Vector{Int}}(length(kept))
    ctr = 0
    for fam in kept
	ctr += 1
	ind1[ctr] = find(x->x==fam, SpecName1)
    end

    ind2 = Vector{Vector{Int}}(length(kept))
    ctr = 0
    for fam in kept
	ctr += 1
	ind2[ctr] = find(x->x==fam, SpecName2)
    end

    lind1 = Int[]
    for i in eachindex(ind1)
	val = findfirst(kept, SpecName1[ind1[i][1]])
	for j in ind1[i]
	    SpecId1[j] = val
	    push!(lind1, j)
	end
    end

    lind2 = Int[]
    for i in eachindex(ind2)
	val = findfirst(kept, SpecName2[ind2[i][1]])
	for j in ind2[i]
	    SpecId2[j] = val
	    push!(lind2, j)
	end
    end

    al1 = Alignment(N1, length(lind1), q1, 0, Z1[lind1, :], Sequence1[lind1], Header1[lind1], SpecName1[lind1], SpecId1[lind1], UniprotId1[lind1])
    al2 = Alignment(N2, length(lind2), q2, 0, Z2[lind2, :], Sequence2[lind2], Header2[lind2], SpecName2[lind2], SpecId2[lind2], UniprotId2[lind2])

    return al1, al2
end

#####################BLOCK FOR WRITING FASTA FROM ALIGNMENTS OBJECTS########

# Writes a given Alignments under "name"
function write_fasta(X1, name)
    g = open(name, "w")
    for i in 1:(X1.M)
	head = join([">", X1.Header[i]])
	println(g, head)
	println(g, X1.Sequence[i])
    end
    close(g)
end

# Rewrites the output as a FASTA for two given Alignments and a given matching, under the "name"
function rewrite_fasta_match(X1, X2, match, name)
    g = open(name, "w")
    for (i,edge) in enumerate(match)
	edge == 0 && continue
	X1.SpecName[i] == X2.SpecName[edge] || error("do you have a well formed match ?")

	head = join([">", X1.UniprotId[i], "with", X2.UniprotId[edge], "/", X1.SpecName[i]])
	println(g, head)
	print(g, X1.Sequence[i])
	println(g, X2.Sequence[edge])
    end
    close(g)
end

######################BLOCK FOR COMPARING TWO FASTA###########################

# Computes the Fasta that is the intersection of two FASTA of names "names1" and "names2"
# If no output, it simply returns the overlap value, with both lengths of the FASTA
function overlap_fasta(name1, name2; opt::AbstractString = "no-out")
    X1 = readdata(name1)
    X2 = readdata(name2)
    pool = unique(X1.SpecName)
    len1 = X1.M
    len2 = X2.M
    over = 0
    if opt == "out"
	str = string("OLP", name1, "_", name2, ".fasta")
	g = open(str, "w")
    end
    for (i,j) in enumerate(X1.Header)
	j ∈ X2.Header || continue
	if opt == "out"
	    println(g, string(">", j))
	    println(g, X1.Sequence[i])
	end
	over += 1
    end
    if opt == "out"
	close(g)
    end
    return over, len1, len2
end

# Computes the FASTA that is the Union of the two FASTA of names "names1" and "names2"
# Same return than Overlap
function union_fasta(name1, name2; opt::AbstractString = "no-out")
    X1 = readdata(name1)
    X2 = readdata(name2)
    pool = unique(X1.SpecName)
    len1 = X1.M
    len2 = X2.M
    over = 0
    if opt == "out"
	str = string("UNION", name1, "_", name2, ".fasta")
	g = open(str, "w")
	for i in 1:len1
	    println(g, string(">", X1.Header[i]))
	    println(g, X1.Sequence[i])
	end
    end
    for (i,j) in enumerate(X2.Header)
	j ∈ X1.Header && continue
	if opt == "out"
	    println(g, string(">", j))
	    println(g, X2.Sequence[i])
	end
	over += 1
    end
    if opt == "out"
	close(g)
    end
    return over, len1, len2
end