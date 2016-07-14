######################BLOCK FOR MANIPULATION OF THE ALIGNMENTS###################

# Orders and removes families that are too large
function order_and_cut(Xtot::Alignment, cutoff::Int64)
    @extract Xtot : spec_name N q split Z sequence uniprot_id spec_id header
    cand = tally(spec_name)
    ncand = [a[1] for a in filter(x->x[2]≤cutoff, cand)]
    kept = Int64[]

    # You may add more criteria for filtering here...

    # Orders so that families have contiguous sequences organized in blocks
    for fam in ncand
	append!(kept, find(x->x==fam, spec_name))
    end

    return Alignment(N, length(kept), q, split, Z[kept,:],
		     sequence[kept], header[kept], spec_name[kept], spec_id[kept], uniprot_id[kept])
end

# Harmonizing consists of creating the following structure for a FASTA:
# Filter the species in common between both FASTA and remove the others, that cannot be matched
# Giving those species a common labeling
# This routine is a bit inefficient
function harmonize_fasta(X1, X2)
    @extract X1 : spec_name1=spec_name N1=N q1=q Z1=Z sequence1=sequence uniprot_id1=uniprot_id spec_id1=spec_id header1=header
    @extract X2 : spec_name2=spec_name N2=N q2=q Z2=Z sequence2=sequence uniprot_id2=uniprot_id spec_id2=spec_id header2=header
    _s1 = unique(spec_name1)
    _s2 = unique(spec_name2)
    kept = intersect(_s1, _s2)

    ind1 = Vector{Vector{Int}}(length(kept))
    ctr = 0
    for fam in kept
	ctr += 1
	ind1[ctr] = find(x->x==fam, spec_name1)
    end

    ind2 = Vector{Vector{Int}}(length(kept))
    ctr = 0
    for fam in kept
	ctr += 1
	ind2[ctr] = find(x->x==fam, spec_name2)
    end

    lind1 = Int[]
    for i in eachindex(ind1)
	val = findfirst(kept, spec_name1[ind1[i][1]])
	for j in ind1[i]
	    spec_id1[j] = val
	    push!(lind1, j)
	end
    end

    lind2 = Int[]
    for i in eachindex(ind2)
	val = findfirst(kept, spec_name2[ind2[i][1]])
	for j in ind2[i]
	    spec_id2[j] = val
	    push!(lind2, j)
	end
    end

    al1 = Alignment(N1, length(lind1), q1, 0, Z1[lind1, :], sequence1[lind1], header1[lind1],
		    spec_name1[lind1], spec_id1[lind1], uniprot_id1[lind1])
    al2 = Alignment(N2, length(lind2), q2, 0, Z2[lind2, :], sequence2[lind2], header2[lind2],
		    spec_name2[lind2], spec_id2[lind2], uniprot_id2[lind2])

    return al1, al2
end

#################### BLOCK FOR WRITING FASTA FROM ALIGNMENTS OBJECTS #######

# Writes a given Alignment under "name"
function write_fasta(X1::Alignment, name::AbstractString)
    @extract X1 : header sequence
    FastaWriter(name) do f
	for (h,s) in zip(header, sequence)
	    writeentry(f, (h, s))
	end
    end
end

# Rewrites the output as a FASTA for two given Alignments and a given matching, under the "name"
function rewrite_fasta_match(X1::Alignment, X2::Alignment, match::Vector{Int}, name::AbstractString)
    FastaWriter(name) do f
	for (i,edge) in enumerate(match)
	    edge == 0 && continue
	    X1.spec_name[i] == X2.spec_name[edge] || error("do you have a well formed match ?")
	    h = string(">", X1.uniprot_id[i], "with", X2.uniprot_id[edge], "/", X1.spec_name[i])
	    write(f, h)
	    write(f, X1.sequence[i])
	    write(f, X2.sequence[edge])
	end
    end
end

##################### BLOCK FOR COMPARING TWO FASTA ##########################

# Computes the intersection of two FASTA files "name1" and "name2", and optionally writes it to a file.
# If printtofile is false, it simply returns the intersection value, with both lengths of the FASTA.
function intersaction_fasta(name1::AbstractString, name2::AbstractString; printtofile::Bool = false)
    X1 = read_fasta_alignment(name1)
    X2 = read_fasta_alignment(name2)
    M1 = X1.M
    M2 = X2.M
    isize = 0
    if printtofile
	f = FastaWriter(string("INTERSECTION_", name1, "-", name2, ".fasta.gz"))
    end
    try
	for (h,s) in zip(X1.header, X1.sequence)
	    h ∈ X2.header || continue
	    printtofile && writeentry(f, (h, s))
	    isize += 1
	end
    finally
	printtofile && close(f)
    end
    return isize, M1, M2
end

# Computes the union of two FASTA files "name1" and "name2", and optionally writes it to a file.
# If printtofile is false, it simply returns the union value, with both lengths of the FASTA.
function union_fasta(name1::AbstractString, name2::AbstractString; printtofile::Bool = false)
    X1 = read_fasta_alignment(name1)
    X2 = read_fasta_alignment(name2)
    M1 = X1.M
    M2 = X2.M
    usize = M1
    if printtofile
	f = FastaWriter(string("UNION_", name1, "-", name2, ".fasta.gz"))
	for (h,s) in zip(X1.header, X1.sequence)
	    writeentry(f, (h, s))
	end
    end
    try
	for (h,s) in zip(X2.header, X2.sequence)
	    h ∈ X1.header && continue
	    printtofile && writeentry(f, (h, s))
	    usize += 1
	end
    finally
	printtofile && close(f)
    end
    return usize, M1, M2
end
