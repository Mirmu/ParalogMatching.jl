# This file is a part of ParalogMatching.jl. License is GPL3+: http://github.com/Mirmu/ParalogMatching.jl/LICENCE.md

##################### BLOCK FOR MANIPULATION OF THE ALIGNMENTS ##################

# Orders and removes families that are too large
function order_and_cut(Xtot::Alignment, cutoff::Integer)
    @extract Xtot : spec_name N q Z sequence uniprot_id spec_id header
    cand = tally_backref(spec_name)
    kept = Int[]
    # Orders so that families have contiguous sequences organized in blocks
    for (fam,n,inds) in cand
        n > cutoff && continue
        append!(kept, inds)
    end

    return Alignment(N, length(kept), q, Z[kept,:],
                     sequence[kept], header[kept], spec_name[kept], spec_id[kept], uniprot_id[kept])
end

# Harmonizing consists of creating the following structure for a FASTA:
# Filter the species in common between both FASTA and remove the others, that cannot be matched
# Giving those species a common labeling
function harmonize_fasta(X1::Alignment, X2::Alignment)
    @extract X1 : spec_name1=spec_name N1=N q1=q Z1=Z sequence1=sequence uniprot_id1=uniprot_id header1=header
    @extract X2 : spec_name2=spec_name N2=N q2=q Z2=Z sequence2=sequence uniprot_id2=uniprot_id header2=header
    us1 = unique(spec_name1)
    us2 = unique(spec_name2)
    kept = intersect(us1, us2)

    isempty(kept) && error("No species in common found in the alignments")

    kd = Dict(s=>i for (i,s) in enumerate(kept)) # associate an index to each kept species

    ind1, sid1 = compute_new_inds(spec_name1, kd)
    ind2, sid2 = compute_new_inds(spec_name2, kd)

    al1 = Alignment(N1, length(ind1), q1, Z1[ind1, :], sequence1[ind1], header1[ind1],
                    spec_name1[ind1], sid1, uniprot_id1[ind1])
    al2 = Alignment(N2, length(ind2), q2, Z2[ind2, :], sequence2[ind2], header2[ind2],
                    spec_name2[ind2], sid2, uniprot_id2[ind2])

    return HarmonizedAlignments(al1, al2)
end

# auxiliary function for harmonize_fasta
function compute_new_inds(spec_name::Vector{String}, kd::Dict{String,Int})
    ind = Int[] # the subset of indices to keep
    sid = Int[] # the new species ids
    for (i,s) in enumerate(spec_name)
        id = get(kd, s, 0) # get the species index in the `kept` vector
        id == 0 && continue # the species is not in `kept`
        push!(sid, id)
        push!(ind, i)
    end
    # here, ind is sorted but sid is not
    # we want to sort everything according to the new species ids
    p = sortperm(sid)

    return ind[p], sid[p]
end


#################### BLOCK FOR WRITING FASTA FROM ALIGNMENTS OBJECTS #######

# Writes a given Alignment under "name"
function write_fasta(X1::Alignment, name::AbstractString)
    @extract X1 : header sequence
    FastaWriter(name) do f
        for (h,s) in zip(header, sequence)
            writeentry(f, h, s)
        end
    end
end

"""
    write_fasta_match(X12::HarmonizedAlignments, match::Vector{Int}, outfile::AbstractString)

Writes a new FASTA file obtained from the two alignments contained in the `X12`
object and the matching contained in `match`.

The sequences headers in the output file have the format `>ID1::ID2/SPECIES`, where `ID1` represents
the ID read from the first alignment, `ID2` that for the second alignment, and `SPECIES` is the species
name.

The new sequences in the output file are the concatenation of the matched sequences.
"""
function write_fasta_match(X12::HarmonizedAlignments, match::Vector{Int}, outfile::AbstractString)
    @extract X12 : X1 X2
    FastaWriter(outfile) do f
        for (i,edge) in enumerate(match)
            edge == 0 && continue
            X1.spec_name[i] == X2.spec_name[edge] || error("do you have a well formed match ?")
            h = string(">", X1.uniprot_id[i], "::", X2.uniprot_id[edge], "/", X1.spec_name[i])
            write(f, h)
            write(f, X1.sequence[i])
            write(f, X2.sequence[edge])
        end
    end
end

##################### BLOCK FOR COMPARING TWO FASTA ##########################

# Computes the intersection of two FASTA files "name1" and "name2", and optionally writes it to a file.
# If printtofile is false, it simply returns the intersection value, with both lengths of the FASTA.
function intersection_fasta(name1::AbstractString, name2::AbstractString; printtofile::Bool = false)
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
            printtofile && writeentry(f, h, s)
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
            writeentry(f, h, s)
        end
    end
    try
        for (h,s) in zip(X2.header, X2.sequence)
            h ∈ X1.header && continue
            printtofile && writeentry(f, h, s)
            usize += 1
        end
    finally
        printtofile && close(f)
    end
    return usize, M1, M2
end

# Special constructor to create an Alignment from one HarmonizedAlignment and a matching
function Alignment(X12::ParalogMatching.HarmonizedAlignments, match::Vector{Int})
    @extract X12 : X1 X2

    N = X1.N + X2.N
    X1.q == X2.q || warn("X1.q = $(X1.q) X2.q = $(X2.q); should be the same!")
    q = X1.q

    header = String[]
    sequence = String[]
    spec_name = String[]
    spec_id = Int[]
    uniprot_id = String[]
    M = sum(match .> 0)
    Z = zeros(Int8, M, N)
    ctr = 0
    for (i,edge) in enumerate(match)
        edge == 0 && continue
        ctr += 1
        X1.spec_name[i] == X2.spec_name[edge] || error("do you have a well formed match?")
        X1.spec_id[i] == X2.spec_id[edge] || error("do you have a well formed match?")
        push!(header, string(X1.uniprot_id[i], "::", X2.uniprot_id[edge], "/", X1.spec_name[i]))
        push!(sequence, X1.sequence[i] * X2.sequence[edge])
        push!(spec_name, X1.spec_name[i])
        push!(spec_id, X1.spec_id[i])
        push!(uniprot_id, X1.uniprot_id[i] * "::" * X2.uniprot_id[edge])
        for j=1:X1.N Z[ctr,j] = X1.Z[i,j] end
        for j=1:X2.N Z[ctr,j+X1.N] = X2.Z[edge,j] end
    end
    return ParalogMatching.Alignment(N, M, q, Z, sequence, header, spec_name, spec_id, uniprot_id)
end
