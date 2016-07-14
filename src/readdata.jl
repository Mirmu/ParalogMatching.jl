abstract DataSet

immutable Alignment <: DataSet
    N::Int
    M::Int
    q::Int
    split::Int
    Z::Matrix{Int8}
    sequence::Vector{ASCIIString}
    header::Vector{ASCIIString}
    spec_name::Vector{ASCIIString}
    spec_id::Vector{Int}
    uniprot_id::Vector{ASCIIString}
end

function read_fasta_alignment(filename::AbstractString, max_gap_fraction::Float64 = 1.0)
    f = FastaReader(filename)

    # pass 1

    seqs = Int[]
    inds = Int[]
    fseqlen = 0

    for (name, seq) in f
        ngaps = 0
        if f.num_parsed == 1
            ls = length(seq)
            resize!(inds, ls)
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    fseqlen += 1
                    inds[fseqlen] = i
                    c == '-' && (ngaps += 1)
                end
            end
        else
            ls = length(seq)
            ls == length(inds) || error("inputs are not aligned")
            tstfseqlen = 0
            for i = 1:ls
                c = seq[i]
                if c != '.' && c == uppercase(c)
                    tstfseqlen += 1
                    inds[tstfseqlen] == i || error("inconsistent inputs")
                    c == '-' && (ngaps += 1)
                end
            end
            tstfseqlen == fseqlen || error("inconsistent inputs")
        end
        ngaps / fseqlen <= max_gap_fraction && push!(seqs, f.num_parsed)
    end

    length(seqs) > 0 || error("Out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

    # pass 2

    Z = Array(Int8, fseqlen, length(seqs))
    header =  Array(ASCIIString, length(seqs));
    sequence =  Array(ASCIIString, length(seqs));
    seqid = 1
    for (name, seq) in f
        #header[seqid] = specname.captures[1];
        header[seqid] = name
        sequence[seqid] = seq
        seqs[end] < f.num_parsed && break
        seqs[seqid] == f.num_parsed || continue
        for i = 1:fseqlen
            c = seq[inds[i]]
            Z[i, seqid] = letter2num(c)
        end
        seqid += 1
    end
    @assert seqid == length(seqs) + 1

    close(f)

    spec_id, spec_name, uniprot_id = compute_spec(header)

    return Alignment(size(Z, 1), size(Z, 2), Int(maximum(Z)), 12,  Z', sequence, header, spec_name, spec_id, uniprot_id)
end

function specname(s::ASCIIString)
    if ismatch(r"\[(.*?)\]", s)
        spec_name = convert(ASCIIString, match(r"\[(.*?)\]", s).captures[1])
        uniprot_id="000000"
        return (uniprot_id, spec_name)
    elseif ismatch(r"^([A-Z,0-9].*?)\_([A-Z,0-9].*?)/", s)
        uniprot_id, spec_name = match(r"^([A-Z,0-9].*?)\_([A-Z,0-9].*?)/", s).captures
        return (convert(ASCIIString, uniprot_id), convert(ASCIIString, spec_name))
    elseif ismatch(r"^(.*?)with(.*?)/(.*?)$", s)
        spec_name = convert(ASCIIString, match(r"^(.*?)with(.*?)/(.*?)$", s).captures[3])
        uniprot_id= "000000"
        return (uniprot_id, spec_name)
    end

    n = length(s)
    i1 = -1
    i2 = -1
    flag1 = 0
    flag2 = 0
    for i = 1:length(s)
        if s[i] == '/'
            i1 = i
            flag1 == 1 && error("badly formed string")
            flag1 = 1
        elseif s[i] == '_'
            i2 = i
            flag2 == 1 && error("badly formed string")
            flag2 = 1
        end
    end
    i1 < 0 || i2 < 0 && error("badly formed string")

    return (s[1:i1-1], s[i2+1:end])
end


function compute_spec(header::Vector{ASCIIString})
    M = length(header)

    spec_name = Array{ASCIIString}(M)
    uniprot_id  = Array{ASCIIString}(M)

    for i = 1:M
        protname, spec = specname(header[i])
        spec_name[i] = spec
        uniprot_id[i] = protname
    end

    specunique = unique(spec_name)
    nspec = length(specunique)
    spec_id = zeros(Int, M)
    idx = zeros(Int, M)

    @inbounds for j = 1:M
        nname = length(spec_name[j])
        for i = 1:nspec
            spec_name[j] == specunique[i] || continue
            spec_id[j] = i
            break
        end
    end

    return spec_id, spec_name, uniprot_id
end

@compat let alphabet = [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
                       # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    global letter2num
    function letter2num(c::Union{Char,UInt8})
        i = UInt(c) - 0x40
        1 <= i <= 25 && return alphabet[i]
        return 21
    end
end
nothing
