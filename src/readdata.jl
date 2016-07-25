function read_fasta_alignment(filename::AbstractString, max_gap_fraction::Float64 = 1.0; header_regex::Union{Void,Regex} = nothing)
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

    spec_id, spec_name, uniprot_id = compute_spec(header, header_regex)

    return Alignment(size(Z, 1), size(Z, 2), Int(maximum(Z)), Z', sequence, header, spec_name, spec_id, uniprot_id)
end

function specname(s::ASCIIString, header_regex::Union{Void,Regex} = nothing)
    if header_regex ≢ nothing
        # user-defined format
        if ismatch(header_regex)
            captures = match(header_regex, s).captures
            length(captures) == 2 || error("Invalid header regex: should always return 2 captured groups if it matches; has returned: $(length(captures))")
            uniprot_id, spec_name = captures
        else
            error("unrecognized spec string: $s")
        end
    else
        regex_uniprot = r"^(?:(?:(?:[^|/_]+?\|){2})|(?:[^|/_]+?/))([^|/_]+?)_([^|/_\s]+)"

        regex_oldskrr = r"\[(.*?)\]"
        regex_joined = r"^(.*?)::(.*?)/(.*)$"

        # standard format
        if ismatch(regex_uniprot, s)
            uniprot_id, spec_name = match(regex_uniprot, s).captures

        # custom internal formats
        elseif ismatch(regex_oldskrr, s)
            spec_name = match(regex_oldskrr, s).captures[1]
            uniprot_id = "000000"
        elseif ismatch(regex_joined, s)
            spec_name = match(regex_joined, s).captures[3]
            uniprot_id = "000000"
        else
            error("unrecognized spec string: $s")
        end
    end
    return convert(ASCIIString, uniprot_id), convert(ASCIIString, spec_name)
end

function compute_spec(header::Vector{ASCIIString}, header_regex::Union{Void,Regex} = nothing)
    M = length(header)

    spec_name = Array{ASCIIString}(M)
    uniprot_id  = Array{ASCIIString}(M)

    for i = 1:M
        uniprot_id[i], spec_name[i] = specname(header[i], header_regex)
    end

    specunique = unique(spec_name)
    sdict = [s=>i for (i,s) in enumerate(specunique)]
    spec_id = Int[sdict[sn] for sn in spec_name]

    return spec_id, spec_name, uniprot_id
end

let alphabet = [ 1,21, 2, 3, 4, 5, 6, 7, 8,21, 9,10,11,12,21,13,14,15,16,17,21,18,19,21,20]
               # A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y
    global letter2num
    letter2num(c) = letter2num(Char(c))
    function letter2num(c::Char)
        i = c - 'A' + 1
        1 <= i <= 25 && return alphabet[i]
        return 21
    end
end
