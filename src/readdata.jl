# This file is a part of ParalogMatching.jl. License is GPL3+: http://github.com/Mirmu/ParalogMatching.jl/LICENCE.md


"""
    read_fasta_alignment(filename, max_gap_fraction=1.0; header_regex=nothing)

Parse a (possibly gzipped) FASTA file, converting it into an alignment which can be used
in the Gaussian model used by ParalogMatching. Optionally, discards the sequences which
have more than a fraction `max_gap_fraction` of missing values (gaps) in the alignment.

The function automatically tries to parse the species names and the Uniprot ID, using
a standard format specification. The species names are then used in the matching procedure;
the Uniprot IDs are only used for writing the output files.

You can parse arbitrary files by providing a custom `header_regex` (see also the
[Julia regex documentation](http://docs.julialang.org/en/stable/manual/strings/#regular-expressions)).
In its simples form, it needs to have at least two capture groups, the first one returning the
Uniprot ID and the second the species name (additional capture groups are ignored).
For example, if the headers in your FASTA file look like this:

```text
>SOMELOCATION/SOMESPECIES/OTHERINFO
```

then you can use a regex like this one: `read_fasta_alignment(..., header_regex=r"^([^/]+)/([^/]+)/")`,
i.e. line start, anything except a slash (captured), followed by a slash, then anything except a slash
(captured), then a slash — the remainder of the line is then simply ignored. In more complicated cases
(e.g. if the UniprotID and the species name are out of order), you can use named capture groups: in this
case the `header_regex` needs to contain both an `id` group and a `species` group (all additional groups are
ignored). For example, if the headers in your FASTA file look like this:

```text
>SOMESPECIES/OTHERINFO/SOMELOCATION
```

then you can use a regex like this one: `header_regex=r"^(?<species>[^/]+)/([^/]+)/(?<id>.+)"`,
i.e. line start, anything except a slash (capture as the species name), followed by a slash, then anything
except a slash (captured but ignored), followed by a slash, then the rest of the line (captured as the ID).
"""
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

    length(seqs) > 0 || error("out of $(f.num_parsed) sequences, none passed the filter (max_gap_fraction=$max_gap_fraction)")

    # pass 2

    Z = Array{Int8}(fseqlen, length(seqs))
    header =  Array{String}(length(seqs));
    sequence =  Array{String}(length(seqs));
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

function specname(s::String, header_regex::Union{Void,Regex}, captureinds::NTuple{2,Integer})
    if header_regex ≢ nothing
        # user-defined format
        if ismatch(header_regex, s)
            captures = match(header_regex, s).captures
            length(captures) ≥ 2 ||
                error("invalid header regex: should always return at least 2 captured groups if it matches; has returned: $(length(captures))")
            uniprot_id, spec_name = captures[captureinds[1]], captures[captureinds[2]]
            isa(uniprot_id, AbstractString) || error("the capture group for `id` did not match the spec string: $s")
            isa(spec_name, AbstractString) || error("the capture group for `species` did not match the spec string: $s")
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
    return convert(String, uniprot_id), convert(String, spec_name)
end

function compute_spec(header::Vector{String}, header_regex::Union{Void,Regex} = nothing)
    M = length(header)

    spec_name = Array{String}(M)
    uniprot_id  = Array{String}(M)

    captureinds = (1,2)
    if header_regex ≢ nothing
        cnamesdict = Base.PCRE.capture_names(header_regex.regex)
        cnames = collect(values(cnamesdict))
        if "id" ∈ cnames && "species" ∈ cnames
            ckeys = collect(keys(cnamesdict))
            captureinds = (ckeys[findfirst(i->cnamesdict[i]=="id", ckeys)],
                           ckeys[findfirst(i->cnamesdict[i]=="species", ckeys)])
        elseif "id" ∉ cnames && "species" ∉ cnames
            nothing
        else
            error("only one of `id` and `species` capture names found in `header_regex`; either none or both should be used")
        end
    end

    for i = 1:M
        uniprot_id[i], spec_name[i] = specname(header[i], header_regex, captureinds)
    end

    specunique = unique(spec_name)
    sdict = Dict{String,Int}(s=>i for (i,s) in enumerate(specunique))
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
