abstract DataSet

immutable Alignment <: DataSet
    N::Int
    M::Int
    q::Int
    split::Int
    Z::Array{Int8,2}
    Sequence::Array{ASCIIString,1}
    Header::Array{ASCIIString,1}
    SpecName::Array{ASCIIString,1}
    SpecId::Array{Int,1}
    UniprotId::Array{ASCIIString,1}
end


readdata(filename::ASCIIString) = read_fasta_alignment(filename::ASCIIString)

read_fasta_alignment(filename::ASCIIString) = read_fasta_alignment(filename::ASCIIString, 1.)

function read_fasta_alignment(filename::ASCIIString, max_gap_fraction::Float64)

    f = FastaReader(filename)

    max_gap_fraction = Float64(max_gap_fraction)

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
    Header =  Array(ASCIIString,length(seqs));
    Sequence =  Array(ASCIIString,length(seqs));
    seqid = 1
    for (name, seq) in f
#        
#        Header[seqid] = specname.captures[1];
        Header[seqid] = name
        Sequence[seqid] = seq
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

    SpecId,SpecName,UniprotId = ComputeSpec(Header)

    Fa = Alignment(size(Z,1),size(Z,2), Int(maximum(Z)), 12,  Z', Sequence, Header, SpecName, SpecId,UniprotId)

    return Fa 
end

function specname(s::ASCIIString)

    if ismatch(r"\[(.*?)\]",s)
        SpecName = convert(ASCIIString,match(r"\[(.*?)\]", s).captures[1])
        UniprotId="000000"
        return (UniprotId,SpecName)
    elseif ismatch(r"^([A-Z,0-9].*?)\_([A-Z,0-9].*?)/",s)
        UniprotId,SpecName = match(r"^([A-Z,0-9].*?)\_([A-Z,0-9].*?)/", s).captures
        return (convert(ASCIIString,UniprotId),convert(ASCIIString,SpecName))

    
    elseif ismatch(r"^(.*?)with(.*?)/(.*?)$",s)
        SpecName = convert(ASCIIString,match(r"^(.*?)with(.*?)/(.*?)$", s).captures[3])
        UniprotId= "000000"
        return (UniprotId,SpecName)
    end

    n = length(s)
    i1=-1
    i2=-1
    flag1 = 0
    flag2 = 0
    for i=1:length(s)
        if s[i] == '/'
            i1 = i
            if flag1 == 1
                error("badly formed string")
            end
            flag1 = 1
        elseif s[i] == '_' 
            i2 = i
            if flag2 == 1
                error("badly formed string")
                end
            flag2 = 1        
        end
    end
    i1 < 0 || i2 < 0 && error("badly formed string")

    return (s[1:i1-1],s[i2+1:end])
end


function ComputeSpec(Header::Array{ASCIIString,1})

    M = length(Header)

    SpecName = Array{ASCIIString}(M)
    UniprotId  = Array{ASCIIString}(M)

    for i=1:M
        protname, spec = specname(Header[i])
        SpecName[i] = spec
        UniprotId[i] = protname
    end

    specunique = unique(SpecName)
    nspec = length(specunique)
    SpecId = zeros(Int,M)
    idx = zeros(Int,M)
  
    @inbounds for j=1:M
        nname = length(SpecName[j])
        for i=1:nspec
            SpecName[j] == specunique[i] || continue 
            SpecId[j] = i
            break
        end
    end

    return SpecId,SpecName,UniprotId
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
