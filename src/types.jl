# This file is a part of ParalogMatching.jl. License is GPL3+: http://github.com/Mirmu/ParalogMatching.jl/LICENCE.md

########################### BLOCK FOR TYPES DEFINITION AND ALLOCATION #####################

immutable Alignment
    N::Int
    M::Int
    q::Int
    Z::Matrix{Int8}
    sequence::Vector{String}
    header::Vector{String}
    spec_name::Vector{String}
    spec_id::Vector{Int}
    uniprot_id::Vector{String}
end

@compat Base.:(==)(X1::Alignment, X2::Alignment) = all(fn->getfield(X1, fn) == getfield(X2, fn), fieldnames(Alignment))

immutable HarmonizedAlignments
    X1::Alignment
    X2::Alignment
end

@compat Base.:(==)(H1::HarmonizedAlignments, H2::HarmonizedAlignments) =
    all(fn->getfield(H1, fn) == getfield(H2, fn), fieldnames(HarmonizedAlignments))

# FreqC is the type of the frequency matrix
# Pij contains the joint (2 points) frequency
# Pi contains the average (1 point) frequency
# M contains the number of sequences and their pseudo count
# matching contains the current matching between the two alignments used to compute FreqC
# q is the number of different amino acids
type FreqC
    Pij::Matrix{Float64}
    Pi::Vector{Float64}
    specs::Vector{Int}
    M::Vector{Int}
    matching::Vector{Int}
    q::Int
    function FreqC(X1::Alignment, X2::Alignment)
	@extract X1 : N1=N M q
	@extract X2 : N2=N
	N = N1 + N2
	s = q - 1
	return new(zeros(s * N, s * N), zeros(N * s), [], [0, 0], zeros(M), q)
    end
end

# FastC is the type of the correlation matrix
# specs is the specs presently matched and used for the computation
# M contains the number of sequences and their pseudo count
type FastC
    Cij::Matrix{Float64}
    specs::Vector{Int}
    M::Vector{Int}
    FastC(N) = new(zeros(N, N), [], [0, 0])
end
FastC(F::FreqC) = FastC(length(F.Pi))
