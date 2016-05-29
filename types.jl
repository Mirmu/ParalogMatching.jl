############################BLOCK FOR TYPES DEFINITION AND ALLOCATION######################

#FastC is the type of the correlation matrix
#specs is the specs presently matched and used for the computation
#M contains the nber of seqs and their pseudo count

type FastC
    Cij::Array{Float64,2}
    specs::Array{Int,1}
    M::Array{Int,1}
end

#FreqC is the type of the frequency matrix
#Pij contains the joint (2 points) frequency
#Pi contains the average (1 point) frequency
#M contains the nber of seqs and their pseudo count
#matching contains the current matching between the two alignment used to compute FreqC
#q is the number of different amino acids

type FreqC
    Pij::Array{Float64,2}
    Pi::Array{Float64,1}
    specs::Array{Int,1}
    M::Array{Int,1}
    matching::Array{Int,1}
    q::Int
end 

#MatchOut is a type containing the matching as an output of the Gurobi optimization module
#N1 is the number of seqs for alignment 1
#N2 is the number of seqs for alignment 2
#fmin is the optimal value
#sol is the solution matching
#cost is the cost matrix given for the matching problem


type MatchOut
    N1::Int
    N2::Int
    fmin::Float64
    sol::Array{Float64,2}
    cost::Array{Float64,2}
    status::Symbol
end

#nullF pre allocates an empty freq matrix for two given alignments

function nullF(X1::Alignment,X2::Alignment)
    @extract X1 q N1=N M
    @extract X2 N2=N
    N=N1+N2
    s=q-1
    return FreqC(zeros(s*N, s*N),zeros(N*s),[],[0,0],zeros(M),q)
end

#Pre allocates an empty corr matrix for two given alignments

function nullC(F::FreqC)
    @extract F Pi
    N=length(Pi)
    return FastC(zeros(N, N),[],[0,0])
end
