using MacroUtils

type FastC
    Cij::Array{Float64,2}
    specs::Array{Int,1}
    M::Array{Int,1}
end

type FreqC
    Pij::Array{Float64,2}
    Pi::Array{Float64,1}
    specs::Array{Int,1}
    M::Array{Int,1}
    matching::Array{Int,1}
    q::Int
end 

type MatchOut
    N1::Int
    N2::Int
    fmin::Float64
    sol::Array{Float64,2}
    cost::Array{Float64,2}
    status::Symbol
end

function nullF(X1::Alignment,X2::Alignment)
    @extract X1 q N1=N M
    @extract X2 N2=N
    N=N1+N2
    s=q-1
    return FreqC(zeros(s*N, s*N),zeros(N*s),[],[0,0],zeros(M),q)
end

function nullC(F::FreqC)
    @extract F Pi
    N=length(Pi)
    return FastC(zeros(N, N),[],[0,0])
end
