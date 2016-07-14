###########################OPTIMIZATION AND MATCHING BLOCK####################################

function mpbmatch{T<:AbstractFloat}(D::DenseMatrix{T}, solver::MathProgBase.AbstractMathProgSolver)
    N1, N2 = size(D)

    A = matching_matrix(N1, N2)
    dir1 = N2 ≤ N1 ? '=' : '<'
    dir2 = N1 ≤ N2 ? '=' : '<'
    dirs = [k ≤ N2 ? dir1 : dir2 for k = 1:(N1+N2)]

    f = vec(D)

    #sol = mixintprog(f, A, dirs, 1, :Int, 0, 1, solver)
    sol = linprog(f, A, dirs, 1.0, 0.0, 1.0, solver)
    sol.status == :Optimal || error("failed: status = $(sol.status)")

    return reshape(sol.sol, N1, N2)::Matrix{Float64}
end

function matching_matrix(N1::Int, N2::Int)
    I = Array{Int}(2 * N1 * N2)
    J = Array{Int}(2 * N1 * N2)
    V = ones(2 * N1 * N2)

    t = 1
    for i = 1:N2, j = 1:N1
        k = (i-1) * N1 + j
        I[t] = i
        J[t] = k
        t += 1
    end
    for i = 1:N1, j = 1:N2
        k = i + (j-1) * N1
        I[t] = i + N2
        J[t] = k
        t += 1
    end

    return sparse(I, J, V, N1 + N2, N1 * N2)
end

######################HELPERS FOR THE MATCHING PROBLEM#################################

# Gets both the lines of alignments matched for a given list of specs
function get_edges(X1, X2, match, specs::Vector{Int64})
    ind = find([a in specs for a in X1.spec_id])
    filter!(x->match[x]!=0, ind)
    return ind, match[ind]
end

# Function that computes the frequency count, for a given matching.
# nspecs corresponds to the specs added at this step, defined by their spec_id
# You have to pass the prior because it is in-place
function unitFC!(X1, X2, match, nspecs, prior)
    @extract X1 : N1=N Z1=Z N1=N q
    @extract X2 : N2=N Z2=Z N2=N
    @extract prior : Pij Pi M

    N = N1 + N2
    s = q - 1
    Ns  = N * s

    # First removing the species to be recomputed
    if !isempty(intersect(prior.specs, nspecs))

        recomp_specs = intersect(prior.specs, nspecs)
        #println("ALERT : removing species :", recomp_specs)
        prior_match = prior.matching
        edg1, edg2 = get_edges(X1, X2, prior_match, recomp_specs)
        len = length(edg1)

        Z1a = Z1[edg1,:]
        Z2a = Z2[edg2,:]
        Zt = hcat(Z1a, Z2a)'

        ZZ = Vector{Int8}[vec(Zt[i,:]) for i = 1:N]

        # first removing
        @inbounds begin
            i0 = 0
            for i = 1:N
                Zi = ZZ[i]
                for k in 1:len
                    a = Zi[k]
                    a == q && continue
                    Pi[i0 + a] -= 1.0
                end
                i0 += s
            end

            i0 = 0
            for i = 1:N
                Zi = ZZ[i]
                j0 = 0
                for j = 1:N
                    Zj = ZZ[j]
                    for k = 1:len
                        a = Zi[k]
                        b = Zj[k]
                        (a == q || b == q) && continue
                        Pij[i0+a, j0+b] -= 1.0
                    end
                    j0 += s
                end
                i0 += s
            end
        end

        ####
        prior.specs = setdiff(prior.specs, recomp_specs)
        M[1] = M[1] - len
    end

    # check the species have been removed
    isempty(intersect(prior.specs, nspecs)) || error("redundant species")

    edg1, edg2 = get_edges(X1, X2, match, nspecs)
    len = length(edg1)

    Z1a = Z1[edg1,:]
    Z2a = Z2[edg2,:]
    Zt = hcat(Z1a, Z2a)'

    ZZ = Vector{Int8}[vec(Zt[i,:]) for i = 1:N]

    @inbounds begin
        i0 = 0
        for i = 1:N
            Zi = ZZ[i]
            for k in 1:len
                a = Zi[k]
                a == q && continue
                Pi[i0 + a] += 1.0
            end
            i0 += s
        end

        i0 = 0
        for i = 1:N
            Zi = ZZ[i]
            j0 = 0
            for j = 1:N
                Zj = ZZ[j]
                for k = 1:len
                    a = Zi[k]
                    b = Zj[k]
                    (a == q || b == q) && continue
                    Pij[i0+a, j0+b] += 1.0
                end
                j0 += s
            end
            i0 += s
        end
    end

    prior.specs = sort([prior.specs;nspecs])
    M[1] = M[1] + len
    prior.matching[:] = match[:]
    return nothing
end

# Computes the correlation matrix knowing the frequency matrix
function full_COD!(prev::FastC, freq::FreqC)
    @extract freq : Pij Pi specs M
    @extract prev : Cij

    m = M[1] + M[2]

    m == 0 && return nothing

    #Cij[:] = 1.0 / m * Pij[:] - 1.0 / m^2 * (Pi * Pi')[:]
    N = length(Pi)
    copy!(Cij, Pij)
    scale!(Cij, 1/m)
    @inbounds for j = 1:N
        Pj = Pi[j]
        @simd for i = 1:N
            Cij[i, j] -= (Pi[i] * Pj) / m^2
        end
    end

    prev.M[1] = M[1]
    prev.M[2] = M[2]
    resize!(prev.specs, length(specs))
    copy!(prev.specs, specs)

    return nothing
end

# Converts the sequences in binary sequences for the Gaussian model
function expand_binary(Z, s)
    a, b = size(Z)
    expZ = zeros(a, b * s)
    for i in 1:b, j in 1:a
        Z[j,i] != s + 1 ? expZ[j,(i-1)*s+Z[j,i]] += 1 : continue
    end
    return expZ
end

########################################BLOCK FOR COMPUTING THE MATCHING OF A GIVEN SPECIES, KNOWING THE PRIORS#################################

# This function takes alignments, priors matrices and ONE spec, and computes the matching following various strategies
# The helpers for the strategies are below
function give_correction(X1, X2, freq::FreqC, invertC::Matrix{Float64}, spec::Int, strat::AbstractString)
    @extract X1 : spec1=spec_id Z1=Z N1=N s=q-1
    @extract X2 : spec2=spec_id Z2=Z N2=N
    @extract freq : M Pi

    r1 = 1:(s*N1)
    r2 = (s*N1+1):(s*(N1+N2))

    #solver = GurobiSolver(OutputFlag=false)
    #solver = ClpSolver()
    #solver = GLPKSolverLP()
    solver = MathProgBase.defaultLPsolver # TODO: allow choosing the solver
    #solver = MathProgBase.defaultMIPsolver # TODO: allow choosing the solver

    # First covariation/greedy strategy
    if strat ∈ ("covariation", "greedy")
        ind1 = find(spec1 .== spec)
        ind2 = find(spec2 .== spec)

        # takes the sequences of the "spec"
        Zb1 = expand_binary(Z1[ind1,:], s)
        Zb2 = expand_binary(Z2[ind2,:], s)

        m = M[1] + M[2]

        # extracts the mean of the prior (one could add the additional contribution brought by the "spec")
        # but it does not change the results and we consider non bijective case here
        vmean1 = repmat(Pi' / m, length(ind1))[:,r1]
        vmean2 = repmat(Pi' / m, length(ind2))[:,r2]

        # the cost function as the matching algorithm wants it
        cost = (Zb1-vmean1) * invertC[r1,r2] * (Zb2-vmean2)'

        if strat == "covariation"
            # Compute the matching via linear programming
            rematch = mpbmatch(cost, solver)
            # and finally the result is converted to a matching
            permres = convert_perm_mat(rematch)
        elseif strat == "greedy"
            # Compute the matching using a greedy heuristic
            permres = greedy_from_cost(cost)
        else
            error("bug")
        end

    # Genetic matching strategy by genetic distance
    elseif strat == "genetic"
        # computes the genetic distance between sequences and returns the matching
        cost = cost_from_annot(X1, X2, spec)
        rematch = mpbmatch(cost, solver)

        permres = filter_gen_dist(convert_perm_mat(rematch), cost)

    # Random matching to test null hypothesis
    elseif strat == "random"
        ind1 = find(spec1 .== spec)
        ind2 = find(spec2 .== spec)
        mini = min(length(ind1), length(ind2))
        permres = (randperm(length(ind1))[1:mini], randperm(length(ind2))[1:mini])
    else
        error("option not known")
    end

    #println("for family $(spec), good matching is : ", permres)
    return permres
end

##################################HELPERS FOR THE MATCHING########################

# Convert the permutation matrix of Gurobi into a matchin array
# The matching format is
# ind1:rows of the first alignment
# ind2:matched rowd of the second alignment
function convert_perm_mat(mat)
    len1, len2 = size(mat)
    ind1 = find([!isempty(find(mat[i,:])) for i in 1:len1])
    ind2 = [findfirst(mat[i,:]) for i in ind1]
    return ind1, ind2
end


# Inverts the matrix of correlation after adding a pseudo count
let Cdict = Dict{Int,Matrix{Float64}}()
    global inverse_with_pseudo!, clear_inverse_mem
    function inverse_with_pseudo!(prev::FastC, freq::FreqC, pc::Float64)
        @extract freq : M Pi
        Mfake = round(Int, ((M[1]+M[2])*pc - M[2]) / (1-pc))

        L = length(Pi)
        CC = Base.@get!(Cdict, L, Array{Float64}(L, L))

        add_pseudocount!(freq, Mfake)
        full_COD!(prev, freq)
        copy!(CC, prev.Cij)
        inter = Base.LinAlg.inv!(cholfact!(CC))
        return inter
    end
    clear_inverse_mem() = empty!(Cdict)
end


# Adds a pseudo count with Mfake sequences
function add_pseudocount!(freq::FreqC, Mfake::Int64)
    Mfake != 0 || return

    @extract freq : Pi Pij M specs q

    s = q - 1
    L = length(Pi)
    @assert L % s == 0
    N = L ÷ s

    pcq = Mfake / q
    pcqq = pcq / q

    @inbounds for j = 1:L, i = 1:L
        Pij[i, j] += pcqq
    end
    @inbounds for i = 1:L
        Pi[i] += pcq
    end

    i0 = 0
    @inbounds for i = 1:N
        xr = i0 + (1:s)
        for x2 in xr, x1 in xr
            Pij[x1, x2] -= pcqq
        end
        i0 += s
    end
    @inbounds for i = 1:L
        Pij[i, i] += pcq
    end

    M[2] += Mfake

    return
end

function annot2num(annotb::AbstractString)
    if length(annotb) < 6
        println("non valid uniprot : ", annotb)
        return 0
    end
    annot = annotb[end-5:end]
    bs = [36, 10, 36, 36, 36, 10]
    dgts = [ parse(Int, string(annot[i]), bs[i]) for i = 1:6 ]
    mbs = cumprod(reverse([10, 36, 36, 36, 10, 1]))
    return dot(reverse(dgts), mbs)
end

function cost_from_annot(X1, X2, spec)
    @extract X1 : spec1=spec_id Z1=Z N1=N
    @extract X2 : spec2=spec_id Z2=Z N2=N

    ind1 = find(spec1 .== spec)
    ind2 = find(spec2 .== spec)
    cost = Float64[abs(annot2num(X1.uniprot_id[i]) - annot2num(X2.uniprot_id[j])) for i in ind1, j in ind2]
    return cost
end

function filter_gen_dist(match, cost)
    res1 = Int64[]
    res2 = Int64[]

    #println(collect(zip(match...)))

    for (i,j) in zip(match...)
        abs(cost[i,j]) < 100 || continue
        push!(res1, i)
        push!(res2, j)
    end
    return res1, res2
end

function greedy_from_cost(mat)
    len1, len2 = size(mat)
    count = 0
    perm1 = Int64[]
    perm2 = Int64[]
    rank = 1

    while count < min(len1, len2)
        a, b = findnth_mat(mat, rank)
        rank += 1
        if !(a in perm1) && !(b in perm2)
            push!(perm1, a)
            push!(perm2, b)
            count += 1
        end
    end
    return perm1, perm2
end

# Finds the MINIMUM of a matrice and returns the index
function findnth_mat(A::Matrix{Float64}, rank::Int64)
    ls = sortperm(A[:])
    i, j = ind2sub(size(A), ls[rank])
    return i, j
end
