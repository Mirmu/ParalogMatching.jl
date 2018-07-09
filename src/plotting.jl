# This file is a part of ParalogMatching.jl. License is GPL3+: http://github.com/Mirmu/ParalogMatching.jl/LICENCE.md

########################### BLOCK FOR RESULTS PLOTTING/SCORING ########################

# Computing the ROC Curve
#   ndata: distance file
#   plm : scoring file for the contact
#   split: separation between the two proteins
function ROC_curve(ndata, plm, split)
    score = readdlm(plm)
    data = readdlm(ndata)
    comp = data[:,1] + im * data[:,2]

    intra = 0
    trueintra = 0
    inter = 0
    trueinter = 0

    curveinter = Float64[]
    curveintra = Float64[]

    for a in 1:size(score)[1]

        # ugliness
        ind = findfirst(comp, (score[a,1]) + (score[a,2]) * im)

        ind != 0 || continue

        if (((score[a,1] <= split) && (score[a,2] <= split)) || ((score[a,1] > split) && (score[a,2] > split))) && (abs(score[a,1] - score[a,2]) > 4)
            data[ind,3] > 8 || (trueintra += 1)
            intra += 1
            push!(curveintra, trueintra / intra)

        elseif (score[a,1] <= split) && (score[a,2] > split)
            data[ind,3] > 8 || (trueinter += 1)
            inter += 1
            push!(curveinter, trueinter / inter)
        end
    end
    return curveinter, curveintra
end

# Plots the best possible ROC curve
function curve_true(ordering)
    trueprop = 0
    count = 0
    sorted = sort([ordering...], by=x->x[2])
    rocc = Float64[]
    for el in sorted
        if el[1][1] == el[1][2]
            trueprop += 1
        end
        count += 1
        push!(rocc, trueprop / count)
    end
    return rocc
end


######################## EDGES SCORING BLOCK #########################################

# This function takes the result and recomputes the model
# then it ranks the edges for the giving model
function post_order_scores(X1, X2, match)
    println("Recomputing the model and sorting the edges... <3")

    freq = nullF(X1, X2)
    corr = nullC(freq)
    single = X1.spec_id[find(match)]
    unitFC!(X1, X2, match, single, freq)
    full_COD!(corr, freq)
    invC = inverse_with_pseudo!(corr, freq, 0.8)

    lines = length(match)

    # Here we consider the pseudo count contribution as well
    Meff = freq.M[1] + freq.M[2]
    aver = freq.Pi / Meff
    flag = 1
    ordering = Tuple{Float64,Float64,Float64}[]

    len1 = 20 * X1.N
    len2 = 20 * X2.N
    nullconst = 1/2 * (logdet(corr.Cij) - logdet(corr.Cij[1:len1,1:len1]) - logdet(corr.Cij[len1+1:end,len1+1:end]))

    for i in 1:lines
        match[i] != 0 || continue

        seq1 = expand_binary(X1.Z[i,:], 20)[:]
        seq2 = expand_binary(X2.Z[match[i],:], 20)[:]
        seq = [seq1;seq2]
        score = -(((seq-aver)[1:len1]') * invC[1:len1,len1+1:end] * ((seq-aver)[len1+1:end]))[1] - nullconst
        println((i, match[i], score))

        push!(ordering, (i, match[i], score))
    end
    sort!(ordering, by=x->x[3], rev=true)
    tab = hcat([[a...] for a in ordering]...)'
    lab = [convert(Float64, tab[i,1] == tab[i,2]) for i in 1:size(tab)[1]]
    return [tab lab]
end

function full_edges_scores(X1, X2, match)
    println("Recomputing the model and sorting the edges... <3")

    freq = nullF(X1, X2)
    corr = nullC(freq)
    single = X1.spec_id[find(match)]
    unitFC!(X1, X2, match, single, freq)
    full_COD!(corr, freq)
    invC = inverse_with_pseudo!(corr, freq, 0.8)

    lines = length(match)

    # Here we consider the pseudo count contribution as well
    Meff = freq.M[1] + freq.M[2]
    aver = freq.Pi / Meff
    flag = 1
    ordering = Vector{Tuple{Tuple{Float64,Float64},Float64}}[]
    suborder = Tuple{Tuple{Float64,Float64},Float64}[]

    # Creates all the possibles within species match
    pairings = [find(x->x==i, X2.spec_id) for i in X1.spec_id]
    len1 = 20 * X1.N
    len2 = 20 * X2.N
    nullconst = 1/2 * (logdet(corr.Cij) - logdet(corr.Cij[1:len1,1:len1]) - logdet(corr.Cij[len1+1:end,len1+1:end]))

    println(nullconst)
    for i in 1:lines
        for j in pairings[i]
            seq1 = expand_binary(X1.Z[i,:], 20)[:]
            seq2 = expand_binary(X2.Z[j,:], 20)[:]

            seq = [seq1;seq2]
            score = -(((seq-aver)[1:len1]') * invC[1:len1,len1+1:end] * ((seq-aver)[len1+1:end]))[1] - nullconst
            println(((i, j), score))
            if (X1.spec_id[i] == flag)
                push!(suborder, ((i, j), score))
            end

            if (X1.spec_id[i] != flag) || (i == lines)
                sort!(suborder, by=x->x[2], rev=true)
                push!(ordering, suborder)
                suborder = Tuple{Tuple{Float64,Float64},Float64}[]
                push!(suborder, ((i, j), score))
                flag = X1.spec_id[i]
            end
        end
    end
    return X1, X2, match, ordering
end

function selecting_edges(ordering, thres)
    order2 = Vector{Tuple{Tuple{Float64,Float64},Float64}}[]
    for el in ordering
        len = length(el)
        extr = sort(el,by=x->x[2])[1:min(len,thres)]
        push!(order2, extr)
    end
    return order2
end

function trim_cov_match(nameX1, nameX2, namematch, prop, nameoutput)
    X1, X2, match, ret = post_order_scores(nameX1, nameX2, namematch)
    println("The distribution of family sizes is the following (size,number):")
    println(tally([length(a) for a in ret]))
    retsor = sort(vcat(ret...), by=x->x[2])
    len = length(retsor)
    propo = round(Int, prop * len)
    edges1 = Int[a[1][1] for a in retsor[1:propo]]
    edges2 = Int[a[1][2] for a in retsor[1:propo]]
    nuovomatch = fill(0.0, X1.M)
    nuovomatch[edges1] = edges2

    write_fasta_match(X1, X2, nuovomatch, nameoutput)
end

#Function to return the Feinauer score of an alignment
function compute_score(score, split)
    counter = 0
    mean = 0
    l = size(score)[1]
    for i in 1:l
        if (score[i,1] <= split) && (score[i,2] > split)
            counter += 1
            mean += score[i, 3]
        end
        counter == 4 && break
    end
    return mean / 4
end

#For computing quickly the number of effective sequences in the Fasta
function my_read_fasta(filename::AbstractString)
    theta = :auto
    max_gap_fraction = 0.9
    remove_dups = true

    Z = GaussDCA.read_fasta_alignment(filename, max_gap_fraction)
    if remove_dups
        Z, _ = GaussDCA.remove_duplicate_seqs(Z)
    end

    N, M = size(Z)
    q = round(Int, maximum(Z))

    q > 32 && error("parameter q=$q is too big (max 31 is allowed)")
    W, Meff = GaussDCA.compute_weights(Z, q, theta)
    scale!(W, 1.0/Meff)
    Zint = round(Int, Z)
    return N, Meff
end
