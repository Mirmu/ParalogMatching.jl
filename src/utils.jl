###########################OPTIMIZATION AND MATCHING BLOCK####################################

#First the Gurobi module that solves the linear matching problem

function gurobimatch{T<:AbstractFloat}(D::DenseArray{T,2})
    
    env = Gurobi.Env()
    setparam!(env, "OutputFlag", false) # set verbose to 0
    N1,N2 = size(D)
    mylb = zeros(T,N1*N2)
    myub = ones(T,N1*N2)
    model = gurobi_model(env; sense=:minimize, lb=mylb, ub=myub, f=D[:])    
    constr1 = zeros(Int,N1)
    myones1 = ones(Int,N1)        
    constr2 = zeros(Int,N2)
    myones2 = ones(T,N2)    
    if N1<N2    
        for i=1:N1
            for j=1:N2
                constr2[j] = sub2ind((N1,N2), i,j)
            end
            add_constr!(model, constr2, myones2,'=',one(T)) 
        end
        for j=1:N2
            for i=1:N1
                constr1[i] = sub2ind((N1,N2),i,j)
            end
            add_constr!(model, constr1, myones1,'<',one(T)) 
        end
    elseif N1 > N2
        for j=1:N2
            for i=1:N1
                constr1[i] = sub2ind((N1,N2),i,j)
            end
            add_constr!(model, constr1, myones1,'=',one(T)) 
        end
        for i=1:N1
            for j=1:N2
                constr2[j] = sub2ind((N1,N2),i,j)
            end
            add_constr!(model, constr2, myones2,'<',one(T)) 
        end
    else
        for i=1:N1
            for j=1:N2
                constr2[j] = sub2ind((N1,N2), i,j)
            end
            add_constr!(model, constr2, myones2,'=',one(T)) 
        end
        for j=1:N2
            for i=1:N1
                constr1[i] = sub2ind((N1,N2),i,j)
            end
            add_constr!(model, constr1, myones1,'=',one(T)) 
        end
    end

    update_model!(model)


    mytime = @elapsed optimize(model)
    #println("elapsed time = $mytime")
    status = Gurobi.get_status(model)
    
    return MatchOut(N1,
                    N2, 
                    get_objval(model),
                    reshape(get_solution(model),(N1,N2)),
                    D,
                    status)

end


######################HELPERS FOR THE MATCHING PROBLEM#################################

#gets both the lines of alignments matched for a given list of specs
function GetEdges(X1,X2,match,specs::Array{Int64,1})
    ind=find([a in specs for a in X1.SpecId])
    filter!(x->match[x]!=0,ind)
    return ind,match[ind]
end

#Function that computes the frequency count, for a given matching.
#nspecs corresponds to the specs added at this step, defined by their SpecId
#You have to pass the prior because it is in-place

function UnitFC!(X1,X2,match,nspecs,prior)
    
    @extract X1 N1=N Z1=Z q N1=N
    @extract X2 N2=N Z2=Z N2=N
    @extract prior Pij Pi M
    

    N=N1+N2
    s = q - 1 
    Ns  = N * s
    #first removing the species to be recomputed

    if intersect(prior.specs,nspecs)!=[]

        recomp_specs=intersect(prior.specs,nspecs)
        #println("ALERT : removing species :",recomp_specs)
        prior_match=prior.matching
        edg1,edg2=GetEdges(X1,X2,prior_match,recomp_specs)
        len=length(edg1)

        Z1a=Z1[edg1,:]
        Z2a=Z2[edg2,:]
        Zt = (hcat(Z1a,Z2a))'
        
        ZZ = Vector{Int8}[vec(Zt[i,:]) for i = 1:N]

        #first removing 
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
        prior.specs=setdiff(prior.specs,recomp_specs)
        M[1] = M[1] - len
    end

    #check the species have been removed
    intersect(prior.specs,nspecs)==[]||error("redundent species")

    edg1,edg2=GetEdges(X1,X2,match,nspecs)
    len=length(edg1)

    Z1a=Z1[edg1,:]
    Z2a=Z2[edg2,:]
    Zt = (hcat(Z1a,Z2a))'
    
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

    prior.specs=sort([prior.specs;nspecs])
    M[1] = M[1] + len
    prior.matching[:]=match[:]
    return nothing  
end

#Computes the correlation matrix knowing the frequency matrix

function FullCOD!(prev::FastC,freq::FreqC)

    @extract freq Pij Pi specs M
    @extract prev Cij

    m=M[1]+M[2]

    if m==0 
        nothing
    else
        Cij[:] = 1.0/m * Pij[:] - 1.0/m^2 * (Pi * Pi')[:]
        prev.M[1]=M[1]
        prev.M[2]=M[2]
        prev.specs=copy(specs)
    end
    return nothing
end

#Converts the sequences in binary sequences for the Gaussian model
function expandBinary(Z,s)
    a,b=size(Z)
    expZ=zeros(a,b*s)
    for i in 1:b, j in 1:a
        Z[j,i]!=s+1?expZ[j,(i-1)*s+Z[j,i]]+=1:continue
    end
    return expZ
end

########################################BLOCK FOR COMPUTING THE MATCHING OF A GIVEN SPECIES, KNOWING THE PRIORS#################################

#This function takes alignments, priors matrices and ONE spec, and computes the matching following various strategies
#The helpers for the strategies are below

function giveCorrection(X1,X2, freq::FreqC, invertC::Array{Float64,2}, spec::Int,strat::AbstractString)
    
    @extract X1 Spec1=SpecId Z1=Z N1=N
    @extract X2 Spec2=SpecId Z2=Z N2=N
    @extract freq M Pi
    
    #first covariation strategy

    if strat=="covariation"

        ind1=find(Spec1.==spec)
        ind2=find(Spec2.==spec)
        #takes the sequences of the "spec"
        Zb1=expandBinary(Z1[ind1,:],20)
        Zb2=expandBinary(Z2[ind2,:],20)
    
        m=M[1]+M[2]

	#extracts the mean of the prior (one could add the additional contribution brought by the "spec")
	#but it does not change the results and we consider non bijective case here

        vmean1=repmat(1.0/m*Pi',length(ind1))[:,1:20*N1]
        vmean2=repmat(1.0/m*Pi',length(ind2))[:,20*N1+1:end]
        
	#The cost function as Gurobi wants it
        cost=(Zb1-vmean1)*invertC[1:20*N1,20*N1+1:end]*(Zb2-vmean2)'
        rematch=gurobimatch(cost)

	#and finally the results is converted to a matching
        permres=convertPermMat(rematch.sol)

    #genetic matching strategy by genetic distance
    elseif strat=="genetic"
	
    	#Computes the genetic distance between sequences and returns the matching
        cost=CostFromAnnot(X1,X2,spec)
        rematch=gurobimatch(cost)

        permres=FilterGenDist(convertPermMat(rematch.sol),cost)

    #Greedy does like covariation, but the returned matching follows a greedy heuristic
    #Works well also in general

    elseif strat=="greedy"

        ind1=find(Spec1.==spec)
        ind2=find(Spec2.==spec)
        Zb1=expandBinary(Z1[ind1,:],20)
        Zb2=expandBinary(Z2[ind2,:],20)
    
        m=M[1]+M[2]
        vmean1=repmat(1.0/m*Pi',length(ind1))[:,1:20*N1]
        vmean2=repmat(1.0/m*Pi',length(ind2))[:,20*N1+1:end]
    
        cost=(Zb1-vmean1)*invertC[1:20*N1,20*N1+1:end]*(Zb2-vmean2)'
        permres=GreedyFromCost(cost)

    #Random matching to test null hypothesis
    elseif strat=="random"
        ind1=find(Spec1.==spec)
        ind2=find(Spec2.==spec)
        mini=min(length(ind1),length(ind2))
        permres=(randperm(length(ind1))[1:mini],randperm(length(ind2))[1:mini])
    
    else 
        error("option not known")

    end

    #println("for family $(spec), good matching is : ", permres)
    return permres
end

##################################HELPERS FOR THE MATCHING########################

#Convert the permutation matrix of Gurobi into a matchin array
#The matching format is 
#ind1:rows of the first alignment
#ind2:matched rowd of the second alignment

function convertPermMat(mat)
    len1,len2=size(mat)
    ind1=find([find(mat[i,:])!=[] for i in 1:len1])
    ind2=[findfirst(mat[i,:]) for i in ind1]
    return ind1,ind2
end


#Inverts the matrix of correlation after adding a pseudo count

function inverseWithPseudo!(prev::FastC,freq::FreqC,pc::Float64)
    @extract freq M Pij Pi
    Mfake=round(Int,((M[1]+M[2])*pc-M[2])/(1-pc))

    add_pseudocount!(freq,Mfake)
    FullCOD!(prev,freq)
    inter=inv(cholfact(prev.Cij))
    return inter
end


#Adds a pseudo count with Mfake sequences


function add_pseudocount!(freq::FreqC, Mfake::Int64)
    
    Mfake!=0||return nothing

    @extract freq Pi Pij M specs q
    
    s = q - 1
    N=round(Int,length(Pi)/s)
    m=M[1]+M[2]

    pc = Mfake/(Mfake+m)
    pcq = pc / q

    Pij_c = 1 / (m+Mfake) * Pij .+ pcq / q
    Pi_c = 1 / (m+Mfake) * Pi .+ pcq

   
    i0 = 0
    for i = 1:N
    xr = i0 + (1:s)
    Pij_c[xr, xr] = 1 / (m+Mfake) * Pij[xr, xr]
        for alpha = 1:s
            x = i0 + alpha
            Pij_c[x, x] += pcq
        end
        i0 += s
    end

    Pi[:]=(Mfake+m)*Pi_c[:]
    Pij[:]=(Mfake+m)*Pij_c[:]
    M[:]=[M[1],M[2]+Mfake][:]
end


function annot2num(annotb::AbstractString)
    if length(annotb) < 6
        println("non valid uniprot : ",annotb)
        return 0
    end
    annot=annotb[end-5:end]
    bs = [36, 10, 36, 36, 36, 10]
    dgts = [ parse(Int,string(annot[i]), bs[i]) for i = 1:6 ]
    mbs = cumprod(reverse([10, 36, 36, 36, 10, 1]))
    return dot(reverse(dgts), mbs)
end

function CostFromAnnot(X1,X2,spec)
    @extract X1 Spec1=SpecId Z1=Z N1=N
    @extract X2 Spec2=SpecId Z2=Z N2=N

    ind1=find(Spec1.==spec)
    ind2=find(Spec2.==spec)
    cost=Float64[abs(annot2num(X1.UniprotId[i])-annot2num(X2.UniprotId[j])) for i in ind1, j in ind2]
    return cost
end

function FilterGenDist(match,cost)
    res1=Int64[]
    res2=Int64[]

    #println(collect(zip(match...)))

    for (i,j) in zip(match...)
        abs(cost[i,j])<100||continue
        push!(res1,i)
        push!(res2,j)
    end
    return res1,res2
end

function GreedyFromCost(mat)
    len1,len2=size(mat)
    count=0
    perm1=Int64[]
    perm2=Int64[]
    rank=1

    while count<min(len1,len2)
        a,b=findnth_mat(mat,rank)
        rank+=1
        if !(a in perm1)&&!(b in perm2)
            push!(perm1,a)
            push!(perm2,b)
            count+=1
        end
    end
    return perm1,perm2
end

#Finds the MINIMUM of a matrice and returns the index
function findnth_mat(A::Array{Float64,2},rank::Int64)
    ls=sortperm(A[:])
    i,j= ind2sub(size(A), ls[rank])
    return i,j
end

