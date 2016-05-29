############################BLOCK FOR RESULTS PLOTTING/SCORING#########################

#Computing the ROC Curve
#ndata: distance file
#plm : scoring file for the contact
#Split: separation between the two proteins

function ROCCurve(ndata,plm,split)
	score=readdlm(plm)
	data=readdlm(ndata)
	comp=data[:,1]+im*data[:,2]
	
	intra=0
	trueintra=0
	inter=0
	trueinter=0
	
	curveinter=Float64[]
	curveintra=Float64[]

	for a in 1:size(score)[1]
	
		#ugliness
		ind=findfirst(comp,(score[a,1])+(score[a,2])*im)
		
		ind!=0||continue

		if (((score[a,1]<=split)&&(score[a,2]<=split))||((score[a,1]>split)&&(score[a,2]>split)))&&(abs(score[a,1]-score[a,2])>4)
			data[ind,3]>8||(trueintra+=1)
			intra+=1
			push!(curveintra,trueintra/intra)

		elseif (score[a,1]<=split)&&(score[a,2]>split)
			data[ind,3]>8||(trueinter+=1)
			inter+=1
			push!(curveinter,trueinter/inter)w		end
	end
	return curveinter,curveintra
end

#Plots the best possible ROC curve

function CurveTrue(ordering)
	trueprop=0
	count=0
	sorted=sort([ordering...],by=x->x[2])
	rocc=Float64[]
	for el in sorted
		if el[1][1]==el[1][2]
			trueprop+=1
		end
		count+=1
		push!(rocc,trueprop/count)
	end
	return rocc
end


#########################EDGES SCORING BLOCK##########################################

#This function takes the result and recomputes the model
#then it ranks the edges for the giving model

function PostOrderScores(X1,X2,match)
	println("Recomputing the model and sorting the edges... <3")

	freq=nullF(X1,X2)
	corr=nullC(freq)
	single=X1.SpecId[find(match)]
	UnitFC!(X1,X2,match,single,freq)
	FullCOD!(corr,freq)
	invC=inverseWithPseudo!(corr,freq,0.8)

	lines=length(match)
	#here we consider the pseudo count contribution as well
	Meff=freq.M[1]+freq.M[2]
	aver=freq.Pi/Meff
	flag=1
	ordering=Tuple{Float64,Float64,Float64}[]

	len1=20*X1.N
	len2=20*X2.N
	nullconst=1/2*(logdet(corr.Cij)-logdet(corr.Cij[1:len1,1:len1])-logdet(corr.Cij[len1+1:end,len1+1:end]))

	for i in 1:lines
		match[i]!=0||continue

		seq1=expandBinary(X1.Z[i,:],20)[:]
		seq2=expandBinary(X2.Z[match[i],:],20)[:]
		seq=[seq1;seq2]
		score=-(((seq-aver)[1:len1]')*invC[1:len1,len1+1:end]*((seq-aver)[len1+1:end]))[1]-nullconst
		println((i,match[i],score))

		push!(ordering,(i,match[i],score))
	end
	sort!(ordering,by=x->x[3],rev=true)
	tab=hcat([[a...] for a in ordering]...)'
	lab=[convert(Float64,tab[i,1]==tab[i,2]) for i in 1:size(tab)[1]]
	return [tab lab]
end

function FullEdgesScores(X1,X2,match)
	println("Recomputing the model and sorting the edges... <3")

	freq=nullF(X1,X2)
	corr=nullC(freq)
	single=X1.SpecId[find(match)]
	UnitFC!(X1,X2,match,single,freq)
	FullCOD!(corr,freq)
	invC=inverseWithPseudo!(corr,freq,0.8)

	lines=length(match)
	#here we consider the pseudo count contribution as well
	Meff=freq.M[1]+freq.M[2]
	aver=freq.Pi/Meff
	flag=1
	ordering=Array{Tuple{Tuple{Float64,Float64},Float64},1}[]
	suborder=Tuple{Tuple{Float64,Float64},Float64}[]
	#creates all the possibles within species match
	pairings=[find(x->x==i,X2.SpecId) for i in X1.SpecId]
	len1=20*X1.N
	len2=20*X2.N
	nullconst=1/2*(logdet(corr.Cij)-logdet(corr.Cij[1:len1,1:len1])-logdet(corr.Cij[len1+1:end,len1+1:end]))
	
	println(nullconst)
	for i in 1:lines
		for j in pairings[i]
			seq1=expandBinary(X1.Z[i,:],20)[:]
			seq2=expandBinary(X2.Z[j,:],20)[:]

			seq=[seq1;seq2]
			score=-(((seq-aver)[1:len1]')*invC[1:len1,len1+1:end]*((seq-aver)[len1+1:end]))[1]-nullconst
			println(((i,j),score))
			if (X1.SpecId[i]==flag)
				push!(suborder,((i,j),score))
			end

			if (X1.SpecId[i]!=flag)||(i==lines)
				sort!(suborder,by=x->x[2],rev=true)
				push!(ordering,suborder)
				suborder=Tuple{Tuple{Float64,Float64},Float64}[]
				push!(suborder,((i,j),score))
				flag=X1.SpecId[i]
			end
		end
	end
	return X1,X2,match,ordering
end

function SelectingEdges(ordering,thres)
	order2=Array{Tuple{Tuple{Float64,Float64},Float64},1}[]
	for el in ordering
		len=length(el)
		extr=sort(el,by=x->x[2])[1:min(len,thres)]
		push!(order2,extr)
	end
	return order2
end

function TrimCovMatch(nameX1,nameX2,namematch,prop,nameoutput)
	X1,X2,match,ret=PostOrderScores(nameX1,nameX2,namematch)
	println("The distribution of family sizes is the following (size,number):")
	println(Tally([length(a) for a in ret]))
	retsor=sort(vcat(ret...),by=x->x[2])
	len=length(retsor)
	propo=round(Int,prop*len)
	edges1=Int64[a[1][1] for a in retsor[1:propo]]
	edges2=Int64[a[1][2] for a in retsor[1:propo]]
	nuovomatch=zeros(X1.M)
	nuovomatch[edges1]=edges2

	RewriteFastaMatch(X1,X2,nuovomatch,nameoutput)

end


