using FastaIO, MacroUtils, Iterators
using HDF5, JLD
using Distributions

include("readdata.jl")
include("types.jl")
include("utils.jl")

function Tally(list)
	slist=sort(list)
	counter=0
	comp=slist[1]
	tallist=Tuple{eltype(list),Int64}[]

	for el in slist
		if el==comp
			counter+=1
		else
			push!(tallist,(comp,counter))
			comp=el
			counter=1
		end
	end
	push!(tallist,(comp,counter))
	return tallist
end

#WARNING !: works only for very specific lists in which the equal elements are clustered together

function IndexOfUnique(list)
	list==sort(list)||error("would you mind sorting the list ?")
	change=list[1]
	ind=Int64[]

	if list[1]!=list[2]
		push!(ind,1)
	end

	for i in 2:length(list)-1
		if list[i]!=list[i-1]&&list[i]!=list[i+1]
			push!(ind,i)
		end
	end

	if list[end]!=list[end-1]
		push!(ind,length(list))
	end
	return ind
end
#This function cleans and orders the Fasta file

function OrderAndCut(Xtot::Alignment,cutoff::Int64)
 	@extract Xtot SpecName N q split Z Sequence UniprotId SpecId Header
 	cand=Tally(SpecName)
 	ncand=[a[1] for a in filter(x->x[2]<=cutoff,cand)]
 	kept=Int64[]
 	#you may add more criteria for filtering here...
 	#Orders so that families have contiguous sequences
 	for fam in ncand
 		kept=[kept;find(x->x==fam,SpecName)]
 	end
 	return Alignment(N,length(kept),q,split,Z[kept,:],Sequence[kept],Header[kept],SpecName[kept],SpecId[kept],UniprotId[kept])
end

#This function prepares both fasta alignment by only keeping species in commom. 
#It also orders the families in the same way for both, giving them similar Id.

function HarmonizeFasta(X1,X2)
	@extract X1 SpecName1=SpecName N1=N q1=q Z1=Z Sequence1=Sequence UniprotId1=UniprotId SpecId1=SpecId Header1=Header
	@extract X2 SpecName2=SpecName N2=N q2=q Z2=Z Sequence2=Sequence UniprotId2=UniprotId SpecId2=SpecId Header2=Header

	kept=unique(intersect(SpecName1,SpecName2))
	ind1=Int64[]
	ind2=Int64[]

	for fam in kept
		ind1=[ind1;find(x->x==fam,SpecName1)]
		ind2=[ind2;find(x->x==fam,SpecName2)]
	end

	#we perform the re labeling of SpecId based over a new, ordered labeling
	spec1=SpecId1[ind1]
	spec2=SpecId2[ind2]
	specu=collect(1:length(unique(spec1)))
	
	change=spec1[1]
	spec1[1]=specu[1]
	pos=1
	#relabeling spec1

	for i in 2:length(spec1)
		if spec1[i]!=change
			change=spec1[i]
			pos+=1
		end
		spec1[i]=specu[pos]
	end
	spec1[end]=specu[end]
	#relabeling spec2
	change=spec2[1]
	spec2[1]=specu[1]
	pos=1

	for i in 2:length(spec2)
		if spec2[i]!=change
			change=spec2[i]
			pos+=1
		end
		spec2[i]=specu[pos]
	end
	spec2[end]=specu[end]

	al1=Alignment(N1,length(ind1),q1,0,Z1[ind1,:],Sequence1[ind1],Header1[ind1],SpecName1[ind1],spec1,UniprotId1[ind1])
	al2=Alignment(N2,length(ind2),q2,0,Z2[ind2,:],Sequence2[ind2],Header2[ind2],SpecName2[ind2],spec2,UniprotId2[ind2])
	return al1,al2
end

#This function rewrites an new fasta file from two and a Matching array

#Works only with harmonized Alignments
#Returns the initial matching between single species

function StartMatching(X1,X2)

	@extract X1 SpecId1=SpecId
	@extract X2 SpecId2=SpecId
	match=zeros(SpecId1)

	ind1=IndexOfUnique(SpecId1)
	ind2=IndexOfUnique(SpecId2)

	candi=intersect(SpecId1[ind1],SpecId2[ind2])

	for el in candi
		a1=findfirst(SpecId1,el)
		a2=findfirst(SpecId2,el)
		match[a1]=a2
	end
	return match
end

#initializes the problem by preparing the fast, constructing frequency matrix
#and correlation matrice, inverting them and returning them
function Initialize(Xi1,Xi2,cut)
	println("initializing the matching, removing families larger than $(cut)... <3")
	X1,X2=HarmonizeFasta(OrderAndCut(Xi1,cut),OrderAndCut(Xi2,cut))
	match=StartMatching(X1,X2)
	freq=nullF(X1,X2)
	corr=nullC(freq)
	#first compute corr from single matched families
	single=X1.SpecId[find(match)]
	UnitFC!(X1,X2,match,single,freq)
	FullCOD!(corr,freq)
	invC=inverseWithPseudo!(corr,freq,0.8)
	println("Setup computed \\o/")
	return X1,X2,match,freq,corr,invC
end

function InitializeWithCorrect(Xi1,Xi2,cut)
	println("initializing the matching, removing families larger than $(cut)... <3")
	X1,X2=HarmonizeFasta(OrderAndCut(Xi1,cut),OrderAndCut(Xi2,cut))
	match=ComputeTrueMatch(X1,X2)
	freq=nullF(X1,X2)
	corr=nullC(freq)
	#first compute corr from single matched families
	correct=unique(X1.SpecId)
	UnitFC!(X1,X2,match,correct,freq)
	FullCOD!(corr,freq)
	invC=inverseWithPseudo!(corr,freq,0.8)
	println("Setup computed \\o/")
	return X1,X2,match,freq,corr,invC
end


function InitializeWithStart(Xi1,Xi2,cut,startmatch)
	println("initializing the matching, removing families larger than $(cut)... <3")
	X1,X2=HarmonizeFasta(OrderAndCut(Xi1,cut),OrderAndCut(Xi2,cut))
	match=startmatch
	freq=nullF(X1,X2)
	corr=nullC(freq)
	#first compute corr from single matched families
	correct=unique(X1.SpecId[find(match)])
	UnitFC!(X1,X2,match,correct,freq)
	FullCOD!(corr,freq)
	invC=inverseWithPseudo!(corr,freq,0.8)
	println("Setup computed \\o/")
	return X1,X2,match,freq,corr,invC
end

#This function rewrites an new fasta file from two and a Matching array
function ComputeTrueMatch(X1,X2)
	return [findfirst(X2.Header,el) for el in X1.Header]
end

function ParCorr(X1,X2, freq::FreqC, invertC, specl::Array{Int,1},strat)
	res=[giveCorrection(X1,X2, freq, invertC, i,strat) for i in specl]
	return res
end

#Works for harmonized Fasta 
#options for the matching are :
# "covariation" or "genetic" or "random" or "greedy"

function RunMatching(a,batch,strat::AbstractString)
	X1,X2,match,freq,corr,invC=a
	spec=Entropy(X1,X2)
	len=length(spec)
	batchl=[spec[i*batch+1:min((i+1)*batch,len)] for i in 0:div(len,batch)]
	savematch=Tuple{Array{Int64,1},Array{Int64,1}}[]
	for el in batchl
		println("optimizing matching for batch :",el)
		res=ParCorr(X1,X2,freq,invC,el,strat)
		println("batch of species")
		println(el)
		println(res)
		ApplyMatching!(X1,X2,match,el,res)


		if strat=="covariation"||strat=="greedy"
			println("updating matrices for batch  :",el)
			# println(freq.specs)
			UnitFC!(X1,X2,match,el,freq)
			FullCOD!(corr,freq)
			invC=inverseWithPseudo!(corr,freq,0.8)
		end
		push!(savematch,(deepcopy(match),deepcopy(freq.specs)))
	end
	return X1,X2,match,freq,corr,invC,savematch
end

#Recomputes several times the model for the species already included
function RunRecompMatching(a,batch,passes)
	X1,X2,match,freq,corr,invC=a
	spec=freq.specs
	len=length(spec)

	#we randomize the order of the species
	for pass in 1:passes
		spec=spec[randperm(len)]
		batchl=[spec[i*batch+1:min((i+1)*batch,len)] for i in 0:div(len,batch)]

		for el in batchl
			el!=[]||continue
			println("optimizing matching for batch :",el)
			res=ParCorr(X1,X2,freq,invC,el,"covariation")
			println("batch of species")
			println(el)
			println(res)
			ApplyMatching!(X1,X2,match,el,res)

			
			println("updating matrices for batch  :",el)
			# println(freq.specs)
			UnitFC!(X1,X2,match,el,freq)
			FullCOD!(corr,freq)
			invC=inverseWithPseudo!(corr,freq,0.8)
		
		end
	end
	return X1,X2,match,freq,corr,invC
end

function RunMatchingWithRecomp(a,batch,strat::AbstractString)
	X1,X2,match,freq,corr,invC=a
	spec=Entropy(X1,X2)
	len=length(spec)
	batchl=[spec[i*batch+1:min((i+1)*batch,len)] for i in 0:div(len,batch)]
	savematch=Tuple{Array{Int64,1},Array{Int64,1}}[]

	for el in batchl
		println("ADDING SPECIES--------------------------------------------")
		println(freq.specs)
		println("optimizing matching for batch :",el)
		res=ParCorr(X1,X2,freq,invC,el,strat)
		println("batch of species")
		println(el)
		println(res)
		ApplyMatching!(X1,X2,match,el,res)


		if strat=="covariation"||strat=="greedy"
			println("updating matrices for batch  :",el)
			# println(freq.specs)
			UnitFC!(X1,X2,match,el,freq)
			FullCOD!(corr,freq)
			invC=inverseWithPseudo!(corr,freq,0.8)
		end
		#recomputation of all the edges
		println("RECOMPUTE ALL THE EDGES WITH 3 PASSES-----------------------------")
		println(freq.specs)
		X1,X2,match,freq,corr,invC=RunRecompMatching((X1,X2,match,freq,corr,invC),batch,3)
		push!(savematch,(deepcopy(match),deepcopy(freq.specs)))
	end
	return X1,X2,match,freq,corr,invC,savematch
end

#we add a random component in that

function ReRunMatching(a,batch,strat::AbstractString)
	X1,X2,match,freq,corr,invC=a
	spec=Entropy(X1,X2)
	len=length(spec)
	#randomization
	batchl=[spec[i*batch+1:min((i+1)*batch,len)] for i in 0:div(len,batch)]
	lenbatch=length(batchl)
	batchl=batchl[randperm(lenbatch)]

	savematch=Tuple{Array{Int64,1},Array{Int64,1}}[]
	for el in batchl
		println("optimizing matching for batch :",el)
		res=ParCorr(X1,X2,freq,invC,el,strat)
		println("batch of species")
		println(el)
		println(res)
		ApplyMatching!(X1,X2,match,el,res)


		if strat=="covariation"||strat=="greedy"
			println("updating matrices for batch  :",el)
			# println(freq.specs)
			UnitFC!(X1,X2,match,el,freq)
			FullCOD!(corr,freq)
			invC=inverseWithPseudo!(corr,freq,0.8)
		end
		push!(savematch,(deepcopy(match),deepcopy(freq.specs)))
	end
	save("HarmX1.jld","X1",X1)
	save("HarmX2.jld","X2",X2)
	save("matchingX1X2.jld","match",match)
	return X1,X2,match,freq,corr,invC,savematch
end


function IterativeMatching(a,batch,stratm,iter)
	res=RunMatching(a,batch,"covariation")
	savematch=typeof(res[7])[res[7]]
	for i in 1:iter
		res=ReRunMatching(res[1:6],batch,stratm)
		push!(savematch,res[7])
	end
	return res,savematch
end

#given a list of specs with their respective matching, updates the match vector
function ApplyMatching!(X1,X2,match,lspec,lmatch)
	@extract X1 SpecId1=SpecId
	@extract X2 SpecId2=SpecId
	length(lspec)==length(lmatch)||error("data non compatible")

	for (i,el) in enumerate(lspec)
	    ind1=find(SpecId1.==el)
	    ind2=find(SpecId2.==el)
	    println("matching edges of X1 :",el)
	    println(ind1[lmatch[i][1]])
	    println("matching edges of X2 :",el)
	    println(ind2[lmatch[i][2]])
	    match[ind1[lmatch[i][1]]]=ind2[lmatch[i][2]]
	end
	return nothing
end

#Works for harmonized Fasta 
function Entropy(X1,X2)
	@extract X1 SpecId1=SpecId
	@extract X2 SpecId2=SpecId
	bib1=Tally(SpecId1)
	bib2=Tally(SpecId2)

	entropy=Tuple{Int64,Float64}[]

	for i in 1:length(bib1)
		bib1[i][1]==bib2[i][1]||error("non harmonized fasta")

		mini=min(bib1[i][2],bib2[i][2])
		maxi=max(bib1[i][2],bib2[i][2])

		ent=sum([log(i) for i in (maxi-mini+1):maxi])
		push!(entropy,(bib1[i][1],ent))
	end
	return [a[1] for a in filter(x->x[2]!=0,sort(entropy,by=x->x[2]))]
end

function Entropy2(X1,X2)
	@extract X1 SpecId1=SpecId
	@extract X2 SpecId2=SpecId
	bib1=Tally(SpecId1)
	bib2=Tally(SpecId2)

	entropy=Tuple{Int64,Float64}[]

	for i in 1:length(bib1)
		bib1[i][1]==bib2[i][1]||error("non harmonized fasta")

		mini=min(bib1[i][2],bib2[i][2])
		maxi=max(bib1[i][2],bib2[i][2])

		ent=sum([log(i) for i in (maxi-mini+1):maxi])
		push!(entropy,(bib1[i][1],ent))
	end
	return [a[1] for a in sort(entropy,by=x->x[2])],[a[2] for a in sort(entropy,by=x->x[2])]
end

function RewriteFastaMatch(X1,X2,match,name)
	g=open(name,"w")
	for (i,edge) in enumerate(match)
		edge!=0.||continue
		X1.SpecName[i]==X2.SpecName[edge]||error("do you have a well formed match ?")

		head=join([">",X1.UniprotId[i],"with",X2.UniprotId[edge],"/",X1.SpecName[i]])
		println(g,head)
		print(g,X1.Sequence[i])
		println(g,X2.Sequence[edge])
	end
	close(g)
end

function WriteFasta(X1,name)
	g=open(name,"w")
	for i in 1:(X1.M)
		head=join([">",X1.Header[i]])
		println(g,head)
		println(g,X1.Sequence[i])
	end
	close(g)
end


function OverlapFasta(name1,name2;opt::AbstractString="no-out")
	X1=readdata(name1)
	X2=readdata(name2)
	pool=unique(X1.SpecName)
	len1=X1.M
	len2=X2.M
	over=0
	if opt=="out"
		str=string("OLP",name1,"_",name2,".fasta")
		g=open(str,"w")
	end
	for (i,j) in enumerate(X1.Header)
		in(j,X2.Header)||continue
		if opt=="out"
			println(g,string(">",j)
				)
			println(g,X1.Sequence[i])
		end
		over+=1
	end
	if opt=="out"
		close(g)
	end
	return over,len1,len2
end

function UnionFasta(name1,name2;opt::AbstractString="no-out")
	X1=readdata(name1)
	X2=readdata(name2)
	pool=unique(X1.SpecName)
	len1=X1.M
	len2=X2.M
	over=0
	if opt=="out"
		str=string("UNION",name1,"_",name2,".fasta")
		g=open(str,"w")
		for i in 1:len1
			println(g,string(">",X1.Header[i]))
			println(g,X1.Sequence[i])
		end
	end
	for (i,j) in enumerate(X2.Header)
		!in(j,X1.Header)||continue
		if opt=="out"
			println(g,string(">",j))
			println(g,X2.Sequence[i])
		end
		over+=1
	end
	if opt=="out"
		close(g)
	end
	return over,len1,len2
end

#ndata: distance file, plm : inference file

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
			push!(curveinter,trueinter/inter)
		end
	end
	return curveinter,curveintra
end


function CutFasta(name,split)
	f=open(readlines,name)
	g=open("SKshortExact.fasta","w")
	h=open("RRshortExact.fasta","w")

	for i in 1:4:length(f)
		print(g,f[i])
		println(h,strip(f[i]))
		print(g,strip(f[i+1]))
		len=length(f[i+1])
		println(g,strip(f[i+2][1:split-len]))
		print(h,strip(f[i+2][split-len+2:end]))
		println(h,strip(f[i+3]))
	end
	
	close(g)
	close(h)
end

function DepleteFamily(X,fam)
	@extract X N M q Z Sequence Header SpecName UniprotId SpecId
	#finds ind that are not belonging to unique families
	ind1=find(x->!(x in fam),X1.SpecId)
	return Alignment(N,length(ind1),q,0,Z[ind1,:],Sequence[ind1],Header[ind1],SpecName[ind1],SpecId[ind1],UniprotId[ind1])
end

function HistoTrue(X1,X2,match)

	truearr=Tuple{Bool,Int,Int,Float64}[]

	for (i,edge) in enumerate(match)
		edge!=0||continue
		card1=length(find(x->x==X1.SpecId[i],X1.SpecId))
		card2=length(find(x->x==X2.SpecId[edge],X2.SpecId))
		mini=min(card1,card2)
		maxi=max(card1,card2)

		ent=sum([log(i) for i in (maxi-mini+1):maxi])
		push!(truearr,(X1.Header[i]==X2.Header[edge],card1,card2,ent))		
	end
	return truearr
end


function EntropyCurve(X1,X2)
	@extract X1 SpecId1=SpecId
	@extract X2 SpecId2=SpecId
	bib1=Tally(SpecId1)
	bib2=Tally(SpecId2)

	entropy=Tuple{Int64,Float64}[]

	for i in 1:length(bib1)
		bib1[i][1]==bib2[i][1]||error("non harmonized fasta")

		mini=min(bib1[i][2],bib2[i][2])
		maxi=max(bib1[i][2],bib2[i][2])

		ent=sum([log(i) for i in (maxi-mini+1):maxi])
		push!(entropy,(bib1[i][1],ent))
	end
	return entropy
end


function InterDistance(distfile,cut)
	data=readdlm(distfile)
	len=size(data)[1]
	contacts=Tuple{Float64,Float64,Float64}[]
	for i in 1:len
		a,b,c=data[i,:]
		if (a<=cut)&&(b>cut)&&(c<8)
			push!(contacts,(a,b,c))
		end
	end
	return contacts
end


function EdgeStat(fast1,fast2)
	dat1=readdlm(fast1)
	dat2=readdlm(fast2)[1:2:end]
	fastacom=AbstractString[]
	score = 0
	for (i,str) in enumerate(dat1[1:2:end])
		if in(str,dat2)
			score+=1
			push!(fastacom,str)
			push!(fastacom,dat1[i*2])
		end
		println(i)
	end
	return score,fastacom
end

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

#gives restrict a specific number of sites

function FullEdgesScores2(nameX1,nameX2,namematch,restrict)
	println("Recomputing the model and sorting the edges... <3")
	X1=load(nameX1,"X1")
	X2=load(nameX2,"X2")
	match=load(namematch,"match")
	siteskept=Float64[]
	for el in restrict
		append!(siteskept,[(el-1)*20+1:el*20])
	end

	nbcontacts1=length(find(restrict.<88))
	freq=nullF(X1,X2)
	corr=nullC(freq)
	single=X1.SpecId[find(match)]
	UnitFC!(X1,X2,match,single,freq)
	FullCOD!(corr,freq)

	
	invC=Restrict_inverseWithPseudo!(corr,freq,0.8,siteskept)
	newCij=corr.Cij[siteskept,siteskept]
	
	lines=length(match)
	#here we consider the pseudo count contribution as well
	Meff=freq.M[1]+freq.M[2]
	aver=(freq.Pi/Meff)[siteskept]
	flag=1
	ordering=Array{Tuple{Tuple{Float64,Float64},Float64},1}[]
	suborder=Tuple{Tuple{Float64,Float64},Float64}[]
	#creates all the possibles within species match
	pairings=[find(x->x==i,X2.SpecId) for i in X1.SpecId]
	len1=20*nbcontacts1
	nullconst=1/2*(logdet(newCij)-logdet(newCij[1:len1,1:len1])-logdet(newCij[(len1+1):end,(len1+1):end]))

	println(nullconst)
	for i in 1:lines
		for j in pairings[i]
			seq1=expandBinary(X1.Z[i,:],20)[:]
			seq2=expandBinary(X2.Z[j,:],20)[:]

			seq=[seq1;seq2][siteskept]
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

function ShuffleSeq(seq)
	len=length(seq)
	newseq=zeros(len)
	for i in 0:len/20-1
		randind=collect(i*20+1:(i+1)*20)[randperm(20)]
		newseq[i*20+1:(i+1)*20]=seq[randind]
	end
	return newseq
end

function matrixComp(nameX1,nameX2,namematch)
	println("Recomputing the model and sorting the edges... <3")
	X1=load(nameX1,"X1")
	X2=load(nameX2,"X2")
	match=load(namematch,"match")

	freq=nullF(X1,X2)
	corr=nullC(freq)
	single=X1.SpecId[find(match)]
	UnitFC!(X1,X2,match,single,freq)
	FullCOD!(corr,freq)
	invC=inverseWithPseudo!(corr,freq,0.8)

	freq2=nullF(X1,X2)
	corr2=nullC(freq)
	single2=X1.SpecId[find(sort(match))]
	UnitFC!(X1,X2,sort(match),single2,freq2)
	FullCOD!(corr2,freq2)
	invC2=inverseWithPseudo!(corr2,freq2,0.8)
	return corr,corr2,invC,invC2
end


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


function CountPercent(listoflist)
	
	len=length(listoflist)
	resu=zeros(len,3)
	count=1
	for (el1,el2) in listoflist
		ind=find((el1.!=0).*(el1.==collect(1:8998)))
		nbtrue=length(ind)
		resu[count,1]=nbtrue
		nbposs=length(find((el1.!=0)))
		resu[count,2]=nbposs
		percent=nbtrue/nbposs
		resu[count,3]=percent
		count+=1
	end
	return resu
end

function ROCedges(matching,scores,X1,X2)
	perm_scores=Tuple{Tuple{Float64,Float64},Float64}[]
	len=length(scores)

	for (i,j) in enumerate(matching)
		
		spec=(X1.SpecId)[i]
		spec<=len||break
		set_edges=scores[spec]
		ind=find(x->(x[1][1]==i)&&(x[1][2]==j),set_edges)
		push!(perm_scores,set_edges[ind][1])
	end

	sort!(perm_scores,by=x->x[2],rev=true)
	count_true=0
	count_tot=0
	roc=Float64[]

	for el in perm_scores
		count_tot+=1
		if el[1][1]==el[1][2]||count_tot==185
			count_true+=1 
		end
		push!(roc,count_true/count_tot)
	end
	return perm_scores,roc
end

function ROCedges2(matching,scores)
	perm_scores=Tuple{Tuple{Float64,Float64},Float64}[]
	len=length(scores)
	s=[scores...][:]

	L=map(x->x[1][1]==x[1][2],s)
	S=[x[2] for x in s]

	p=sortperm(S,rev=true)
	L=L[p]
	S=S[p]


	sl=sum(L)
	nl=sum(L.==false)
	TPRN=cumsum(L)/sl
	FPRN=cumsum(!L)/nl

	return TPRN,FPRN
end

function trapz{Tx<:Number, Ty<:Number}(x::Vector{Tx}, y::Vector{Ty})
    # Trapezoidal integration rule
    local n = length(x)
    if (length(y) != n)
        error("Vectors 'x', 'y' must be of same length")
    end
    r = zero(zero(Tx) + zero(Ty))
    if n == 1; return r; end
    for i in 2:n
        r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    #= correction -h^2/12 * (f'(b) - f'(a))
    ha = x[2] - x[1]
    he = x[end] - x[end-1]
    ra = (y[2] - y[1]) / ha
    re = (y[end] - y[end-1]) / he
    r/2 - ha*he/12 * (re - ra)
    =#
    return r/2
end


function DistanceMatrx(listofedge)
	len1=length(unique([a[1][1] for a in listofedge]))
	len2=length(unique([a[1][2] for a in listofedge]))
	base1=minimum(unique([a[1][1] for a in listofedge]))
	base2=minimum(unique([a[1][2] for a in listofedge]))

	mat=zeros(len1,len2)

	for el in listofedge
		mat[el[1][1]+1-base1,el[1][2]+1-base2]=1/el[2]
	end
	labels1=["prot1_p$(i)" for i in 1:len1]
	labels2=["prot2_p$(i)" for i in 1:len2]'
	labels=["labels";labels1]
	final_mat=hcat(labels,vcat(labels2,mat))
	return final_mat
end


function cycleDecomp(permlist)
    l=length(permlist)
    totcycle=Array{Int64,1}[]
    seed=1
    res=copy(permlist)
    
    while true
        cycle=[seed]
        accu=permlist[seed]

            while accu!=seed
                push!(cycle,accu)
                accu=permlist[accu]
            end
        push!(totcycle,cycle)

        res=sort(setdiff(res,cycle))      
        
        length(res)!=0||break
        
        seed=res[1]
        
    end
    return totcycle
end
