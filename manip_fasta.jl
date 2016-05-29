######################BLOCK FOR MANIPULATION OF THE ALIGNMENTS###################
#Orders and removes families that are too large

function OrderAndCut(Xtot::Alignment,cutoff::Int64)
 	@extract Xtot SpecName N q split Z Sequence UniprotId SpecId Header
 	cand=Tally(SpecName)
 	ncand=[a[1] for a in filter(x->x[2]<=cutoff,cand)]
 	kept=Int64[]

 	#you may add more criteria for filtering here...
 	#Orders so that families have contiguous sequences organized in blocks
 	for fam in ncand
 		kept=[kept;find(x->x==fam,SpecName)]
 	end

 	return Alignment(N,length(kept),q,split,Z[kept,:],Sequence[kept],Header[kept],SpecName[kept],SpecId[kept],UniprotId[kept])

end

#Harmonizing consists of creating the following structure for a FASTA:
#Filter the species in common between both FASTA and remove the others, that cannot be matched
#Giving those species a common labeling
#This routine is a bit inefficient

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
	#the new labeling is given by kept
	for i in ind1
		SpecId1[i]=findfirst(kept,SpecName1[i])
	end
	
	for i in ind2
		SpecId2[i]=findfirst(kept,SpecName2[i])
	end


	al1=Alignment(N1,length(ind1),q1,0,Z1[ind1,:],Sequence1[ind1],Header1[ind1],SpecName1[ind1],SpecId1[ind1],UniprotId1[ind1])
	al2=Alignment(N2,length(ind2),q2,0,Z2[ind2,:],Sequence2[ind2],Header2[ind2],SpecName2[ind2],SpecId2[ind2],UniprotId2[ind2])
	return al1,al2
end


#####################BLOCK FOR WRITING FASTA FROM ALIGNMENTS OBJECTS########
#Writes a given Alignments under "name"

function WriteFasta(X1,name)
	g=open(name,"w")
	for i in 1:(X1.M)
		head=join([">",X1.Header[i]])
		println(g,head)
		println(g,X1.Sequence[i])
	end
	close(g)
end


#Rewrites the output as a FASTA for two given Alignments and a given matching, under the "name"
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

######################BLOCK FOR COMPARING TWO FASTA###########################
#Computes the Fasta that is the intersection of two FASTA of names "names1" and "names2"
#If no output, it simply returns the overlap value, with both lengths of the FASTA

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

#Computes the FASTA that is the Union of the two FASTA of names "names1" and "names2"
#Same return than Overlap

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

