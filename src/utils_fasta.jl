#very useful function that classifies elements along with their counts in a given set

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

#Returns the indexes of the list that have unique labels
#WARNING !: works only for sorted lists 

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


