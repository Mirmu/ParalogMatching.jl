# Very useful function that classifies elements along with their counts in a given set
function tally(list)
    slist = sort(list)
    counter = 0
    comp = slist[1]
    tallist = Tuple{eltype(list),Int64}[]

    for el in slist
	if el == comp
	    counter += 1
	else
	    push!(tallist, (comp, counter))
	    comp = el
	    counter = 1
	end
    end
    push!(tallist, (comp, counter))
    return tallist
end

# Returns the indexes of the list that have unique labels
# WARNING!: works only for sorted lists
function index_of_unique(list)
    ind = Int64[]
    L = length(list)
    L == 0 && return ind

    chgd = true
    prev = list[1]
    for i = 2:L
	curr = list[i]
	curr ≥ prev || error("list is not sorted")
	if curr ≠ prev
	    chgd && push!(ind, i-1)
	    chgd = true
	else
	    chgd = false
	end
	prev = curr
    end
    chgd && push!(ind, L)
    return ind
end
