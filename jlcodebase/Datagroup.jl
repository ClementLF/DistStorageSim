module Datagroup

	using Costs.totalcost

	export onecoaliform, coaliform, datagroupCF, datagroupES

	function onecoaliform(datasize, popul, transmat, 
	                      pararetr, parasurv,
	                      cratio, pfail, ingroup, init)
	    N = length(datasize)
	    group = [init]
	    ingroup[init] = 1
	    cost =  totalcost(group, datasize, popul, transmat, 
	                      pararetr, parasurv, pfail, cratio)
	    ind = (init + 1) % N 
	    ind = ind == 0 ? N : ind
	    while ind != init
	        if ingroup[ind] == 1
	            ind = (ind + 1) % N 
	            ind = ind == 0 ? N : ind
	            continue
	        end
	        join = copy(group)
	        push!(join, ind)
	        cost_join = totalcost(join, datasize, popul, transmat, 
	                              pararetr, parasurv, pfail, cratio)
	        alone = [ind]
	        cost_alone = totalcost(alone, datasize, popul, transmat, 
	                               pararetr, parasurv, pfail, cratio)
	        if cost_join < cost + cost_alone
	            cost = cost_join
	            group = join
	            ingroup[ind] = 1
	        end
	        ind = (ind + 1) % N 
	        ind = ind == 0 ? N : ind
	     end
	    return group, cost
	end


	function coaliform(datasize, popul, transmat, 
	                   pararetr, parasurv, 
	                   cratio, pfail, init)
	    N = length(datasize)
	    ingroup = zeros(Int, N)
	    groups = Array{Int, 1}[]
	    costs = Float64[]
	    ind = init
	    for i = 1:N
	        if ingroup[ind] == 1
	            ind = (ind + 1) % N 
	            ind = ind == 0 ? N : ind
	            continue
	        end
	        group, cost = onecoaliform(datasize, popul, transmat, 
	                                   pararetr, parasurv,
	                                   cratio, pfail, ingroup, ind)
	        push!(groups, group)
	        push!(costs, cost)  
	        ind = (ind + 1) % N 
	        ind = ind == 0 ? N : ind
	    end
	    return groups, costs
	end

	function datagroupCF(datasize, popul, transmat, 
	                     pararetr, parasurv, cratio, pfail)
	    groupsfinal = nothing
	    costsfinal = nothing
	    costsum = Inf 
	    N = length(datasize)
	    for i = 1:N
	        groups, costs = coaliform(datasize, popul, transmat, 
	                                  pararetr, parasurv, cratio, pfail, i)
	        if sum(costs) < costsum
	            costsum = sum(costs)
	            groupsfinal = groups
	            costsfinal = costs
	        end
	    end
	    return groupsfinal, costsfinal
	end

	function datagroupES(datasize, popul, transmat, 
	                     pararetr, parasurv, 
	                     cratio, pfail, allpartitions)
	    groupsfinal = nothing
	    costsfinal = nothing
	    costsum = Inf 
	    for partitions in allpartitions
	        part = partitions
	        costs = Float64[]
	        for group in part
	            temp = totalcost(group, datasize, popul, transmat11[:, :, 1], 
	                             pararetr, parasurv, pfail, cratio)
	            push!(costs, temp)
	        end
	        if sum(costs) < costsum
	            costsum = sum(costs)
	            groupsfinal = part
	            costsfinal = costs
	        end
	    end
	    return groupsfinal, costsfinal
	end
end