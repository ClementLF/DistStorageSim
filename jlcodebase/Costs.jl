module Costs

    using Utilities.binomcdf

	export costretrival, costsurvival, totalcost

	function costretrival(group, datasize, popul, transmat)
    mem = length(group)
    groupsize = datasize[group]
    groupind = sum(groupsize) ./ groupsize 
    sum(groupind)
	end

	function costsurvival(group, datasize, popul, transmat, cratio, pfail)
	    groupsize = datasize[group]
	    k = sum(groupsize)
	    n = ceil(k * cratio)
	    p = 1 - pfail
	    return -log(binomcdf(n, p, k))
	end

	function totalcost(group, datasize, popul, transmat, 
	                   pararetr, parasurv, pfail, cratio)
	    costretr = costretrival(group, datasize, popul, transmat)
	    costsurv = costsurvival(group, datasize, popul, transmat, cratio, pfail)
	    return pararetr * costretr + parasurv * costsurv
	end
end