module Datagroup

    using Costs.totalcost2

    export onecoaliform, coaliform, datagroupCF, datagroupES

    function onecoaliform(datasize, codesize, 
                          popul, transmat, 
                          pararetr, parasurv,
                          paccess, psurv,
                          ingroup, init)
        N = length(datasize)
        group = [init]
        ingroup[init] = 1
        cost =  totalcost2(group,
                           datasize, codesize,
                           popul, transmat,
                           pararetr, parasurv,
                           paccess, psurv)
        ind = (init + 1) % N 
        ind = ind == 0 ? N : ind
        while ind != init
            if ingroup[ind] == 1
                ind = (ind + 1) % N 
                ind = ind == 0 ? N : ind
                continue
            end
            join = deepcopy(group)
            push!(join, ind)
            cost_join = totalcost2(join,
                                   datasize, codesize,
                                   popul, transmat,
                                   pararetr, parasurv,
                                   paccess, psurv)
            alone = [ind]
            cost_alone = totalcost2(alone, 
                                    datasize, codesize,
                                    popul, transmat,
                                    pararetr, parasurv,
                                    paccess, psurv)
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


    function coaliform(datasize, codesize, 
                       popul, transmat, 
                       pararetr, parasurv,
                       paccess, psurv,
                       init)
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
            group, cost = onecoaliform(datasize, codesize, 
                                       popul, transmat, 
                                       pararetr, parasurv,
                                       paccess, psurv,
                                       ingroup, ind)
            push!(groups, group)
            push!(costs, cost)  
            ind = (ind + 1) % N 
            ind = ind == 0 ? N : ind
        end
        return groups, costs
    end

    function datagroupCF(datasize, codesize, 
                         popul, transmat, 
                         pararetr, parasurv,
                         paccess, psurv)
        groupsfinal = nothing
        costsfinal = nothing
        costsum = Inf 
        N = length(datasize)
        for i = 1:N
            groups, costs = coaliform(datasize, codesize, 
                                      popul, transmat, 
                                      pararetr, parasurv,
                                      paccess, psurv,
                                      i)
            if sum(costs) < costsum
                costsum = sum(costs)
                groupsfinal = groups
                costsfinal = costs
            end
        end
        return groupsfinal, costsfinal
    end

    function datagroupES(datasize, codesize, 
                         popul, transmat, 
                         pararetr, parasurv,
                         paccess, psurv, 
                         allpartitions)
        groupsfinal = nothing
        costsfinal = nothing
        costsum = Inf 
        for partitions in allpartitions
            part = partitions
            costs = Float64[]
            for group in part
                temp = totalcost2(group,
                                  datasize, codesize,
                                  popul, transmat,
                                  pararetr, parasurv,
                                  paccess, psurv)
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