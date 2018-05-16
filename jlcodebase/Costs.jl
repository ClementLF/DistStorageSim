module Costs

    using Utilities.binomcdf

    export costretrival2, costsurvival2, totalcost2

    # function costretrival(group, datasize, popul, transmat)
    #     mem = length(group);
    #     groupsize = datasize[group];
    #     groupind = sum(groupsize) - groupsize ;
    #     mean(groupind);
    # end

    function costsurvival2(group, datasize, codesize, psurv)
        k = sum(datasize[group]);
        n = sum(codesize[group]);
        return -log(binomcdf(n, psurv, k));
    end

    # function totalcost(group, datasize, codesize, popul, transmat, 
    #                    pararetr, parasurv, psurv, paccess)
    #     costretr = costretrival(group, datasize, popul, transmat)
    #     costsurv = costsurvival(group, datasize, codesize, psurv)
    #     return pararetr * costretr + parasurv * costsurv
    # end


    function costretrival2(group, datasize, codesize, popul, transmat, paccess)
        d, r = getGroupRetrieval(group, transmat, datasize, codesize, paccess);
        return sum(popul[group] .* (r - d));
    end

    function totalcost2(group,
                        datasize, codesize,
                        popul, transmat,
                        pararetr, parasurv,
                        paccess, psurv)
        costretr = costretrival2(group, datasize, codesize, popul, transmat, paccess)
        costsurv = costsurvival2(group, datasize, codesize, psurv)
        return pararetr * costretr + parasurv * costsurv
    end

    function getRetrieval(kr, K, N, pa, Gp)
        # direct 
        eta_d = pa^kr * kr;
        # indirect
        eta_id = 0;
        for i in K:N
            cnt = binomial(BigInt(N), BigInt(i)) - binomial(BigInt(N-kr), BigInt(i-kr))
            prb = pa^(i) * (1 - pa)^(N-i)
            eta_id += cnt * prb * i
        end
        # penalty
        prb = binomcdf(N, pa, K, "lower")
        eta_pnt = prb * Gp
        return  eta_d + eta_id + eta_pnt
    end

    function getAllPathesU(s, t, V, G, visited, path, pathes)
        visited[s] = 1
        push!(path, s)
        if s == t
            push!(pathes, deepcopy(path))
        else
            for i in 1:V
                if G[s, i] == 0
                    continue
                elseif visited[i] == 1
                    continue
                else
                    getAllPathesU(i, t, V, G, visited, path, pathes)
                end
            end
        end
        pop!(path)
        visited[s] = 0
    end

    function getDemand(s, t, G, datasize)
        (V, _) = size(G)
        path = Array{Int, 1}()
        pathes = Array{Int, 1}[]
        visited = zeros(V, 1)
        getAllPathesU(s, t, V, G, visited, path, pathes)
        pii = zeros(V, 1)
        for path in pathes
            tempPrb = 1
            if length(path) == 2
                continue
            end
            temp = 1
            prev = s
            for i in 2:length(path)-1
                temp *= G[prev, path[i]]
                prev = path[i]
            end
            pii[path[end-1]] += temp
        end
        return sum(pii[1:end-1] .* datasize) + datasize[s]
    end

    function getSubGraph(group, G)
        Gsub = G[group, group]
        lastC =  1 -  sum(Gsub, 2)
        lastR = zeros(1, length(group) + 1)
        lastR[end] = 1
        return vcat(hcat(Gsub, lastC), lastR)
    end
    
    function getGroupRetrieval(group, G, datasize, codesize, pa)
        subgraph = getSubGraph(group, G)
        sink = size(subgraph)[1]
        K = sum(datasize[group])
        N = sum(codesize[group])
        demands = zeros(1, length(group))
        retrievals = zeros(1, length(group))
        for i in 1:length(group)
            temp = getDemand(i, sink, subgraph, datasize[group])
            demands[i] = temp
            kr = datasize[group[i]]
            retrievals[i] = getRetrieval(kr, K, N, pa, N)
        end
        return demands, retrievals
    end

end