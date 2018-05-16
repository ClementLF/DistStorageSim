module Utilities

    export binomcdf, zipfpdf

    function binomcdf(n, p, k, cdftype="upper")
        cdf = 0
        for i = k:n
            cdf += binomial(BigInt(n), BigInt(i)) * p^i * (1-p)^(n-i)
        end
        if cdftype == "lower"
            cdf = 1 - cdf
        end
        return cdf
    end

    function zipfpdf(n, g)
        rankinfo = collect(1:n)
        numera  = 1 ./ rankinfo.^g
        denomi = sum(numera)
        freq = numera / denomi
        return freq
    end

end
