push!(LOAD_PATH, ".\\jlcodebase")
using Utilities.zipfpdf

function gettransmat(filepath)
    mat = readdlm(filepath, ',')
    (n, m) = size(mat)
    transmat = reshape(mat, (n, n, Int(m/n)))
    return transmat
end

function genpopul(n, g)
    freq = zipfpdf(n, g)
    randrank = randperm(n)
    freq = freq[randrank]
    return freq, randrank
end

function gendatasize(n, minval, maxval)
    rand(minval:maxval, n)
end

# setup
n = 10
minval = 3
maxval = 12

datasize = gendatasize(N, minsize, maxsize)
popul, rankinfo = genpopul(N, para_popular)

datasize = gendatasize(N, minsize, maxsize)
popul, rankinfo = genpopul(N, para_popular)

path = joinpath(outputdir, "datasize.csv");
writedlm(path, datasize, ",")

path = joinpath(outputdir, "popul.csv");
writedlm(path, popul, ",")

path = joinpath(outputdir, "rankinfo.csv");
writedlm(path, rankinfo, ",")
