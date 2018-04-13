push!(LOAD_PATH, "./jlcodebase")
using Utilities
using Costs
using Datagroup

function gettransmat(filepath)
    mat = readdlm(filepath, ',')
    (n, m) = size(mat)
    transmat = reshape(mat, (n, n, Int(m/n)))
    return transmat
end

function getdatasize(filepath)
    mat = readdlm(filepath, ',')
    return mat
end

function getpopul(filepath)
    mat = readdlm(filepath, ',')
    return mat
end

function getrankinfo(filepath)
    mat = readdlm(filepath, ',')
    return mat
end



