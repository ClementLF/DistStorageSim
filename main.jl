push!(LOAD_PATH, "./jlcodebase")

using Utilities
using Costs
using Datagroup
using Combinatorics

function gettransmat(filepath)
    mat = readdlm(filepath, ',')
    (n, m) = size(mat)
    transmat = reshape(mat, (n, n, Int(m/n)))
    return transmat
end;

function getdatasize(filepath)
    mat = readdlm(filepath, ',')
    return mat
end;

function evalJointP(groups, datasize, codesize, p)
    P = 1
    for group in groups
        P *= binomcdf(sum(codesize[group]), p, sum(datasize[group]))
    end
    return P
end

function evalRetrieval(groups, datasize, codesize,
                       popul, transmat, paccess)
    oh = 0
    for group in groups
        oh += costretrival2(group, datasize, codesize, popul, transmat, paccess)
    end
    return oh                   
end

transmat11 = gettransmat("transmat11.csv");
datasize = getdatasize("datasize1.csv");
NN = length(datasize)
allpartitions = ncpartitions(NN)
lambda_pop = 0.5;


# different similarity cases
G1 = transmat11[:, :, 1];
G2 = transmat11[:, :, 2];
G3 = transmat11[:, :, 3];
G4 = transmat11[:, :, 4];
G5 = transmat11[:, :, 5];

# different popularity cases
popul1 = zipfpdf(NN, lambda_pop);
popul2 = ones(1, NN) / NN;

# system setup
psurv = 0.51;
paccess = 0.5;
re_ratio = 2;
codesize =  ceil.(datasize * 2);

# utility setup
pararetr = 0.05;
parasurv = 10;

psurv_r = 0.6:0.02:1;
pacc_r = [0.5, 0.6, 0.7];

# CF_JointP = zeros(length(pacc_r), length(psurv_r));
# CF_OH = zeros(length(pacc_r), length(psurv_r));
# CF_TotalU = zeros(length(pacc_r), length(psurv_r));

ES_JointP = zeros(length(pacc_r), length(psurv_r));
ES_OH = zeros(length(pacc_r), length(psurv_r));
ES_TotalU = zeros(length(pacc_r), length(psurv_r));

for (inda, pa) in enumerate(pacc_r)
    for (inds, ps) in enumerate(psurv_r)
        # CF_groups, CF_costs = datagroupCF(datasize, codesize, 
        #                                   popul1, G5, 
        #                                   pararetr, parasurv,
        #                                   pa, ps)
        # CF_TotalU[inda, inds] = sum(CF_costs)
        # CF_JointP[inda, inds] = evalJointP(CF_groups, datasize, codesize, ps)
        # CF_OH[inda, inds] = evalRetrieval(CF_groups, datasize, codesize,
        #                       popul2, G1, pa)


        ES_groups, ES_costs = datagroupES(datasize, codesize, 
                                          popul2, G1, 
                                          pararetr, parasurv,
                                          pa, ps,
                                          allpartitions)

        ES_TotalU[inda, inds] = sum(ES_costs)
        ES_JointP[inda, inds] = evalJointP(ES_groups, datasize, codesize, ps)
        ES_OH[inda, inds] = evalRetrieval(ES_groups, datasize, codesize,
                              popul2, G1, pa)
    end
end
# writedlm("datasize1_G5_pop0.5_CF_JointP.csv", CF_JointP)
# writedlm("datasize1_G5_pop0.5_CF_OH.csv", CF_OH)
# writedlm("datasize1_G5_pop0.5_CF_TotalU.csv", CF_TotalU)

writedlm("datasize1_G1_pop0_ES_JointP.csv", ES_JointP)
writedlm("datasize1_G1_pop0_ES_OH.csv", ES_OH)
writedlm("datasize1_G1_pop0_ES_TotalU.csv", ES_TotalU)

