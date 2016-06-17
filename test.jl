include("ParaMatch.jl")

X1=ParaMatch.readdata("test_set/SKshortExact.fasta")
X2=ParaMatch.readdata("test_set/RRshortExact.fasta")

init=ParaMatch.initialize(X1,X2,4)
resu=ParaMatch.run_matching(init,5,"covariation");

