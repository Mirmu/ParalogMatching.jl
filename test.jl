include("ParaMatch.jl")

X1=ParaMatch.readdata("test_set/SKshortExact.fasta")
X2=ParaMatch.readdata("test_set/RRshortExact.fasta")

init=ParaMatch.Initialize(X1,X2,4)
resu=ParaMatch.RunMatching(init,5,"covariation");

