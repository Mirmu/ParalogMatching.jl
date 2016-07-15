include("src/ParaMatch.jl")

X1 = ParaMatch.read_fasta_alignment("test_set/SKshortExact.fasta")
X2 = ParaMatch.read_fasta_alignment("test_set/RRshortExact.fasta")

init = ParaMatch.initialize_matching(X1, X2, 4)
resu = ParaMatch.run_matching(init, 5, "covariation");

