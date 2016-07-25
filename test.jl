include("src/ParaMatch.jl")

X1 = ParaMatch.read_fasta_alignment("test_set/SKshortExact.fasta")
X2 = ParaMatch.read_fasta_alignment("test_set/RRshortExact.fasta")

X1new, X2new = ParaMatch.prepare_alignments(X1, X2, 4)
resu = ParaMatch.run_matching(X1new, X2new, 5, strat="covariation");

