include("src/ParaMatch.jl")

X1 = ParaMatch.read_fasta_alignment("test_set/SKshortExact.fasta")
X2 = ParaMatch.read_fasta_alignment("test_set/RRshortExact.fasta")

X12 = ParaMatch.prepare_alignments(X1, X2, 4)
match = ParaMatch.run_matching(X12, 5, strat="covariation");

