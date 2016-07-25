include("src/ParalogMatching.jl")

using Gurobi

X1 = ParaMatch.read_fasta_alignment("test_set/SKshortExact.fasta")
X2 = ParaMatch.read_fasta_alignment("test_set/RRshortExact.fasta")

X12 = ParaMatch.prepare_alignments(X1, X2, cutoff=4)
match = ParaMatch.run_matching(X12, batch=5, strategy="covariation",
                               lpsolver = GurobiSolver(OutputFlag=false));

