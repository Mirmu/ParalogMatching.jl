module ParalogMatchingTests

using ParalogMatching
using FastaIO
using Base.Test

const testfile1 = "X1.fasta.gz"
const testfile2 = "X2.fasta.gz"
const outfile_ok = "Xmatched.fasta.gz"

function test()
    outfile_tst = tempname()
    try
        paralog_matching(testfile1, testfile2, outfile_tst, cutoff=4, pseudo_count=0.8, batch=5)

        M1 = readfasta(outfile_tst)
        M2 = readfasta(outfile_ok)
        @test M1 == M2
    finally
        rm(outfile_tst)
    end
end

test()

end # module
