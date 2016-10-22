module ParalogMatchingTests

using ParalogMatching
using FastaIO
using Base.Test
using Compat: String

const testfile1 = "X1.fasta.gz"
const testfile2 = "X2.fasta.gz"
const outfile_ok = "Xmatched.fasta.gz"

function test()
    outfile_tst = tempname()
    outfile_harmonized = tempname()
    try
        X12, match = paralog_matching(testfile1, testfile2, outfile_tst, cutoff=4, pseudo_count=0.8, batch=5)
        X = ParalogMatching.Alignment(X12, match)
        ParalogMatching.write_fasta(X, outfile_harmonized)
        M1 = readfasta(outfile_tst)
        M2 = readfasta(outfile_ok)
        Mharmo = readfasta(outfile_harmonized)
        @test M1 == M2
        @test M1 == Mharmo
    finally
        isfile(outfile_tst) && rm(outfile_tst)
        isfile(outfile_harmonized) && rm(outfile_harmonized)
    end
end

function regex_tests()
    header = String["SOMELOCATION/SOMESPECIES/OTHERINFO"]
    answer = [1], ["SOMESPECIES"], ["SOMELOCATION"]

    @test ParalogMatching.compute_spec(header, r"^([^/]+)/([^/]+)/") == answer
    @test ParalogMatching.compute_spec(header, r"^([^/]+)/([^/]+)/(.*)") == answer
    @test ParalogMatching.compute_spec(header, r"^(?<id>[^/]+)/(?<species>[^/]+)/") == answer
    @test ParalogMatching.compute_spec(header, r"^(?<id>[^/]+)(/)(?<species>[^/]+)/") == answer
    @test ParalogMatching.compute_spec(header, r"^([^/]+)/([^/]+)/(?<rest>.*)") == answer
    @test ParalogMatching.compute_spec(header, r"^(?<id>[^/]+)(/)(?<species>[^/]+)/(?<rest>.*)") == answer

    @test_throws ErrorException ParalogMatching.compute_spec(header, r"^([^/]+)/[^/]+/")
    @test_throws ErrorException ParalogMatching.compute_spec(header, r"^(?<id>[^/]+)/([^/]+)/")
    @test_throws ErrorException ParalogMatching.compute_spec(header, r"^([^/]+)/(?<species>[^/]+)/")

    header = String["SOMESPECIES/OTHERINFO/SOMELOCATION"]
    @test ParalogMatching.compute_spec(header, r"^(?<species>[^/]+)/([^/]+)/(?<id>.+)") == answer

    header = String["SOMELOCATION/SOMESPECIES/OTHERINFO", "OTHERLOC/DOESNTMATCH"]
    @test_throws ErrorException ParalogMatching.compute_spec(header, r"^([^/]+)/([^/]+)/")
    @test_throws ErrorException ParalogMatching.compute_spec(header, r"^(?<id>[^/]+)/(?<species>[^/]+)/")
end

test()
regex_tests()

end # module
