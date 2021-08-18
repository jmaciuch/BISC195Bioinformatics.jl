using BISC195Bioinformatics
using Test

@testset "BISC195Bioinformatics" begin
    
@testset "Using Strings" begin
    
    @testset "normalizeDNA" begin
        @test normalizeDNA("aatgn") == "AATGN"
        @test_throws Exception normalizeDNA("ZCA")
        @test_throws Exception normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == "C"
        @test typeof(c) == String
    end # normalizeDNA

    @testset "composition" begin
        seq = rand(['A','T','G','C','N'], 20) |> join
        bc = composition(seq)
        @test bc isa Dict

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)

        bc = composition(lowercase(seq))

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)
    end # composition

    @testset "gc_content" begin
        @test gc_content("ANTG") == 0.25
        @test gc_content("cccggg") * 100 == 100.0
        @test gc_content("ATta") == 0.0
        @test_throws Exception gc_content("ATty")
    end # gc_content

    @testset "complement" begin
        @test complement("ATTAN") == "TAATN"
        @test complement("gcta") == "CGAT"
        @test complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception complement("ABC")
    end # complement

    @testset "reverse_complement" begin
        @test reverse_complement("ATTAN") == "NTAAT"
        @test reverse_complement("gcta") == "TAGC"
        @test reverse_complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception reverse_complement("ABC")
    end # reverse_complement

    @testset "parse_fasta" begin
        testpath = normpath(joinpath(@__DIR__, "..", "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")

        ex1 = parse_fasta(ex1_path)
        @test ex1 isa Tuple
        @test all(x-> x isa AbstractString, ex1[1])
        @test all(x-> x isa String, ex1[2])

        @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
        @test ex1[2] == ["AATTATAGC", "CGCCCCCCAGTCGGATT"]

        cov2 = parse_fasta(genomes)
        @test length(cov2[1]) == 8
        @test length(cov2[2]) == 8
    end #parse_fasta

    @test_set "splice_fasta" begin
        variant_dict = Dict("Alpha" => "Alpha",
                "B.1.1.7" => "Alpha",
                "Beta" => "Beta",
                "B.1.351" => "Beta",
                "B.1.351.2" => "Beta",
                "B.1.351.3" => "Beta",
                "Gamma" => "Gamma",
                "P.1" => "Gamma",
                "P.1.1" => "Gamma",
                "P.1.2" => "Gamma",
                "Delta" => "Delta",
                "B.1.617.2" => "Delta",
                "AY.1" => "Delta",
                "AY.2" => "Delta",
                "AY.3" => "Delta",
                "AY.3.1" => "Delta")

        testpath = normpath(joinpath(@__DIR__, "..", "data"))
        Analysis1_path = joinpath(testpath, "Analysis1_test.fasta")

        Analysis1 = splice_fasta(Analysis1_path, variant_dict, 2)
        @test all(x -> x isa Tuple)

        @test Analysis1[1] == ["Beta", "TTTGCGTTTTTAAAGCGCCCCGATAAGCTAGATCGATCGCGTAGCGCTCAGCTAGCTTAG"]

        @test_throws Exception splice_fasta(Analysis1_path, variant_dict, 100)
        
end # strings

# @testset "Using BioSequences" begin
    
#     @testset "normalizeDNA" begin
#         @test normalizeDNA("aatgn") == dna"AATGN"
#         @test_throws Exception normalizeDNA("ZCA")
#         @test_throws Exception normalizeDNA(42)
#         c = normalizeDNA('C') 
#         @test c == dna"c"
#         @test c isa LongSequence
#     end #  normalizeDNA

#     @testset "gc_content" begin
#         @test gc_content(dna"ANTG") == 0.25
#         @test gc_content(dna"cccggg") * 100 == 100.0
#         @test gc_content(dna"ATta") == 0.0
#     end #  composition

#     @testset "composition" begin
#         seq = rand(['A','T','G','C','N'], 20) |> join |> LongDNASeq
#         bc = composition(seq)

#         @test bc[DNA_A] == count(==(DNA_A), collect(seq))
#         @test bc[DNA_C] == count(==(DNA_C), collect(seq))
#         @test bc[DNA_G] == count(==(DNA_G), collect(seq))
#         @test bc[DNA_T] == count(==(DNA_T), collect(seq))
#         @test bc[DNA_N] == count(==(DNA_N), collect(seq))
#     end #  gc_content

#     @testset "complement" begin
#         @test complement(dna"ATTAN") == dna"TAATN"
#         @test complement(dna"gcta") == dna"CGAT"
#         @test complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  complement

#     @testset "reverse_complement" begin
#         @test reverse_complement(dna"ATTAN") == dna"NTAAT"
#         @test reverse_complement(dna"gcta") == dna"TAGC"
#         @test reverse_complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  reverse_complement

#     @testset "parse_fasta" begin
#         testpath = normpath(joinpath(@__DIR__, "..", "data"))
#         genomes = joinpath(testpath, "cov2_genomes.fasta")
#         ex1_path = joinpath(testpath, "ex1.fasta")
#         ex2_path = joinpath(testpath, "ex2.fasta")

#         ex1 = parse_fasta(ex1_path)
#         @test ex1 isa Tuple
#         @test all(x-> x isa String, ex1[1])
#         @test all(x-> x isa LongSequence, ex1[2])

#         @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
#         @test ex1[2] == [dna"AATTATAGC", dna"CGCCCCCCAGTCGGATT"]

#         @test_throws Exception parse_fasta(ex2_path)

#         cov2 = parse_fasta(genomes)
#         @test length(cov2[1]) == 8
#         @test length(cov2[2]) == 8
#     end # parse_fasta

# end # BioSequences

end # Assignment07
