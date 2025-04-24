
@testset "Proportion distances" begin
    @testset "Biological Sequences" begin
        dnaA = dna"ATCGCCA-"
        dnaB = dna"ATCCCCTA"
        rnaA = rna"AUCGCCA-"
        rnaB = rna"AUCCCCUA"

        # Test Mismatch distance
        @test isapprox(pdistance(dnaA, dnaB), 0.375; atol=1e-3)
        @test isapprox(pdistance(rnaA, rnaB), 0.375; atol=1e-3)

        # Test Mutated distance
        @test isapprox(pdistance_mutated(dnaA, dnaB), 0.2857142857142857; atol=1e-3)
        @test isapprox(pdistance_mutated(rnaA, rnaB), 0.2857142857142857; atol=1e-3)
    end
end