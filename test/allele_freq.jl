@testset "Allele frequenies" begin
    @testset "Gene frequencies" begin

        sequences = ["ATCGATCG", "AGGGGG", "CCCCCCCCCCCCCCC", "TTTTTCCCC",
                    "ATCGATCG", "AGGGGG", "ATCGATCG", "ATCGATCG"]

        answer = Dict{String, Float64}("TTTTTCCCC" => 0.125,
                                        "CCCCCCCCCCCCCCC" => 0.125,
                                        "ATCGATCG" => 0.5,
                                        "AGGGGG" => 0.25)

        function test_gene_frequencies(genes, answer)
            test_genes = [bioseq(sq) for sq in sequences]
            test_answer = Dict{LongDNA, Float64}(bioseq(key) => val for (key, val) in answer)
            @test gene_frequencies(test_genes) == test_answer
        end

        for n in (2, 4)
            test_gene_frequencies(sequences, answer)
        end
    end
end