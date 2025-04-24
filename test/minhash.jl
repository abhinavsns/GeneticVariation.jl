@testset "MASH distances" begin
    # Define sequences
    seq_a = "ATCGCCA-"
    seq_b = "ATCGCCTA"

    # Set parameters
    k = 4  # k-mer size
    s = 100  # sketch size

    # Create sketches
    sketch_a = create_sketch(seq_a, k, s)
    sketch_b = create_sketch(seq_b, k, s)

    # Compute Mash distance
    mash_dist = mash(sketch_a, sketch_b, k)

    # Perform tests
    @test isapprox(mash_dist, 0.1277, atol=0.001)
    @test mash(sketch_a, sketch_a, k) == 0.0
    @test issorted(sketch_a.hashes)
end