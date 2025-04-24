@testset "Site counting" begin
    @testset "Site Counting" begin
        using BioSequences

        # Define test sequences
        seq1 = dna"ATCG"
        seq2 = dna"ATGG"

        # Test matches and mismatches
        @test matches(seq1, seq2) == 3
        @test mismatches(seq1, seq2) == 1

        # Test iscertain count
        @test count(iscertain, seq1) == 4
        @test count(((x, y),) -> iscertain(x) && iscertain(y), zip(seq1, seq2)) == 4
    end

    @testset "Randomized Tests" begin
        # Generate random sequences
        seq_length = 100
        seq1 = randdnaseq(seq_length)
        seq2 = randdnaseq(seq_length)

        # Count matches and mismatches
        match_count = matches(seq1, seq2)
        mismatch_count = mismatches(seq1, seq2)

        # Total should equal sequence length
        @test match_count + mismatch_count == seq_length
    end


        @testset "Pairwise Methods" begin
        using BioSequences

        sequences = [
            dna"ATCG",
            dna"ATGG",
            dna"TTGG"
        ]

        # Compute pairwise mismatches
        for i in eachindex(sequences)
            for j in i:length(sequences)
                seq1 = sequences[i]
                seq2 = sequences[j]
                mismatch = mismatches(seq1, seq2)
                match = matches(seq1, seq2)
                @test mismatch + match == length(seq1)
            end
        end
    end


    @testset "Windowed Methods" begin
        using BioSequences

        seq1 = dna"ATCGATCG"
        seq2 = dna"ATGGATCC"
        window_size = 4

        for i in 1:(length(seq1)-window_size+1)
            window1 = seq1[i:i+window_size-1]
            window2 = seq2[i:i+window_size-1]
            match = matches(window1, window2)
            mismatch = mismatches(window1, window2)
            @test match + mismatch == window_size
        end
    end
end