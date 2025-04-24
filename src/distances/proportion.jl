using BioSequences

"""
    pdistance(seq1::BioSequence, seq2::BioSequence) -> Float64

Calculate the proportion of sites that differ between two aligned sequences.
This is defined as the number of differing positions divided by the total length.
"""
function pdistance(seq1::BioSequence, seq2::BioSequence)
    n = length(seq1)
    if n != length(seq2)
        error("Sequences must be aligned and of equal length")
    end
    diff = sum(seq1[i] != seq2[i] for i in 1:n)
    return diff / n
end

"""
    pdistance(seq1::BioSequence, seq2::BioSequence) -> Float64

Calculate the proportion of mutated sites between two aligned sequences.
This is defined as the number of differing positions (excluding ambiguous symbols and gaps)
divided by the number of valid (certain) positions.
"""
function pdistance_mutated(seq1::BioSequence, seq2::BioSequence)
    if length(seq1) != length(seq2)
        error("Sequences must be aligned and of equal length")
    end
    mutated_sites = 0
    total_certain_sites = 0
    for i in eachindex(seq1)
        base1 = seq1[i]
        base2 = seq2[i]
        if iscertain(base1) && iscertain(base2)
            total_certain_sites += 1
            if base1 != base2
                mutated_sites += 1
            end
        end
    end
    if total_certain_sites == 0
        return 0.0  # or throw an error if appropriate
    end
    return mutated_sites / total_certain_sites
end
