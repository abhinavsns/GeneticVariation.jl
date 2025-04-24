"""
    jukes_cantor(seq1::BioSequence, seq2::BioSequence) -> Float64

Calculate the Jukes–Cantor corrected distance between two aligned DNA sequences.
The formula is: d = -3/4 * log(1 - (4/3) * p)
where p is the observed proportion of differences.
Returns Inf if p ≥ 0.75.
"""
function jukes_cantor(seq1::BioSequence, seq2::BioSequence)
    n = length(seq1)
    if n != length(seq2)
        error("Sequences must be the same length")
    end
    p = sum(seq1[i] != seq2[i] for i in 1:n) / n
    if p >= 0.75
        return Inf
    end
    return -0.75 * log(1 - (4/3) * p)
end

"""
    kimura_distance(seq1::BioSequence, seq2::BioSequence) -> Float64

Calculate the Kimura 2–parameter distance between two aligned DNA sequences.
This method takes into account both transitions (ts) and transversions (tv)
using the formula:
    
    d = -0.5*log(1 - 2P - Q) - 0.25*log(1 - 2Q)

where P is the proportion of transitions and Q is the proportion of transversions.
"""
function kimura_distance(seq1::BioSequence, seq2::BioSequence)
    n = length(seq1)
    if n != length(seq2)
        error("Sequences must be the same length")
    end
    transitions = 0
    transversions = 0
    # Define transitions: A<->G and C<->T.
    transitions_set = Set([(dna"A", dna"G"), (dna"G", dna"A"), (dna"C", dna"T"), (dna"T", dna"C")])
    for i in 1:n
        a, b = seq1[i], seq2[i]
        if a != b
            if (a, b) in transitions_set
                transitions += 1
            else
                transversions += 1
            end
        end
    end
    P = transitions / n
    Q = transversions / n
    # Avoid log of a non-positive number:
    if 1 - 2*P - Q <= 0 || 1 - 2*Q <= 0
        return Inf
    end
    d = -0.5 * log(1 - 2*P - Q) - 0.25 * log(1 - 2*Q)
    return d
end
