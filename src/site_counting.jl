"""
    segregating_sites(seqs::Vector{<:BioSequence}) -> BitVector

Given an alignment of sequences (each assumed to have the same length),
returns a BitVector with `true` at positions that are segregating (i.e., variable),
and `false` otherwise.
"""
function segregating_sites(seqs::Vector{<:BioSequence})
    n = isempty(seqs) ? 0 : minimum(length.(seqs))
    segregating = falses(n)
    for pos in 1:n
        bases = Set(seq[pos] for seq in seqs)
        if length(bases) > 1
            segregating[pos] = true
        end
    end
    return segregating
end

"""
    count_segregating_sites(seqs::Vector{<:BioSequence}) -> Int

Counts and returns the number of segregating sites (i.e., variable sites)
in the given alignment.
"""
function count_segregating_sites(seqs::Vector{<:BioSequence})
    n = isempty(seqs) ? 0 : minimum(map(length, seqs))
    segregating = falses(n)
    for pos in 1:n
        bases = Set()
        for seq in seqs
            push!(bases, seq[pos])
        end
        if length(bases) > 1
            segregating[pos] = true
        end
    end
    return count(segregating), n
end
