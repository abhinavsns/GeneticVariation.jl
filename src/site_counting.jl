module SiteCounting

export segregating_sites, count_segregating_sites

using BioSequences: BioSequence

"""
    segregating_sites(seqs::Vector{<:BioSequence}) -> BitVector

Given an alignment of sequences (each assumed to have the same length),
returns a BitVector with `true` at positions that are segregating (i.e. variable),
and `false` otherwise.
"""
function segregating_sites(seqs::Vector{<:BioSequence})
    n = isempty(seqs) ? 0 : length(seqs[1])
    segregating = falses(n)
    # Iterate over each column (site)
    for pos in 1:n
        bases = Set{UInt8}()
        for seq in seqs
            push!(bases, seq[pos])
        end
        # Mark site as segregating if more than one unique base exists.
        # (You can add additional logic to ignore gaps or missing symbols if needed.)
        if length(bases) > 1
            segregating[pos] = true
        end
    end
    return segregating
end

"""
    count_segregating_sites(seqs::Vector{<:BioSequence}) -> Int

Counts and returns the number of segregating sites (i.e. variable sites)
in the given alignment.
"""
function count_segregating_sites(seqs::Vector{<:BioSequence})
    seg_sites = segregating_sites(seqs)
    return count(seg_sites)
end

end # module SiteCounting
