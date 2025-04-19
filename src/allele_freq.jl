module AlleleFreq

export gene_frequencies

using BioSequences 

"""
    gene_frequencies(iterable)

Compute allele frequencies from any iterable whose element type is a subtype of BioSequence.
This function iterates over the input, counts occurrences of unique sequences,
and then returns a dictionary mapping each sequence to its relative frequency.
Example:
    freqs = gene_frequencies(["ATGC", "ATGC", "ATGT"])
"""
function gene_frequencies(iterable)
    counts = Dict{eltype(iterable),Int}()
    total = 0
    for seq in iterable
        counts[seq] = get(counts, seq, 0) + 1
        total += 1
    end
    return Dict(k => v / total for (k, v) in counts)
end

end  # module AlleleFreq
