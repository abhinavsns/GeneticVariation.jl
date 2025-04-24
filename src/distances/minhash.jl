using MinHash
using Statistics
import Base: log
"""
    create_sketch(seq::AbstractString, k::Int, s::Int) -> Vector{UInt64}

Computes a MinHash sketch for the given sequence `seq` using k-mer length `k`
and sketch size `s`. This implementation uses a simple sliding-window
approach to extract k-mers. For more advanced k-mer iteration, consider using Kmers.jl.
""" 
function create_sketch(seq::AbstractString, k::Int, s::Int)
    n = length(seq)
    if n < k
        error("Sequence length ($(n)) is shorter than k ($(k)).")
    end
    kmers = [seq[i:i+k-1] for i in 1:(n-k+1)]
    # Compute the sketch using MinHash.jl
    return MinHash.sketch(kmers, s)
end

"""
    jaccard(sketch1, sketch2) -> Float64

Returns the approximate Jaccard similarity between two MinHash sketches.
It is computed as the ratio of shared hashes to the total number of unique hashes.
"""
function jaccard(sketch1::MinHash.MinHashSketch, sketch2::MinHash.MinHashSketch)
    s1 = Set(sketch1.hashes)
    s2 = Set(sketch2.hashes)
    return length(intersect(s1, s2)) / length(union(s1, s2))
end

"""
    mash(sketch1, sketch2, k::Int) -> Float64

Computes the Mash distance between two MinHash sketches using the formula:

    D = -1/k * log( (2*J)/(1+J) )

where J is the Jaccard similarity. Returns `Inf` if J is 0.
"""
function mash(sketch1::MinHash.MinHashSketch, sketch2::MinHash.MinHashSketch, k::Int)
    J = jaccard(sketch1, sketch2)
    if J == 0.0
        return Inf
    end
    return -1.0 / k * log((2 * J) / (1 + J))
end

"""
    distance(metric::Symbol, sketch1, sketch2, k::Int) -> Float64

Dispatches to the appropriate distance metric. Currently supported metrics:
 - `:jaccard`: returns (1 - Jaccard similarity).
 - `:mash`: returns the Mash distance.
"""
function distance(metric::Symbol, sketch1::MinHash.MinHashSketch, sketch2::MinHash.MinHashSketch, k::Int)
    if metric == :jaccard
        return 1.0 - jaccard(sketch1, sketch2)
    elseif metric == :mash
        return mash(sketch1, sketch2, k)
    else
        error("Unsupported metric: $metric")
    end
end
