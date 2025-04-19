module ProportionDistance

export pdistance

"""
    pdistance(seq1::AbstractString, seq2::AbstractString) -> Float64

Calculate the proportion of sites that differ between two aligned sequences.
This is defined as the number of differing positions divided by the total length.
"""
function pdistance(seq1::AbstractString, seq2::AbstractString)
    n = length(seq1)
    if n != length(seq2)
        error("Sequences must be aligned and of equal length")
    end
    diff = sum(seq1[i] != seq2[i] for i in 1:n)
    return diff / n
end

end  # module ProportionDistance
