module TestGeneticVariation

using Test

using BioSequences, GeneticVariation
using TranscodingStreams
import IntervalTrees: IntervalValue
import YAML
using FormatSpecimens

import GeneticVariation.VCF: isfilled, metainfotag, metainfoval, VCF, VCF.Reader, VCF.Writer, VCF.Record
import BioGenerics.Exceptions: MissingFieldException

function random_seq(::Type{A}, n::Integer) where A <: Alphabet
    nts = alphabet(A)
    probs = Vector{Float64}(undef, length(nts))
    fill!(probs, 1 / length(nts))
    return BioSequence{A}(random_seq(n, nts, probs))
end

include("vcf.jl")
include("bcf.jl")
include("site_counting.jl")
include("minhash.jl")
include("allele_freq.jl")
include("diversity_measures.jl")
include("seg_sites.jl")

end # Module TestGeneticVariation
