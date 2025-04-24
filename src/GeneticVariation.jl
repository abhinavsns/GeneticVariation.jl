# GeneticVariation.jl
# ===================
#
# A julia package for the representation and analysis of genetic variation.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

__precompile__()

module GeneticVariation

export
    segregating_sites,
    count_segregating_sites,

    # Distances
    create_sketch, 
    jaccard, 
    mash, 
    distance,
    pdistance,
    pdistance_mutated,
    jukes_cantor,
    kimura_distance,
    # Allele frequencies
    gene_frequencies,

    # Diversity measures
    avg_mut,
    NL79,

    # VCF and BCF
    VCF,
    BCF,
    header

# Import only the necessary symbols from BioSequences v3.
import BioSequences:
    BioSequences,
    Alphabet,
    AA_Term,
    BioSequence,
    DNAAlphabet,
    GeneticCode,
    ispurine,
    RNAAlphabet

# Import metadata functions from BioGenerics
import BioGenerics: header, metainfotag, metainfoval, isfilled
using Indexes
using TranscodingStreams
#import BioGenerics.Exceptions: MissingFieldException, missingerror

import Combinatorics.permutations
import IntervalTrees: Interval, IntervalValue
import Twiddle:
    enumerate_nibbles,
    nibble_mask,
    count_0000_nibbles,
    count_nonzero_nibbles,
    count_1111_nibbles

include("vcf/vcf.jl")
include("bcf/bcf.jl")
include("site_counting.jl")
include("distances/minhash.jl")
include("distances/proportion.jl")
include("distances/evodistances.jl")
include("allele_freq.jl")
include("diversity_measures.jl")

end # Module GeneticVariation
