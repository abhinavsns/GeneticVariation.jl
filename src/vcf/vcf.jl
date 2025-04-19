# vcf.jl
# ======
#
# A submodule for reading, writing, and working with VCF formatted files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE
#const MissingFieldException = BioGenerics.Exceptions.MissingFieldException

module VCF
import Automa
import Automa:@re_str, onenter!, onexit!, RegExp
import BioGenerics
import BioGenerics.Exceptions: MissingFieldException

using Indexes
using TranscodingStreams 


include("record.jl")
include("metainfo.jl")
include("header.jl")
include("reader.jl")
include("writer.jl")

end
