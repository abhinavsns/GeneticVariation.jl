# vcf.jl
# ======
#
# A submodule for reading, writing, and working with VCF formatted files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE
#const MissingFieldException = BioGenerics.Exceptions.MissingFieldException

module VCF
using BioGenerics
import BioGenerics: BioGenerics, isfilled, header
import BioGenerics.Exceptions: MissingFieldException, missingerror
import BioGenerics.Automa: State
import TranscodingStreams: TranscodingStreams, TranscodingStream

#using Printf: @sprintf


function unsafe_parse_decimal(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int}) where {T<:Unsigned}
    x = zero(T)
    @inbounds for i in range
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return x
end

# r"[-+]?[0-9]+" must match `data[range]`.
function unsafe_parse_decimal(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int}) where {T<:Signed}
    lo = first(range)
    if data[lo] == UInt8('-')
        sign = T(-1)
        lo += 1
    elseif data[lo] == UInt8('+')
        sign = T(+1)
        lo += 1
    else
        sign = T(+1)
    end
    x = zero(T)
    @inbounds for i in lo:last(range)
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return sign * x
end

include("record.jl")
include("metainfo.jl")
include("header.jl")
include("reader.jl")
include("readrecord.jl")
include("writer.jl")

end
