# writer.jl
# =========
#
# A writer for VCF formatted files.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

# --- VCF Writer Type ---
"""
    mutable struct VCF.Writer{T<:IO} <: BioGenerics.IO.AbstractWriter

A VCF writer that sends data to a wrapped IO stream.
"""
mutable struct Writer <: BioGenerics.IO.AbstractWriter
    stream::IO
end

# Convenience constructor: if given a plain IO, wrap it in a NoopStream.
#Writer(io::IO) = Writer(NoopStream(io))

# --- VCF Writer: Header Writing ---
"""
    VCF.Writer(output::IO, header::VCF.Header)

Create a VCF writer by writing the header to the provided output.
# Arguments
* `output`: data sink (IO)
* `header`: a VCF.Header object that contains metadata (e.g. metainfo, sampleID, etc.)
"""
function Writer(output::IO, header::Header)
    writer = Writer(output)
    write(writer, header)
    return writer
end

# Expose the underlying stream.
function BioGenerics.IO.stream(writer::Writer)
    return writer.stream
end

# Write the VCF header.
function Base.write(writer::Writer, header::Header)
    n = 0
    # Write each meta-information line (assumed to be preformatted strings)
    for metainfo in header.metainfo
        n += write(writer.stream, metainfo, '\n')
    end
    # Write the VCF header line with columns.
    n += write(writer.stream, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    if !isempty(header.sampleID)
        n += write(writer.stream, "\tFORMAT")
    end
    for id in header.sampleID
        n += write(writer.stream, '\t', id)
    end
    n += write(writer.stream, '\n')
    return n
end

# Write an individual VCF record.
function Base.write(writer::Writer, record::Record)
    # Here we assume that the record object implements its own write method,
    # or overload Base.write to be able to convert it to a line of text.
    return write(writer.stream, record, '\n')
end
