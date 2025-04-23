# metainfo.jl
# ===========
#
# A representation of metainformation in a VCF file.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

mutable struct MetaInfo
    data::Vector{UInt8}             # Raw header data as bytes
    filled::UnitRange{Int}          # Range of data that has been "filled"
    dict::Bool                      # True if value is a dictionary (<...>)
    tag::UnitRange{Int}             # Range in data for the tag (after "##" until '=')
    val::UnitRange{Int}             # Range in data for the value (after '=' until end)
    dictkey::Vector{UnitRange{Int}} # For dict values: ranges for each inner key
    dictval::Vector{UnitRange{Int}} # For dict values: ranges for each corresponding value
end

"""
    MetaInfo()

Create an unfilled VCF metainfo object.
"""
function MetaInfo()
    return MetaInfo(UInt8[], 1:0, false, 1:0, 1:0, UnitRange{Int}[], UnitRange{Int}[])
end

"""
    MetaInfo(data::Vector{UInt8})

Create a VCF metainfo object from raw data. Ownership of `data` is transferred.
The header line is validated and indexed.
"""
function MetaInfo(data::Vector{UInt8})
    return convert(MetaInfo, data)
end

function Base.convert(::Type{MetaInfo}, data::Vector{UInt8})
    mi = MetaInfo(data, 1:0, false, 1:0, 1:0, UnitRange{Int}[], UnitRange{Int}[])
    index!(mi)
    return mi
end

"""
    MetaInfo(str::AbstractString)

Create a VCF metainfo object from a header string.
"""
function MetaInfo(str::AbstractString)
    return convert(MetaInfo, str)
end

function Base.convert(::Type{MetaInfo}, str::AbstractString)
    return MetaInfo(Vector{UInt8}(str))
end

# =============================================================================
# Indexing Functions
# =============================================================================

"""
    index!(mi::MetaInfo)

Index a MetaInfo object by splitting the header line into tag and value portions.
If the value is delimited by '<' and '>', the record is treated as a dictionary;
in that case, `index_dict!` is called.
"""
function index!(mi::MetaInfo)
    # Verify the header begins with "##"
    if length(mi.data) < 3 || mi.data[1] != UInt8('#') || mi.data[2] != UInt8('#')
        throw(ArgumentError("Invalid meta-information format: must start with ##"))
    end

    # Mark the header as completely filled.
    mi.filled = 1:length(mi.data)

    # Find first '=' (separator between tag and value)
    eq_index = findfirst(x -> x == UInt8('='), mi.data)
    if eq_index === nothing
        throw(ArgumentError("No '=' found in meta-information"))
    end

    # The tag is from after "##" until just before '='.
    mi.tag = 3:(eq_index-1)
    # The value starts immediately after '=' and goes until the end.
    mi.val = (eq_index+1):length(mi.data)

    # Check if the value appears to be a dictionary (enclosed in < ... >)
    if mi.val.start <= mi.val.stop &&
        mi.data[mi.val.start] == UInt8('<') && mi.data[mi.val.stop] == UInt8('>')
        mi.dict = true
        index_dict!(mi)
    else
        mi.dict = false
    end

    return mi
end

"""
    index_dict!(mi::MetaInfo)

Parse the inner content of a dictionary metainfo (between '<' and '>').
Splits on unquoted commas and then on the first '=' in each token to compute
byte ranges for the keys and values.
"""
function index_dict!(mi::MetaInfo)
    inner_start = mi.val.start + 1
    inner_end = mi.val.stop - 1
    inner = mi.data[inner_start:inner_end]
    empty!(mi.dictkey)
    empty!(mi.dictval)
    ranges = split_unquoted_commas(inner)
    for r in ranges
        part = inner[r]
        eq_rel = findfirst(x -> x == UInt8('='), part)
        if eq_rel === nothing
            throw(ArgumentError("Malformed dictionary entry in metainfo"))
        end
        key_start = inner_start + first(r) - 1
        key_end = key_start + eq_rel - 2
        val_start = key_end + 2
        val_end = inner_start + last(r) - 1
        push!(mi.dictkey, key_start:key_end)
        push!(mi.dictval, val_start:val_end)
    end
end

"""
    initialize!(mi::MetaInfo)

Reset the MetaInfo object to an unfilled state.
"""
function initialize!(mi::MetaInfo)
    mi.filled = 1:0
    mi.dict = false
    mi.tag = 1:0
    mi.val = 1:0
    empty!(mi.dictkey)
    empty!(mi.dictval)
    return mi
end

# =============================================================================
# Basic Properties and Equality
# =============================================================================

datarange(mi::MetaInfo) = mi.filled
isfilled(mi::MetaInfo) = !isempty(mi.filled)

function Base.:(==)(mi1::MetaInfo, mi2::MetaInfo)
    if isfilled(mi1) && isfilled(mi2)
        r1 = datarange(mi1)
        r2 = datarange(mi2)
        return length(r1) == length(r2) &&
            ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
            pointer(mi1.data, first(r1)),
            pointer(mi2.data, first(r2)),
            length(r1)) == 0
    else
        return !isfilled(mi1) && !isfilled(mi2)
    end
end

function checkfilled(mi::MetaInfo)
    if !isfilled(mi)
        throw(ArgumentError("unfilled VCF metainfo"))
    end
end

# =============================================================================
# Editing/Reconstructing MetaInfo
# =============================================================================

"""
    MetaInfo(base::MetaInfo; tag=nothing, value=nothing)

Create a new VCF metainfo based on an existing one. Optionally replace the tag or value.
If `value` is a dictionary (as AbstractDict or a vector of key-value pairs), the new meta line
will enclose the content in angle brackets and quote parts if necessary.
"""
function MetaInfo(base::MetaInfo; tag=nothing, value=nothing)
    checkfilled(base)
    buf = IOBuffer()
    print(buf, "##")
    if tag === nothing
        write(buf, base.data[base.tag])
    else
        print(buf, tag)
    end
    print(buf, '=')
    if value === nothing
        write(buf, base.data[base.val])
    elseif isa(value, String)
        print(buf, value)
    elseif isa(value, AbstractDict) || isa(value, Vector)
        print(buf, '<')
        # Iterate over key–value pairs.
        firstpair = true
        for (k, v) in value
            if !firstpair
                print(buf, ',')
            end
            print(buf, k, '=')
            if isa(v, String) && needs_quote(v)
                print(buf, '"', escape_string(v), '"')
            else
                print(buf, v)
            end
            firstpair = false
        end
        print(buf, '>')
    else
        throw(ArgumentError("value must be String, AbstractDict, or Vector"))
    end
    return MetaInfo(take!(buf))
end

# A simple helper that returns true if the string should be quoted.
needs_quote(val::String) = occursin(r"[ ,\"\\]", val)

# Dummy escape function; in practice, implement proper escaping as needed.
escape_string(val::String) = replace(val, "\"" => "\\\"")

# =============================================================================
# Equality for the Tag Field
# =============================================================================

function isequaltag(mi::MetaInfo, tag::AbstractString)
    checkfilled(mi)
    return length(mi.tag) == sizeof(tag) &&
            ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
        pointer(mi.data, first(mi.tag)),
        pointer(tag), length(mi.tag)) == 0
end

# =============================================================================
# Accessor Functions
# =============================================================================

function metainfotag(mi::MetaInfo)
    checkfilled(mi)
    return String(mi.data[mi.tag])
end

function metainfoval(mi::MetaInfo)
    checkfilled(mi)
    return String(mi.data[mi.val])
end

# =============================================================================
# Dictionary Interface
# =============================================================================

function Base.keys(mi::MetaInfo)
    checkfilled(mi)
    if !mi.dict
        throw(ArgumentError("not a dictionary"))
    end
    return [String(mi.data[r]) for r in mi.dictkey]
end

function Base.values(mi::MetaInfo)
    checkfilled(mi)
    if !mi.dict
        throw(ArgumentError("not a dictionary"))
    end
    vals = String[]
    for r in mi.dictval
        push!(vals, extract_metainfo_value(mi.data, r))
    end
    return vals
end

function Base.getindex(mi::MetaInfo, key::String)
    checkfilled(mi)
    if !mi.dict
        throw(ArgumentError("not a dictionary"))
    end
    for (i, r) in enumerate(mi.dictkey)
        n = length(r)
        if n == sizeof(key) &&
           ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
            pointer(mi.data, first(r)),
            pointer(key), n) == 0
            return extract_metainfo_value(mi.data, mi.dictval[i])
        end
    end
    throw(KeyError(key))
end

# Helper to extract a value from the raw data range.
function extract_metainfo_value(data::Vector{UInt8}, range::UnitRange{Int})
    lo = first(range)
    hi = last(range)
    # If the string is enclosed in quotes then remove them and unescape.
    if length(range) ≥ 2 && data[lo] == UInt8('"') && data[hi] == UInt8('"')
        return unescape_string(String(data[lo+1:hi-1]))
    else
        return String(data[lo:hi])
    end
end

# Dummy unescape (for demonstration, you may want a more robust version)
function unescape_string(s::String)
    return replace(s, "\\\"" => "\"")
end

# =============================================================================
# Display and Write
# =============================================================================

function Base.show(io::IO, mi::MetaInfo)
    print(io, "VCF MetaInfo (", summary(mi), "):")
    if isfilled(mi)
        println(io)
        println(io, "    tag: ", metainfotag(mi))
        print(io, "  value:")
        if mi.dict
            for (k, v) in zip(keys(mi), values(mi))
                print(io, ' ', k, "=\"", v, "\"")
            end
        else
            print(io, ' ', metainfoval(mi))
        end
    else
        print(io, " <not filled>")
    end
end

function Base.write(io::IO, mi::MetaInfo)
    checkfilled(mi)
    return write(io, mi.data)
end

function split_unquoted_commas(data::Vector{UInt8})
    parts = UnitRange{Int}[]
    in_quotes = false
    start = 1
    for i in 1:length(data)
        b = data[i]
        if b == UInt8('"')
            in_quotes = !in_quotes
        elseif b == UInt8(',') && !in_quotes
            push!(parts, start:(i-1))
            start = i + 1
        end
    end
    push!(parts, start:length(data))
    return parts
end