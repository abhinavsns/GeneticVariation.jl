# record.jl
# =========
#
# Representation of a record of a VCF file.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/GeneticVariation.jl/blob/master/LICENSE

# -- Data structure --------------------------------------------------

mutable struct Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    ncols::Int
    # indexes into data vector for fields (all fields are ranges or arrays of ranges)
    chrom::UnitRange{Int}
    pos::UnitRange{Int}
    id::Vector{UnitRange{Int}}
    ref::UnitRange{Int}
    alt::Vector{UnitRange{Int}}
    qual::UnitRange{Int}
    filter::Vector{UnitRange{Int}}
    infokey::Vector{UnitRange{Int}}
    format::Vector{UnitRange{Int}}
    genotype::Vector{Vector{UnitRange{Int}}}
end

# Default (empty/unfilled) constructor
"""
    VCF.Record()

Create an unfilled VCF record.
"""
function Record()
    return Record(
        UInt8[], 1:0, 0,
        1:0, 1:0, UnitRange{Int}[], 1:0, UnitRange{Int}[],
        1:0, UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[]
    )
end

#copy constructor
"""
    Record(record::Record)
Create a copy of the given VCF record.
"""
function Record(base::Record;
    chrom=nothing, pos=nothing, id=nothing,
    ref=nothing, alt=nothing, qual=nothing,
    filter=nothing, info=nothing, genotype=nothing)
    checkfilled(base)
    buf = IOBuffer()

    if chrom == nothing
        write(buf, base.data[base.chrom])
    else
        print(buf, string(chrom))
    end

    print(buf, '\t')
    if pos == nothing
        write(buf, base.data[base.pos])
    else
        print(buf, convert(Int, pos))
    end

    print(buf, '\t')
    if id == nothing
        if isempty(base.id)
            print(buf, '.')
        else
            for (i, r) in enumerate(base.id)
                if i != 1
                    print(buf, ';')
                end
                write(buf, base.data[r])
            end
        end
    else
        if !isa(id, Vector)
            id = [id]
        end
        if isempty(id)
            print(buf, '.')
        else
            for (i, x) in enumerate(id)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if ref == nothing
        write(buf, base.data[base.ref])
    else
        print(buf, string(ref))
    end

    print(buf, '\t')
    if alt == nothing
        if isempty(base.alt)
            print(buf, '.')
        else
            for (i, r) in enumerate(base.alt)
                if i != 1
                    print(buf, ';')
                end
                write(buf, base.data[r])
            end
        end
    else
        if !isa(alt, Vector)
            alt = [alt]
        end
        if isempty(alt)
            print(buf, '.')
        else
            for (i, x) in enumerate(alt)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if qual == nothing
        write(buf, base.data[base.qual])
    else
        print(buf, convert(Float64, qual))
    end

    print(buf, '\t')
    if filter == nothing
        if isempty(base.filter)
            print(buf, '.')
        else
            for (i, r) in enumerate(base.filter)
                if i != 1
                    print(buf, ';')
                end
                write(buf, base.data[r])
            end
        end
    else
        if !isa(filter, Vector)
            filter = [filter]
        end
        if isempty(filter)
            print(buf, '.')
        else
            for (i, x) in enumerate(filter)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(x))
            end
        end
    end

    print(buf, '\t')
    if info == nothing
        if isempty(base.infokey)
            print(buf, '.')
        else
            write(buf, base.data[first(base.infokey[1]):last(infovalrange(base, lastindex(base.infokey)))])
        end
    else
        if !isa(info, AbstractDict)
            throw(ArgumentError("info must be an AbstractDict"))
        elseif isempty(info)
            print(buf, '.')
        else
            for (i, (key, val)) in enumerate(info)
                if i != 1
                    print(buf, ';')
                end
                print(buf, string(key))
                if val != nothing
                    print(buf, '=', vcfformat(val))
                end
            end
        end
    end

    print(buf, '\t')
    if genotype == nothing
        if isempty(base.format)
            print(buf, '.')
        else
            write(buf, base.data[first(base.format[1]):last(base.format[end])])
        end
        if !isempty(base.genotype)
            for indiv in base.genotype
                print(buf, '\t')
                for (i, r) in enumerate(indiv)
                    if i != 1
                        print(buf, ':')
                    end
                    write(buf, base.data[r])
                end
            end
        end
    else
        if !isa(genotype, Vector)
            genotype = [genotype]
        end
        if isempty(genotype)
            print(buf, '.')
        else
            allkeys = String[]
            for indiv in genotype
                if !isa(indiv, AbstractDict)
                    throw(ArgumentError("individual must be anabstract dictionary"))
                end
                append!(allkeys, keys(indiv))
            end
            allkeys = sort!(unique(allkeys))
            if !isempty(allkeys)
                join(buf, allkeys, ':')
                for indiv in genotype
                    print(buf, '\t')
                    for (i, key) in enumerate(allkeys)
                        if i != 1
                            print(buf, ':')
                        end
                        print(buf, vcfformat(get(indiv, key, '.')))
                    end
                end
            end
        end
    end

    return Record(take!(buf))
end

# Constructor from a data vector, delegating to convert
"""
    Record(data::Vector{UInt8})

Create a VCF record from the given data vector.
This will index the record and validate its structure.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

# Conversion from raw data vector to a Record
function Base.convert(::Type{Record}, data::Vector{UInt8})
    rec = Record(
        data, 1:0, 0,
        1:0, 1:0, UnitRange{Int}[], 1:0, UnitRange{Int}[],
        1:0, UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[], UnitRange{Int}[]
    )
    index!(rec)  # note: index! must update rec.filled and all field ranges
    return rec
end

# Constructor from a string (by converting to Vector{UInt8})
"""
    Record(str::AbstractString)

Create a VCF record from the given string.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return convert(Record, Vector{UInt8}(str))
end

function Base.empty!(record::Record)
    record.filled = 1:0
    record.chrom = 1:0
    record.pos = 1:0
    empty!(record.id)
    record.ref = 1:0
    empty!(record.alt)
    record.qual = 1:0
    empty!(record.filter)
    empty!(record.infokey)
    empty!(record.format)
    empty!(record.genotype)
    return record
end

# Check whether record has been filled (i.e. indexed)
function isfilled(record::Record)
    return !isempty(record.filled)
end

# Return the range of valid data within the record
function datarange(record::Record)
    return record.filled
end

# Equality: compare filled regions
function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = record1.filled
        r2 = record2.filled
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    end

    return isfilled(record1) == isfilled(record2) == false
end

# A helper for formatting INFO field values (simply converts to string)
#vcfformat(val) = string(val)
#vcfformat(val::Vector) = join(map(vcfformat, val), ",")

function Base.copy(rec::Record)
    return Record(
        rec.data[rec.filled],
        rec.filled,
        rec.ncols,
        rec.chrom,
        rec.pos,
        copy(rec.id),
        rec.ref,
        copy(rec.alt),
        rec.qual,
        copy(rec.filter),
        copy(rec.infokey),
        copy(rec.format),
        deepcopy(rec.genotype)
    )
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled))
end
function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, "VCF Record")
    if isfilled(record)
        println(io)
        println(io, "  chromosome: ", haschrom(record) ? chrom(record) : "<missing>")
        println(io, "      position: ", haspos(record) ? pos(record) : "<missing>")
        println(io, "  identifier(s): ", hasid(record) ? join(id(record), " ") : "<missing>")
        println(io, "   reference: ", hasref(record) ? ref(record) : "<missing>")
        println(io, "   alternate: ", hasalt(record) ? join(alt(record), " ") : "<missing>")
        println(io, "     quality: ", hasqual(record) ? qual(record) : "<missing>")
        println(io, "      filter: ", hasfilter(record) ? join(filter(record), " ") : "<missing>")
        print(io, "     information: ")
        if hasinfo(record)
            for (key, val) in info(record)
                print(io, key)
                if !isempty(val)
                    print(io, '=', val, ' ')
                else
                    print(io, ' ')
                end
            end
        else
            print(io, "<missing>")
        end
        println(io)
        print(io, "       format: ", hasformat(record) ? join(format(record), " ") : "<missing>")
        if hasformat(record)
            println(io)
            print(io, "     genotype:")
            for i in 1:length(record.genotype)
                print(io, " [$(i)] ", begin
                    g = genotype(record, i)
                    isempty(g) ? "." : join(g, " ")
                end)
            end
        end
    else
        print(io, " <not filled>")
    end
end

# Accessor functions
# ------------------

"""
    chrom(record::Record)::String

Return the chromosome as a String. Throws MissingFieldException if missing.
"""
function chrom(record::Record)::String
    checkfilled(record)
    return String(record.data[record.chrom])
end

function haschrom(record::Record)
    return isfilled(record) && !ismissing(record, record.chrom)
end

"""
    pos(record::Record)::Int

Return the reference position as an Int. Throws MissingFieldException if missing.
"""
function pos(record::Record)::Int
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.pos)
end

function haspos(record::Record)
    return isfilled(record) && !ismissing(record, record.pos)
end

"""
    id(record::Record)::Vector{String}

Return the identifiers from the record as a Vector of Strings.
Throws MissingFieldException if missing.
"""
function id(record::Record)::Vector{String}
    checkfilled(record)
    if isempty(record.id) || ismissing(record, record.id[1])
        missingerror(:id)
    end
    return [String(record.data[r]) for r in record.id]
end

function hasid(record::Record)
    return record.ncols ≥ 3 && !isempty(record.id) && !ismissing(record, record.id[1])
end

"""
    ref(record::Record)::String

Return the reference bases as a String. Throws MissingFieldException if missing.
"""
function ref(record::Record)::String
    checkfilled(record)
    if ismissing(record, record.ref)
        missingerror(:ref)
    end
    return String(record.data[record.ref])
end

function hasref(record::Record)
    return record.ncols ≥ 4 && !ismissing(record, record.ref)
end

"""
    alt(record::Record)::Vector{String}

Return the alternate alleles as a Vector of Strings.
Throws MissingFieldException if missing.
"""
function alt(record::Record)::Vector{String}
    checkfilled(record)
    if isempty(record.alt)
        missingerror(:alt)
    end
    return [String(record.data[r]) for r in record.alt]
end

function hasalt(record::Record)
    return record.ncols ≥ 5 && !ismissing(record, record.alt[1])
end

"""
    qual(record::Record)::Float64

Return the quality score as a Float64.
Throws MissingFieldException if missing.
"""
function qual(record::Record)::Float64
    checkfilled(record)
    if ismissing(record, record.qual)
        missingerror(:qual)
    end
    # TODO: no-copy parse
    return parse(Float64, String(record.data[record.qual]))
end

function hasqual(record::Record)
    return record.ncols ≥ 6 && !ismissing(record, record.qual)
end

"""
    filter(record::Record)::Vector{String}

Return the filter field as a Vector of Strings.
Throws MissingFieldException if missing.
"""
function filter(record::Record)::Vector{String}
    checkfilled(record)
    if ismissing(record,record.filter[1])
        missingerror(:filter)
    end
    return [String(record.data[r]) for r in record.filter]
end

function hasfilter(record::Record)
    return record.ncols ≥ 7 && !ismissing(record, record.filter[1])
end 

"""
    info(record::Record)::Vector{Pair{String,String}}

Return the INFO field as a vector of key=>value pairs.
Throws MissingFieldException if missing.
"""
function info(record::Record)::Vector{Pair{String,String}}
    checkfilled(record)
    if isempty(record.infokey)
        missingerror(:info)
    end
    ret = Pair{String,String}[]
    for i in eachindex(record.infokey)
        key = record.infokey[i]
        val = infovalrange(record, i)
        push!(ret, String(record.data[key]) => String(record.data[val]))
    end
    return ret
end

function hasinfo(record::Record)
    return record.ncols ≥ 7 && !isempty( record.infokey)
end

"""
    info(record::Record, key::String)::String

Return the INFO field value for key.
Throws KeyError if key is not present.
"""
function info(record::Record, key::String)::String
    checkfilled(record)
    i = findinfokey(record, key)
    if i == 0
        throw(KeyError(key))
    end
    val = infovalrange(record, i)
    return isempty(val) ? "" : String(record.data[val])
end

function hasinfo(record::Record, key::String)
    return record.ncols ≥ 7 && findinfokey(key) > 0
end

# Find the index of a key in the INFO keys
function findinfokey(record::Record, key::String)
    for i in eachindex(record.infokey)
        if isequaldata(key, record.data, record.infokey[i])
            return i
        end
    end
    return 0
end

"""
    infokeys(record::Record)::Vector{String}

Return all keys in the INFO field.
"""
function infokeys(record::Record)::Vector{String}
    checkfilled(record)
    return [String(record.data[r]) for r in record.infokey]
end

# Helper: given the index of an INFO key, return the range in data for its value.
function infovalrange(record::Record, i::Int)
    checkfilled(record)
    data = record.data
    key_range = record.infokey[i]
    if last(key_range) + 1 ≤ lastindex(data) && data[last(key_range)+1] == UInt8('=')
        endpos = findnext(x -> x == UInt8(';'), data, last(key_range) + 1)
        if endpos === nothing
            endpos = findnext(x -> x == UInt8('\t'), data, last(key_range) + 1)
            @assert endpos !== nothing "Could not determine end position in INFO value"
        end
        return (last(key_range)+2):(endpos-1)
    else
        return (last(key_range)+1):last(key_range)
    end
end

"""
    format(record::Record)::Vector{String}

Return the genotype format fields as a vector of Strings.
Throws MissingFieldException if missing.
"""
function format(record::Record)::Vector{String}
    checkfilled(record)
    if isempty(record.format)
        missingerror(:format)
    end
    return [String(record.data[r]) for r in record.format]
end

function hasformat(record::Record)
    return record.ncols ≥ 7 && !isempty(record.format) && !ismissing(record, record.format[1])
end

"""
    genotype(record::Record)::Vector{Vector{String}}

Return all genotype fields as an array of arrays of Strings.
"""
function genotype(record::Record)
    checkfilled(record)
    ret = Vector{Vector{String}}()
    for i in 1:length(record.genotype)
        push!(ret, genotype_impl(record, i, 1:length(record.format)))
    end
    return ret
end

"""
    genotype(record::Record, index::Integer)::Vector{String}

Return genotype fields for a specific individual.
"""
function genotype(record::Record, index::Integer)
    checkfilled(record)
    return genotype_impl(record, index, 1:length(record.format))
end

"""
    genotype(record::Record, index::Integer, key::String)::String

Return the genotype field for a specific individual and key.
"""
function genotype(record::Record, index::Integer, key::String)::String
    checkfilled(record)
    k = findgenokey(record, key)
    if k === nothing
        throw(KeyError(key))
    end
    return genotype_impl(record, index, k)
end

function genotype(record::Record, index::Integer, keys::AbstractVector{String})::Vector{String}
    checkfilled(record)
    return [genotype(record, index, key) for key in keys]
end

function genotype(record::Record, indexes::AbstractVector{T}, key::String) where T<:Integer
    checkfilled(record)
    k = findgenokey(record, key)
    if k === nothing
        throw(KeyError(key))
    end
    return [genotype_impl(record, i, k) for i in indexes]
end

function genotype(record::Record, indexes::AbstractVector{T}, keys::AbstractVector{String}) where T<:Integer
    checkfilled(record)
    ks = map(key -> begin
            k = findgenokey(record, key)
            if k === nothing
                throw(KeyError(key))
            end
            k
        end, keys)
    return [genotype_impl(record, i, ks) for i in indexes]
end

function genotype(record::Record, ::Colon, key::String)::Vector{String}
    return genotype(record, 1:length(record.genotype), key)
end

# Find the index of a genotype key in the format field ranges
function findgenokey(record::Record, key::String)
    return findfirst(r -> isequaldata(key, record.data, r), record.format)
end

# Implementation for genotype accessor for a given genotype (index) and key(s)
function genotype_impl(record::Record, index::Int, keys::Union{Int,AbstractVector{Int}})
    if isa(keys, Int)
        keys = [keys]
    end
    return [genotype_impl(record, index, k) for k in keys]
end

function genotype_impl(record::Record, index::Int, key::Int)
    geno = record.genotype[index]
    if key > length(geno)
        return "."
    else
        return String(record.data[geno[key]])
    end
end

# Check if the data in a given range equals the provided string.
function isequaldata(str::String, data::Vector{UInt8}, range::UnitRange{Int})
    return length(range) == sizeof(str) && ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t),
        pointer(data, first(range)), pointer(str), length(range)) == 0
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled VCF record"))
    end
end

function ismissing(record::Record, range::UnitRange{Int})
    return length(range) == 1 && record.data[first(range)] == UInt8('.')
end



function vcfformat(val)
    return string(val)
end

function vcfformat(val::Vector)
    return join(map(vcfformat, val), ',')
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end
