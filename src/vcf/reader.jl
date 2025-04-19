function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

mutable struct Reader <: BioGenerics.IO.AbstractReader
    state::BioGenerics.Automa.State
    index::Union{Indexes.Tabix,Nothing}
    header
    # The constructor now creates a Reader with an empty header.
    function Reader(stream::TranscodingStreams.TranscodingStream, index=nothing)
        # Header is assumed to be created by Header type constructor.
        hdr = Header()
        return new(BioGenerics.Automa.State(stream, 1, 1, false), index, hdr)
    end
end

# --- Constructors ---
"""
    VCF.Reader(input::IO; index=nothing)
    VCF.Reader(input::AbstractString; index=:auto)

Create a data reader of the VCF file format.

The first argument specifies the data source.

Arguments
---------
- `input`: data source
- `index`: path to a tabix file
"""
function Reader(input::IO; index=nothing)
    if isa(index, AbstractString)
        index = Indexes.Tabix(index)
    end
    # Wrap the raw IO in a NoopStream from TranscodingStreams.
    stream = TranscodingStreams.NoopStream(input)
    return Reader(stream, index)
end

function Reader(filepath::AbstractString; index=:auto)
    if isa(index, Symbol) && index != :auto
        throw(ArgumentError("invalid index argument: ':$(index)'"))
    end
    input = open(filepath)
    return Reader(input, index=index)
end

# --- Iterator Interface & Stream Getter ---
function Base.eltype(::Type{Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.iterate(reader::Reader, nextone=Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    # return a copy of the record and empty the record (to reuse its buffers)
    return copy(nextone), empty!(nextone)
end

# The header accessor – the Header type must provide the proper interface.
header(reader::Reader) = reader.header

# --- Automata Machine Constants and Actions ---
# (The full automa machinery is defined elsewhere.
#  Here we simply refer to the compiled machines.)
# VCF Reader Implementation
# -------------------------
# This complete VCF reader handles the meta-information header lines, the main
# header, and the record lines (data lines). We use a consistent marker (:pos)
# for record fields and separate markers (:mpos_key, :mpos_val) for meta-information
# to capture both key and value parts. This improves compactness and consistency.

# --- Automata Machine Constants and Regex Definitions ---
# VCF Reader Implementation
# -------------------------
# This complete VCF reader handles the meta-information header lines, the main
# header line, and the data records. We use a consistent marker (:pos) for capturing
# record fields and distinct markers for meta–information (e.g. :mpos_key, :mpos_val).
# The ALT field has been re‐defined to remove ambiguity by requiring that a missing ALT
# (a literal period) be immediately followed by a tab, newline, or end–of–input.

# --- Automata Machine Constants and Regex Definitions ---
const metainfo_machine, record_machine, header_machine, body_machine = let
    # Aliases for regex primitives from Automa.RegExp.
    ralt = Automa.RegExp.alt
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    opt = Automa.RegExp.opt
    delim(x, sep) = opt(cat(x, rep(cat(sep, x))))

    #### 1. FILEFORMAT Line ####
    fileformat = let
        key = onexit!(onenter!(cat("fileformat"), :mpos_key), :metainfo_tag)
        version = onexit!(onenter!(re"[!-~]+", :mpos_val), :metainfo_val)
        cat("##", key, '=', version)
    end
    onenter!(fileformat, :anchor)
    onexit!(fileformat, :metainfo)

    #### 2. META-INFORMATION Lines ####
    metainfo = let
        tag = onexit!(onenter!(re"[0-9A-Za-z_\.]+", :mpos_key), :metainfo_tag)
        simple_val = re"[ -;=-~][ -~]*"
        dict_val = let
            dictkey = onexit!(onenter!(re"[0-9A-Za-z_]+", :mpos_key), :metainfo_dict_key)
            dictval = onexit!(onenter!(ralt(cat('"', rep(ralt(re"[ !#-[\]-~]", "\\\"", "\\\\")), '"'),
                        rep(re"[ -~]" \ re"[\",>]")), :mpos_val), :metainfo_dict_val)
            cat('<', delim(cat(dictkey, '=', dictval), ','), '>')
        end
        local val = ralt(simple_val, dict_val)
        val = onexit!(onenter!(val, :mpos_val), :metainfo_val)
        cat("##", tag, '=', val)
    end
    onenter!(metainfo, :anchor)
    onexit!(metainfo, :metainfo)

    #### 3. HEADER Line ####
    header_line = let
        sampleID = onexit!(onenter!(re"[ -~]+", :mpos_key), :header_sampleID)
        cat("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", opt(re"\tFORMAT" * rep(re"\t" * sampleID)))
    end
    onenter!(header_line, :anchor)

    #### 4. RECORD (Data Line) Definition ####
    record = let
        # CHROM field: allowed printable characters (except those that conflict with delimiters).
        chrom = onexit!(onenter!(re"" | re"[^# \t\v\r\n\f][ -~]*", :pos), :record_chrom)
        # POS field: a number or the missing field (“.”).
        pos_field = onexit!(onenter!(re"[0-9]+|\.", :pos), :record_pos)
        # ID field: either missing “.” or a nonmissing id consisting of allowed characters plus a dot.
        id_field = let
            missing_id = onexit!(onenter!(re"\.", :pos), :record_id)
            nonmissing_id = onexit!(onenter!(cat(re"[!-:<-~]+", "."), :pos), :record_id)
            ralt(missing_id, nonmissing_id)
        end
        # REF field.
        ref_field = onexit!(onenter!(re"[!-~]+", :pos), :record_ref)
        # ALT field: either missing “.” or nonmissing allele.
        # To remove ambiguity, the branch for missing ALT (a literal period)
        # requires that it be immediately followed by a tab, newline, or end–of–input.
        alt_field = let
            alt_is_missing = onexit!(onenter!(re"\.", :pos), :record_alt)
            nonmissing_alt = onexit!(onenter!(re"[!-+--~]+", :pos), :record_alt)
            ralt(alt_is_missing, nonmissing_alt)
        end
        # QUAL field: numeric (including scientific), "NaN", ±Inf, or “.”
        qual = onexit!(onenter!(re"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?|NaN|[-+]Inf|\.", :pos), :record_qual)
        # FILTER field.
        filter_field = let
            missing_filter = onexit!(onenter!(re"\.", :pos), :record_filter)
            nonmissing_filter = onexit!(onenter!(cat(re"[!-:<-~]+", "."), :pos), :record_filter)
            ralt(missing_filter, nonmissing_filter)
        end
        # INFO field: a list of key[=value] pairs delimited by ";" or the missing field “.”
        info_field = let
            key = onexit!(onenter!(re"[A-Za-z_][0-9A-Za-z_.]*", :pos), :record_info_key)
            val = opt(cat('=', re"[ -:<-~]+"))
            ralt(delim(cat(key, val), ';'), re"\.")
        end
        # FORMAT field: similar structure with ':' as a delimiter.
        format_field = let
            elm = onexit!(onenter!(re"[A-Za-z_][0-9A-Za-z_.]*", :pos), :record_format)
            ralt(delim(elm, ':'), re"\.")
        end
        # GENOTYPE fields: sample genotype entries separated by ':'.
        genotype = let
            elm = onexit!(onenter!(re"[ -9;-~]+", :pos), :record_genotype_elm)
            delim(elm, ':')
        end
        onenter!(genotype, :record_genotype)
        # Compose the full record line: mandatory fields followed by optional FORMAT/GENOTYPE.
        cat(
            chrom, '\t',
            pos_field, '\t',
            id_field, '\t',
            ref_field, '\t',
            alt_field, '\t',
            qual, '\t',
            filter_field, '\t',
            info_field,
            opt(cat('\t', format_field, rep(cat('\t', genotype))))
        )
    end
    onenter!(record, :anchor)
    onexit!(record, :record)

    #### 5. Newline Pattern (supporting LF with an optional CR) ####
    newline = let
        lf = onenter!(re"\n", :countline)
        cat(opt('\r'), lf)
    end

    #### 6. Full Header and Body Definitions ####
    vcfheader = cat(
        fileformat, newline,
        rep(cat(metainfo, newline)),
        header_line, newline
    )
    onexit!(vcfheader, :vcfheader)
    vcfbody = rep(cat(record, newline))

    map(Automa.compile, (metainfo, record, vcfheader, vcfbody))
end


# --- Automata Actions for VCF Meta-information ---
const vcf_metainfo_actions = Dict(
    :mpos_key => :(@mark),
    :mpos_val => :(@mark),
    :anchor => :(),
    :metainfo_tag => :(record.tag = data[mpos_key:p-1]),
    :metainfo_val => :(record.val = data[mpos_val:p-1]; record.dict = (data[mpos_key] == UInt8('<'))),
    :metainfo_dict_key => :(push!(record.dictkey, data[mpos_key:p-1])),
    :metainfo_dict_val => :(push!(record.dictval, data[mpos_key:p-1])),
    :metainfo => quote
        copyto!(record.data, offset + 1, data, 1, p - offset - 1)
        record.filled = (offset+1):(p-1)
        @assert isfilled(record)
        push!(reader.header.metainfo, record)
        record = MetaInfo()
    end
)

# --- Automata Actions for VCF Data Record ---
const vcf_record_actions = Dict(
    #:mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :anchor => :(),
    :record_chrom => :(record.chrom = (1:@relpos(p - 1)-1); record.ncols += 1),
    :record_pos => :(record.pos = (pos-1:@relpos(p - 1)); record.ncols += 1),
    :record_id => :(push!(record.id, (pos:@relpos(p - 1))); record.ncols += 1),
    :record_ref => :(record.ref = (pos:@relpos(p - 1)); record.ncols += 1),
    :record_alt => :(push!(record.alt, (pos:@relpos(p - 1))); record.ncols += 1),
    :record_qual => :(record.qual = (pos:@relpos(p - 1)); record.ncols += 1),
    :record_filter => :(push!(record.filter, (pos:@relpos(p - 1))); record.ncols += 1),
    :record_info_key => :(push!(record.infokey, (pos:@relpos(p - 1))); record.ncols += 1),
    :record_format => :(push!(record.format, (pos:@relpos(p - 1))); record.ncols += 1),
    :record_genotype => :(push!(record.genotype, UnitRange{Int}[]); record.ncols += 1),
    :record_genotype_elm => :(push!(record.genotype[end], (pos:@relpos(p - 1)))),
    :record => quote
        copyto!(record.data, offset + 1, data, 1, p - offset - 1)
        record.filled = (offset+1):(p-1)
    end,
)

# --- Generate Index Readers ---
Automa.generate_reader(
    :index!,
    record_machine,
    arguments=(:(record::Record),),
    actions=vcf_record_actions,
    initcode=:(pos = 0; offset = 0)
) |> eval

Automa.generate_reader(
    :index!,
    metainfo_machine,
    arguments=(:(mi::MetaInfo),),
    actions=vcf_metainfo_actions,
    initcode=:(pos = 0; offset = 0)
) |> eval

# --- Full Record Reader Setup ---
const initcode = quote
    pos = 0            # initial input position
    linenum = 0        # initial line number
    found_record = false
    cs, linenum = state   # extract state from the stream state
end

const loopcode = quote
    if found_record
        @goto __return__
    end
end

Automa.generate_reader(
    :readrecord!,
    body_machine,
    arguments=(:(record::Record), :(state::Tuple{Int,Int})),
    actions=merge(vcf_record_actions, Dict(
        :record => quote
            appendfrom!(record.data, 1, data, @markpos, p - @markpos)
            record.filled = (1):(p-@markpos)
            found_record = true
            @escape
        end,
        :countline => :(linenum += 1),
        :anchor => :(anchor!(stream, p); offset = p - 1)
    )),
    initcode=initcode,
    loopcode=loopcode,
    returncode=:(return cs, linenum, found_record)
) |> eval

# --- VCF Record Indexing Function ---
function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    cs = index!(stream, record)
    if cs != 0
        throw(ArgumentError("Invalid VCF record. Machine failed to transition from state $(cs)."))
    end
    return record
end



"""
    read!(rdr::Reader, rec::Record)

Read a `Record` into `rec`, overwriting or adding to existing field values.
It is assumed that `rec` is already initialized or empty.
"""
function Base.read!(rdr::Reader, record::Record)
    cs, ln, found = readrecord!(rdr.state.stream, record, (rdr.state.state, rdr.state.linenum))
    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = found

    if found
        return record
    end

    if cs == 0 || eof(rdr.state.stream)
        throw(EOFError())
    end

    throw(ArgumentError("Malformed VCF record at line $(ln). Machine failed to transition from state $(cs)."))
end