# Automa.jl generated readrecord! and readmetainfo! functions
# ========================================

import Automa: @re_str, rep, onexit!, onenter!, CodeGenContext, generate_reader, compile, RegExp

# --- Automata Machine Constants and Regex Definitions ---
const vcf_machine_metainfo, vcf_machine_record, vcf_machine_header, vcf_machine_body, vcf_machine = let
    # Aliases for regex primitives from Automa.RegExp.
    ralt = RegExp.alt
    cat = RegExp.cat
    rep = RegExp.rep
    opt = RegExp.opt
    delim(x, sep) = opt(cat(x, rep(cat(sep, x))))
    newline = "\r?" * onenter!(re"\n", :countline)


    #### 1. FILEFORMAT Line ####
    fileformat = let
        key = onexit!(onenter!(cat("fileformat"), :pos1), :metainfo_tag)
        version = onexit!(onenter!(re"[!-~]+", :pos2), :metainfo_val)
        cat("##", key, '=', version)
    end
    onexit!(onenter!(fileformat, :mark), :metainfo)

    #### 2. META-INFORMATION Lines ####
    metainfo = let
        tag = onexit!(onenter!(re"[0-9A-Za-z_\.]+", :pos1), :metainfo_tag)
        simple_val = re"[ -;=-~][ -~]*"
        dict_val = let
            dictkey = onexit!(onenter!(re"[0-9A-Za-z_]+", :pos1), :metainfo_dict_key)
            # Removed the extra onenter! here to avoid overwriting pos2.
            dictval = onexit!(ralt(
                    cat('"', rep(ralt(re"[ !#-[\]-~]", "\\\"", "\\\\")), '"'),
                    rep(re"[ -~]" \ re"[\",>]")
                ), :metainfo_dict_val)
            cat('<', delim(cat(dictkey, '=', dictval), ','), '>')
        end
        local val = ralt(simple_val, dict_val)
        # Wrap the combined value with onenter!(…, :pos2)
        val = onexit!(onenter!(val, :pos2), :metainfo_val)
        cat("##", tag, '=', val)
    end
    onexit!(onenter!(metainfo, :mark), :metainfo)

    # 3. Header line with sample IDs
    header_line = let
        sampleID = onexit!(onenter!(re"[^\t\r\n]+", :mark_sampleid), :header_sampleID)
        cat(
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            opt(cat("\tFORMAT", rep(cat('\t', sampleID))))
        )
    end
    onenter!(header_line, :mark)

    record = let
        # CHROM field: allowed printable characters (except those that conflict with delimiters).
        chrom = onexit!(onenter!(re"[!-9;-~]+", :pos), :record_chrom)
        # POS field: a number or the missing field (“.”).
        pos_field = onexit!(onenter!(re"[0-9]+|\.", :pos), :record_pos)
        # ID field: either missing “.” or a nonmissing id consisting of allowed characters plus a dot.
        id_field = let
            missing_id = onexit!(onenter!(re"\.", :pos), :record_id)
            nonmissing_id = onexit!(onenter!(re"[!-:<-~]+", :pos), :record_id)
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
            nonmissing_filter = onexit!(onenter!(re"[!-:<-~]+", :pos), :record_filter)
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
        # GENOTYPE fields: vcfple genotype entries separated by ':'.
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
    onexit!(onenter!(record, :mark), :record)

    header = onexit!(
        cat(fileformat, newline, rep(metainfo * newline), header_line, newline),
        :header
    )

    body = onexit!(rep(record * newline), :body)
    vcf = header * body

    map(compile, (metainfo, record, header, body, vcf))
end

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

# --- Automata Actions for VCF Meta-information ---
const vcf_actions_metainfo = Dict(
    :mark => :(@mark),
    :pos1 => :(pos1 = @relpos(p)),
    :pos2 => :(pos2 = @relpos(p)),
    :metainfo_tag => :(metainfo.tag = pos1:@relpos(p - 1)),
    :metainfo_val => :(metainfo.val = pos2:@relpos(p - 1); metainfo.dict = data[pos2] == UInt8('<')),
    :metainfo_dict_key => :(push!(metainfo.dictkey, pos1:@relpos(p - 1))),
    :metainfo_dict_val => :(push!(metainfo.dictval, pos1:@relpos(p - 1))),
    :metainfo => quote
        appendfrom!(metainfo.data, 1, data, @markpos, p - @markpos)
        metainfo.filled = 1:(p-@markpos)
    end
)


const vcf_actions_header = merge(
    vcf_actions_metainfo,
    Dict(
        :countline => :(linenum += 1),
        :metainfo => quote
            $(vcf_actions_metainfo[:metainfo])
            push!(header, metainfo)
            metainfo = MetaInfo()
        end,
        :mark_sampleid => :(@mark),  # Mark position where sample ID actually starts
        :header_sampleID => :(push!(header.sampleID, String(data[@markpos():p-1]))),
        :header => :(@escape)
    )
)
# --- Automata Actions for VCF Data Record ---
const vcf_actions_record = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :record_chrom => :(record.chrom = (pos:@relpos(p - 1)); record.ncols += 1),
    :record_pos => :(record.pos = (pos:@relpos(p - 1)); record.ncols += 1),
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
        appendfrom!(record.data, 1, data, @markpos, p - @markpos)
        record.filled = 1:(p-@markpos)
    end
)

const vcf_actions_body = merge(
    vcf_actions_record,
    Dict(
        :countline => :(linenum += 1),
        :record => quote
            found_record = true
            $(vcf_actions_record[:record])
            @escape
        end,
        :body => :(@escape)
    )
)

const vcf_context = CodeGenContext(generator=:goto)

const vcf_initcode_metainfo = quote
    pos1 = 0
    pos2 = 0
end

const vcf_initcode_header = quote
    $(vcf_initcode_metainfo)
    metainfo = MetaInfo()
    cs, linenum = state
end

const vcf_initcode_record = quote
    pos = 0
end

const vcf_initcode_body = quote
    $(vcf_initcode_record)
    found_record = false
    cs, linenum = state
end

generate_reader(
    :index!,
    vcf_machine_metainfo,
    arguments=(:(metainfo::MetaInfo),),
    actions=vcf_actions_metainfo,
    context=vcf_context,
    initcode=vcf_initcode_metainfo,
) |> eval

const vcf_returncode_header = quote
    return cs, linenum
end

generate_reader(
    :readheader!,
    vcf_machine_header,
    arguments=(:(header::VCF.Header), :(state::Tuple{Int,Int})),
    actions=vcf_actions_header,
    context=vcf_context,
    initcode=vcf_initcode_header,  # using the updated init code without offset
    returncode=vcf_returncode_header,
    errorcode=quote
        # Accept states -2, -1, or 0 as a signal for proper header termination.
        if cs in (-58,-2, -1, 0)

            @goto __return__
        else
            error("Expected input byte after vcf header, got state $(cs)")
        end
    end
) |> eval

const vcf_loopcode_body = quote
    if found_record
        @goto __return__
    end
end

generate_reader(
    :index!,
    vcf_machine_record,
    arguments=(:(record::Record),),
    actions=vcf_actions_record,
    context=vcf_context,
    initcode=vcf_initcode_record,
) |> eval


const vcf_returncode_body = quote
    return cs, linenum, found_record
end

generate_reader(
    :readrecord!,
    vcf_machine_body,
    arguments=(:(record::Record), :(state::Tuple{Int,Int})),
    actions=vcf_actions_body,
    context=vcf_context,
    initcode=vcf_initcode_body,
    loopcode=vcf_loopcode_body,
    returncode=vcf_returncode_body,
    errorcode=quote
        if cs in (-11, -2, -1, 0)
            @goto __return__
        else
            error("Malformed VCF file body. Machine failed to transition from state $(cs).")
        end
    end
) |> eval