#!/usr/bin/env julia


function readbam(bam::String)

    reader = open(BAM.Reader, bam)
    record = BAM.Record()
    treads = Int64(0)
    passreads = Int64(0)
    pareads = Int64(0)
    pclus = Dict{String,Array{Int16,1}}()

    while !eof(reader)
        read!(reader, record)
        if BAM.ismapped(record)
            rfname::String = BAM.refname(record)
            flag = BAM.flag(record)
            treads += 1

            # skip not primary aligned, supplementary alignment
            if (flag&256 == 0) && (flag&2048 == 0)
                spl = split(BAM.tempname(record),":")

                if flag&16 == 0
                    strness = "+"
                    pos = BAM.rightposition(record)
                elseif flag&16 == 16
                    strness = "-"
                    pos = BAM.position(record)
                else
                    strness = "."
                end

                if spl[2] == "A" && strness != "."
                    pareads += 1
                    cps = rfname * "::" * string(pos) * "::" * strness
                    pclus = clust!(pclus, cps, parse(Int16, spl[1]))
                end

                passreads += 1
            end
        end
    end

    return treads, passreads, pareads, pclus
end


function clust!(d::Dict, p::Any, l::Any)::Dict

    if haskey(d, p)
        v = d[p]
        d[p] = push!(v, l)
    else
        d[p] = [l]
    end

    return d
end


"""
    stats_poly_a(d::Dict{String,Array{Int16,1}})

    Collects simple statistics on polyA sites.
    Returns a dataframe.
"""
function stats_poly_a(d::Dict{String,Array{Int16,1}})::DataFrame
    dframe = DataFrame(Chrmosome=String[], Position=Int32[], Strand=String[], Median=Float16[], Minimum=Int16[], Maximum=Int16[], Counts=Int32[])

    for k in keys(d)
        md = median(d[k])
        mn = minimum(d[k])
        mx = maximum(d[k])
        l = length(d[k])
        push!(dframe, [split(k,"::")[1] parse(Int32,split(k,"::")[2]) split(k,"::")[3] md mn mx l])
    end

    return dframe
end


function wrframe(fname::String, data::DataFrame, delim::Char)
    CSV.write(fname, data, delim=delim)
end


function parseGFF3(gff3file::String)::Dict{String, IntervalCollection{String}}
    intcol = Dict{String, IntervalCollection{String}}()

    for record in open(GFF3.Reader, gff3file)
        chr::String = GFF3.seqid(record)
        start::Int = GFF3.seqstart(record)
        seqend::Int = GFF3.seqend(record)
        strand::Strand = GFF3.strand(record)
        split_k = get_split_key(chr, start, seqend)
        mdstr = parserecord(record)

        for i in split_k
            if haskey(intcol, i)
                push!(intcol[i], Interval(chr, start, seqend, strand, mdstr))
            else
                intcol[i] = IntervalCollection{String}()
                push!(intcol[i], Interval(chr, start, seqend, strand, mdstr))
            end
        end
    end

    return intcol
end

"""
    get_split_key(chr::String, x::Int64, y::Int64; step::Int64=10000)

    Calculates ranges depending on the step and gff3 coordinates.
"""
function get_split_key(chr::String, x::Int64, y::Int64; step::Int64=10000)::Array{String,1}

    a = Array{String,1}()

    function get_interval(x::Int64,y::Int64)
        if x >= y
            if x%y == 0
                ix1 = x - y + 1
                ix2 = x
            else
                ix1 = x - x%y + 1
                ix2 = x - x%y + y
            end
        else
            ix1 = 1
            ix2 = y
        end

        return ix1:ix2
    end

    interv = get_interval(x, step)

    if y in interv
        push!(a, chr*":"*string(interv))
    else
        push!(a, chr*":"*string(interv))
        interv2 = get_interval(y, step)

        if interv[end] == interv2[1]
            push!(a, chr*":"*string(interv2))
        else
            for i in interv[end]:step:interv2[1]-1
                ix1=i+1
                ix2=i+step
                push!(a, chr*":"*string(ix1:ix2))
            end
        end
    end

    return a
end


function parserecord(r::GFF3.Record)::String

    md::String = ""

    try s = GFF3.attributes(r, "ID")
        md = string(md,"ID=")
        for x in s
            if x == s[end]
                md = string(md,x)
            else
                md = string(md,x,",")
            end
        end
        md = string(md,";")
    catch   md = string(md,"ID=;")
    end

    try s = GFF3.attributes(r, "Parent")
        md = string(md,"Parent=")
        for x in s
            if x == s[end]
                md = string(md,x)
            else
                md = string(md,x,",")
            end
        end
        md = string(md,";")
    catch   md = string(md,"Parent=;")
    end

    try s = GFF3.attributes(r, "gene_id")
        md = string(md,"gene_id=")
        for x in s
            if x == s[end]
                md = string(md,x)
            else
                md = string(md,x,",")
            end
        end
        md = string(md,";")
    catch   md = string(md,"gene_id=;")
    end

    try s = GFF3.attributes(r, "transcript_id")
        md = string(md,"transcript_id=")
        for x in s
            if x == s[end]
                md = string(md,x)
            else
                md = string(md,x,",")
            end
        end
        md = string(md,";")
    catch   md = string(md,"transcript_id=;")
    end

    try s = GFF3.attributes(r, "gene_type")
        md = string(md,"gene_type=")
        for x in s
            if x == s[end]
                md = string(md,x)
            else
                md = string(md,x,",")
            end
        end
        md = string(md,";")
    catch   md = string(md,"gene_type=;")
    end

    try s = GFF3.attributes(r, "gene_name")
        md = string(md,"gene_name=")
        for x in s
            if x == s[end]
                md = string(md,x)
            else
                md = string(md,x,",")
            end
        end
        md = string(md,";")
    catch   md = string(md,"gene_name=;")
    end

    md = string(md,"ftype=",GFF3.featuretype(r))

    return md
end


"""
    getintervals(dframe::DataFrame)

    Returns a collection of intervals from a polyA sites dataframe.
"""
function getintervals(dframe::DataFrame)::IntervalCollection{String}

    intcol = IntervalCollection{String}()

    for i in eachrow(dframe)
        if i[3] == "+"
            s = Strand('+')
        elseif i[3] == "-"
            s = Strand('-')
        elseif i[3] == "."
            s = Strand('.')
        else
            s = Strand('?')
        end

        push!(intcol, Interval(i[1], i[2], i[2], s,
              "Median="*string(i[4])*";"*"Min="*string(i[5])*
              ";"*"Max="*string(i[6])*";"*"Counts="*string(i[7])))
    end

    return intcol
end


"""
    annotate_polya_sites(a::IntervalCollection, b::Dict{String, IntervalCollection{String}})

    Function searches for overlaping intervals. If overlapped - adds annotation
    to the dataframe. If not - add NA.

    # Arguments
    - `a::IntervalCollection`: a collection of intervals of PolyA sites.
    - `b::Dict{String, IntervalCollection{String}}`: a dictionary with a sorted collections of intervals from gff3.
"""
function annotate_polya_sites(a::IntervalCollection, b::Dict{String, IntervalCollection{String}})::DataFrame

    # Two dataframes for '+' and '-' strands.
    dfp = DataFrame(Chr=String[], Start=Int64[], End=Int64[], Name=String[],
                    Counts=Int32[], Strand=String[], Feature=String[],
                    Median=Float32[], Min=Int16[], Max=Int16[], Biotype=String[])
    dfn = DataFrame(Chr=String[], Start=Int64[], End=Int64[], Name=String[],
                    Counts=Int32[], Strand=String[], Feature=String[],
                    Median=Float32[], Min=Int16[], Max=Int16[], Biotype=String[])
    # Dataframe for not annotated sites.
    dfna = DataFrame(Chr=String[], Start=Int64[], End=Int64[], Name=String[],
                    Counts=Int32[], Strand=String[], Feature=String[],
                    Median=Float32[], Min=Int16[], Max=Int16[], Biotype=String[])

    for i1 in a
        ct = Int64(0)
        # Parse metadata from polyA site.
        splt1::Array{SubString{String},1} = split(metadata(i1),";")
        c = parse(Int32,split(splt1[4],"=")[2])
        str = string(strand(i1))
        me = parse(Float32,split(splt1[1],"=")[2])
        mi = parse(Int16,split(splt1[2],"=")[2])
        mx = parse(Int16,split(splt1[3],"=")[2])
        chr = seqname(i1)
        it1 = first(i1)
        it2 = last(i1)

        # Assuming that start and end is the same of the polyA site.
        split_k = get_split_key(chr, it1, it2)[1]

        try
            for i2 in b[split_k]
                if chechoverlap(i1, i2)
                    # parse metadata from GFF3 anotation
                    splt::Array{SubString{String},1} = split(metadata(i2),";")
                    ft = split(splt[7],"=")[2]
                    gn = split(splt[6],"=")[2]
                    bt = split(splt[5],"=")[2]
                    ct += 1
                    if strand(i1) == Strand('-')
                        push!(dfn, [seqname(i1) it1 it2 gn c str ft me mi mx bt])
                    else
                        push!(dfp, [seqname(i1) it1 it2 gn c str ft me mi mx bt])
                    end
                end
            end
        catch
            KeyError(split_k)
        end

        # if no overlaping of polyA coordinates with GFF3 annotated coordinates - NA is added.
        if ct == 0
            ft = "NA"
            gn = "NA"
            bt = "NA"

            push!(dfna, [seqname(i1) first(i1) last(i1) gn c str ft me mi mx bt])
        end
    end

    sort!(dfn, [:Chr, :Start, :End, :Name], rev=true)
    sort!(dfp, [:Chr, :Start, :End, :Name])
    sort!(dfna, [:Chr, :Start, :End, :Name])

    dfp = rmdups(dfp)
    dfn = rmdups(dfn)

    enumeratenames!(dfn)
    enumeratenames!(dfp)
    enumeratenames!(dfna)

    df = vcat(dfp, dfn, dfna)

    return df
end


"""
    chechoverlap(i1::Interval, i2::Interval)

    Checks if intervals overlap.
"""
function chechoverlap(i1::Interval, i2::Interval)::Bool

    # only if chr match and strand match
    if cmp(seqname(i1), seqname(i2)) == 0 && strand(i1) == strand(i2)
        if i2.first <= i1.last && i2.first >= i1.first
            return true
        elseif i2.last <= i1.last && i2.last >= i1.first
            return true
        elseif i1.first <= i2.last && i1.first >= i2.first
            return true
        elseif i1.last >= i2.first && i1.last <= i2.last
            return true
        else
            return false
        end
    else
        return false
    end

end


"""
    enumeratenames!(df::DataFrame)

    Enumerates gene names in the data frame by adding incrementing number.
"""
function enumeratenames!(df::DataFrame)::DataFrame

    ct = Int64(0)
    l = size(df)[1]
    na_ct = Int64(0)

    for (i, x) in enumerate(eachrow(df))
        ct += 1

        if x[4] == "NA" && i+1 <= l
            na_ct += 1
            x[4] = x[4]*".$na_ct"
        elseif i+1 <= l
            if x[4] == eachrow(df)[i+1][4]
                x[4] = x[4]*".$ct"
            else
                x[4] = x[4]*".$ct"
                ct = 0
            end
        else
            x[4] = x[4]*".$ct"
            ct = 0
        end
    end

    return df
end


"""
    rmdups(dframe::DataFrame)

    Removes duplicated entries differentiating by feature.
    DataFrame should be sorted by Chr, Start, End and gene name.
    Featues are rated by their annotation hierarchy.
    Keeps most accurate feature (gene has exons and introns, etc.).
"""
function rmdups(dframe::DataFrame)::DataFrame

    df = DataFrame(Chr=String[], Start=Int64[], End=Int64[], Name=String[],
                   Counts=Int32[], Strand=String[], Feature=String[],
                   Median=Float32[], Min=Int16[], Max=Int16[], Biotype=String[])

    d = Dict("gene"=>13, "transcript"=>12,"exon"=>11, "intron"=>10, "CDS"=>9,
             "UTR"=>8, "three_prime_UTR"=>7, "five_prime_UTR"=>6,
             "polyA_sequence"=>5, "polyA_site"=>4, "start_codon"=>3,
             "stop_codon"=>2, "Selenocysteine"=>1, "NA"=>0)

    ct = Int64(0)
    ct2 = Int64(0)

    for y in eachrow(dframe[:,1:6])
        ct2 += 1

        if ct == 0
            push!(df, [dframe[1][ct2] dframe[2][ct2] dframe[3][ct2] dframe[4][ct2] dframe[5][ct2] dframe[6][ct2] dframe[7][ct2] dframe[8][ct2] dframe[9][ct2] dframe[10][ct2] dframe[11][ct2]])
            ct += 1
        end

        if ct > 0
            if df[end,:][1] == [y[1]] && df[end,:][2] == [y[2]] && df[end,:][3] == [y[3]] && df[end,:][4] == [y[4]] && df[end,:][5] == [y[5]] && df[end,:][6] == [y[6]]

                if d[dframe[7][ct2]] < d[df[end,:][7][1]]
                    df[7][ct] = dframe[7][ct2]
                end

            else
                ct += 1
                push!(df, [dframe[1][ct2] dframe[2][ct2] dframe[3][ct2] dframe[4][ct2] dframe[5][ct2] dframe[6][ct2] dframe[7][ct2] dframe[8][ct2] dframe[9][ct2] dframe[10][ct2] dframe[11][ct2]])
            end
        end
    end

    return df
end
