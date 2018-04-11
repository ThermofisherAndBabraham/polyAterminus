#!/usr/bin/env julia

using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args
using Bio.Align
using DataFrames
using DataStructures
using CSV
using Bio.Intervals


function BamRead(bam::String)

    reader = open(Bio.Align.BAMReader, bam)
    record = Bio.Align.BAMRecord()
    treads = Int64(0)
    passreads = Int64(0)
    pareads = Int64(0)
    pclus = Dict{String,Array{Int16,1}}()

    while !eof(reader)
        read!(reader, record)
        rfname::String = refname(record)
        flag = Bio.Align.flag(record)
        treads += 1

        if (flag&4 == 0) && (flag&256 == 0) && (flag&2048 == 0)
            spl = split(seqname(record),":")

            if flag&16 == 0
                strness = "+"
            elseif flag&16 == 16
                strness = "-"
            else
                strness = "."
            end

            if spl[2] == "A"
                pareads += 1
                rpos = rightposition(record)
                cps = rfname * "::" * string(rpos) * "::" * strness
                pclus = Clust(pclus, cps, parse(Int16, spl[1]))
            end
            passreads += 1
        end
    end

    return treads, passreads, pareads, pclus
end


function Clust(d::Dict, p::Any, l::Any)::Dict

    if haskey(d, p)
        v = d[p]
        d[p] = push!(v, l)
    else
        d[p] = [l]
    end

    return d
end


function PolACalculus(d::Dict{String,Array{Int16,1}})::DataFrame
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


function WrFrame(fname::String, data::DataFrame, delim::Char)
    if !isfile(fname)
        CSV.write(fname, data, delim=delim)
    else
        CSV.write(fname, data, delim=delim, append=true, header=false)
    end
end


function ParseGFF3(gff3file::String)::IntervalCollection{String}
    intcol = IntervalCollection{String}()

    for record in open(GFF3.Reader, gff3file)
        chr::String = GFF3.seqid(record)
        start::Int = GFF3.seqstart(record)
        seqend::Int = GFF3.seqend(record)
        strand::Bio.Intervals.Strand = GFF3.strand(record)

        mdstr = ParseRecord(record)
        push!(intcol, Interval(chr, start, seqend, strand, mdstr))
    end

    return intcol
end


function ParseRecord(r::GFF3.Record)::String

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


function GetIntervalSet(dframe::DataFrame)::IntervalCollection{String}

    intcol = IntervalCollection{String}()

    for i in eachrow(dframe)
        if i[3] == "+"
            s = Strand('+')
        elseif i[3] == '-'
            s = Strand('-')
        elseif i[3] == '.'
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


function ItsectCollection(a::IntervalCollection, b::IntervalCollection)::DataFrame

    dfp = DataFrame(Chr=String[], Start=Int64[], End=Int64[], Name=String[], Counts=Int32[], Strand=String[], Feature=String[], Median=Float32[], Min=Int16[], Max=Int16[], Biotype=String[])
    dfn = DataFrame(Chr=String[], Start=Int64[], End=Int64[], Name=String[], Counts=Int32[], Strand=String[], Feature=String[], Median=Float32[], Min=Int16[], Max=Int16[], Biotype=String[])

    for x in a, y in b
        if isoverlapping(x, y)
            if strand(x) == strand(y)
                geneID = split(split(metadata(y),";")[3],"=")[2]
                gn = split(split(metadata(y),";")[6],"=")[2]
                c = parse(Int32,split(split(metadata(x),";")[4],"=")[2])
                str = string(strand(x))
                ft = split(split(metadata(y),";")[7],"=")[2]
                me = parse(Float32,split(split(metadata(x),";")[1],"=")[2])
                mi = parse(Int16,split(split(metadata(x),";")[2],"=")[2])
                mx = parse(Int16,split(split(metadata(x),";")[3],"=")[2])
                bt = split(split(metadata(y),";")[5],"=")[2]
                if strand(x) == Strand('-')
                    push!(dfn, [seqname(x) first(x) last(x) gn c str ft me mi mx bt])
                else
                    push!(dfp, [seqname(x) first(x) last(x) gn c str ft me mi mx bt])
                end
            end
        end
    end

    sort!(dfn, rev=true, cols = [:Start, :End, :Name])
    sort!(dfp, cols = [:Start, :End, :Name])

    dfp = rmdups!(dfp)
    dfn = rmdups!(dfn)

    ct = Int64(0)
    l = size(dfn)[1]

    for (i, x) in enumerate(eachrow(dfn))
        ct += 1
        if i+1 <= l
            if x[4] == eachrow(dfn)[i+1][4]
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

    ct = Int64(0)
    l = size(dfp)[1]

    for (i, x) in enumerate(eachrow(dfp))
        ct += 1
        if i+1 <= l
            if x[4] == eachrow(dfp)[i+1][4]
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

    df = vcat(dfp, dfn)

    return df
end


function rmdups!(dframe::DataFrame)::DataFrame
    allrows = eachrow(dframe)
    ct = Int32(0)

    for (i, x) in enumerate(eachrow(dframe))
        ct += 1
        while i+ct <= size(dframe)[1]
            while i+ct <= size(dframe)[1] && x[1] == allrows[i+ct][1] && x[2] == allrows[i+ct][2] && x[3] == allrows[i+ct][3] && x[4] == allrows[i+ct][4]
                if x[7] == "stop_codon"
                    deleterows!(dframe, i+ct)
                elseif x[7] == "start_codon"
                    deleterows!(dframe, i+ct)
                elseif x[7] == "CDS"
                    deleterows!(dframe, i+ct)
                elseif x[7] == "exon"
                    deleterows!(dframe, i+ct)
                elseif x[7] == "transcript"
                    deleterows!(dframe, i+ct)
                elseif x[7] == "gene"
                    deleterows!(dframe, i+ct)
                end
            end

            ct += 1
        end
        ct = 0
    end

    return dframe
end

function main(args)
    arg_parse_settings = ArgParseSettings(description="Calculates uncovered windows number and GC bias.")
    @add_arg_table arg_parse_settings begin
        "--bam", "-b"
            nargs = '+'
            arg_type = String
            help = "Bam files separated by spaces."
            required = true
            dest_name = "bams"
        "--out", "-o"
            arg_type = String
            help = "Prefix of output file."
            required = false
            dest_name = "outfile"
            default = "output"
        "--gff", "-g"
            arg_type = String
            help = "GFF3 file"
            required = true
            dest_name = "gff3file"
    end

    parsed_args = parse_args(arg_parse_settings)

    gffcollection = ParseGFF3(parsed_args["gff3file"])

    for bam in parsed_args["bams"]
        treads, passreads, pareads, pclus = BamRead(bam)
        statdframe = PolACalculus(pclus)
        statdframe[:Smaple] = basename(bam)
        statdframe[:TotalReads] = treads
        statdframe[:PassedReads] = passreads
        statdframe[:PolyAReads] = pareads
        WrFrame(parsed_args["outfile"]*"detected_polyA.tsv", statdframe, '\t')
        intervalcolection = GetIntervalSet(statdframe)
        joinedcollection = ItsectCollection(intervalcolection, gffcollection)
        WrFrame(parsed_args["outfile"]*"_mapped_polyA.bed", joinedcollection, '\t')
    end
end

main(ARGS)
