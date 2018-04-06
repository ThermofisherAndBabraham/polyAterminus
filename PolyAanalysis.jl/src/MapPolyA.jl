#!/usr/bin/env julia

using ArgParse.ArgParseSettings
using ArgParse.@add_arg_table
using ArgParse.parse_args
using Bio.Align
using DataFrames
using DataStructures
using CSV
using Bio.Intervals

"""
"Turi buti: PolyA parsing is read name -> clustering -> statistics -> gtf, arba bed."
"Count total reads, aligned reads in intervals and general coverage."
"""
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

        if (flag&4 == 0) && (flag&256 == 0) && (flag&2048 == 0)  # filter out only mapped and primary alignments
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
                pclus = PolAClus(pclus, cps, parse(Int16, spl[1]))
            end
            passreads += 1
        end
    end

    return treads, passreads, pareads, pclus
end


function PolAClus(d::Dict{String,Array{Int16,1}}, p::String, l::Int16)

    if haskey(d, p)
        v = d[p]
        d[p] = push!(v, l)
    else
        d[p] = [l]
    end

    return d
end


function PolACalculus(d::Dict{String,Array{Int16,1}})
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


function WrFrame(fname::String, data::DataFrame)
    if !isfile(fname)
        CSV.write(fname, data)
    else
        CSV.write(fname, data, append=true, header=false)
    end
end


function ParseGFF3(gff3file::String)
    gffcol = IntervalCollection{Bio.Intervals.GFF3.Record}()
    # id = ""
    # geneid = ""
    # tic()
    # # attr = GetAttributes(gff3file)
    # toc()
    # geneatr = ""
    # gpatr = ""
    #
    # for i in attr
    #     if match(r"gene.i|Id|D", i) !== nothing
    #         geneatr = i
    #     end
    #
    #     if match(r"P|parent", i) !== nothing
    #         gpatr = i
    #     end
    # end

    # println(attr)
    #
    # if geneatr == "" || gpatr == ""
    #     println("Gene ID or Parent was not found in gff3 attributes, recomendation is to download gff3 from: https://www.gencodegenes.org/releases")
    #     exit(1)
    # end

    for record in open(GFF3.Reader, gff3file)
        id = GFF3.attributes(record, "ID")[1]
        geneid = GFF3.attributes(record, "gene_id")[1]
        push!(gffcol, Interval(record))
    end

    return gffcol
end


function GetAttributes(gff3file::String)

    s = Array{String,1}()

    for record in open(GFF3.Reader, gff3file)
        for i in GFF3.attributes(record)
            if first(i) in s
                continue
            else
                push!(s, first(i))
            end
        end
    end

    return s
end


function GetIntervalSet(dframe::DataFrame)
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
              string([4])*","*string(i[5])*","*string(i[6])*","*string(i[7])))
    end

    return intcol
end


function IntersectCollections(a::IntervalCollection, b::IntervalCollection)

    overlaps = Dict{String,Array{String,1}}()

    for x in a, y in b
        if isoverlapping(x, y)

            println(x)
            println(y)
            exit()
        end
    end

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

    for bam in parsed_args["bams"]
        treads, passreads, pareads, pclus = BamRead(bam)
        statdframe = PolACalculus(pclus)
        statdframe[:Smaple] = basename(bam)
        statdframe[:TotalReads] = treads
        statdframe[:PassedReads] = passreads
        statdframe[:PolyAReads] = pareads
        # println(statdframe)
        WrFrame(parsed_args["outfile"]*".csv", statdframe)
        intervalcolection = GetIntervalSet(statdframe)
        gffcollection = ParseGFF3(parsed_args["gff3file"])
        IntersectCollections(intervalcolection, gffcollection)
        # genes = ParseGFF3(parsed_args["gff3file"])

    end
end

main(ARGS)
