#!/usr/bin/env julia

"""
    get_transcripts_from_gff(fa::String, gff::String)

    Computes and returns array of transcripts from reference genome and gff3 file.

    # Arguments
    - `fa::String`: path to fasta file of reference genome
    - `gff::String`: path to gff3 transcripts file
"""
function get_transcripts_from_gff(fa::String, gff::String)::Array{FASTA.Record,1}

    # Reads Gff file generates dict chr => [transcriptsPositions]
    println("Reading GFF")
    reader = open(GFF3.Reader, gff)
    transdict = Dict{String,Array{String,1}}()
    ct = 0
    num = 0
    currtrans = ""
    coordex = ""
    strnd = ""
    chr = ""
    tic()

    for record in reader
        featureType = GFF3.featuretype(record)
        # Check if new feature found, then push exons
        if featureType == "transcript" || featureType == "mRNA"
            if currtrans != string(GFF3.attributes(record, "ID")[1]) && currtrans != ""
                if coordex == ""
                    println(STDERR,"ERROR! Exon is empty!")
                else
                    coordex = coordex * strnd
                end

                if !(in(coordex, transdict[chr]))
                    push!(transdict[chr],coordex)
                    num += 1
                else
                    ct += 1
                end
            end
            chr = string(GFF3.seqid(record))
            currtrans = string(GFF3.attributes(record, "ID")[1])
            strnd = string(GFF3.strand(record))
            # Check if chr exist in dict
            if !(haskey(transdict,chr))
                transdict[chr] = []
            end
            coordex = ""

        elseif featureType == "exon"
            strpos = string(GFF3.seqstart(record))
            endpos = string(GFF3.seqend(record))

            if currtrans == GFF3.attributes(record, "Parent")[1]
                coordex = coordex * strpos * "-" * endpos * ","
            else
                println(STDERR,"ERROR! Exon before transcript.")
                println(record)
            end
        end
    end

    # Adds If exon was last record in gff3 file.
    if coordex != ""
        coordex = coordex * strnd

        if !(in(coordex, transdict[chr]))
            push!(transdict[chr],coordex)
            num += 1
        else
            ct += 1
        end
    end

    toc()
    println(STDERR,"Dublicated transcripts count in gff = ",ct)
    println(STDERR,"Number of transcripts in gff = ",num)
    close(reader)

    println(STDERR,"Collecting sequences from ref genome using data from GFF")
    tic()
    # Read ref genome FASTA File and optains transcripts sequences
    file_stream = open(fa,"r")
    result = @parallel (vcat) for record in collect(FASTA.Reader(file_stream))
        get_transcripts_from_dict(record, transdict)
    end
    toc()
    close(file_stream)
    return result
end


"""
    get_transcripts_from_dict(record::BioSequences.FASTA.Record, transdict::Dict)

    Function returns transcripts list as an array of FASTA records.

    # Arguments
    - `record::BioSequences.FASTA.Record`: fasta record of sequence
    - `transdict::Dict`: dictionary with transcripts coordinates
"""
function get_transcripts_from_dict(record::FASTA.Record,
    transdict::Dict)::Array{FASTA.Record,1}

    outrecords = Array{FASTA.Record,1}()
    chr = FASTA.identifier(record)

    if haskey(transdict,chr)
        # Reads transcripts in Chromosome from transcripts dict
        d = transdict[chr]

        for trans in d
            transseq = dna""
            transcoordsandstrand = split(trans,",")
            transcoords = transcoordsandstrand[1:end-1]
            strand = transcoordsandstrand[end]

            for coords in transcoords
                coords = split(coords, "-")
                coor1 = parse(Int, coords[1])
                coor2 = parse(Int, coords[2])
                transseq = transseq * FASTA.sequence(record,coor1:coor2)
            end
            # If strand - makes reverse complement
            if strand == "-"
                transseq = reverse_complement!(transseq)
            end
            push!(outrecords, FASTA.Record("Seq", transseq))
        end
    else
        println(STDERR,"ERROR! Chromosome '",chr, "' is not in Gff file.")
    end
    return outrecords
end
