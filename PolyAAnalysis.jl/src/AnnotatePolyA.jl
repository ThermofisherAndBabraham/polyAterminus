#!/usr/bin/env julia


function readbam(bam::String, strandness::String, k::Int64, cluster::Bool; verbose=false)

    reader = open(BAM.Reader, bam)
    record = BAM.Record()
    treads = Int64(0)
    passreads = Int64(0)
    pareads = Int64(0)
    pclus = Dict{String,Array{Int16,1}}()
    pclus2 = Dict()

    if !(strandness in ["+" "-"])
        println(STDERR, "ERROR: provided srandness is not correct: $strandness, [+/-].")
        exit(1)
    end

    while !eof(reader)
        read!(reader, record)

        if BAM.ismapped(record)
            rfname::String = BAM.refname(record)
            flag = BAM.flag(record)
            treads += 1

            iftrue::Bool = false

            # if pair end
            if flag&1 == 1
                # if strandness +, take first read
                if strandness == "+" && (flag&64==64)
                    iftrue = true
                # if strandness -, take second read
                elseif strandness == "-" && (flag&128==128)
                    iftrue = true
                else
                    iftrue = false
                end
            else
                iftrue = true
            end

            # skip not primary aligned, supplementary alignment.
            if (flag&256 == 0) && (flag&2048 == 0) && iftrue
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

                    if cluster
                        pclus2 = cluster!(pclus2, rfname, pos, [parse(Int64, spl[1])],
                                          strness, k; verbose=verbose)
                    end
                end

                passreads += 1
            end
        end
    end

    return treads, passreads, pareads, pclus, pclus2
end


"""
    cluster!(d::Dict, chr::String, pos::Int64, l::Array{Int64,1},
             str::String, k::Int64; verbose::Bool=false)::Dict

    Function adds new termination site (TS) to a given dataset by searching most
    adjacent clusters and if needed correcting clusters to meat definition of
    +/-k from center.
    Algorithm:
        1. Search for a cluster center.
        2. If !found - add new cluster of one.
        3. Else
            3.1. Add value to the cluster
            3.2. Check center
                3.2.1 If !changed - return
                    3.2.1.1 Check adjacent clusters, comapare values.
                    3.2.1.2 Recenter recursively if needed.
                3.2.2. If changed
                    3.2.2.1. Recenter
                    3.2.2.2. Check if all values holds in
                        3.2.2.2.1. If holds in - return
                        3.2.2.2.2. Else drop those values
                                   and call function recursively.
                    3.2.2.3 Check adjacent clusters, comapare values.
                    3.2.2.4 Recenter recursively if needed.

    # Arguments
    - `d::Dict`: a collection of clusters where key is combination of chromosome
       mapping position.
    - `chr::String`: Chromosome.
    - `pos::Int64`: Mapping position of TS.
    - `l::Array{Int64,1}`: Lengths of detected PolyA. First values is count.
    - `str::String`: Strand
    - `k::Int64`: k - distance from center allowed.
    - `verbose::Bool=false`: print clustering proceeding to STDOUT.
"""
function cluster!(d::Dict, chr::String, pos::Int64, l::Array{Int64,1},
                  str::String, k::Int64; verbose::Bool=false)::Dict

    if verbose
        println(STDOUT, "SEARCHING FOR: ", chr*"::"*string(pos)*"::"*str, " TS in +/- $k")
    end

    # search if TS fall into any clusters in range of +/- k, rember keys if so.
    a = Array{String,1}()
    v = Array{Int,1}()
    found::Bool = false

    for i in pos-k:pos+k
        if i >= 0
            key = chr*"::"*string(i)*"::"str
            if key in keys(d)
                push!(a, key)
                push!(v, abs(pos-i))
                found = true
                if verbose
                    println(STDOUT, "FOUND EXISTING IN +/- $k RANGE CLUSTER: ", key)
                end
            end
        end
    end

    if !found
        if verbose
            println(STDOUT, "NO EXISTING CLUSTERS FOUND: ", chr*"::"*string(pos)*"::"*str)
            println(STDOUT, "CREATE NEW. CHECK ADJ.")
        end

        # in some cases with recursion l might be already with frequencies.
        if length(l) > 1
            d[chr*"::"*string(pos)*"::"*str] = Dict(pos => l)
            search_adj_clust!(d, pos, chr*"::"*string(pos)*"::"*str,
                              chr, str, k; verbose=verbose)
            @goto stop_cluster!
        else
            d[chr*"::"*string(pos)*"::"*str] = Dict(pos => vcat([1], l))
            search_adj_clust!(d, pos, chr*"::"*string(pos)*"::"*str,
                              chr, str, k; verbose=verbose)
            @goto stop_cluster!
        end
    end

    minv = k + 1
    near_center = ""
    ct = 0

    # find which is the nearest cluster if it is and remember its key.
    for i in v
        ct += 1
        if i < minv
            minv = i
            near_center = a[ct]
            if verbose
                println(STDOUT, "NEAREST CLUSTER FOUND TO BE: ", a[ct])
            end

        elseif i == minv
            # means that this TS is in the midle between two clusters
            if verbose
                println(STDOUT, "POSITION TO BE IN THE MIDLE OF 2 CLUSTERS: ",
                        chr*"::"*string(pos - i)*"::"*"str",
                        " --> ",
                        chr*"::"*string(pos)*"::"*"str",
                        " <-- ",
                        chr*"::"*string(i + pos)*"::"*"str"
                        )
            end

            # try to decide by weight of cluster.
            weight1 = 0
            for x in keys(d[chr*"::"*string(pos-i)*"::"*str])
                weight1 += d[chr*"::"*string(pos-i)*"::"*str][x][1]
            end

            weight2 = 0
            for x in keys(d[chr*"::"*string(pos+i)*"::"*str])
                weight2 += d[chr*"::"*string(pos+i)*"::"*str][x][1]
            end

            if weight1 > weight2
                near_center = chr*"::"*string(pos-i)*"::"*str
            elseif weight2 > weight1
                near_center = chr*"::"*string(pos+i)*"::"*str

            else
                # in some cases with recursion l might be already with frequencies.
                if length(l) > 1
                    d[chr*"::"*string(pos)*"::"*str] = Dict(pos => l)
                    search_adj_clust!(d, pos, chr*"::"*string(pos)*"::"*str,
                                      chr, str, k; verbose=verbose)
                    @goto stop_cluster!

                else
                    d[chr*"::"*string(pos)*"::"*str] = Dict(pos => vcat([1], l))
                    search_adj_clust!(d, pos, chr*"::"*string(pos)*"::"*str,
                                      chr, str, k; verbose=verbose)
                    @goto stop_cluster!
                end
            end
        end
    end

    # search in a given cluster
    # check if TS is in this cluster
    # else add new
    # check center and recenter if needed.

    if pos in keys(d[near_center])
        if verbose
            println(STDOUT, "FOUND SAME TS IN THE CLUSTER: ", pos, " in ", near_center)
        end

        d[near_center][pos][1] += 1
        d[near_center][pos] = vcat(d[near_center][pos], l)
        recenter!(d, d[near_center], near_center, k, chr, str; verbose=verbose)

    else
        if verbose
            println(STDOUT, "NO TS FOUND IN THE CLUSTER, ADDING NEW: ", pos, " in ", near_center)
        end

        # in some cases with recursion l might be already with frequencies.
        if length(l) > 1
            d[near_center][pos] = l
            recenter!(d, d[near_center], near_center, k, chr, str; verbose=verbose)
        else
            d[near_center][pos] = vcat([1], l)
            recenter!(d, d[near_center], near_center, k, chr, str; verbose=verbose)
        end
    end

    @label stop_cluster!

    return d
end


function recenter!(d1::Dict, d2::Dict, near_center::String, k::Int64,
                   chr::String, str::String; verbose::Bool=false)::Dict

    # find new most frequent element
    center_pos = parse(Int64, split(near_center, "::")[2])

    if center_pos in keys(d2)
        center_freq = d2[center_pos][1]
    else
        if verbose
            println(STDOUT, "CENTER NOT FOUND, MUST HAVE CHANGED IN RECURSION: $center_pos")
        end
        @goto stop_recenter!
    end

    new_center = 0
    new_center_freq = 0
    found_new_center::Bool = false
    equals = Int64[]
    all_positions = Int64[]

    if verbose
        println(STDOUT, "CHECK CENTER IF PERSISTS.")
    end

    for position in keys(d2)
        push!(all_positions, position)

        if position != center_pos
            if !found_new_center && d2[position][1] > center_freq
                if verbose
                    println(STDOUT, "PREVIOUS POSITION: $center_pos with FREQUENCY: $center_freq")
                    println(STDOUT, "NEW POSITION: $position with FREQUENCY: ", d2[position][1])
                end
                new_center = position
                new_center_freq = d2[position][1]
                found_new_center = true
                equals = []

            elseif !found_new_center && d2[position][1] == center_freq
                push!(equals, position)
                if verbose
                    println(STDOUT, "FOUND EQUAL FREQUENCIES: ", equals)
                end

            elseif found_new_center && d2[position][1] > new_center_freq
                new_center = position
                new_center_freq = d2[position][1]
                equals = []
                if verbose
                    println(STDOUT, "PREVIOUS POSITION: $center_pos with FREQUENCY: $center_freq")
                    println(STDOUT, "NEW POSITION: $position with FREQUENCY: ", d2[position][1])
                end

            elseif found_new_center && d2[position][1] == new_center_freq
                push!(equals, position)
            end
        end
    end

    # find witch one is nearest to median if found > 1 most frequent
    # and select it as center.
    if length(equals) > 0
        med = round(median(all_positions))

        if found_new_center
            diff = abs(new_center - med)
        else
            diff = abs(center_pos - med)
        end

        for i in equals
            if abs(i - med) < diff
                new_center = i
                found_new_center = true
                if verbose
                    println(STDOUT, "FOUND NEW CENTER MOST ADJ TO MEDIAN: ", i)
                end
            end
        end
    end

    if found_new_center
        new_cluster_id = chr*"::"*string(new_center)*"::"*str
    else
        # for checking adjecent clusters even if new center not found.
        new_cluster_id = near_center
        new_center = center_pos
    end

    # != arrises from recursive call when two adj clusters are merged.
    # in this case we want to cat all values of TS.
    if d1[near_center] != d2
        if verbose
            println(STDOUT, "CAT VALUES of $d2 to ", d1[near_center])
        end

        for position in keys(d2)
            if position in keys(d1[near_center])
                d1[near_center][position] = vcat(d2[position], d1[near_center][position][2:end])
                d1[near_center][position][1] += d2[position][1]
            else
                d1[near_center][position] = d2[position]
            end
        end

        if verbose
            println(STDOUT, "CATTED ARRAY: ", d1[near_center])
        end

        if found_new_center
            d1[new_cluster_id] = d1[near_center]
            delete!(d1, near_center)
            recenter!(d1, d1[new_cluster_id], new_cluster_id, k, chr, str; verbose=verbose)
            @goto stop_recenter!
        else
            recenter(d1, d1[near_center], near_center, k, chr, str; verbose=verbose)
            @goto stop_recenter!
        end
    end

    # check if values of a cluster is still holds +/- k
    # if not - reject those TS to recursion.
    # evaluate near stand alone TS and join if holds.
    rejected_TS1 = Dict{Int64,Array{Int64,1}}()
    if found_new_center
        for i in all_positions
            if i > new_center + k || i < new_center - k

                if verbose
                    println(STDOUT, "REJECTING POSITION FOR RECURSION: ", i)
                end

                rejected_TS1[i] = d2[i]
                if verbose
                    println(STDOUT, "!!! DELETING: $i --> ", d2[i])
                end

                delete!(d2, i)
            end
        end

        # correct main collection of cluster.
        if verbose
            println(STDOUT, "!!! DELETING FROM MAIN DATASET: $near_center --> ", d1[near_center])
            println(STDOUT, "RECREATE IN MAIN DATASET: $new_cluster_id --> ", d2)
        end

        delete!(d1, near_center)
        d1[new_cluster_id] = d2
    end

    search_adj_clust!(d1, new_center, new_cluster_id, chr, str, k; verbose=verbose)

    # recursively correct rejected
    if found_new_center
        for position in keys(rejected_TS1)

            if verbose
                println(STDOUT, "CLUSTER RECURSION OF REJECTED VALUES:")
                println(STDOUT, "NEW CLUSTER: $new_cluster_id, POSITION FOR RECURSION: $position")
            end

            cluster!(d1, chr, position, rejected_TS1[position], str, k; verbose=verbose)
            if verbose
                println(STDOUT, "RETURNED FROM CLUSTER RECURSION2.")
            end
        end
    end

    if verbose && !found_new_center
        println(STDOUT, "NEW CENTER WAS NOT FOUND.")
    end

    @label stop_recenter!

    return d1
end


function search_adj_clust!(d1, new_center::Int64, new_cluster_id::String,
                           chr::String, str::String, k::Int64; verbose=false)

        # search for adjacent clusters
        for adj_center in new_center-k*2:new_center+k*2
            adj_cluster_id = chr*"::"*string(adj_center)*"::"str

            # new cluster id might change from recursion.
            if !(new_cluster_id in keys(d1))
                break
            end

            if adj_cluster_id != new_cluster_id && adj_cluster_id in keys(d1)
                if verbose
                    println(STDOUT, "FOUND ADJ CLUSTER IN +/- 2*$k RANGE: $adj_cluster_id near $new_cluster_id")
                    println(STDOUT, "$adj_cluster_id:")
                end

                ck = collect(keys(d1[adj_cluster_id]))
                ck2 = collect(keys(d1[new_cluster_id]))
                min = minimum(ck)
                max = maximum(ck)

                if verbose
                    println(STDOUT, d1[adj_cluster_id])
                    println(STDOUT, ck)
                end

                # if found cluster falls in new cluster - add and check center
                # with recursion.
                if min in new_center-k:new_center+k && max in new_center-k:new_center+k
                    if verbose
                        println(STDOUT, "FOUND CLUSTER FALLS IN: $adj_cluster_id --> $new_cluster_id")
                        println(STDOUT, "RECENTERING WITH NEW CLUSTER ID: $new_cluster_id")
                    end

                    recenter!(d1, d1[adj_cluster_id], new_cluster_id, k, chr, str; verbose=verbose)

                    if verbose
                        println(STDOUT, "RETURNED FROM RECENTER RECURSION: $new_cluster_id")
                    end

                # if whole cluster does not fall, still some
                # values might be more adjacent to either cluster center
                else
                    if verbose
                        println(STDOUT, "FOUND CLUSTER DID NOT FALL IN: $adj_cluster_id --> $new_cluster_id")
                        println(STDOUT, "CHECKING IF ANY VALUES ARE MORE ADJACENT.")
                    end

                    check_adj_values!(d1, ck, new_center, adj_center,
                                      adj_cluster_id, chr, str,
                                      k; verbose=verbose)

                    check_adj_values!(d1, ck2, adj_center, new_center,
                                      new_cluster_id, chr, str,
                                      k; verbose=verbose)

                    if verbose
                        println(STDOUT, "FINISHED CHECKING ADJ VALUES.")
                    end
                end
            end
        end

    return d1
end


function check_adj_values!(d1::Dict, ck::Array{Int64,1}, new_center::Int64,
                           adj_center::Int64, adj_cluster_id::String,
                           chr::String, str::String, k::Int64; verbose=false)

    for position in ck
        if position != new_center && abs(adj_center - position) > abs(new_center - position)
            if verbose
                println(STDOUT, "FOUND VALUE TO BE MORE ADJ TO $new_center RATHER THEN $adj_cluster_id: $position")
                println(STDOUT, ck)
                println(STDOUT, d1[adj_cluster_id])
            end

            # sometimes keys are deleted in recursion
            if position in keys(d1[adj_cluster_id])
                val = values(d1[adj_cluster_id][position])
                if verbose
                    println(STDOUT, "!!! DELETE $position in $adj_cluster_id")
                    println(STDOUT, "AND PUSH IT FOR RECURSION.")
                end

                delete!(d1[adj_cluster_id], position)
                cluster!(d1, chr, position, val, str, k; verbose=verbose)

                if verbose
                    println(STDOUT, "RETURNED FROM CLUSTER RECURSION1")
                end

            elseif verbose
                println(STDOUT, "POSITION $position WAS ALREADY REMOVED.")
            end
        end
    end

    return d1
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


function stats_clusters(d::Dict)::DataFrame
    dframe = DataFrame(Chrmosome=String[], Start=Int64[], End=Int64[],
                       Strand=String[], ClusterCenter=Int64[], ClusterSize=Int32[],
                       ClusterMedian=Float64[], ClusterMean=Float64[],
                       ClusterMin=Int64[], ClusterMax=Int64[],
                       Cluster1stQuartile=Float64[], Cluster3rdQuartile=Float64[],
                       TSMedian=Float64[], TSMean=Float64[],
                       TSMin=Int64[], TSMax=Int64[],
                       TS1stQuartile=Float64[], TS3rdQuartile=Float64[])

    for k in keys(d)
        id = split(k,"::")
        chr = id[1]
        center = parse(Int64, id[2])
        str = id[3]

        all_positions = Int64[]
        frequencies = Int64[]
        ts_lengths = Int64[]
        for position in keys(d[k])
            push!(all_positions, position)
            push!(frequencies, d[k][position][1])
            ts_lengths = vcat(ts_lengths, d[k][position][2:end])
         end

        start = minimum(all_positions)
        end_  = maximum(all_positions)
        md = median(frequencies)
        avg = mean(frequencies)
        mn = minimum(frequencies)
        mx = maximum(frequencies)
        fq = nquantile(frequencies,4)[2]
        tq = nquantile(frequencies,4)[4]
        md2 = median(ts_lengths)
        avg2 = mean(ts_lengths)
        mn2 = minimum(ts_lengths)
        mx2 = maximum(ts_lengths)
        fq2 = nquantile(ts_lengths,4)[2]
        tq2 = nquantile(ts_lengths,4)[4]
        push!(dframe, [chr start end_ str center sum(frequencies) md avg mn mx fq tq md2 avg2 mn2 mx2 fq2 tq2])
    end

    return dframe
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
    catch md = string(md,"ID=;")
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
    catch md = string(md,"Parent=;")
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
    catch md = string(md,"gene_id=;")
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
    catch md = string(md,"transcript_id=;")
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
    catch md = string(md,"gene_type=;")
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
    catch md = string(md,"gene_name=;")
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
function annotate_polya_sites(a::IntervalCollection,
                              b::Dict{String, IntervalCollection{String}};
                              verbose::Bool=false)::DataFrame

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

            if verbose
                println(STDERR, "KeyError: $split_k")
            end
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

    # remove duplicated gene features
    dfp = rmdups(dfp)
    dfn = rmdups(dfn)

    enumeratenames!(dfn)
    enumeratenames!(dfp)
    enumeratenames!(dfna)

    df = vcat(dfp, dfn, dfna)

    return df
end


function annotate_clusters(a::DataFrame, b::Dict{String, IntervalCollection{String}};
                           verbose::Bool=false)::DataFrame

    # Two dataframes for '+' and '-' strands.
    dfp = DataFrame(Chr=String[], Start=Int64[], End=Int64[], GeneName=String[],
                    ClusterSize=Int32[], Strand=String[], Feature=String[],
                    ClusterCenter=Int64[], Biotype=String[],
                    ClusterMedian=Float64[], ClusterMean=Float64[],
                    ClusterMin=Int64[], ClusterMax=Int64[],
                    Cluster1stQuartile=Float64[], Cluster3rdQuartile=Float64[],
                    TSMedian=Float64[], TSMean=Float64[],
                    TSMin=Int64[], TSMax=Int64[],
                    TS1stQuartile=Float64[], TS3rdQuartile=Float64[])

    dfn = DataFrame(Chr=String[], Start=Int64[], End=Int64[], GeneName=String[],
                    ClusterSize=Int32[], Strand=String[], Feature=String[],
                    ClusterCenter=Int64[], Biotype=String[],
                    ClusterMedian=Float64[], ClusterMean=Float64[],
                    ClusterMin=Int64[], ClusterMax=Int64[],
                    Cluster1stQuartile=Float64[], Cluster3rdQuartile=Float64[],
                    TSMedian=Float64[], TSMean=Float64[],
                    TSMin=Int64[], TSMax=Int64[],
                    TS1stQuartile=Float64[], TS3rdQuartile=Float64[])

    # Dataframe for not annotated sites.
    dfna = DataFrame(Chr=String[], Start=Int64[], End=Int64[], GeneName=String[],
                    ClusterSize=Int32[], Strand=String[], Feature=String[],
                    ClusterCenter=Int64[], Biotype=String[],
                    ClusterMedian=Float64[], ClusterMean=Float64[],
                    ClusterMin=Int64[], ClusterMax=Int64[],
                    Cluster1stQuartile=Float64[], Cluster3rdQuartile=Float64[],
                    TSMedian=Float64[], TSMean=Float64[],
                    TSMin=Int64[], TSMax=Int64[],
                    TS1stQuartile=Float64[], TS3rdQuartile=Float64[])
    multiple_features = Int64(0)

    for i1 in eachrow(a)
        ct = Int64(0)
        # Parse metadata from polyA site.
        chr = i1[1]
        it1 = i1[2]
        it2 = i1[3]
        str = i1[4]
        c = i1[6]
        center = i1[5]
        cluster_multiple_features = Int64(0)

        split_keys = get_split_key(chr, it1, it2)

        for split_k in split_keys
            try
                for i2 in b[split_k]
                    if chechoverlap(Interval(chr, it1, it2, str[1]), i2)
                        cluster_multiple_features += 1
                        # parse metadata from GFF3 anotation
                        splt::Array{SubString{String},1} = split(metadata(i2),";")
                        ft = split(splt[7],"=")[2]
                        gn = split(splt[6],"=")[2]
                        bt = split(splt[5],"=")[2]
                        ct += 1

                        if str == "-"
                            push!(dfn, hcat([chr it1 it2 gn c str ft center bt],
                                            [i1[7] i1[8] i1[9] i1[10] i1[11]],
                                            [i1[12] i1[13] i1[14] i1[15] i1[16]],
                                            [i1[17] i1[18]]))

                        else
                            push!(dfp, hcat([chr it1 it2 gn c str ft center bt],
                                            [i1[7] i1[8] i1[9] i1[10] i1[11]],
                                            [i1[12] i1[13] i1[14] i1[15] i1[16]],
                                            [i1[17] i1[18]]))
                        end
                    end
                end

            catch
                KeyError(split_k)
                if verbose
                    println(STDERR, "KeyError: $split_k")
                end
            end
        end

        if cluster_multiple_features > 1
            multiple_features += 1
        end

        # if no overlaping of polyA coordinates with GFF3 annotated coordinates - NA is added.
        if ct == 0
            ft = "NA"
            gn = "NA"
            bt = "NA"

            push!(dfna, hcat([chr it1 it2 gn c str ft center bt],
                            [i1[7] i1[8] i1[9] i1[10] i1[11]],
                            [i1[12] i1[13] i1[14] i1[15] i1[16]],
                            [i1[17] i1[18]]))
        end
    end

    sort!(dfn, [:Chr, :Start, :End, :GeneName], rev=true)
    sort!(dfp, [:Chr, :Start, :End, :GeneName])
    sort!(dfna, [:Chr, :Start, :End, :GeneName])

    # remove duplicated gene features
    dfp = rmdups(dfp; cluster=true, verbose=verbose)
    dfn = rmdups(dfn; cluster=true, verbose=verbose)

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

    # Arguments
    - `dframe::DataFrame`: a sorted dataframe.
"""
function rmdups(dframe::DataFrame; cluster::Bool=false, verbose::Bool=false)::DataFrame

    if cluster
        df = DataFrame(Chr=String[], Start=Int64[], End=Int64[], GeneName=String[],
                        ClusterSize=Int32[], Strand=String[], Feature=String[],
                        ClusterCenter=Int64[], Biotype=String[],
                        ClusterMedian=Float64[], ClusterMean=Float64[],
                        ClusterMin=Int64[], ClusterMax=Int64[],
                        Cluster1stQuartile=Float64[], Cluster3rdQuartile=Float64[],
                        TSMedian=Float64[], TSMean=Float64[],
                        TSMin=Int64[], TSMax=Int64[],
                        TS1stQuartile=Float64[], TS3rdQuartile=Float64[])

    else
        df = DataFrame(Chr=String[], Start=Int64[], End=Int64[], Name=String[],
                       Counts=Int32[], Strand=String[], Feature=String[],
                       Median=Float32[], Min=Int16[], Max=Int16[], Biotype=String[])
    end

    d = Dict("gene"=>13, "transcript"=>12, "mRNA"=>12, "exon"=>11, "intron"=>10, "CDS"=>9,
             "UTR"=>8, "three_prime_UTR"=>7, "five_prime_UTR"=>6,
             "polyA_sequence"=>5, "polyA_site"=>4, "start_codon"=>3,
             "stop_codon"=>2, "Selenocysteine"=>1)

    ct = Int64(0)
    ct2 = Int64(0)

    for y in eachrow(dframe[:,1:6])
        ct2 += 1

        if ct == 0
            if cluster
                push!(df, hcat([dframe[1][ct2] dframe[2][ct2] dframe[3][ct2]],
                               [dframe[4][ct2] dframe[5][ct2] dframe[6][ct2]],
                               [dframe[7][ct2] dframe[8][ct2] dframe[9][ct2]],
                               [dframe[10][ct2] dframe[11][ct2]],
                               [dframe[12][ct2] dframe[13][ct2]],
                               [dframe[14][ct2] dframe[15][ct2]],
                               [dframe[16][ct2] dframe[17][ct2]],
                               [dframe[18][ct2] dframe[19][ct2]],
                               [dframe[20][ct2] dframe[21][ct2]]))

            else
                push!(df, hcat([dframe[1][ct2] dframe[2][ct2] dframe[3][ct2]],
                               [dframe[4][ct2] dframe[5][ct2] dframe[6][ct2]],
                               [dframe[7][ct2] dframe[8][ct2] dframe[9][ct2]],
                               [dframe[10][ct2] dframe[11][ct2]]))
            end
            ct += 1
        end

        if ct > 0
            if (df[end,:][1] == [y[1]] && df[end,:][2] == [y[2]]
                && df[end,:][3] == [y[3]] && df[end,:][4] == [y[4]]
                && df[end,:][5] == [y[5]] && df[end,:][6] == [y[6]])

                try
                    if d[dframe[7][ct2]] < d[df[end,:][7][1]]
                        df[7][ct] = dframe[7][ct2]
                    end

                catch
                    KeyError(dframe[7][ct2])
                    if verbose
                        println(STDERR, "KeyError: $ct2 ")
                    end
                end

            else
                ct += 1

                if cluster
                    push!(df, hcat([dframe[1][ct2] dframe[2][ct2] dframe[3][ct2]],
                                   [dframe[4][ct2] dframe[5][ct2] dframe[6][ct2]],
                                   [dframe[7][ct2] dframe[8][ct2] dframe[9][ct2]],
                                   [dframe[10][ct2] dframe[11][ct2]],
                                   [dframe[12][ct2] dframe[13][ct2]],
                                   [dframe[14][ct2] dframe[15][ct2]],
                                   [dframe[16][ct2] dframe[17][ct2]],
                                   [dframe[18][ct2] dframe[19][ct2]],
                                   [dframe[20][ct2] dframe[21][ct2]]))

                else
                    push!(df, hcat([dframe[1][ct2] dframe[2][ct2] dframe[3][ct2]],
                                   [dframe[4][ct2] dframe[5][ct2] dframe[6][ct2]],
                                   [dframe[7][ct2] dframe[8][ct2] dframe[9][ct2]],
                                   [dframe[10][ct2] dframe[11][ct2]]))

                end

            end
        end
    end

    return df
end
