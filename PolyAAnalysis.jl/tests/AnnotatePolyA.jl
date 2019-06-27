
# function clust!(d::Dict, p::Any, l::Any)::Dict

dc = Dict("chr1::150::+" => [Int16(100),Int16(200)], "chr2::160::-" => [Int16(300),Int16(400)])
k1 = "chr1::150::+"
k2 = "chr2::160::-"
k3 = "chr3::170::."
vl1 = Int16(300)
vl2 = Int16(400)
vl3 = Int16(600)

@test Dict("chr1::150::+" => [100,200,300], "chr2::160::-" => [300,400]) == clust!(dc, k1, vl1)
@test Dict("chr1::150::+" => [100,200,300], "chr2::160::-" => [300,400,400]) == clust!(dc, k2, vl2)
@test Dict("chr1::150::+" => [100,200,300], "chr2::160::-" => [300,400,400], "chr3::170::." => [600] ) == clust!(dc, k3, vl3)


# function cluster!(d::Dict, chr::String, pos::Int, l::Array{Int16,1}, str::String, k::Int)::Dict

dict = Dict("chr1::10::+" => Dict(10 => [2,400,120]))
chr1 = "chr1"
pos1 = 8
l1 = [Int64(80)]
str1 = "+"
k1 = 5

@test Dict("chr1::10::+" => Dict(10=>[2,400,120], 8=>[1,80])) == cluster!(dict, chr1, pos1, l1, str1, k1)

cluster!(dict, chr1, 20, l1, str1, k1)
test_dict1 = Dict("chr1::10::+" => Dict(10=>[2,400,120], 8=>[1,80]), "chr1::20::+" => Dict(20=>[1,80]))

for key in sort(collect(keys(dict)))
    @test dict[key] == test_dict1[key]
end

cluster!(dict, chr1, 20, l1, str1, k1)
test_dict2 = Dict("chr1::10::+" => Dict(10=>[2,400,120], 8=>[1,80]), "chr1::20::+" => Dict(20=>[2,80,80]))

for key in sort(collect(keys(dict)))
    @test dict[key] == test_dict2[key]
end

# when center changes with a median.
for i in 1:10
    cluster!(dict, chr1, 20+i, l1, str1, k1)
    cluster!(dict, chr1, 20-i, l1, str1, k1)
end

test_dict3 = Dict(14=>[1,80], 10=>[3,400,120,80], 13=>[1,80], 11=>[1,80],
                  8=>[1,80], 12=>[1,80])
for i in keys(dict["chr1::10::+"])
    @test dict["chr1::10::+"][i] == test_dict3[i]
end

test_dict4 = Dict(16=>[1,80], 21=>[1,80], 25=>[1,80], 19=>[1,80], 17=>[1,80],
                  22=>[1,80], 20=>[2,80,80], 15=>[1,80], 18=>[1,80])
for i in keys(dict["chr1::20::+"])
    @test dict["chr1::20::+"][i] == test_dict4[i]
end

test_dict5 = Dict(23=>[1,80], 24=>[1,80], 25=>[1,80], 26=>[1,80], 29=>[1,80], 27=>[1,80],
                  28=>[1,80], 30=>[1,80])
for i in keys(dict["chr1::26::+"])
    @test dict["chr1::26::+"][i] == test_dict5[i]
end

# when center changes with frequency.
for i in 1:10
    cluster!(dict, chr1, 27, l1, str1, k1)
end

test_dict6 = Dict(14=>[1,80], 10=>[3, 400, 120,80], 13=>[1,80], 11=>[1,80],
                  8=>[1,80], 12=>[1,80])
for i in keys(dict["chr1::10::+"])
    @test dict["chr1::10::+"][i] == test_dict6[i]
end

test_dict7 = Dict(16=>[1,80], 21=>[1,80], 25=>[1,80], 19=>[1,80], 17=>[1,80],
                  22=>[1,80], 24=>[1,80], 20=>[2,80,80], 23=>[1,80], 15=>[1,80],
                  18=>[1,80])
for i in keys(dict["chr1::20::+"])
    @test dict["chr1::20::+"][i] == test_dict7[i]
end

test_dict8 = Dict(24=>[1,80], 25=>[1,80], 27=>[11,80,80,80,80,80,80,80,80,80,80,80],
                  29=>[1,80], 26=>[1,80], 28=>[1,80], 30=>[1,80])
for i in keys(dict["chr1::27::+"])
    @test dict["chr1::27::+"][i] == test_dict8[i]
end

# when center changes and recursion follows for discarded TS on the left.
for i in 1:5
    cluster!(dict, chr1, 40+i, l1, str1, k1)
end

cluster!(dict, chr1, 46, l1, str1, k1)

for i in 1:5
    cluster!(dict, chr1, 47, l1, str1, k1)
end

for i in keys(dict["chr1::10::+"])
    @test dict["chr1::10::+"][i] == test_dict6[i]
end

for i in keys(dict["chr1::20::+"])
    @test dict["chr1::20::+"][i] == test_dict7[i]
end

for i in keys(dict["chr1::27::+"])
    @test dict["chr1::27::+"][i] == test_dict8[i]
end

test_dict9 = Dict(43=>[1,80], 41=>[1,80], 42=>[1,80], 44=>[1,80])
for i in keys(dict["chr1::42::+"])
    @test dict["chr1::42::+"][i] == test_dict9[i]
end

test_dict10 = Dict(46=>[1,80], 47=>[5,80,80,80,80,80], 45=>[1,80])
for i in keys(dict["chr1::47::+"])
    @test dict["chr1::47::+"][i] == test_dict10[i]
end

# when center changes and recursion follows for discarded TS on the rigth.

for i in 1:5
    cluster!(dict, chr1, 71-i, l1, str1, k1)
end

for i in 1:2
    cluster!(dict, chr1, 65, l1, str1, k1)
end

cluster!(dict, chr1, 72, l1, str1, k1)

for i in 1:5
    cluster!(dict, chr1, 61, l1, str1, k1)
end

for i in keys(dict["chr1::10::+"])
    @test dict["chr1::10::+"][i] == test_dict6[i]
end

for i in keys(dict["chr1::20::+"])
    @test dict["chr1::20::+"][i] == test_dict7[i]
end

for i in keys(dict["chr1::27::+"])
    @test dict["chr1::27::+"][i] == test_dict8[i]
end

for i in keys(dict["chr1::42::+"])
    @test dict["chr1::42::+"][i] == test_dict9[i]
end

for i in keys(dict["chr1::47::+"])
    @test dict["chr1::47::+"][i] == test_dict10[i]
end

test_dict11 = Dict(61=>[5,80,80,80,80,80], 65=>[2,80,80], 66=>[1,80], 67=>[1,80])

for i in keys(dict["chr1::61::+"])
    @test dict["chr1::61::+"][i] == test_dict11[i]
end

test_dict12 = Dict(68=>[1,80], 69=>[1,80], 70=>[1,80], 72=>[1,80])
for i in keys(dict["chr1::70::+"])
    @test dict["chr1::70::+"][i] == test_dict12[i]
end

# when center changes and adjacent left cluster falls in.
for i in 1:5
    cluster!(dict, chr1, 99+i, l1, str1, k1)
end

for i in 1:2
    cluster!(dict, chr1, 105, l1, str1, k1)
end

cluster!(dict, chr1, 98, l1, str1, k1)

for i in 1:5
    cluster!(dict, chr1, 101, l1, str1, k1)
end

for i in keys(dict["chr1::10::+"])
    @test dict["chr1::10::+"][i] == test_dict6[i]
end

for i in keys(dict["chr1::20::+"])
    @test dict["chr1::20::+"][i] == test_dict7[i]
end

for i in keys(dict["chr1::27::+"])
    @test dict["chr1::27::+"][i] == test_dict8[i]
end

for i in keys(dict["chr1::42::+"])
    @test dict["chr1::42::+"][i] == test_dict9[i]
end

for i in keys(dict["chr1::47::+"])
    @test dict["chr1::47::+"][i] == test_dict10[i]
end

for i in keys(dict["chr1::61::+"])
    @test dict["chr1::61::+"][i] == test_dict11[i]
end

for i in keys(dict["chr1::70::+"])
    @test dict["chr1::70::+"][i] == test_dict12[i]
end

test_dict13 = Dict(100=>[1,80],102=>[1,80],98=>[1,80],101=>[6,80,80,80,80,80,80],
                   103=>[1,80],104=>[1,80],105=>[2,80,80])
for i in keys(dict["chr1::101::+"])
    @test dict["chr1::101::+"][i] == test_dict13[i]
end

# when center changes and adjacent right cluster falls in.
for i in 1:5
    cluster!(dict, chr1, 199+i, l1, str1, k1)
end

for i in 1:2
    cluster!(dict, chr1, 199, l1, str1, k1)
end

cluster!(dict, chr1, 205, l1, str1, k1)

for i in 1:5
    cluster!(dict, chr1, 205, l1, str1, k1)
end

for i in keys(dict["chr1::10::+"])
    @test dict["chr1::10::+"][i] == test_dict6[i]
end

for i in keys(dict["chr1::20::+"])
    @test dict["chr1::20::+"][i] == test_dict7[i]
end

for i in keys(dict["chr1::27::+"])
    @test dict["chr1::27::+"][i] == test_dict8[i]
end

for i in keys(dict["chr1::42::+"])
    @test dict["chr1::42::+"][i] == test_dict9[i]
end

for i in keys(dict["chr1::47::+"])
    @test dict["chr1::47::+"][i] == test_dict10[i]
end

for i in keys(dict["chr1::61::+"])
    @test dict["chr1::61::+"][i] == test_dict11[i]
end

for i in keys(dict["chr1::70::+"])
    @test dict["chr1::70::+"][i] == test_dict12[i]
end

for i in keys(dict["chr1::101::+"])
    @test dict["chr1::101::+"][i] == test_dict13[i]
end

test_dict14 = Dict(102=>[1,80], 103=>[1,80], 104=>[1,80], 105=>[2,80,80])
for i in keys(dict["chr1::105::+"])
    @test dict["chr1::105::+"][i] == test_dict14[i]
end

test_dict15 = Dict(199=>[2,80,80], 200=>[1,80], 201=>[1,80], 202=>[1,80])
for i in keys(dict["chr1::199::+"])
    @test dict["chr1::199::+"][i] == test_dict15[i]
end

test_dict16 = Dict(203=>[1,80], 204=>[1,80], 205=>[6,80,80,80,80,80,80])
for i in keys(dict["chr1::205::+"])
    @test dict["chr1::205::+"][i] == test_dict16[i]
end

# when center changes and adjecent cluster falls in with greater frequency of the center.
for i in 1:5
    cluster!(dict, chr1, 300+i, l1, str1, k1)
end

for i in 1:3
    cluster!(dict, chr1, 302, l1, str1, k1)
end

cluster!(dict, chr1, 307, l1, str1, k1)

for i in 1:5
    cluster!(dict, chr1, 306, l1, str1, k1)
end

test_dict17 = Dict(301=>[1,80], 302=>[4,80,80,80,80], 303=>[1,80], 304=>[1,80],
                   305=>[1,80], 306=>[5,80,80,80,80,80], 307=>[1,80])
for i in keys(dict["chr1::306::+"])
    @test dict["chr1::306::+"][i] == test_dict17[i]
end

# test all positions.

collection_of_positions = Int64[]
for i in keys(dict)
    for position in keys(dict[i])
        push!(collection_of_positions, position)
    end
end

test_collection_of_positions = [10, 8, 20, 21, 19, 22, 18, 23, 17, 24, 16, 25,
                                15, 26, 14, 27, 13, 28, 12, 29, 11, 30, 41, 42,
                                43, 44, 45, 46, 47, 70, 69, 68, 67, 66, 65, 61,
                                100, 101, 102, 103, 104, 105, 98, 199, 200, 201,
                                202, 203, 204, 205, 72, 301, 302, 303, 304, 305,
                                306, 307]

@test sort(test_collection_of_positions) == sort(collection_of_positions)

# function merge_adj_clusters!(d::Dict, m=Int64; verbose::Bool=false)::Dict

merge_adj_clusters!(dict, 1)

test_dict18 = Dict("chr1::205::+"=>Dict(200=>[1,80], 203=>[1,80], 199=>[2,80,80],
                                        204=>[1,80], 201=>[1,80], 205=>[6,80,80,80,80,80,80],
                                        202=>[1,80]),
                   "chr1::27::+"=>Dict(30=>[1,80], 16=>[1,80], 11=>[1,80], 21=>[1,80],
                                       26=>[1,80], 25=>[1,80], 10=>[3,400,120,80],
                                       29=>[1,80], 19=>[1,80], 17=>[1,80], 8=>[1,80],
                                       22=>[1,80], 24=>[1,80], 28=>[1,80], 20=>[2,80,80],
                                       23=>[1,80], 14=>[1,80], 13=>[1,80],
                                       27=>[11,80,80,80,80,80,80,80,80,80,80,80],
                                       15=>[1,80], 12=>[1,80], 18=>[1,80]),
                   "chr1::65::+"=>Dict(68=>[1,80], 66=>[1,80], 69=>[1,80], 67=>[1,80],
                                       72=>[1,80], 65=>[2,80,80], 70=>[1,80]),
                   "chr1::61::+"=>Dict(61=>[5,80,80,80,80,80]),
                   "chr1::101::+"=>Dict(100=>[1,80], 102=>[1,80], 98=>[1,80],
                                        101=>[6,80,80,80,80,80,80], 103=>[1,80],
                                        104=>[1,80], 105=>[2,80,80]),
                   "chr1::47::+"=>Dict(43=>[1,80], 47=>[5,80,80,80,80,80], 45=>[1,80],
                                       42=>[1,80], 41=>[1,80], 44=>[1,80], 46=>[1,80]),
                   "chr1::306::+"=>Dict(305=>[1,80], 306=>[5,80,80,80,80,80],
                                        304=>[1,80], 307=>[1,80], 301=>[1,80],
                                        303=>[1,80], 302=>[4,80,80,80,80]))

for k in keys(dict)
    @test test_dict18[k] == sort(dict[k])
end


# function stats_poly_a(d::Dict{String,Array{Int16,1}})::DataFrame

@test DataFrame(Chrmosome=String["chr1","chr2","chr3"],
             Position=Int32[150,160,170],
             Strand=String["+","-","."],
             Median=Float16[200,400,600],
             Minimum=Int16[100,300,600],
             Maximum=Int16[300,400,600],
             Counts=Int32[3,3,1]) == sort!(stats_poly_a(dc))


# function parserecord(r::GFF3.Record)::String

record = GFF3.Record(b"chr1\tHAVANA\texon\t43930132\t43931077\t.\t+\t.\tID=exon:ENST00000643283.1:12;Parent=ENST00000643283.1;gene_id=ENSG00000126091.20;transcript_id=ENST00000643283.1;gene_type=protein_coding;gene_name=ST3GAL3;transcript_type=nonsense_mediated_decay;transcript_name=RP11-7O11.1-079;exon_number=12;exon_id=ENSE00003818117.1;level=2;protein_id=ENSP00000494746.1;tag=RNA_Seq_supported_only;havana_gene=OTTHUMG00000007561.18;havana_transcript=OTTHUMT00000494570.1")
record2 = GFF3.Record(b"chr1\tHAVANA\texon\t43930132\t43931077\t.\t+\t.\tParent=ENST00000643283.1;gene_id=ENSG00000126091.20;transcript_id=ENST00000643283.1;gene_type=protein_coding;gene_name=ST3GAL3;transcript_type=nonsense_mediated_decay;transcript_name=RP11-7O11.1-079;exon_number=12;exon_id=ENSE00003818117.1;level=2;protein_id=ENSP00000494746.1;tag=RNA_Seq_supported_only;havana_gene=OTTHUMG00000007561.18;havana_transcript=OTTHUMT00000494570.1")
record3 = GFF3.Record(b"chr1\tHAVANA\texon\t43930132\t43931077\t.\t+\t.\tgene_id=ENSG00000126091.20;transcript_id=ENST00000643283.1;gene_type=protein_coding;gene_name=ST3GAL3;transcript_type=nonsense_mediated_decay;transcript_name=RP11-7O11.1-079;exon_number=12;exon_id=ENSE00003818117.1;level=2;protein_id=ENSP00000494746.1;tag=RNA_Seq_supported_only;havana_gene=OTTHUMG00000007561.18;havana_transcript=OTTHUMT00000494570.1")
record4 = GFF3.Record(b"chr1\tHAVANA\texon\t43930132\t43931077\t.\t+\t.\t")
record5 = GFF3.Record(b"chr1\tHAVANA\texon\t43930132\t43931077\t.\t+\t.\tgene_id=ENSG00000126091.20")

@test String("ID=exon:ENST00000643283.1:12;"*
             "Parent=ENST00000643283.1;"*
             "gene_id=ENSG00000126091.20;"*
             "transcript_id=ENST00000643283.1;"*
             "gene_type=protein_coding;"*
             "gene_name=ST3GAL3;"*
             "ftype=exon") == parserecord(record)
@test String("ID=;"*
             "Parent=ENST00000643283.1;"*
             "gene_id=ENSG00000126091.20;"*
             "transcript_id=ENST00000643283.1;"*
             "gene_type=protein_coding;"*
             "gene_name=ST3GAL3;"*
             "ftype=exon") == parserecord(record2)
@test String("ID=;"*
             "Parent=;"*
             "gene_id=ENSG00000126091.20;"*
             "transcript_id=ENST00000643283.1;"*
             "gene_type=protein_coding;"*
             "gene_name=ST3GAL3;"*
             "ftype=exon") == parserecord(record3)
@test String("ID=;"*
             "Parent=;"*
             "gene_id=;"*
             "transcript_id=;"*
             "gene_type=;"*
             "gene_name=;"*
             "ftype=exon") == parserecord(record4)
@test String("ID=;"*
             "Parent=;"*
             "gene_id=ENSG00000126091.20;"*
             "transcript_id=;"*
             "gene_type=;"*
             "gene_name=;"*
             "ftype=exon") == parserecord(record5)

# function getintervals(dframe::DataFrame)::IntervalCollection{String}

tdf = DataFrame(Chrmosome=String["chr1","chr2","chr3"],
                Position=Int32[150,160,170],
                Strand=String["+","-","."],
                Median=Float16[200,400,600],
                Minimum=Int16[100,300,600],
                Maximum=Int16[300,400,600],
                Counts=Int32[3,3,1],
                Smaple=String["A","A","A"],
                TotalReads=Int64[10000,10000,10000],
                PassedReads=Int64[1000,1000,1000],
                PolyAReads=Int64[500,500,500]
                )

@test IntervalCollection([Interval("chr1",150,150,Strand('+'),"Median=200.0;Min=100;Max=300;Counts=3")
                          Interval("chr2",160,160,Strand('-'),"Median=400.0;Min=300;Max=400;Counts=3")
                          Interval("chr3",170,170,Strand('.'),"Median=600.0;Min=600;Max=600;Counts=1")
                          ]) == getintervals(tdf)

# function annotate_polya_sites(a::IntervalCollection,  b::Dict{String, IntervalCollection{String}})::DataFrame

intcol1 = IntervalCollection([Interval("chr1",150,150,Strand('+'),"Median=200.0;Min=100;Max=300;Counts=3")
                          Interval("chr2",160,160,Strand('-'),"Median=400.0;Min=300;Max=400;Counts=3")
                          Interval("chr3",170,170,Strand('.'),"Median=600.0;Min=600;Max=600;Counts=1")
                          Interval("chr3",3000,3000,Strand('.'),"Median=600.0;Min=600;Max=600;Counts=1")
                          ])

intcol2 = Dict("chr1:1:10000" => IntervalCollection([Interval("chr1",150,150,Strand('+'),
                                                          String("ID=exon:ENST00000643283.1:12;"*
                                                                 "Parent=ENST00000643283.1;"*
                                                                 "gene_id=ENSG00000126091.20;"*
                                                                 "transcript_id=ENST00000643283.1;"*
                                                                 "gene_type=protein_coding;"*
                                                                 "gene_name=ST3GAL3;"*
                                                                 "ftype=exon"))
                                                         Interval("chr1",160,160,Strand('-'),
                                                          String("ID=exon:ENST00000643282.1:12;"*
                                                                 "Parent=ENST00000643282.1;"*
                                                                 "gene_id=ENSG00000126092.20;"*
                                                                 "transcript_id=ENST00000643282.1;"*
                                                                 "gene_type=protein_coding;"*
                                                                 "gene_name=ST3GAL32;"*
                                                                 "ftype=intron"))
                                                         Interval("chr1",170,170,Strand('.'),
                                                          String("ID=exon:ENST00000643283.4:12;"*
                                                                 "Parent=ENST00000643283.4;"*
                                                                 "gene_id=ENSG00000126094.20;"*
                                                                 "transcript_id=ENST00000643284.1;"*
                                                                 "gene_type=protein_coding;"*
                                                                 "gene_name=ST3GAL43;"*
                                                                 "ftype=exon"))
                                                         ]),
               "chr2:1:10000" => IntervalCollection([Interval("chr2",150,150,Strand('+'),
                                                          String("ID=exon:ENST00000643283.1:12;"*
                                                                 "Parent=ENST00000643283.1;"*
                                                                 "gene_id=ENSG00000126091.20;"*
                                                                 "transcript_id=ENST00000643283.1;"*
                                                                 "gene_type=protein_coding;"*
                                                                 "gene_name=ST3GAL3;"*
                                                                 "ftype=exon"))
                                                         Interval("chr2",160,160,Strand('-'),
                                                          String("ID=exon:ENST00000643282.1:12;"*
                                                                 "Parent=ENST00000643282.1;"*
                                                                 "gene_id=ENSG00000126092.20;"*
                                                                 "transcript_id=ENST00000643282.1;"*
                                                                 "gene_type=protein_coding;"*
                                                                 "gene_name=ST3GAL32;"*
                                                                 "ftype=intron"))
                                                         Interval("chr2",170,170,Strand('.'),
                                                          String("ID=exon:ENST00000643283.4:12;"*
                                                                 "Parent=ENST00000643283.4;"*
                                                                 "gene_id=ENSG00000126094.20;"*
                                                                 "transcript_id=ENST00000643284.1;"*
                                                                 "gene_type=protein_coding;"*
                                                                 "gene_name=ST3GAL43;"*
                                                                 "ftype=exon"))
                                                         ]),
               "chr3:1:10000" => IntervalCollection([Interval("chr3",150,150,Strand('+'),
                                                      String("ID=exon:ENST00000643283.1:12;"*
                                                             "Parent=ENST00000643283.1;"*
                                                             "gene_id=ENSG00000126091.20;"*
                                                             "transcript_id=ENST00000643283.1;"*
                                                             "gene_type=protein_coding;"*
                                                             "gene_name=ST3GAL3;"*
                                                             "ftype=exon"))
                                                     Interval("chr3",160,160,Strand('-'),
                                                      String("ID=exon:ENST00000643282.1:12;"*
                                                             "Parent=ENST00000643282.1;"*
                                                             "gene_id=ENSG00000126092.20;"*
                                                             "transcript_id=ENST00000643282.1;"*
                                                             "gene_type=protein_coding;"*
                                                             "gene_name=ST3GAL32;"*
                                                             "ftype=intron"))
                                                     Interval("chr3",170,170,Strand('.'),
                                                      String("ID=exon:ENST00000643283.4:12;"*
                                                             "Parent=ENST00000643283.4;"*
                                                             "gene_id=ENSG00000126094.20;"*
                                                             "transcript_id=ENST00000643284.1;"*
                                                             "gene_type=protein_coding;"*
                                                             "gene_name=ST3GAL43;"*
                                                             "ftype=exon"))
                                                     ])
               )

@test DataFrame(Chr=String["chr1","chr2","chr3","chr3"],
                Start=Int64[149,159,169,2999],
                End=Int64[150,160,170,3000],
                Name=String["ST3GAL3.1","ST3GAL32.1","ST3GAL43.1","NA.1"],
                Counts=Int32[3,3,1,1],
                Strand=String["+","-",".","."],
                Feature=String["exon","intron","exon","NA"],
                Median=Float32[200.0,400.0,600.0,600.0],
                Min=Int16[100,300,600,600],
                Max=Int16[300,400,600,600],
                Biotype=String["protein_coding","protein_coding","protein_coding","NA"]
                ) == sort!(annotate_polya_sites(intcol1, intcol2), [:Chr, :Start])

intcol1 = IntervalCollection([Interval("chr1",1500,1500,Strand('+'),"Median=200.0;Min=100;Max=300;Counts=3")
                              Interval("chr1",3000,3000,Strand('.'),"Median=600.0;Min=600;Max=600;Counts=1")
                              Interval("chr1",10000,10000,Strand('-'),"Median=400.0;Min=300;Max=400;Counts=3")
                              Interval("chr1",10001,10001,Strand('.'),"Median=600.0;Min=600;Max=600;Counts=1")
                              Interval("chr1",17000,17000,Strand('.'),"Median=600.0;Min=600;Max=600;Counts=1")
                              ])

intcol2 = Dict("chr1:1:10000" => IntervalCollection([Interval("chr1",1,1000,Strand('+'),
                                                              String("ID=exon:ENST00000643283.1:12;"*
                                                                     "Parent=ENST00000643283.1;"*
                                                                     "gene_id=ENSG00000126091.20;"*
                                                                     "transcript_id=ENST00000643283.1;"*
                                                                     "gene_type=protein_coding;"*
                                                                     "gene_name=1;"*
                                                                     "ftype=exon"))
                                                     Interval("chr1",1,20000,Strand('.'),
                                                              String("ID=exon:ENST00000643283.4:12;"*
                                                                     "Parent=ENST00000643283.4;"*
                                                                     "gene_id=ENSG00000126094.20;"*
                                                                     "transcript_id=ENST00000643284.1;"*
                                                                     "gene_type=protein_coding;"*
                                                                     "gene_name=2;"*
                                                                     "ftype=exon"))
                                                     Interval("chr1",1000,10000,Strand('-'),
                                                              String("ID=exon:ENST00000643282.1:12;"*
                                                                     "Parent=ENST00000643282.1;"*
                                                                     "gene_id=ENSG00000126092.20;"*
                                                                     "transcript_id=ENST00000643282.1;"*
                                                                     "gene_type=protein_coding;"*
                                                                     "gene_name=3;"*
                                                                     "ftype=intron"))
                                                         ]),
                "chr1:10001:20000" => IntervalCollection([Interval("chr1",1,20000,Strand('.'),
                                                                   String("ID=exon:ENST00000643283.4:12;"*
                                                                          "Parent=ENST00000643283.4;"*
                                                                          "gene_id=ENSG00000126094.20;"*
                                                                          "transcript_id=ENST00000643284.1;"*
                                                                          "gene_type=three_prime_UTR;"*
                                                                          "gene_name=4;"*
                                                                          "ftype=exon"))
                                                          Interval("chr1",150,150,Strand('+'),
                                                                   String("ID=exon:ENST00000643283.1:12;"*
                                                                          "Parent=ENST00000643283.1;"*
                                                                          "gene_id=ENSG00000126091.20;"*
                                                                          "transcript_id=ENST00000643283.1;"*
                                                                          "gene_type=protein_coding;"*
                                                                          "gene_name=5;"*
                                                                          "ftype=exon"))
                                                          Interval("chr1",160,160,Strand('-'),
                                                                   String("ID=exon:ENST00000643282.1:12;"*
                                                                          "Parent=ENST00000643282.1;"*
                                                                          "gene_id=ENSG00000126092.20;"*
                                                                          "transcript_id=ENST00000643282.1;"*
                                                                          "gene_type=protein_coding;"*
                                                                          "gene_name=6;"*
                                                                          "ftype=intron"))
                                                      ])
                  )

@test DataFrame(Chr=String["chr1","chr1","chr1","chr1","chr1"],
              Start=Int64[1499,2999,9999,10000,16999],
              End=Int64[1500,3000,10000,10001,17000],
              Name=String["NA.1","2.1","3.1","4.1","4.2"],
              Counts=Int32[3,1,3,1,1],
              Strand=String["+",".","-",".","."],
              Feature=String["NA","exon","intron","exon","exon"],
              Median=Float32[200.0,600.0,400.0,600.0,600.0],
              Min=Int16[100,600,300,600,600],
              Max=Int16[300,600,400,600,600],
              Biotype=String["NA","protein_coding","protein_coding","three_prime_UTR","three_prime_UTR"]
              ) == sort!(annotate_polya_sites(intcol1, intcol2), [:Chr, :Start])

# function rmdups(dframe::DataFrame)::DataFrame

tdf = DataFrame(Chr=String["chr1","chr1","chr2","chr2","chr2","chr2"],
                Start=Int64[150,160,120,120,120,180],
                End=Int64[150,160,120,120,120,180],
                Name=String["A.1","A.2","B.1","B.1","B.1","C.1"],
                Counts=Int32[1,2,4,4,4,5],
                Strand=String["+","-","+","+","+","+"],
                Feature=String["exon","transcript","gene","transcript","CDS","gene"],
                Median=Float32[200.0,400.0,600.0,600.0,600.0,100.0],
                Min=Int16[100,300,100,100,100,200],
                Max=Int16[300,400,800,800,800,400],
                Biotype=String["AC","AC","BC","BC","BC","GC"]
                )

tdf2 = DataFrame(Chr=String["chr1","chr1","chr2","chr2"],
                 Start=Int64[150,160,120,180],
                 End=Int64[150,160,120,180],
                 Name=String["A.1","A.2","B.1","C.1"],
                 Counts=Int32[1,2,4,5],
                 Strand=String["+","-","+","+"],
                 Feature=String["exon","transcript","CDS","gene"],
                 Median=Float32[200.0,400.0,600.0,100.0],
                 Min=Int16[100,300,100,200],
                 Max=Int16[300,400,800,400],
                 Biotype=String["AC","AC","BC","GC"]
                 )

tdf4 = DataFrame(Chr=String["chr1","chr1","chr2","chr2","chr2","chr2"],
                 Start=Int64[150,160,120,120,120,180],
                 End=Int64[150,160,120,120,120,180],
                 Name=String["A.1","A.2","B.1","B.2.1","B.3.1","C.1"],
                 Counts=Int32[1,2,4,4,4,5],
                 Strand=String["+","-","+","+","+","+"],
                 Feature=String["exon","transcript","gene","transcript","CDS","gene"],
                 Median=Float32[200.0,400.0,600.0,600.0,600.0,100.0],
                 Min=Int16[100,300,100,100,100,200],
                 Max=Int16[300,400,800,800,800,400],
                 Biotype=String["AC","AC","BC","BC","BC","GC"]
                 )

tdf5 = DataFrame(Chr=String["chr1","chr1","chr2","chr2","chr2","chr2"],
                 Start=Int64[150,160,120,120,120,180],
                 End=Int64[150,160,120,120,120,180],
                 Name=String["A.1","A.2","B.1","B.2.1","B.3.1","C.1"],
                 Counts=Int32[1,2,4,4,4,5],
                 Strand=String["+","-","+","+","+","+"],
                 Feature=String["exon","transcript","gene","transcript","CDS","gene"],
                 Median=Float32[200.0,400.0,600.0,600.0,600.0,100.0],
                 Min=Int16[100,300,100,100,100,200],
                 Max=Int16[300,400,800,800,800,400],
                 Biotype=String["AC","AC","BC","BC","BC","GC"]
                 )

@test tdf2 == rmdups(tdf)
@test tdf5 == rmdups(tdf4)

# function enumeratenames!(df::DataFrame)::DataFrame

tdf3 = DataFrame(Chr=String["chr1","chr1","chr2","chr2","chr2","chr2"],
                 Start=Int64[150,160,120,120,120,180],
                 End=Int64[150,160,120,120,120,180],
                 Name=String["A.1.1","A.2.1","B.1.1","B.1.2","B.1.3","C.1.1"],
                 Counts=Int32[1,2,4,4,4,5],
                 Strand=String["+","-","+","+","+","+"],
                 Feature=String["exon","transcript","gene","transcript","CDS","gene"],
                 Median=Float32[200.0,400.0,600.0,600.0,600.0,100.0],
                 Min=Int16[100,300,100,100,100,200],
                 Max=Int16[300,400,800,800,800,400],
                 Biotype=String["AC","AC","BC","BC","BC","GC"]
                 )

@test tdf3 == enumeratenames!(tdf)

# function get_split_key(chr::String, x::Int; step::Int64=10000)::String

@test ["chr1:1:10"] == get_split_key("chr1",2,4; step=10)
@test ["chr2:1291:1300",
       "chr2:1301:1310",
       "chr2:1311:1320"] == get_split_key("chr2",1300,1315; step=10)
@test ["chr2:1291:1300",
      "chr2:1301:1310",
      "chr2:1311:1320",
      "chr2:1321:1330"] == get_split_key("chr2",1300,1325; step=10)
@test ["chr4:1:10000"] == get_split_key("chr4",1300,1300; step=10000)
@test ["chr4:1:10000"] == get_split_key("chr4",10000,10000; step=10000)
@test ["chr4:10001:20000"] == get_split_key("chr4",10290,10290; step=10000)
@test ["chr4:10001:20000"] == get_split_key("chr4",10590,10590; step=10000)
@test ["chr4:500001:510000"] == get_split_key("chr4",500590,500590; step=10000)
@test ["chr21:45070001:45080000"] == get_split_key("chr21",45073853,45073853)
@test ["chr21:46680001:46690000"] == get_split_key("chr21", 46689753,46689753)
@test ["chr1:11:20",
       "chr1:21:30",
       "chr1:31:40"] == get_split_key("chr1",12,35; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",12,25; step=10)
@test ["chr1:11:20"] == get_split_key("chr1",11,20; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",11,21; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",20,30; step=10)
@test ["chr1:11:20",
       "chr1:21:30",
       "chr1:31:40"] == get_split_key("chr1",20,31; step=10)
@test ["chr1:11:20"] == get_split_key("chr1",12,20; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",12,21; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",11,22; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",11,30; step=10)
@test ["chr1:11:20",
       "chr1:21:30",
       "chr1:31:40"] == get_split_key("chr1",11,31; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",12,30; step=10)
@test ["chr1:11:20",
       "chr1:21:30",
       "chr1:31:40"] == get_split_key("chr1",11,32; step=10)
@test ["chr1:11:20",
       "chr1:21:30",
       "chr1:31:40"] == get_split_key("chr1",12,32; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",11,29; step=10)
@test ["chr1:11:20",
       "chr1:21:30"] == get_split_key("chr1",12,29; step=10)
@test ["chr1:11:20"] == get_split_key("chr1",11,19; step=10)
@test ["chr1:11:20"] == get_split_key("chr1",12,19; step=10)
