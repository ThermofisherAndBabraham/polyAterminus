
# function Clust!(d::Dict, p::Any, l::Any)::Dict

dc = Dict("chr1::150::+" => [Int16(100),Int16(200)], "chr2::160::-" => [Int16(300),Int16(400)])
k1 = "chr1::150::+"
k2 = "chr2::160::-"
k3 = "chr3::170::."
vl1 = Int16(300)
vl2 = Int16(400)
vl3 = Int16(600)

@test Dict("chr1::150::+" => [100,200,300], "chr2::160::-" => [300,400]) == Clust!(dc, k1, vl1)
@test Dict("chr1::150::+" => [100,200,300], "chr2::160::-" => [300,400,400]) == Clust!(dc, k2, vl2)
@test Dict("chr1::150::+" => [100,200,300], "chr2::160::-" => [300,400,400], "chr3::170::." => [600] ) == Clust!(dc, k3, vl3)


# function PolACalculus(d::Dict{String,Array{Int16,1}})::DataFrame

@test DataFrame(Chrmosome=String["chr1","chr2","chr3"],
             Position=Int32[150,160,170],
             Strand=String["+","-","."],
             Median=Float16[200,400,600],
             Minimum=Int16[100,300,600],
             Maximum=Int16[300,400,600],
             Counts=Int32[3,3,1]) == sort!(PolACalculus(dc))


# function ParseRecord(r::GFF3.Record)::String

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
             "ftype=exon") == ParseRecord(record)
@test String("ID=;"*
             "Parent=ENST00000643283.1;"*
             "gene_id=ENSG00000126091.20;"*
             "transcript_id=ENST00000643283.1;"*
             "gene_type=protein_coding;"*
             "gene_name=ST3GAL3;"*
             "ftype=exon") == ParseRecord(record2)
@test String("ID=;"*
             "Parent=;"*
             "gene_id=ENSG00000126091.20;"*
             "transcript_id=ENST00000643283.1;"*
             "gene_type=protein_coding;"*
             "gene_name=ST3GAL3;"*
             "ftype=exon") == ParseRecord(record3)
@test String("ID=;"*
             "Parent=;"*
             "gene_id=;"*
             "transcript_id=;"*
             "gene_type=;"*
             "gene_name=;"*
             "ftype=exon") == ParseRecord(record4)
@test String("ID=;"*
             "Parent=;"*
             "gene_id=ENSG00000126091.20;"*
             "transcript_id=;"*
             "gene_type=;"*
             "gene_name=;"*
             "ftype=exon") == ParseRecord(record5)

# function GetIntervalSet(dframe::DataFrame)::IntervalCollection{String}

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
                PolyAReads=Int64[500,500,500])

@test IntervalCollection([Interval("chr1",150,150,Strand('+'),"Median=200.0;Min=100;Max=300;Counts=3")
                          Interval("chr2",160,160,Strand('-'),"Median=400.0;Min=300;Max=400;Counts=3")
                          Interval("chr3",170,170,Strand('.'),"Median=600.0;Min=600;Max=600;Counts=1")
                          ]) == GetIntervalSet(tdf)
