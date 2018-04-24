
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
