
dc = Dict("Key1" => [1,2], "Key2" => [1,2])
k1 = "Key1"
k2 = "Key2"
k3 = "Key3"
vl1 = Int16(1)
vl2 = Int16(2)
vl3 = Int16(3)

@test Dict("Key1" => [1,2,1], "Key2" => [1,2]) == Clust!(dc, k1, vl1)
@test Dict("Key1" => [1,2,1], "Key2" => [1,2,2]) == Clust!(dc, k2, vl2)
@test Dict("Key1" => [1,2,1], "Key2" => [1,2,2], "Key3" => [3]) == Clust!(dc, k3, vl3)
