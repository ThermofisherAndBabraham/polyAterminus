#!/usr/bin/env julia

try using ArgParse
catch
    println(STDERR, "No ArgParse package installed, installing...")
    Pkg.add("ArgParse")
end

try using DataStructures
catch
    println(STDERR, "No DataStructures package installed, installing...")
    Pkg.add("DataStructures")
end

try using DataFrames
catch
    println(STDERR, "no DataFrames package installed, installing.")
    Pkg.add("DataFrames")
end

try using CSV
catch
    println(STDERR, "no CSV package installed, installing.")
    Pkg.add("CSV")
end

try using BioSequences
catch
    println(STDERR, "no BioSequences package installed, installing.")
    Pkg.add("BioSequences")
end

try using GenomicFeatures
catch
    println(STDERR, "no GenomicFeatures package installed, installing.")
    Pkg.add("GenomicFeatures")
end

try using BioAlignments
catch
    println(STDERR, "no BioAlignments package installed, installing.")
    Pkg.add("BioAlignments")
end
