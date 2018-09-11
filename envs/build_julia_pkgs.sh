JULIA=julia
if ! $JULIA -e 'using ArgParse' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("ArgParse", v"0.5.0")'
    $JULIA -e 'using ArgParse'
fi

if ! $JULIA -e 'using AutoHashEquals' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("AutoHashEquals", v"0.2.0")'
    $JULIA -e 'using AutoHashEquals'
fi

if ! $JULIA -e 'using JLD' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("JLD", v"0.8.3")'
    $JULIA -e 'using JLD'
fi

if ! $JULIA -e 'using BioSequences' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("BioSequences", v"0.8.3")'
    $JULIA -e 'using BioSequences'
fi

if ! $JULIA -e 'using CodecZlib' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("CodecZlib", v"0.4.3")'
    $JULIA -e 'using CodecZlib'
fi

if ! $JULIA -e 'using BioAlignments' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("BioAlignments", v"0.3.0")'
    $JULIA -e 'using BioAlignments'
fi

if ! $JULIA -e 'using DataFrames' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("DataFrames", v"0.11.6")'
    $JULIA -e 'using DataFrames'
fi

if ! $JULIA -e 'using CSV' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("CSV", v"0.2.5")'
    $JULIA -e 'using CSV'
fi

if ! $JULIA -e 'using StringDistances' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("StringDistances", v"0.2.1")'
    $JULIA -e 'using StringDistances'
fi

if ! $JULIA -e 'using FastaIO' > /dev/null 2>&1; then
    $JULIA -e 'Pkg.add("FastaIO", v"0.4.1")'
    $JULIA -e 'using FastaIO'
fi
