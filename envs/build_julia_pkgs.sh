# Installs Julia packages.
JULIA=julia
$JULIA -e 'Pkg.init()'

if ! $JULIA -e 'using ArgParse' > /dev/null 2>&1; then
    echo 'Julia pkg ArgParse is beeing installed ...'
    $JULIA -e 'Pkg.add("ArgParse", v"0.5.0")'
    $JULIA -e 'using ArgParse'
else echo 'Julia pkg ArgParse is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using AutoHashEquals' > /dev/null 2>&1; then
    echo 'Julia pkg AutoHashEquals is beeing installed ...'
    $JULIA -e 'Pkg.add("AutoHashEquals", v"0.2.0")'
    $JULIA -e 'using AutoHashEquals'
else echo 'Julia pkg AutoHashEquals is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using JLD' > /dev/null 2>&1; then
    echo 'Julia pkg JLD is beeing installed ...'
    $JULIA -e 'Pkg.add("JLD", v"0.8.3")'
    $JULIA -e 'using JLD'
else echo 'Julia pkg JLD is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using BioSequences' > /dev/null 2>&1; then
    echo 'Julia pkg BioSequences is beeing installed ...'
    $JULIA -e 'Pkg.add("BioSequences", v"0.8.3")'
    $JULIA -e 'using BioSequences'
else echo 'Julia pkg BioSequences is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using CodecZlib' > /dev/null 2>&1; then
    echo 'Julia pkg CodecZlib is beeing installed ...'
    $JULIA -e 'Pkg.add("CodecZlib", v"0.4.3")'
    $JULIA -e 'using CodecZlib'
else echo 'Julia pkg CodecZlib is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using BioAlignments' > /dev/null 2>&1; then
    echo 'Julia pkg BioAlignments is beeing installed ...'
    $JULIA -e 'Pkg.add("BioAlignments", v"0.3.0")'
    $JULIA -e 'using BioAlignments'
else echo 'Julia pkg BioAlignments is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using DataFrames' > /dev/null 2>&1; then
    echo 'Julia pkg DataFrames is beeing installed ...'
    $JULIA -e 'Pkg.add("DataFrames", v"0.11.6")'
    $JULIA -e 'using DataFrames'
else echo 'Julia pkg DataFrames is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using CSV' > /dev/null 2>&1; then
    echo 'Julia pkg CSV is beeing installed ...'
    $JULIA -e 'Pkg.add("CSV", v"0.2.5")'
    $JULIA -e 'using CSV'
else echo 'Julia pkg CSV is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using StringDistances' > /dev/null 2>&1; then
    echo 'Julia pkg StringDistances is beeing installed ...'
    $JULIA -e 'Pkg.add("StringDistances", v"0.2.1")'
    $JULIA -e 'using StringDistances'
else echo 'Julia pkg StringDistances is installed ... Nothing to be done.'
fi

if ! $JULIA -e 'using FMIndexes' > /dev/null 2>&1; then
    echo 'Julia pkg FMIndexes is beeing installed ...'
    $JULIA -e 'Pkg.add("FMIndexes", v"0.1.0")'
    $JULIA -e 'using FMIndexes'
else echo 'Julia pkg FMIndexes is installed ... Nothing to be done.'
fi
