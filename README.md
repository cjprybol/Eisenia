# Eisenia

A package for kmer-based probablistic graph genome assembly and analysis

WIP

How to turn the package into a runnable executable from the command line
```bash
echo 'export PATH="$HOME/Desktop/Microbes/Eisenia/bin:$PATH"'
```

```bash
mkdir data
wget --directory-prefix data/ ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Escherichia_virus_phiX174/latest_assembly_versions/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
```

Making a package
```julia
julia> using Pkg

julia> Pkg.generate("Eisenia")
Generating project Eisenia:
    Eisenia/Project.toml
    Eisenia/src/Eisenia.jl
Dict{String,Base.UUID} with 1 entry:
  "Eisenia" => UUID("ff91eec8-cf20-4569-aeea-763e189e8d12")
# move all of my code into it
(v1.2) pkg> activate .
Activating environment at `~/Repos/Eisenia/Project.toml`

# add things
```

developing packages outside of the standard Julia locations
```bash
mkdir ~/.julia/config
echo 'push!(LOAD_PATH, expanduser("~/Repos/Eisenia/src"))' > ~/.julia/config/startup.jl
```