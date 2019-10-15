# Eisenia

A package for probablistic graph genome assembly, analysis, and identification using kmers

WIP

How to turn the package into a runnable executable from the command line
```bash
mkdir ~/.julia/config
echo 'push!(LOAD_PATH, expanduser("~/Desktop/Microbes/Eisenia/src"))' > ~/.julia/config/startup.jl
echo 'export PATH="$HOME/Desktop/Microbes/Eisenia/bin:$PATH"'
```

```bash
mkdir data
wget --directory-prefix data/ ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/Escherichia_virus_phiX174/latest_assembly_versions/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz
```
