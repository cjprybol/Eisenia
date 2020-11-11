# build distance matrix should be split by sequence type RNA & DNA or protein

# protein should require a scoring model

# dna/rna should require a kmer

# fasta -> distance matrix

# list of fastas -> distance matrix

# identifiers + distance matrix -> tree

# tree + identifiers -> newick

score_model = BioAlignments.AffineGapScoreModel(BioAlignments.BLOSUM45, gap_open=0, gap_extend=0)
distance_matrix = zeros(Float64, DataFrames.nrow(phage_table), DataFrames.nrow(phage_table))
ProgressMeter.@showprogress for i in 1:DataFrames.nrow(phage_table)
    for j in i+1:DataFrames.nrow(phage_table)
        s1 = phage_table[i, "protein_sequence"]
        s2 = phage_table[j, "protein_sequence"]
        d = 1/BioAlignments.pairalign(BioAlignments.GlobalAlignment(), s1, s2, score_model).value
        distance_matrix[i, j] = distance_matrix[j, i] = d
    end
end
distance_matrix

sequence_names = map(row -> join([replace(row["protein_name"], r"[^\w\d\._-]" => ""), row["phage"], row["genus"]], "-"), collect(DataFrames.eachrow(phage_table)))

# this is equivalent to UPGMA
tree = Clustering.hclust(distance_matrix, linkage=:average, branchorder=:optimal)
return identifiers, tree, distance matrix


newick = Dict()
for row in 1:size(tree.merges, 1)
    left, right = tree.merges[row, :]
    if left < 0
        l = sequence_names[abs(left)]
    else
        l = newick[left]
    end
    if right < 0
        r = sequence_names[abs(right)]
    else
        r = newick[right]
    end
    height = tree.heights[row]
    newick[row] = "($l:$height,$r:$height)"
end
open("$(DIR)/proteins.myoviridae.newick", "w") do io
    println(io, newick[size(tree.merges, 1)] * ";")
end