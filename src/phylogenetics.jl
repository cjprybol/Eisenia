# build distance matrix should be split by sequence type RNA & DNA or protein

# protein should require a scoring model

# dna/rna should require a kmer

# fasta -> distance matrix
# list of fastas -> distance matrix

# sanitize protein identifiers
# should grab this from fasta & join id + description
# may make sanitize a true false option, I want it for compatability with dendroscope
sequence_names = map(row -> join([replace(row["protein_name"], r"[^\w\d\._-]" => ""), row["phage"], row["genus"]], "-"), collect(DataFrames.eachrow(phage_table)))

function determine_distance(s1, s2, score_model)
    alignment = BioAlignments.pairalign(BioAlignments.GlobalAlignment(), s1, s2, score_model)
    mean_length = Statistics.mean([length(s1), length(s2)])
    length_normalized_score = alignment.value/mean_length
    # negative scores get higher distances
    # possible distances range from (0 - Inf)
    d = 1/MathConstants.e^(length_normalized_score)
end

# allow user to pass score model | just distance matrix 
score_model = 
BioAlignments.AffineGapScoreModel(
    BioAlignments.BLOSUM45, 
    gap_open=minimum(BioAlignments.BLOSUM45.data), 
    gap_extend=Statistics.median(BioAlignments.BLOSUM45.data))


distance_matrix = zeros(Float64, DataFrames.nrow(phage_table), DataFrames.nrow(phage_table))
ProgressMeter.@showprogress for i in 1:DataFrames.nrow(phage_table)
    for j in i+1:DataFrames.nrow(phage_table)
        s1 = phage_table[i, "protein_sequence"]
        s2 = phage_table[j, "protein_sequence"]
        d = determine_distance(s1, s2, score_model)
        distance_matrix[i, j] = distance_matrix[j, i] = d
    end
end
distance_matrix

# distance matrix -> tree
tree = Clustering.hclust(distance_matrix, linkage=:average, branchorder=:optimal)

# tree + identifiers -> newick
# take identifiers, tree, make newick
# use fasta + newick as default name
# allow user passable name
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



return identifiers, tree, distance matrix


