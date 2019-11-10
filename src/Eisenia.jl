module Eisenia

# import Reexport
# import ArgParse
# for finding the distances between kmers
import BioAlignments
# @reexport import BioSequences
# for all sequence useage
import BioSequences
# for reading and writing gzipped fasta/fastq
import CodecZlib
# for managing colors in graph plots
import Colors
# for colorschemes used in graph plots
# import ColorSchemes
import DataFrames
# import DataStructures
# import Dates
# import Distances
# import Distributions
import GLM
# import IterTools
# import JLD2
# @reexport import LightGraphs
# for representing genome graphs
import LightGraphs
# for drawing graph plots
import Luxor
# import Mmap
# for handling metadata with genome graphs
import MetaGraphs
import Plots
# import ProgressMeter
import Random
# import StaticArrays
import Statistics
# import StatsBase

const NUCLEOTIDES =
    [BioSequences.DNA_A,
     BioSequences.DNA_C,
     BioSequences.DNA_G, 
     BioSequences.DNA_T]

"""
document me
"""
function sequence_to_stranded_path(stranded_kmers, sequence)
    k = length(first(stranded_kmers))
    path = Vector{Pair{Int, Bool}}()
    for (i, kmer) in BioSequences.each(BioSequences.DNAKmer{k}, sequence)
        kmer_index = findfirst(stranded_kmer -> kmer == stranded_kmer, stranded_kmers)
        orientation = true
        push!(path, kmer_index => orientation)
    end
    return path
end

"""
document me
"""
function sequence_to_canonical_path(canonical_kmers, sequence)
    k = length(first(canonical_kmers))
    path = Vector{Pair{Int, Bool}}()
    for (i, kmer) in BioSequences.each(BioSequences.DNAKmer{k}, sequence)
        canonical_kmer = BioSequences.canonical(kmer)
        kmer_index = findfirst(_canonical_kmer -> _canonical_kmer == canonical_kmer, canonical_kmers)
        orientation = kmer == canonical_kmer
        push!(path, kmer_index => orientation)
    end
    return path
end

"""
document me
"""
function path_to_sequence(kmers, path)
    sequence = copy(kmers[first(path)])
    for i in 2:length(path)
        push!(sequence, kmers[path[i]][end])
    end
    return sequence
end

# """
#     edge_probability(stranded_kmer_graph, edge)
#
# Compute the probability of an edge relative to all edges of a given source vertex
# """
function edge_probability(stranded_kmer_graph, edge)
    neighbors = LightGraphs.outneighbors(stranded_kmer_graph, edge.src)
    neighbor_under_consideration = findfirst(neighbor -> neighbor == edge.dst, neighbors)
    edge_weights = [length(stranded_kmer_graph.eprops[LightGraphs.Edge(edge.src, neighbor)][:coverage]) for neighbor in neighbors]
    if sum(edge_weights) == 0
        p = 0.0
    else
        edge_probabilities = edge_weights ./ sum(edge_weights)
        p = edge_probabilities[neighbor_under_consideration]
    end
    return p
end

# """
#     viterbi_maximum_likelihood_traversal(stranded_kmer_graph; error_rate=1-1/k + 1))
#
#
#
# ```jldoctest
# import Primes
# julia> for prime in Primes.primes(3, 101)
#            println("$prime assumes accuracy of $(1-1/(prime + 1))")
#        end
# 3 assumes accuracy of 0.6666666666666667
# 5 assumes accuracy of 0.8
# 7 assumes accuracy of 0.8571428571428572
# 11 assumes accuracy of 0.9090909090909091
# 13 assumes accuracy of 0.9230769230769231
# 17 assumes accuracy of 0.9411764705882353
# 19 assumes accuracy of 0.9473684210526316
# 23 assumes accuracy of 0.9565217391304348
# 29 assumes accuracy of 0.9655172413793104
# 31 assumes accuracy of 0.967741935483871
# 37 assumes accuracy of 0.972972972972973
# 41 assumes accuracy of 0.975609756097561
# 43 assumes accuracy of 0.9767441860465116
# 47 assumes accuracy of 0.9787234042553191
# 53 assumes accuracy of 0.9811320754716981
# 59 assumes accuracy of 0.9830508474576272
# 61 assumes accuracy of 0.9836065573770492
# 67 assumes accuracy of 0.9850746268656716
# 71 assumes accuracy of 0.9859154929577465
# 73 assumes accuracy of 0.9863013698630136
# 79 assumes accuracy of 0.9873417721518988
# 83 assumes accuracy of 0.9879518072289156
# 89 assumes accuracy of 0.9887640449438202
# 97 assumes accuracy of 0.9896907216494846
# 101 assumes accuracy of 0.9900990099009901
# ```
# """
function viterbi_maximum_likelihood_traversals(stranded_kmer_graph;
                                               error_rate::Float64=1/(stranded_kmer_graph.gprops[:k] + 1),
                                               verbosity::String="dataset")
    @assert verbosity in ["debug", "reads", "dataset"]
    if error_rate >= .5
        error("Error rate >= 50%. Did you enter the accuracy by mistake?")
    end

    if verbosity in ["debug", "reads", "dataset"]
        println("computing kmer counts...")
    end
    stranded_kmer_counts = [length(stranded_kmer_graph.vprops[vertex][:coverage]) for vertex in LightGraphs.vertices(stranded_kmer_graph)]
    if verbosity in ["debug", "reads", "dataset"]
        println("computing kmer state likelihoods...")
    end
    stranded_kmer_likelihoods = stranded_kmer_counts ./ sum(stranded_kmer_counts)
    accuracy = 1 - error_rate

    if verbosity in ["debug"]
        println("STATE LIKELIHOODS:")
        println("\tkmer\tcount\tlikelihood")
        for vertex in LightGraphs.vertices(stranded_kmer_graph)
            kmer = stranded_kmer_graph.gprops[:stranded_kmers][vertex]
            count = stranded_kmer_counts[vertex]
            likelihood = stranded_kmer_likelihoods[vertex]
            println("\t$kmer\t$count\t$likelihood")
        end
    end
    if verbosity in ["debug", "reads", "dataset"]
        println("finding shortest paths between kmers...")
    end
    shortest_paths = enumerate_paths(floyd_warshall_shortest_paths(stranded_kmer_graph))
    K = stranded_kmer_graph.gprops[:K]
    for K1 in 1:K
        for K2 in 1:K
            if K1 != K2
                shortest_path = shortest_paths[K1][K2]
                path_likelihood = 1.0
                for ui in 1:length(shortest_path)-1
                    u = shortest_path[ui]
                    v = shortest_path[ui + 1]
                    # likelihood of the transition
                    path_likelihood *= edge_probability(stranded_kmer_graph, LightGraphs.Edge(u, v))
                end
                if path_likelihood == 0.0
                    shortest_paths[K1][K2] = Vector{Int}()
                end
            elseif K1 == K2
                # the shortest path from a kmer to itself is an insertion (no edge)
                # so need to manually check for self loops
                if LightGraphs.has_edge(stranded_kmer_graph, LightGraphs.Edge(K1, K2))
                    if edge_probability(stranded_kmer_graph, LightGraphs.Edge(K1, K2)) != 0.0
                        shortest_paths[K1][K2] = [K1, K2]
                    else
                        shortest_paths[K1][K2] = Vector{Int}()
                    end
                # otherwise, check to see if any outneighbors connect back to the kmer
                else
                    connected_outneighbors = filter(outneighbor -> has_path(stranded_kmer_graph, outneighbor, K2), LightGraphs.outneighbors(stranded_kmer_graph, K1))
                    if !isempty(connected_outneighbors)
                        outneighbor_cycles = [[K1, shortest_paths[outneighbor][K2]...] for outneighbor in connected_outneighbors]
                        cycle_likelihoods = ones(length(outneighbor_cycles))
                        for (i, cycle) in enumerate(outneighbor_cycles)
                            for ui in 1:length(cycle)-1
                                u = cycle[ui]
                                v = cycle[ui + 1]
                                # likelihood of the transition
                                cycle_likelihoods[i] *= edge_probability(stranded_kmer_graph, LightGraphs.Edge(u, v))
                            end
                            # include likelihoods of states
                            for vertex in cycle[2:end-1]
                                cycle_likelihoods[i] *= stranded_kmer_likelihoods[vertex]
                            end
                        end
                        path_likelihood = maximum(cycle_likelihoods)
                        max_likelihood_cycle_indices = findall(cycle_likelihoods .== path_likelihood)
                        shortest_paths[K1][K2] = outneighbor_cycles[first(max_likelihood_cycle_indices)]
                    else
                        shortest_paths[K1][K2] = Vector{Int}()
                    end
                end
            end
            if length(shortest_paths[K1][K2]) == 1
                shortest_paths[K1][K2] = Vector{Int}()
            end
        end
    end

    if verbosity in ["debug"]
        for K1 in 1:K
            for K2 in 1:K
                println("\t$K1\t$K2\t$(shortest_paths[K1][K2])")
            end
        end
    end

    total_bases_observed = 0
    total_edits_accepted = 0

    corrected_observations = BioSequences.FASTA.Record[]
    if verbosity in ["debug", "reads", "dataset"]
        println("finding viterbi maximum likelihood paths for observed sequences...")
    end
    # p = Progress(length(stranded_kmer_graph.gprops[:observed_paths]))
    for (observation_index, observed_path) in enumerate(stranded_kmer_graph.gprops[:observed_paths])
        if verbosity in ["debug", "reads"]
            println("\nevaluating sequence $observation_index of $(length(stranded_kmer_graph.gprops[:observed_paths]))")
        end
        # consider switching to log transform
        kmer_likelihoods = zeros(LightGraphs.nv(stranded_kmer_graph), length(observed_path))
        kmer_arrival_paths = Array{Vector{Int}}(undef, LightGraphs.nv(stranded_kmer_graph), length(observed_path))
        edit_distances = zeros(Int, LightGraphs.nv(stranded_kmer_graph), length(observed_path))
        observed_kmer_index = observed_path[1]
        observed_kmer_sequence = stranded_kmer_graph.gprops[:stranded_kmers][observed_kmer_index]
        for hidden_kmer_index in LightGraphs.vertices(stranded_kmer_graph)
            hidden_kmer_sequence = stranded_kmer_graph.gprops[:stranded_kmers][hidden_kmer_index]
            alignment_result = pairalign(LevenshteinDistance(), observed_kmer_sequence, hidden_kmer_sequence)
            number_of_matches = count_matches(alignment(alignment_result))
            number_of_edits = stranded_kmer_graph.gprops[:k] - number_of_matches
            kmer_likelihoods[hidden_kmer_index, 1] = stranded_kmer_likelihoods[hidden_kmer_index]
            for match in 1:number_of_matches
                kmer_likelihoods[hidden_kmer_index, 1] *= accuracy
            end
            for edit in 1:number_of_edits
                kmer_likelihoods[hidden_kmer_index, 1] *= error_rate
            end
            kmer_arrival_paths[hidden_kmer_index, 1] = Vector{Int}()
            edit_distances[hidden_kmer_index, 1] = number_of_edits
        end
        kmer_likelihoods[:, 1] ./= sum(kmer_likelihoods[:, 1])
        # from here on, all probabilities are log transformed
        kmer_likelihoods[:, 1] .= log.(kmer_likelihoods[:, 1])
        if verbosity in ["debug"]
            println("\tconsidering path state 1")
            println("\t\tobserved kmer $observed_kmer_sequence")
            println("\t\tInitial state log likelihoods:")
            for line in split(repr(MIME("text/plain"), kmer_likelihoods[:, 1]), '\n')
                println("\t\t\t$line")
            end
        end
        for observed_path_index in 2:length(observed_path)
            observed_kmer_index = observed_path[observed_path_index]
            observed_base = stranded_kmer_graph.gprops[:stranded_kmers][observed_kmer_index][end]

            if verbosity in ["debug"]
                println("\tconsidering path state $observed_path_index")
                println("\t\tobserved base $observed_base")
            end

            MATCH = 1
            MISMATCH = 2
            DELETION = 3
            INSERTION = 4
            arrival_likelihoods = ones(K, 4)
            arrival_paths = fill(Vector{Int}(), K, 4)

            for K2 in 1:K
                kmer_base = stranded_kmer_graph.gprops[:stranded_kmers][K2][end]
                base_is_match = kmer_base == observed_base

                maximum_likelihood = log(0.0)
                maximum_likelihood_path = Vector{Int}()
                maximum_likelihood_edit_distance = 0

                for K1 in 1:K
                    shortest_path = shortest_paths[K1][K2]
                    if length(shortest_path) >= 2
                        edit_distance = Int(!base_is_match) + length(shortest_path) - 2
                        if edit_distance == 0
                            p = kmer_likelihoods[K1, observed_path_index-1] +
                                log(accuracy) + log(stranded_kmer_likelihoods[K2])
                        else
                            p = kmer_likelihoods[K1, observed_path_index-1] +
                                log(error_rate^edit_distance) + log(stranded_kmer_likelihoods[K2])
                        end
                        edit_distance += edit_distances[K1, observed_path_index-1]
                    else
                        p = log(0.0)
                    end
                    if K1 == K2 # consider insertion
                        # in theory, I don't think we should care if the base
                        # matches or not because it's an inserted & erroneous
                        # base, but in practice it's necessary to balance
                        # insertion probabilities with deletion probabilities
                        insertion_p = kmer_likelihoods[K1, observed_path_index-1] +
                                      log(error_rate^(1 + Int(!base_is_match))) + log(stranded_kmer_likelihoods[K2])
                        if insertion_p > p
                            p = insertion_p
                            edit_distance = edit_distances[K1, observed_path_index-1] + 1
                            shortest_path = [K2]
                        end
                    end
                    if p > maximum_likelihood
                        maximum_likelihood = p
                        maximum_likelihood_path = shortest_path
                        maximum_likelihood_edit_distance = edit_distance
                    end
                end
                kmer_likelihoods[K2, observed_path_index] = maximum_likelihood
                kmer_arrival_paths[K2, observed_path_index] = maximum_likelihood_path
                edit_distances[K2, observed_path_index] = maximum_likelihood_edit_distance
            end

            if verbosity in ["debug"]
                println("\t\tkmer log likelihoods")
                for line in split(repr(MIME("text/plain"), kmer_likelihoods), '\n')
                    println("\t\t\t$line")
                end
                println("\t\tarrival paths")
                for line in split(repr(MIME("text/plain"), kmer_arrival_paths), '\n')
                    println("\t\t\t$line")
                end
            end
        end

        if verbosity in ["debug"]
            println("\n\tInputs for viterbi maximum likelihood traversal evaluation:")
            println("\t\tkmer log likelihoods")
            for line in split(repr(MIME("text/plain"), kmer_likelihoods), '\n')
                println("\t\t\t$line")
            end
            println("\t\tkmer arrival paths")
            for line in split(repr(MIME("text/plain"), kmer_arrival_paths), '\n')
                println("\t\t\t$line")
            end
            println("\t\tedit distances")
            for line in split(repr(MIME("text/plain"), edit_distances), '\n')
                println("\t\t\t$line")
            end
        end

        ## backtrack
        maximum_likelihood_path_value = maximum(kmer_likelihoods[:, end])
        maximum_likelihood_path_indices = findall(kmer_likelihoods[:, end] .== maximum_likelihood_path_value)
        # if multiple paths are tied, randomly choose one
        maximum_likelihood_path_index = rand(maximum_likelihood_path_indices)
        maximum_likelihood_edit_distance = edit_distances[maximum_likelihood_path_index, end]

        if length(kmer_arrival_paths[maximum_likelihood_path_index, end]) > 0
            maximum_likelihood_path = last(kmer_arrival_paths[maximum_likelihood_path_index, end])
            for observed_path_index in length(observed_path):-1:1
                maximum_likelihood_arrival_path = kmer_arrival_paths[maximum_likelihood_path_index, observed_path_index]
                maximum_likelihood_path = vcat(maximum_likelihood_arrival_path[1:end-1], maximum_likelihood_path)
                maximum_likelihood_path_index = first(maximum_likelihood_path)
            end
        else
            maximum_likelihood_path = [maximum_likelihood_path_index]
        end
        observed_sequence = path_to_sequence(stranded_kmer_graph.gprops[:stranded_kmers], observed_path)
        maximum_likelihood_sequence = path_to_sequence(stranded_kmer_graph.gprops[:stranded_kmers], maximum_likelihood_path)
        if verbosity in ["debug", "reads"]
            println("\tobserved sequence                 $observed_sequence")
            println("\tmaximum likelihood sequence       $maximum_likelihood_sequence")
            println("\tmaximum likelihood edit distance  $maximum_likelihood_edit_distance")
        end
        total_bases_observed += length(observed_sequence)
        total_edits_accepted += maximum_likelihood_edit_distance
        id = stranded_kmer_graph.gprops[:observation_ids][observation_index]
        kmer_stamped_id = id * "_" * string(stranded_kmer_graph.gprops[:k])
        push!(corrected_observations, BioSequences.FASTA.Record(kmer_stamped_id, maximum_likelihood_sequence))
        # progress meter
        # next!(p)
    end
    if verbosity in ["debug", "reads", "dataset"]
        println("\nDATASET STATISTICS:")
        println("\tassumed error rate    $(error_rate * 100)%")
        println("\ttotal bases observed  $total_bases_observed")
        println("\ttotal edits accepted  $total_edits_accepted")
        println("\tinferred error rate   $((total_edits_accepted/total_bases_observed) * 100)%")
    end
    return corrected_observations
end

function plot_stranded_kmer_graph(stranded_kmer_graph; filename=Random.randstring(Random.MersenneTwister(Int(round(time()))), 3) * ".svg")
    canonical_vertices = [i for (i, kmer) in enumerate(stranded_kmer_graph.gprops[:stranded_kmers]) if stranded_kmer_graph.gprops[:reverse_complement_map][i] > i]
    vertex_ordered_pairs = [(plus => minus) for (plus, minus) in zip(canonical_vertices, stranded_kmer_graph.gprops[:reverse_complement_map][canonical_vertices])]
    stranded_vertex_to_unstranded_vertex_map = [findfirst(vertex_ordered_pair -> vertex in vertex_ordered_pair, vertex_ordered_pairs) for vertex in LightGraphs.vertices(stranded_kmer_graph)]

    canonical_kmer_graph = LightGraphs.SimpleGraph(length(vertex_ordered_pairs))

    for vertex_ordered_pair in vertex_ordered_pairs
        plus = first(vertex_ordered_pair)
        minus = last(vertex_ordered_pair)
        for plus_strand_neighbor in [LightGraphs.inneighbors(stranded_kmer_graph, plus)..., LightGraphs.outneighbors(stranded_kmer_graph, plus)...]
            LightGraphs.add_edge!(canonical_kmer_graph, stranded_vertex_to_unstranded_vertex_map[plus], stranded_vertex_to_unstranded_vertex_map[plus_strand_neighbor])
        end
        for minus_strand_neighbor in [LightGraphs.inneighbors(stranded_kmer_graph, minus)..., LightGraphs.outneighbors(stranded_kmer_graph, minus)...]
            LightGraphs.add_edge!(canonical_kmer_graph, stranded_vertex_to_unstranded_vertex_map[minus], stranded_vertex_to_unstranded_vertex_map[minus_strand_neighbor])
        end
    end

    textsize = 12
    radius = textsize/2

    current_global_x_min = 0
    current_global_y_min = 0
    current_global_x_max = 0
    current_global_y_max = 0

    stranded_vertex_coordinates = [Dict{Symbol, Any}() for vertex in LightGraphs.vertices(stranded_kmer_graph)]

    for connected_component in LightGraphs.connected_components(canonical_kmer_graph)
        # reset xmin to left-align new contigs
        current_contig_x_min = current_global_x_min
        current_contig_x_max = current_global_x_min
        current_contig_y_min = current_global_y_max
        current_contig_y_max = current_global_y_max
        ## TODO could also achor by vertex with the smallest indegree
        anchor = minimum(connected_component)
        canonical_bfs_tree = LightGraphs.bfs_tree(canonical_kmer_graph, anchor)
        vertices_in_current_depth_of_field = [anchor]
        while !isempty(vertices_in_current_depth_of_field)
            vertices_in_next_depth_of_field = Vector{Int}()
            current_vertex_x_min = current_contig_x_max
            current_vertex_y_min = current_contig_y_min
            max_coverage = 0
            for vertex in vertices_in_current_depth_of_field
                vertex_ordered_pair = vertex_ordered_pairs[vertex]
                plus = first(vertex_ordered_pair)
                minus = last(vertex_ordered_pair)
                plus_coverage = length(stranded_kmer_graph.vprops[plus][:coverage])
                minus_coverage = length(stranded_kmer_graph.vprops[minus][:coverage])
                vertex_max_coverage = max(plus_coverage, minus_coverage)
                if vertex_max_coverage > max_coverage
                    max_coverage = vertex_max_coverage
                end
            end
            for vertex in vertices_in_current_depth_of_field
                vertex_ordered_pair = vertex_ordered_pairs[vertex]
                plus = first(vertex_ordered_pair)
                minus = last(vertex_ordered_pair)
                plus_coverage = length(stranded_kmer_graph.vprops[plus][:coverage])
                minus_coverage = length(stranded_kmer_graph.vprops[minus][:coverage])
                above_y_offset = 10 + plus_coverage + radius
                below_y_offset = 10 + minus_coverage + radius
                before_x_offset = 10 + max_coverage
                after_x_offset = 10 + max_coverage

                svcp = stranded_vertex_coordinates[plus]
                svcp[:height] = textsize + plus_coverage
                svcp[:width] = (length(stranded_kmer_graph.gprops[:stranded_kmers][plus]) + 1) * textsize
                svcp[:xmin] = current_vertex_x_min + before_x_offset
                svcp[:xmax] = svcp[:xmin] + svcp[:width]
                svcp[:ymin] = current_vertex_y_min + above_y_offset
                svcp[:ymax] = svcp[:ymin] + svcp[:height]
                svcp[:xin] = svcp[:xmin]
                svcp[:xout] = svcp[:xmax]
                svcp[:coverage_to_y_map] = [svcp[:ymin] + radius + i for i in 1:plus_coverage]
                svcp[:center] = Luxor.Point(Statistics.mean((svcp[:xmin], svcp[:xmax])), Statistics.mean((svcp[:ymin], svcp[:ymax])))

                svcm = stranded_vertex_coordinates[minus]
                svcm[:height] = textsize + minus_coverage
                svcm[:width] = (length(stranded_kmer_graph.gprops[:stranded_kmers][minus]) + 1) * textsize
                svcm[:xmin] = current_vertex_x_min + before_x_offset
                svcm[:xmax] = svcm[:xmin] + svcm[:width]
                svcm[:ymin] = svcp[:ymax] + 1
                svcm[:ymax] = svcm[:ymin] + svcm[:height]
                svcm[:xin] = svcm[:xmax]
                svcm[:xout] = svcm[:xmin]
                svcm[:coverage_to_y_map] = [svcm[:ymin] + radius + i for i in 1:minus_coverage]
                svcm[:center] = Luxor.Point(Statistics.mean((svcm[:xmin], svcm[:xmax])), Statistics.mean((svcm[:ymin], svcm[:ymax])))

                @assert svcp[:width] == svcm[:width]
                @assert svcp[:xmin] == svcm[:xmin]
                if svcm[:xmax] + after_x_offset + radius > current_contig_x_max
                    current_contig_x_max = svcm[:xmax] + after_x_offset + radius
                end
                if svcm[:ymax] + minus_coverage + radius > current_contig_y_max
                    current_contig_y_max = svcm[:ymax] + minus_coverage + radius
                end
                current_vertex_y_min = svcm[:ymax] + minus_coverage
                for neighbor in LightGraphs.outneighbors(canonical_bfs_tree, vertex)
                    push!(vertices_in_next_depth_of_field, neighbor)
                end
            end
            current_vertex_x_min = current_contig_x_max
            vertices_in_current_depth_of_field = sort!(unique(vertices_in_next_depth_of_field))
        end
        if current_contig_x_max > current_global_x_max
            current_global_x_max = current_contig_x_max
        end
        if current_contig_y_max > current_global_y_max
            current_global_y_max = current_contig_y_max
        end
    end

    Luxor.Drawing(current_global_x_max, current_global_y_max, filename)
    Luxor.background("white")
    Luxor.setline(1)
    ncolors = length(stranded_kmer_graph.gprops[:observation_ids])
    color_list = Colors.distinguishable_colors(ncolors+1, Colors.RGB(1,1,1))[2:end]
    for sequence_record_index in 1:length(stranded_kmer_graph.gprops[:observed_paths])
        sequence_path = stranded_kmer_graph.gprops[:observed_paths][sequence_record_index]
        sequence_color = color_list[stranded_kmer_graph.gprops[:observation_color_map][sequence_record_index]]
        i = 1
        ui = first(sequence_path[i])
        uiis = findall(observation -> observation == (sequence_record_index => (i => true)), stranded_kmer_graph.vprops[ui][:coverage])
        @assert length(uiis) == 1
        uii = first(uiis)
        ux = stranded_vertex_coordinates[ui][:xout]
        uy = stranded_vertex_coordinates[ui][:coverage_to_y_map][uii]
        if ui == first(vertex_ordered_pairs[stranded_vertex_to_unstranded_vertex_map[ui]])
            source_strand = :plus
        else
            source_strand = :minus
        end
        Luxor.setcolor(sequence_color)
        for i in 2:length(sequence_path)
            vi = first(sequence_path[i])
            stranded_kmer_graph.vprops[vi][:coverage]
            viis = findall(observation -> observation == (sequence_record_index => (i => true)), stranded_kmer_graph.vprops[vi][:coverage])
            @assert length(viis) == 1
            vii = first(viis)
            vy = stranded_vertex_coordinates[vi][:coverage_to_y_map][vii]
            vx = stranded_vertex_coordinates[vi][:xin]
            if ux == vx
                x_direction = :even
            elseif ux < vx
                x_direction = :forward
            else
                x_direction = :backward
            end
            if uy == vy
                y_direction = :even
            elseif uy < vy
                y_direction = :down
            else
                y_direction = :up
            end
            if vi == first(vertex_ordered_pairs[stranded_vertex_to_unstranded_vertex_map[vi]])
                destination_strand = :plus
            else
                destination_strand = :minus
            end
            if x_direction == :even && source_strand == :plus && destination_strand == :minus && y_direction in (:even, :down)
                u_radius = radius + (length(stranded_vertex_coordinates[ui][:coverage_to_y_map]) - uii)
                v_radius = radius + vii - 1
                a = Luxor.Point(ux, uy)
                b = Luxor.Point(ux + u_radius, uy)
                c = Luxor.Point(vx + v_radius, vy)
                d = Luxor.Point(vx, vy)
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            elseif x_direction == :even && source_strand == :plus && destination_strand == :minus && y_direction == :up
                u_radius = radius + uii - 1
                v_radius = radius + (length(stranded_vertex_coordinates[vi][:coverage_to_y_map]) - vii)
                a = Luxor.Point(ux, uy)
                b = Luxor.Point(ux + u_radius, uy)
                c = Luxor.Point(vx + v_radius, vy)
                d = Luxor.Point(vx, vy)
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            elseif x_direction == :even && source_strand == :minus && destination_strand == :plus && y_direction in (:even, :up)
                u_radius = radius + uii - 1
                v_radius = radius + (length(stranded_vertex_coordinates[vi][:coverage_to_y_map]) - vii)
                a = Luxor.Point(ux, uy)
                b = Luxor.Point(ux - u_radius, uy)
                c = Luxor.Point(vx - v_radius, vy)
                d = Luxor.Point(vx, vy)
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            elseif x_direction == :even && source_strand == :minus && destination_strand == :plus && y_direction == :down
                u_radius = radius + (length(stranded_vertex_coordinates[ui][:coverage_to_y_map]) - uii)
                v_radius = radius + vii - 1
                a = Luxor.Point(ux, uy)
                b = Luxor.Point(ux - u_radius, uy)
                c = Luxor.Point(vx - v_radius, vy)
                d = Luxor.Point(vx, vy)
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            elseif x_direction == :forward && source_strand == :plus && destination_strand == :plus
                u_radius = radius + uii - 1
                v_radius = radius + vii - 1
                a = Luxor.Point(ux, uy)
                b = Luxor.Point(Statistics.mean((ux, vx)), uy)
                c = Luxor.Point(Statistics.mean((ux, vx)), vy)
                d = Luxor.Point(vx, vy)
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            elseif x_direction == :forward && source_strand == :plus && destination_strand == :minus
                u_radius = radius + uii - 1
                v_radius = radius + (length(stranded_vertex_coordinates[vi][:coverage_to_y_map]) - vii)
                va = Luxor.Point(vx, vy)
                vb = Luxor.Point(vx + v_radius, vy)
                vc = Luxor.Point(vx + v_radius, vy + 2v_radius)
                vd = Luxor.Point(vx, vy + 2v_radius)
                Luxor.move(va); Luxor.curve(vb, vc, vd); Luxor.strokepath()
                v_extension = Point(stranded_vertex_coordinates[vi][:xout], vd.y)
                move(vd); curve(vd, v_extension, v_extension); Luxor.strokepath()
                a = v_extension
                d = Luxor.Point(ux, uy)
                b = Luxor.Point(Statistics.mean((a.x, d.x)), a.y)
                c = Luxor.Point(Statistics.mean((a.x, d.x)), d.y)
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            elseif x_direction == :forward && source_strand == :minus && destination_strand == :plus
                u_radius = radius + (length(stranded_vertex_coordinates[ui][:coverage_to_y_map]) - uii)
                v_radius = radius
                # downward, c - loop from the source strand
                ua = Luxor.Point(ux, uy)
                ub = ua + Luxor.Point(-u_radius, 0)
                uc = ua + Luxor.Point(-u_radius, 2u_radius)
                ud = ua + Luxor.Point(0, 2u_radius)
                Luxor.move(ua); Luxor.curve(ub, uc, ud); Luxor.strokepath()
                va = Luxor.Point(vx, vy)
                vb = Luxor.Point((stranded_vertex_coordinates[ui][:xmax] + vx)/2, vy)
                vc = Luxor.Point((stranded_vertex_coordinates[ui][:xmax] + vx)/2, ud.y)
                vd = Luxor.Point(stranded_vertex_coordinates[ui][:xmax], ud.y)
                Luxor.move(va); Luxor.curve(vb, vc, vd); Luxor.strokepath()
                Luxor.move(ud); Luxor.curve(ud, vd, vd); Luxor.strokepath()
            elseif x_direction == :forward && source_strand == :minus && destination_strand == :minus
                u_radius = radius + (length(stranded_vertex_coordinates[ui][:coverage_to_y_map]) - uii)
                v_radius = radius + (length(stranded_vertex_coordinates[vi][:coverage_to_y_map]) - vii)
                if stranded_vertex_coordinates[ui][:xin] == stranded_vertex_coordinates[vi][:xin] && ui != vi
                    va = Luxor.Point(vx, vy)
                    vb = va + Luxor.Point(v_radius, 0)
                    vc = va + Luxor.Point(v_radius, 2v_radius)
                    vd = va + Luxor.Point(0, 2v_radius)
                    Luxor.move(va); Luxor.curve(vb, vc, vd); Luxor.strokepath()
                    v_extension = Luxor.Point(stranded_vertex_coordinates[vi][:xout], vd.y)
                    Luxor.move(vd); Luxor.curve(vd, v_extension, v_extension); Luxor.strokepath()
                    a = Luxor.Point(ux, uy)
                    b = Luxor.Point(ux - u_radius, uy)
                    c = Luxor.Point(v_extension.x - v_radius, v_extension.y)
                    d = v_extension
                    Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                else
                    ua = Luxor.Point(ux, uy)
                    ub = ua + Luxor.Point(-u_radius, 0)
                    uc = ua + Luxor.Point(-u_radius, 2u_radius)
                    ud = ua + Luxor.Point(0, 2u_radius)
                    Luxor.move(ua); Luxor.curve(ub, uc, ud); Luxor.strokepath()
                    va = Luxor.Point(vx, vy)
                    vb = va + Luxor.Point(v_radius, 0)
                    vc = va + Luxor.Point(v_radius, 2v_radius)
                    vd = va + Luxor.Point(0, 2v_radius)
                    Luxor.move(va); Luxor.curve(vb, vc, vd); Luxor.strokepath()
                    if ui == vi
                        a = ud
                        b = Luxor.Point(Statistics.mean((ud.x, vd.x)), ud.y)
                        c = Luxor.Point(Statistics.mean((ud.x, vd.x)), vd.y)
                        d = vd
                        Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                    else
                        u_extension = Luxor.Point(stranded_vertex_coordinates[ui][:xin], ud.y)
                        Luxor.move(ud); Luxor.curve(ud, u_extension, u_extension); Luxor.strokepath()
                        v_extension = Luxor.Point(stranded_vertex_coordinates[vi][:xout], vd.y)
                        Luxor.move(vd); Luxor.curve(vd, v_extension, v_extension); Luxor.strokepath()
                        a = u_extension
                        d = v_extension
                        b = Luxor.Point(Statistics.mean((a.x, d.x)), a.y)
                        c = Luxor.Point(Statistics.mean((a.x, d.x)), d.y)
                        Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                    end
                end
            elseif x_direction == :backward && source_strand == :plus && destination_strand == :plus
                u_radius = radius + uii - 1
                v_radius = radius + vii - 1
                ua = Luxor.Point(ux, uy)
                ub = ua + Luxor.Point(u_radius, 0)
                uc = ua + Luxor.Point(u_radius, -2u_radius)
                ud = ua + Luxor.Point(0, -2u_radius)
                Luxor.move(ua); Luxor.curve(ub, uc, ud); Luxor.strokepath()
                if stranded_vertex_coordinates[ui][:xin] == stranded_vertex_coordinates[vi][:xin] && ui != vi
                    u_extension = Luxor.Point(stranded_vertex_coordinates[ui][:xin], ud.y)
                    Luxor.move(ud); Luxor.curve(ud, u_extension, u_extension); Luxor.strokepath()
                    a = u_extension
                    b = Luxor.Point(u_extension.x - u_radius, u_extension.y)
                    c = Luxor.Point(vx - v_radius, vy)
                    d = Luxor.Point(vx, vy)
                    Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                else
                    va = Luxor.Point(vx, vy)
                    vb = va + Luxor.Point(-v_radius, 0)
                    vc = va + Luxor.Point(-v_radius, -2v_radius)
                    vd = va + Luxor.Point(0, -2v_radius)
                    Luxor.move(va); Luxor.curve(vb, vc, vd); Luxor.strokepath()
                    if ui == vi
                        a = ud
                        b = Luxor.Point(Statistics.mean((ud.x, vd.x)), ud.y)
                        c = Luxor.Point(Statistics.mean((ud.x, vd.x)), vd.y)
                        d = vd
                        Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                    else
                        u_extension = Luxor.Point(stranded_vertex_coordinates[ui][:xin], ud.y)
                        Luxor.move(ud); Luxor.curve(ud, u_extension, u_extension); Luxor.strokepath()
                        v_extension = Luxor.Point(stranded_vertex_coordinates[vi][:xout], vd.y)
                        Luxor.move(vd); Luxor.curve(vd, v_extension, v_extension); Luxor.strokepath()
                        a = u_extension
                        d = v_extension
                        b = Luxor.Point(Statistics.mean((a.x, d.x)), a.y)
                        c = Luxor.Point(Statistics.mean((a.x, d.x)), d.y)
                        Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                    end
                end
            elseif x_direction == :backward && source_strand == :plus && destination_strand == :minus
                u_radius = radius + uii - 1
                v_radius = radius + vii - 1
                ua = Luxor.Point(ux, uy)
                ub = Luxor.Point(ux + u_radius, uy)
                uc = Luxor.Point(ux + u_radius, uy - 2u_radius)
                ud = Luxor.Point(ux, uy - 2u_radius)
                Luxor.move(ua); Luxor.curve(ub, uc, ud); Luxor.strokepath()
                a = ud
                d = Luxor.Point(stranded_vertex_coordinates[ui][:xin], ud.y)
                b = a
                c = d
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                va = Luxor.Point(vx, vy)
                vb = Luxor.Point(Statistics.mean((vx, d.x)), vy)
                vc = Luxor.Point(Statistics.mean((vx, d.x)), d.y)
                vd = d
                Luxor.move(va); Luxor.curve(vb, vc, vd); Luxor.strokepath()
            elseif x_direction == :backward && source_strand == :minus && destination_strand == :plus
                u_radius = radius + uii - 1
                v_radius = radius + vii - 1
                va = Luxor.Point(vx, vy)
                vb = Luxor.Point(vx - v_radius, vy)
                vc = Luxor.Point(vx - v_radius, vy - 2v_radius)
                vd = Luxor.Point(vx, vy - 2v_radius)
                Luxor.move(va); Luxor.curve(vb, vc, vd); Luxor.strokepath()
                a = vd
                d = Luxor.Point(stranded_vertex_coordinates[vi][:xout], vy - 2v_radius)
                b = a
                c = d
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
                ua = Luxor.Point(ux, uy)
                ub = Luxor.Point(Statistics.mean((ux, d.x)), uy)
                ud = d
                uc = Luxor.Point(Statistics.mean((ux, d.x)), d.y)
                Luxor.move(ua); Luxor.curve(ub, uc, ud); Luxor.strokepath()
            elseif x_direction == :backward && source_strand == :minus && destination_strand == :minus
                u_radius = radius + uii - 1
                v_radius = radius + vii - 1
                a = Luxor.Point(ux, uy)
                b = Luxor.Point(Statistics.mean((ux, vx)), uy)
                c = Luxor.Point(Statistics.mean((ux, vx)), vy)
                d = Luxor.Point(vx, vy)
                Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            else
                error("how'd we get here?")
            end
            ui = vi
            uii = vii
            ux = stranded_vertex_coordinates[ui][:xout]
            uy = stranded_vertex_coordinates[ui][:coverage_to_y_map][uii]
            if ui == first(vertex_ordered_pairs[stranded_vertex_to_unstranded_vertex_map[ui]])
                source_strand = :plus
            else
                source_strand = :minus
            end
        end
    end
    for vertex_ordered_pair in vertex_ordered_pairs
        plus = first(vertex_ordered_pair)
        minus = last(vertex_ordered_pair)
        stranded_vertex_coordinates[plus][:ymin]
        stranded_vertex_coordinates[minus][:ymax]
        @assert stranded_vertex_coordinates[plus][:xmin] == stranded_vertex_coordinates[minus][:xmin]
        @assert stranded_vertex_coordinates[plus][:xmax] == stranded_vertex_coordinates[minus][:xmax]
        stranded_vertex_coordinates[plus][:xmin]
        stranded_vertex_coordinates[plus][:xmax]
        center = Luxor.Point(Statistics.mean((stranded_vertex_coordinates[plus][:xmin], stranded_vertex_coordinates[plus][:xmax])),
                             Statistics.mean((stranded_vertex_coordinates[plus][:ymin], stranded_vertex_coordinates[minus][:ymax])))
        height = stranded_vertex_coordinates[minus][:ymax] - stranded_vertex_coordinates[plus][:ymin]
        width = stranded_vertex_coordinates[minus][:xmax] - stranded_vertex_coordinates[plus][:xmin]
        Luxor.setcolor("gray")
        Luxor.box(center, width, height, :fill)
        for (i, y) in enumerate(stranded_vertex_coordinates[plus][:coverage_to_y_map])
            record_index, (sequence_index, orientation) = stranded_kmer_graph.vprops[plus][:coverage][i]
            color_index = stranded_kmer_graph.gprops[:observation_color_map][record_index]
            Luxor.setcolor(color_list[color_index])
            Luxor.line(Luxor.Point(stranded_vertex_coordinates[plus][:xmin], y), Luxor.Point(stranded_vertex_coordinates[plus][:xmax], y), :stroke)
            Luxor.setcolor("white")
            Luxor.fontsize(1)
            Luxor.fontface("Menlo Bold")
            Luxor.text(string(sequence_index), Luxor.Point(stranded_vertex_coordinates[plus][:xin] + 1, y), valign=:middle, halign=:center)
            if sequence_index == 1
                Luxor.circle(Luxor.Point(stranded_vertex_coordinates[plus][:xin] + 2, y), 0.5, :fill)
            end
        end
        for (i, y) in enumerate(stranded_vertex_coordinates[minus][:coverage_to_y_map])
            record_index, (sequence_index, orientation) = stranded_kmer_graph.vprops[minus][:coverage][i]
            color_index = stranded_kmer_graph.gprops[:observation_color_map][record_index]
            Luxor.setcolor(color_list[color_index])
            Luxor.line(Luxor.Point(stranded_vertex_coordinates[minus][:xmin], y), Luxor.Point(stranded_vertex_coordinates[minus][:xmax], y), :stroke)
            Luxor.setcolor("white")
            Luxor.fontsize(1)
            Luxor.fontface("Menlo Bold")
            Luxor.text(string(sequence_index), Luxor.Point(stranded_vertex_coordinates[minus][:xin] - 1, y), valign=:middle, halign=:center)
            if sequence_index == 1
                Luxor.circle(Luxor.Point(stranded_vertex_coordinates[minus][:xin] - 2, y), 0.5, :fill)
            end
        end
        plus_sequence = stranded_kmer_graph.gprops[:stranded_kmers][plus]
        Luxor.setcolor("white")
        Luxor.fontsize(textsize)
        Luxor.fontface("Menlo")
        Luxor.text(string(plus_sequence), stranded_vertex_coordinates[plus][:center], valign=:middle, halign=:center)
        minus_sequence = stranded_kmer_graph.gprops[:stranded_kmers][minus]
        Luxor.text(string(minus_sequence), stranded_vertex_coordinates[minus][:center], valign=:middle, halign=:center)
        Luxor.setcolor("black")

        a = Luxor.Point(stranded_vertex_coordinates[plus][:xin], stranded_vertex_coordinates[plus][:ymin])
        b = Luxor.Point(stranded_vertex_coordinates[plus][:xin], stranded_vertex_coordinates[plus][:ymin] + radius + 0.5)
        c = Luxor.Point(stranded_vertex_coordinates[plus][:xin] + radius, stranded_vertex_coordinates[plus][:ymin] + radius + 0.5)
        Luxor.poly([a, b, c], :fill)
        a = Luxor.Point(stranded_vertex_coordinates[plus][:xin], stranded_vertex_coordinates[plus][:ymax])
        b = Luxor.Point(stranded_vertex_coordinates[plus][:xin], stranded_vertex_coordinates[plus][:ymax] - radius + 0.5)
        c = Luxor.Point(stranded_vertex_coordinates[plus][:xin] + radius, stranded_vertex_coordinates[plus][:ymax] - radius + 0.5)
        Luxor.poly([a, b, c], :fill)
        a = Luxor.Point(stranded_vertex_coordinates[minus][:xin], stranded_vertex_coordinates[minus][:ymin])
        b = Luxor.Point(stranded_vertex_coordinates[minus][:xin], stranded_vertex_coordinates[minus][:ymin] + radius + 0.5)
        c = Luxor.Point(stranded_vertex_coordinates[minus][:xin] - radius, stranded_vertex_coordinates[minus][:ymin] + radius + 0.5)
        Luxor.poly([a, b, c], :fill)
        a = Luxor.Point(stranded_vertex_coordinates[minus][:xin], stranded_vertex_coordinates[minus][:ymax])
        b = Luxor.Point(stranded_vertex_coordinates[minus][:xin], stranded_vertex_coordinates[minus][:ymax] - radius + 0.5)
        c = Luxor.Point(stranded_vertex_coordinates[minus][:xin] - radius, stranded_vertex_coordinates[minus][:ymax] - radius + 0.5)
        Luxor.poly([a, b, c], :fill)
    end
    Luxor.finish()
end


function plot_canonical_kmer_graph(canonical_kmer_graph; filename=Random.randstring(Random.MersenneTwister(Int(round(time()))), 3) * ".svg")
    textsize = 12
    radius = textsize/2

    current_global_x_min = 0
    current_global_y_min = 0
    current_global_x_max = 0
    current_global_y_max = 0

    canonical_vertex_coordinates = [Dict{Symbol, Any}() for vertex in LightGraphs.vertices(canonical_kmer_graph)]

    for connected_component in LightGraphs.connected_components(canonical_kmer_graph)
        # reset xmin to left-align new contigs
        current_contig_x_min = current_global_x_min
        current_contig_x_max = current_global_x_min
        current_contig_y_min = current_global_y_max
        current_contig_y_max = current_global_y_max
        ## TODO could also achor by vertex with the smallest indegree
        anchor = minimum(connected_component)
        canonical_bfs_tree = LightGraphs.bfs_tree(canonical_kmer_graph, anchor)
        vertices_in_current_depth_of_field = [anchor]
        while !isempty(vertices_in_current_depth_of_field)
            vertices_in_next_depth_of_field = Vector{Int}()
            current_vertex_x_min = current_contig_x_max
            current_vertex_y_min = current_contig_y_min
            max_coverage = 0
            for vertex in vertices_in_current_depth_of_field
                coverage = length(canonical_kmer_graph.vprops[vertex][:coverage])
                if coverage > max_coverage
                    max_coverage = coverage
                end
            end
            for vertex in vertices_in_current_depth_of_field
                coverage = length(canonical_kmer_graph.vprops[vertex][:coverage])
                above_y_offset = 10 + coverage + radius
                below_y_offset = 10 + coverage + radius
                before_x_offset = 10 + max_coverage
                after_x_offset = 10 + max_coverage
                cvc = canonical_vertex_coordinates[vertex]
                cvc[:height] = textsize + coverage + textsize + 2
                cvc[:width] = (length(canonical_kmer_graph.gprops[:canonical_kmers][vertex]) + 1) * textsize
                cvc[:xmin] = current_vertex_x_min + before_x_offset
                cvc[:xmax] = cvc[:xmin] + cvc[:width]
                cvc[:ymin] = current_vertex_y_min + above_y_offset
                cvc[:ymax] = cvc[:ymin] + cvc[:height]
                cvc[:coverage_to_y_map] = [cvc[:ymin] + textsize + i for i in 1:coverage]
                cvc[:center] = Luxor.Point(Statistics.mean((cvc[:xmin], cvc[:xmax])), Statistics.mean((cvc[:ymin], cvc[:ymax])))

                if cvc[:xmax] + after_x_offset + radius > current_contig_x_max
                    current_contig_x_max = cvc[:xmax] + after_x_offset + radius
                end
                if cvc[:ymax] + coverage + radius > current_contig_y_max
                    current_contig_y_max = cvc[:ymax] + coverage + radius
                end
                current_vertex_y_min = cvc[:ymax] + coverage
                for neighbor in LightGraphs.outneighbors(canonical_bfs_tree, vertex)
                    push!(vertices_in_next_depth_of_field, neighbor)
                end
            end
            current_vertex_x_min = current_contig_x_max
            vertices_in_current_depth_of_field = sort!(unique(vertices_in_next_depth_of_field))
        end
        if current_contig_x_max > current_global_x_max
            current_global_x_max = current_contig_x_max
        end
        if current_contig_y_max > current_global_y_max
            current_global_y_max = current_contig_y_max
        end
    end

    Luxor.Drawing(current_global_x_max, current_global_y_max, filename)
    Luxor.background("white")
    Luxor.setline(1)
    ncolors = length(canonical_kmer_graph.gprops[:observation_ids])
    color_list = Colors.distinguishable_colors(ncolors+1, Colors.RGB(1,1,1))[2:end]
    for sequence_record_index in 1:length(canonical_kmer_graph.gprops[:observed_paths])
        sequence_path = canonical_kmer_graph.gprops[:observed_paths][sequence_record_index]
        sequence_color = color_list[canonical_kmer_graph.gprops[:observation_color_map][sequence_record_index]]
        i = 1
        ui = first(sequence_path[i])
        uiis = findall(observation -> first(observation) == sequence_record_index, canonical_kmer_graph.vprops[ui][:coverage])
        @assert length(uiis) == 1
        uii = first(uiis)
        ux = canonical_vertex_coordinates[ui][:xmax]
        uy = canonical_vertex_coordinates[ui][:coverage_to_y_map][uii]
        Luxor.setcolor(sequence_color)
        for i in 2:length(sequence_path)
            vi = first(sequence_path[i])
            canonical_kmer_graph.vprops[vi][:coverage]
            viis = findall(observation -> first(observation) == sequence_record_index, canonical_kmer_graph.vprops[vi][:coverage])
            @assert length(viis) == 1
            vii = first(viis)
            vy = canonical_vertex_coordinates[vi][:coverage_to_y_map][vii]
            vx = canonical_vertex_coordinates[vi][:xmin]
            u_radius = radius + uii - 1
            v_radius = radius + vii - 1
            a = Luxor.Point(ux, uy)
            b = Luxor.Point(Statistics.mean((ux, vx)), uy)
            c = Luxor.Point(Statistics.mean((ux, vx)), vy)
            d = Luxor.Point(vx, vy)
            Luxor.move(a); Luxor.curve(b, c, d); Luxor.strokepath()
            ui = vi
            uii = vii
            ux = canonical_vertex_coordinates[ui][:xmax]
            uy = canonical_vertex_coordinates[ui][:coverage_to_y_map][uii]
        end
    end
    for vertex in LightGraphs.vertices(canonical_kmer_graph)
        height = canonical_vertex_coordinates[vertex][:height]
        width = canonical_vertex_coordinates[vertex][:width]
        center = canonical_vertex_coordinates[vertex][:center]
        Luxor.setcolor("gray")
        Luxor.box(center, width, height, :fill)
        for (i, y) in enumerate(canonical_vertex_coordinates[vertex][:coverage_to_y_map])
            record_index, (sequence_index, orientation) = canonical_kmer_graph.vprops[vertex][:coverage][i]
            color_index = canonical_kmer_graph.gprops[:observation_color_map][record_index]
            Luxor.setcolor(color_list[color_index])
            Luxor.line(Luxor.Point(canonical_vertex_coordinates[vertex][:xmin], y), Luxor.Point(canonical_vertex_coordinates[vertex][:xmax], y), :stroke)
            Luxor.setcolor("white")
            Luxor.fontsize(1)
            Luxor.fontface("Menlo Bold")
            if orientation == true
                Luxor.text(string(sequence_index), Luxor.Point(canonical_vertex_coordinates[vertex][:xmin] + 1, y), valign=:middle, halign=:center)
            else
                Luxor.text(string(sequence_index), Luxor.Point(canonical_vertex_coordinates[vertex][:xmax] - 1, y), valign=:middle, halign=:center)
            end
            if sequence_index == 1
                if orientation == true
                    Luxor.circle(Luxor.Point(canonical_vertex_coordinates[vertex][:xmin] + 2, y), 0.5, :fill)
                else
                    Luxor.circle(Luxor.Point(canonical_vertex_coordinates[vertex][:xmax] - 2, y), 0.5, :fill)
                end
            end
        end
        plus_sequence = canonical_kmer_graph.gprops[:canonical_kmers][vertex]
        Luxor.setcolor("white")
        Luxor.fontsize(textsize)
        Luxor.fontface("Menlo")
        x = canonical_vertex_coordinates[vertex][:center].x
        y = Statistics.mean((canonical_vertex_coordinates[vertex][:ymin], canonical_vertex_coordinates[vertex][:center].y))
        Luxor.text(string(plus_sequence), Luxor.Point(x, y), valign=:middle, halign=:center)
        
        minus_sequence = BioSequences.reverse_complement(plus_sequence)
        y = Statistics.mean((canonical_vertex_coordinates[vertex][:ymax], canonical_vertex_coordinates[vertex][:center].y))
        Luxor.text(string(minus_sequence), Luxor.Point(x, y), valign=:middle, halign=:center)
    end
    Luxor.finish()
end

"""
document me
"""
function UPGMA(distance_matrix)
    N = size(distance_matrix, 1)
    tree = MetaGraphs.MetaDiGraph(N)
    labels = collect(1:N)
    c_vertex = N
    while N > 1
        c_vertex += 1
        add_vertex!(tree)
        closest_distance = minimum(d for d in distance_matrix if d > 0)
        b_index, a_index = Tuple(findfirst(distance_matrix .== closest_distance))
        a_vertex, b_vertex = labels[a_index], labels[b_index]
        c2a_distance = c2b_distance = distance_matrix[a_index, b_index]/2
        a_paths = filter(path -> !isempty(path), enumerate_paths(dijkstra_shortest_paths(tree, a_vertex)))
        a_weight = length(a_paths) > 1 ? length(a_paths) : 1
        b_paths = filter(path -> !isempty(path), enumerate_paths(dijkstra_shortest_paths(tree, b_vertex)))
        b_weight = length(b_paths) > 1 ? length(b_paths) : 1
        if a_weight > 1
            for (u,v) in zip(first(a_paths)[1:end-1], first(a_paths)[2:end])
                c2a_distance -= LightGraphs.get_prop(tree, LightGraphs.Edge(u, v), :length)
            end
        end
        if b_weight > 1
            for (u,v) in zip(first(b_paths)[1:end-1], first(b_paths)[2:end])
                c2b_distance -= LightGraphs.get_prop(tree, LightGraphs.Edge(u, v), :length)
            end
        end
        c2a_distance
        c2b_distance
        LightGraphs.add_edge!(tree, c_vertex, a_vertex)
        LightGraphs.set_prop!(tree, LightGraphs.Edge(c_vertex, a_vertex), :length, c2a_distance)
        LightGraphs.set_prop!(tree, LightGraphs.Edge(c_vertex, a_vertex), :weight, 1)
        LightGraphs.add_edge!(tree, c_vertex, b_vertex)
        LightGraphs.set_prop!(tree, LightGraphs.Edge(c_vertex, b_vertex), :length, c2b_distance)
        LightGraphs.set_prop!(tree, LightGraphs.Edge(c_vertex, b_vertex), :weight, 1)
        labels = deleteat!(labels, b_index)
        c_index = a_index
        labels[c_index] = c_vertex
        N-=1
        distance_matrix′ = zeros(N, N)
        for row in 1:N
            for column in 1:N
                row′ = row
                column′ = column
                if row′ >= b_index
                    row′ += 1
                end
                if column′ >= b_index
                    column′ += 1
                end
                if row != column
                    if row == c_index
                        a = distance_matrix[a_index, column′]
                        b = distance_matrix[b_index, column′]
                        value = (a * a_weight + b * b_weight) / (a_weight + b_weight)
                    elseif column == c_index
                        a = distance_matrix[row′, a_index]
                        b = distance_matrix[row′, b_index]
                        value = (a * a_weight + b * b_weight) / (a_weight + b_weight)
                    else
                        value = distance_matrix[row′, column′]
                    end
                    distance_matrix′[row, column] = value
                end
            end
        end
        distance_matrix = distance_matrix′
        # display(distance_matrix)
    end
    return tree
end

"""
    build_stranded_kmer_graph(canonical_kmers, observations)

Create a weighted, strand-specific kmer (de bruijn) graph from a set of kmers
and a series of sequence observations in FASTA format.
"""
function build_stranded_kmer_graph(canonical_kmers, observations)
    stranded_kmers = sort!(vcat(canonical_kmers, [BioSequences.reverse_complement(kmer) for kmer in canonical_kmers]))
    stranded_kmer_to_reverse_complement_map = [
        findfirst(stranded_kmer -> BioSequences.reverse_complement(stranded_kmer) == kmer, stranded_kmers) for kmer in stranded_kmers
    ]
    stranded_kmer_graph = MetaGraphs.MetaDiGraph(length(stranded_kmers))
    stranded_kmer_graph.gprops[:stranded_kmers] = stranded_kmers
    stranded_kmer_graph.gprops[:reverse_complement_map] = stranded_kmer_to_reverse_complement_map
    stranded_kmer_graph.gprops[:k] = length(first(stranded_kmers))
    stranded_kmer_graph.gprops[:K] = length(stranded_kmers)
    stranded_kmer_graph.gprops[:observation_color_map] = Vector{Int}()
    stranded_kmer_graph.gprops[:observation_ids] = Vector{String}()
    stranded_kmer_graph.gprops[:observed_paths] = Vector{Vector{Pair{Int, Bool}}}()
    for vertex in 1:LightGraphs.nv(stranded_kmer_graph)
        stranded_kmer_graph.vprops[vertex] = Dict(:coverage => Vector{Pair{Int, Pair{Int, Bool}}}())
    end
    for (observation_index, observation) in enumerate(observations)
        observation_id = BioSequences.FASTA.identifier(observation)
        observed_sequence = BioSequences.FASTA.sequence(observation)
        if length(observed_sequence) < stranded_kmer_graph.gprops[:k]
            @error "skipping sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        else
            observed_path = sequence_to_stranded_path(stranded_kmer_graph.gprops[:stranded_kmers], observed_sequence)
            i = 1
            ui, ui_orientation = observed_path[i]
            ui_coverage = (observation_index => (i => ui_orientation ))
            push!(stranded_kmer_graph.vprops[ui][:coverage], ui_coverage)
            for i in 2:length(observed_path)
                vi, vi_orientation = observed_path[i]
                vi_coverage = (observation_index => (i => vi_orientation))
                push!(stranded_kmer_graph.vprops[vi][:coverage], vi_coverage)
                edge_coverage = ui_coverage => vi_coverage
                if LightGraphs.has_edge(stranded_kmer_graph, ui, vi)
                    push!(stranded_kmer_graph.eprops[LightGraphs.Edge(ui, vi)][:coverage], edge_coverage)
                else
                    LightGraphs.add_edge!(stranded_kmer_graph, ui, vi, Dict(:coverage => [edge_coverage]))
                end
                # not sure this is necessary
#                 ui′ = stranded_kmer_graph.gprops[:reverse_complement_map][ui]
#                 vi′ = stranded_kmer_graph.gprops[:reverse_complement_map][vi]
#                 if !LightGraphs.has_edge(stranded_kmer_graph, vi′, ui′)
#                     LightGraphs.add_edge!(stranded_kmer_graph, vi′, ui′, Dict(:coverage => Vector{typeof(edge_coverage)}()))
#                 end
                ui, ui_orientation = vi, vi_orientation
                ui_coverage = vi_coverage
            end
            push!(stranded_kmer_graph.gprops[:observed_paths], observed_path)
            push!(stranded_kmer_graph.gprops[:observation_ids], observation_id)
            push!(stranded_kmer_graph.gprops[:observation_color_map], observation_index)
        end
    end
    return stranded_kmer_graph
end
    
"""
    build_canonical_kmer_graph(canonical_kmers, observations)

Create a weighted kmer (de bruijn) graph from a set of kmers
and a series of sequence observations in FASTA format.
"""
function build_canonical_kmer_graph(canonical_kmers, observations)
    canonical_kmer_graph = MetaGraphs.MetaGraph(length(canonical_kmers))
    canonical_kmer_graph.gprops[:canonical_kmers] = canonical_kmers
    canonical_kmer_graph.gprops[:k] = length(first(canonical_kmers))
    canonical_kmer_graph.gprops[:K] = length(canonical_kmers)
    canonical_kmer_graph.gprops[:observation_color_map] = Vector{Int}()
    canonical_kmer_graph.gprops[:observation_ids] = Vector{String}()
    canonical_kmer_graph.gprops[:observed_paths] = Vector{Vector{Pair{Int, Bool}}}()
    for vertex in 1:LightGraphs.nv(canonical_kmer_graph)
        canonical_kmer_graph.vprops[vertex] = Dict(:coverage => Vector{Pair{Int, Pair{Int, Bool}}}())
    end
    for (observation_index, observation) in enumerate(observations)
        observation_id = BioSequences.FASTA.identifier(observation)
        observed_sequence = BioSequences.FASTA.sequence(observation)
        if length(observed_sequence) < canonical_kmer_graph.gprops[:k]
            @error "sequence shorter than k with id $observation_id & length $(length(observed_sequence))"
        else
            observed_path = sequence_to_canonical_path(canonical_kmer_graph.gprops[:canonical_kmers], observed_sequence)
            i = 1
            ui, ui_orientation = observed_path[i]
            ui_coverage = (observation_index => (i => ui_orientation))
            push!(canonical_kmer_graph.vprops[ui][:coverage], ui_coverage)
            for i in 2:length(observed_path)
                vi, vi_orientation = observed_path[i]
                vi_coverage = (observation_index => (i => vi_orientation))
                push!(canonical_kmer_graph.vprops[vi][:coverage], vi_coverage)
                edge_coverage = ui_coverage => vi_coverage
                if LightGraphs.has_edge(canonical_kmer_graph, ui, vi)
                    push!(canonical_kmer_graph.eprops[LightGraphs.Edge(ui, vi)][:coverage], edge_coverage)
                else
                    LightGraphs.add_edge!(canonical_kmer_graph, ui, vi, Dict(:coverage => [edge_coverage]))
                end
                ui, ui_orientation = vi, vi_orientation
                ui_coverage = vi_coverage
            end
        end
        push!(canonical_kmer_graph.gprops[:observed_paths], observed_path)
        push!(canonical_kmer_graph.gprops[:observation_ids], observation_id)
        push!(canonical_kmer_graph.gprops[:observation_color_map], observation_index)
    end
    return canonical_kmer_graph
end

function determine_file_type(file)
    if endswith(file, ".gz")
        file = file[1:end-3]
    end
    if endswith(file, ".fa") || endswith(file, ".fasta") || endswith(file, ".fna")
        return BioSequences.FASTA
    elseif endswith(file, ".fq") || endswith(file, ".fastq")
        return FASTQ
    end
end

function open_file(file)
    if endswith(file, ".gz")
        return CodecZlib.GzipDecompressorStream(open(file))
    else
        return open(file)
    end
end

function plot_histogram(parsed_args)
    lines = strip.(readlines(open_file(parsed_args["histogram"])))
    number_of_kmers = Vector{Int}(undef, length(lines))
    depth_of_coverage = Vector{Int}(undef, length(lines))
    for (i, line) in enumerate(lines)
        number_of_kmers[i], depth_of_coverage[i] = parse.(Int, split(line))
    end

    X = log2.(depth_of_coverage)
    Y = log2.(number_of_kmers)
    p = scatter(X,
                Y,
                xlims = (-0.5, maximum(X)+0.5),
                ylims = (-0.5, maximum(Y)+0.5),
                title=basename(parsed_args["histogram"]),
                titlefontsize=11,
                xlabel="log2(depth of coverage)",
                ylabel="log2(# kmers)",
                legend=false);
    if length(X) > 1
        X2 = LinRange(minimum(X), maximum(X), 100)
        Y2 = predict(lm(GLM.@formula(Y ~ X), DataFrame(X = X, Y = Y)), DataFrame(X = X2))
        plot!(p, X2, Y2)
        annotate!(p, 0, 1, Plots.text("r = $(round(cor(X, Y), digits=2))", 10, :left))
    end
    if endswith(parsed_args["histogram"], ".gz")
        savefig(p, join(split(parsed_args["histogram"], '.')[1:end-1], '.') * ".svg");
    else
        savefig(p, parsed_args["histogram"] * ".svg");
    end
end

function plot_rank_frequency(parsed_args)
    kmer_counts = sort(parse.(Int, first.(split.(strip.(readlines(open_file(parsed_args["kmer-counts"])))))), rev=true)
    rank = 1:length(kmer_counts)
    X = log2.(rank)
    Y = log2.(kmer_counts)
    X2 = LinRange(minimum(X), maximum(X), 100)
    Y2 = predict(lm(GLM.@formula(Y ~ X), DataFrame(X = X, Y = Y)), DataFrame(X = X2))

    if endswith(parsed_args["kmer-counts"], ".gz")
        filename = join(split(parsed_args["kmer-counts"], '.')[1:end-1], '.')
    else
        filename = parsed_args["kmer-counts"]
    end

    savefig(
    annotate!(
        plot!(
            scatter(X,
                    Y,
                    xlims = (-0.5, maximum(X)+0.5),
                    ylims = (-0.5, maximum(Y)+0.5),
                    title=basename(parsed_args["kmer-counts"]),
                    titlefontsize=11,
                    xlabel="log2(rank)",
                    ylabel="log2(count)",
                    legend=false),
            X2, Y2),
            0, 1, Plots.text("r = $(round(cor(X, Y), digits=2))", 10, :left)),
            filename * ".png");
end

end
