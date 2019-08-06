 

# == title
# The ConsensusPartitionList class
#
# == alias
# ConsensusPartitionList
#
# == details
# The object contains results from all combinations of top-value methods
# and partition methods.
#
# == Methods
# The `ConsensusPartitionList-class` provides following methods:
#
# -`run_all_consensus_partition_methods`: constructor method.
# -`top_rows_overlap,ConsensusPartitionList-method`: plot the overlaps of top rows under different top-value methods.
# -`top_rows_heatmap,ConsensusPartitionList-method`: plot the heatmap of top rows under different top-value methods.
# -`get_classes,ConsensusPartitionList-method`: get consensus class IDs merging from all methods.
# -`get_matrix,ConsensusPartition-method`: get the original matrix.
# -`get_stats,ConsensusPartitionList-method`: get metrics for a specified k.
# -`get_membership,ConsensusPartitionList-method`: get consensus membership matrix summarized from all methods.
# -`suggest_best_k,ConsensusPartitionList-method`: guess the best number of partitions for all methods.
# -`collect_plots,ConsensusPartitionList-method`: collect plots from all combinations of top-value methods and partition methods with choosing a plotting function.
# -`collect_classes,ConsensusPartitionList-method`: make a plot which contains predicted classes from all combinations of top-value methods and partition methods.
# -`test_to_known_factors,ConsensusPartitionList-method`: test correlation between predicted classes and known annotations, if provided.
# -`cola_report,ConsensusPartitionList-method`: generate a HTML report for the whole analysis.
#
# == seealso
# The `ConsensusPartition-class`.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
ConsensusPartitionList = setClass("ConsensusPartitionList",
	slots = list(
		list = "list",
		top_value_method = "character",
		partition_method = "character",
		consensus_class = "ANY",
		comb = "ANY",
		call = "ANY",
		.env = "environment"
))


# == title
# The ConsensusPartition class
#
# == alias
# ConsensusPartition
#
# == methods
# The `ConsensusPartition-class` has following methods:
#
# -`consensus_partition`: constructor method, run consensus partition with a specified top-value method and a partition method.
# -`select_partition_number,ConsensusPartition-method`: make a list of plots to select optimized number of partitions.
# -`consensus_heatmap,ConsensusPartition-method`: make heatmap of the consensus matrix.
# -`membership_heatmap,ConsensusPartition-method`: make heatmap of the membership in every random sampling.
# -`get_signatures,ConsensusPartition-method`: get the signature rows and make heatmap.
# -`dimension_reduction,ConsensusPartition-method`: make dimension reduction plots.
# -`collect_plots,ConsensusPartition-method`: make heatmaps for consensus matrix and membership matrix with different number of partitions.
# -`collect_classes,ConsensusPartition-method`: make heatmap of classes with different numbers of partitions.
# -`get_param,ConsensusPartition-method`: get parameters for the consensus clustering.
# -`get_matrix,ConsensusPartition-method`: get the original matrix.
# -`get_consensus,ConsensusPartition-method`: get the consensus matrix.
# -`get_membership,ConsensusPartition-method`: get the membership in random samplings.
# -`get_stats,ConsensusPartition-method`: get metrics for the consensus clustering.
# -`get_classes,ConsensusPartition-method`: get the consensus class IDs and other columns.
# -`suggest_best_k,ConsensusPartition-method`: guess the best number of partitions.
# -`test_to_known_factors,ConsensusPartition-method`: test correlation between predicted classes and known factors, if available.
# -`GO_enrichment,ConsensusPartition-method`: perform GO enrichment analysis on significant genes if rows in the matrix can be corresponded to genes.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
ConsensusPartition = setClass("ConsensusPartition",
	slots = list(
		object_list = "list", 
		k = "numeric", 
		n_partition = "numeric",  
		partition_method = "character", 
		top_value_method = "character",
		top_n = "numeric",
		anno = "ANY", 
		anno_col = "ANY", 
		scale_rows = "logical",
		column_index = "numeric",
		running_time = "ANY",
		cache = "list",
		hash = "character",
		.env = "environment"
	),
	prototype = list(hash = "")
)



# == title
# The HierarchicalPartition class
#
# == alias
# HierarchicalPartition
#
# == methods
# The `HierarchicalPartition-class` has following methods:
#
# -`hierarchical_partition`: constructor method.
# -`collect_classes,HierarchicalPartition-method`: plot the hierarchy of subgroups predicted.
# -`get_classes,HierarchicalPartition-method`: get the class IDs of subgroups.
# -`suggest_best_k,HierarchicalPartition-method`: guess the best number of partitions for each node.
# -`get_matrix,HierarchicalPartition-method`: get the original matrix.
# -`get_signatures,HierarchicalPartition-method`: get the signatures for each subgroup.
# -`dimension_reduction,HierarchicalPartition-method`: make dimension reduction plots.
# -`test_to_known_factors,HierarchicalPartition-method`: test correlation between predicted subgrouping and known annotations, if available.
# -`cola_report,HierarchicalPartition-method`: generate a HTML report for the whole analysis.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
HierarchicalPartition = setClass("HierarchicalPartition",
    slots = list(
        list = "list",
        best_k = "numeric",
        hierarchy = "matrix",
        subgroup = "character",
        subgroup_col = "character",
        call = "ANY",
        .env = "environment"
    )
)