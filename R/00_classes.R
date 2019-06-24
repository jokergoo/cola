 

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
		.env = "environment"
))
