 

# == title
# The ConsensusPartitionList class
#
# == alias
# ConsensusPartitionList
#
# == details
# The object contains results from all combinations of top methods
# and partition methods.
#
# == Methods
# The `ConsensusPartitionList-class` provides following methods:
#
# -`run_all_consensus_partition_methods`: constructor method.
# -`get_single_run,ConsensusPartitionList-method`: get a single `ConsensusPartition-class` object from the list.
# -`top_rows_overlap,ConsensusPartitionList-method`: plot the overlaps of top rows under different top methods.
# -`top_rows_heatmap,ConsensusPartitionList-method`: plot the heatmap of top rows under different top methods.
# -`get_class,ConsensusPartitionList-method`: get a consensus class IDs merging from all methods.
# -`get_stat,ConsensusPartitionList-method`: get statistics for a specified k.
# -`collect_plots,ConsensusPartitionList-method`: collect plots from all combinations of top methods and partition methods with choosing a plotting function.
# -`collect_classes,ConsensusPartitionList-method`: make a plot which contains predicted classes from all combination of top methods and partition methods.
# -`test_to_known_factors,ConsensusPartitionList-method`: test correlation between predicted subgrouping and known annotations, if available.
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
		top_method = "character",
		partition_method = "character",
		consensus_class = "ANY",
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
# -`consensus_partition`: constructor method, run consensus partition with a specified top method and a partition method.
# -`select_partition_number,ConsensusPartition-method`: make a list of plots used to select optimized number of partitions.
# -`consensus_heatmap,ConsensusPartition-method`: Heatmap of the consensus matrix.
# -`membership_heatmap,ConsensusPartition-method`: Heatmap of the membership in each random sampling.
# -`get_signatures,ConsensusPartition-method`: get the signature rows and make heatmaps.
# -`dimension_reduction,ConsensusPartition-method`: dimension reduction plots.
# -`collect_plots,ConsensusPartition-method`: Heatmaps for consensus matrix and membership matrix with different number of partitions.
# -`collect_classes,ConsensusPartition-method`: Heatmap of classes with different number of partitions.
# -`get_param,ConsensusPartition-method`: get parameters for the consensus clustering.
# -`get_consensus,ConsensusPartition-method`: get the consensus matrix.
# -`get_membership,ConsensusPartition-method`: get the membership in random samplings.
# -`get_stat,ConsensusPartition-method`: get statistics for the consensus clustering.
# -`get_class,ConsensusPartition-method`: get the consensus class IDs and other columns.
# -`get_best_k,ConsensusPartition-method`: guess the best number of partitions.
# -`signature_density,ConsensusPartition-method`: plot the density distribution of signatures in different groups.
# -`test_to_known_factors,ConsensusPartition-method`: test correlation between predicted subgrouping and known annotations, if available.
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
		top_method = "character",
		top_n = "numeric",
		known_anno = "ANY", 
		known_col = "ANY", 
		scale_rows = "logical",
		column_index = "numeric",
		.env = "environment"
))
