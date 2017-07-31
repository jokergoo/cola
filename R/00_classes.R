 

# == title
# The ConsensusPartitionList class
#
# == alias
# ConsensusPartition
#
# == details
# The object contains results from all combinations of top methods
# and partition methods.
#
# == Methods
# The `ConsensusPartitionList-class` provides following methods:
#
# -`run_all_consensus_partition_methods`: constructor method;
# -`get_single_run,ConsensusPartitionList-method`: get a single `ConsensusPartition-class` object from the list;
# -`top_rows_overlap,ConsensusPartitionList-method`: plot the overlaps of top rows under different top methods;
# -`top_rows_heatmap,ConsensusPartitionList-method`: plot the heatmap of top rows under different top methods;
# -`get_class,ConsensusPartitionList-method`: get a consensus class IDs
# -`collect_plots,ConsensusPartitionList-method`: collect plots from all combination of top methods and partition methods;
# -`collect_classes,ConsensusPartitionList-method`: make a plot which contains predicted classes from all combination of top methods and partition methods.
#
# == seealso
# The `ConsensusPartition-class`.
#
ConsensusPartitionList = setClass("ConsensusPartitionList",
	slots = list(
		list = "list",
		top_method = "character",
		partition_method = "character",
		.env = "environment"
))


# == title
# The ConsensusPartition class
#
# == alias
# ConsensusPartitionList
#
# == methods
# The `ConsensusPartition-class` has following methods:
#
# -`consensus_partition`: constructor method;
# -`plot_ecdf,ConsensusPartition-method`: plot the ecdf of the consensus matrix;
# -`select_partition_number,ConsensusPartition-method`: a list of plots used to pick optimized number of partitions;
# -`consensus_heatmap,ConsensusPartition-method`: Heatmap of the consensus matrix;
# -`membership_heatmap,ConsensusPartition-method`: Heatmap of the membership in each randomization;
# -`get_signatures,ConsensusPartition-method`: Heatmap of signature rows;
# -`dimension_reduction,ConsensusPartition-method`: dimension reduction plots;
# -`collect_plots,ConsensusPartition-method`: Heatmaps for consensus matrix and membership matrix with different number of partitions;
# -`collect_classes,ConsensusPartition-method`: Heatmap of classes with different number of partitions;
# -`get_param,ConsensusPartition-method`: get parameters for the consensus clustering;
# -`get_consensus,ConsensusPartition-method`: get the consensus matrix;
# -`get_membership,ConsensusPartition-method`: get the membership in randomizations;
# -`get_stat,ConsensusPartition-method`: get statistics for the consensus clustering;
# -`get_class,ConsensusPartition-method`: get the consensus class IDs and other columns.
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
