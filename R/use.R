
use_scRNASeq = function() {
	cat(
"Following configurations are suggested for scRNASeq datasets:
- use `ATC_kNN()` for extracting top genes for clustering.
    first call `register_ATC_kNN()` to register ATC_kNN method,
    and later method 'ATC_kNN' can be used for clustering.
- hierarchical_partition():
    min_n_signature = 100
- get_signatures():
    group_diff = 0.5 or higher
")
}

use_methylation = function() {
	cat(
"Following configurations are suggested for methylation datasets:
- consensus_partition(), hierarchical_partition(), consensus_partition_by_down_sampling(),
  run_all_consensus_partition_methods():
    top_value_method = 'MAD' or 'SD'
    scale_rows = FALSE
- get_signatures():
    group_diff = 0.2 or higher
    fdr_cutoff = 0.00001
")
}
