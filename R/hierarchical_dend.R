
has_hierarchy = function(object) {
	nrow(object@hierarchy) > 0
}

subgroup_dend = function(object, merge_node = merge_node_param()) {

	check_pkg("data.tree", bioc = FALSE)

	hierarchy = get_hierarchy_table(object, merge_node)
	
	lt = list()
	lt[["0"]] = data.tree::Node$new("0")
	cn = colnames(object@list[["0"]]@.env$data)
	max_depth = max(nchar(hierarchy))
	lt[["0"]]$node_height = max_depth - 1
	for(i in seq_len(nrow(hierarchy))) {
		lt[[ hierarchy[i, 2] ]] = lt[[ hierarchy[i, 1] ]]$AddChildNode({
			node = data.tree::Node$new(hierarchy[i, 2])
			node$max_height = max_depth - nchar(hierarchy[i, 2])
			node
		})
		l = hierarchy[, 1] == hierarchy[i, 2]
	}
	dend = as.dendrogram(lt[["0"]], heightAttribute = "node_height", edgetext = TRUE)

	dend = dendextend::`order.dendrogram<-`(dend, value = 1:nobs(dend))

	dend = edit_node(dend, function(d, index) {
		if(is.leaf(d)) {
			attr(d, "node_id") = attr(d, "label")
		} else {
			attr(d, "node_id") = attr(d, "edgetext")
		}
		d
	})

	# make sure all nodes have a node_id attribute
	# depth first
	.get_node_id = function(d) {
		node_id = attr(d, "node_id")
		if(is.null(node_id)) {
			n_children = length(d)
			for(i in seq_len(n_children)) {
				child_node_id = .get_node_id(d[[i]])
				if(!is.null(child_node_id)) break
			}
			node_id = gsub("\\d$", "", child_node_id)
		}
		return(node_id)
	}

	dend = edit_node(dend, function(d, index) {
		node_id = attr(d, "node_id")
		if(is.null(node_id)) {
			node_id = .get_node_id(d)
			attr(d, "node_id") = node_id
		}
		d
	})

	oe = try(dend_tmp <- as.dendrogram(as.hclust(dend)), silent = TRUE)

	if(!inherits(oe, "try-error")) {
		dend = edit_node(dend, function(d, ind) {
			if(length(ind) == 0) {
				attr(d, "midpoint") = attr(dend_tmp, "midpoint")
			} else {
				attr(d, "midpoint") = attr(dend_tmp[[ind]], "midpoint")
			}
			d
		})
	}

	max_height = max(nchar(all_nodes(object)))
	edit_node(dend, function(d) {
		attr(d, "height") = (max_height - nchar(attr(d, "node_id")))
		d
	})

}

get_hierarchy_dend = function(object, merge_node = merge_node_param()) {

	dend = subgroup_dend(object, merge_node)
	dend
}

random_dend = function(n, label) {
	x = rnorm(n)
	if(n == 1) {
		dend = structure(1, members = 1, height = 0, label = label[ind], leaf = TRUE, class = "dendrogram")
	} else {
		dend = as.dendrogram(hclust(dist(1:n)))
	}
	# set height to zero

	dendrapply(dend, function(x) {attr(x, "height") = 0; x})
}

zero_height_dend = function(n) {
	check_pkg("data.tree", bioc = FALSE)

	lt = data.tree::Node$new("foo")
	lt$node_height = 0
	for(i in 1:n) {
		lt$AddChildNode({
			node = data.tree::Node$new(paste0("foo", i))
			node$node_height = 0
			node
		})
	}
	dend = as.dendrogram(lt, heightAttribute = "node_height")
	
}

calc_dend = function(object, merge_node = merge_node_param(), mat = NULL) {

	pd = get_hierarchy_dend(object, merge_node)
	classes = get_classes(object, merge_node)
	if(is.null(names(classes))) names(classes) = seq_along(classes)

	if(is.null(mat)) {
		if(!is.null(object[["0"]]@anno)) {
			if(inherits(object[["0"]], "ConsensusPartition")) {
				mat = apply(object[["0"]]@anno, 2, rank)
			} else {
				if(!is.null(object[["0"]]@full_anno)) {
					mat = apply(object[["0"]]@full_anno, 2, rank)
				} else if(!is.null(object[["0"]]@anno)) {
					mat = apply(object[["0"]]@anno, 2, rank)
				} 
			}
			mat = t(mat)
		}
	}

	if(is.null(mat)) {
		cd_list = lapply(tapply(names(classes), classes, function(x) x), function(x) {
			d = random_dend(length(x), colnames(object))
			d = dendextend::`labels<-`(d, value = x)
			d
		})
	} else {
		colnames(mat) = names(classes)
		cd_list = tapply(seq_along(classes), classes, function(ind) {
			if(length(ind) > 1) {
				d = as.dendrogram(hclust(dist(t(mat[, ind, drop = FALSE]))))
			} else {
				d = structure(1, members = 1, height = 0, label = colnames(mat)[ind], leaf = TRUE, class = "dendrogram")
			}
			d = edit_node(d, function(x) {attr(x, "height") = 0; x})
			d
		})
	}
	cd_list = cd_list[labels(pd)]

	dend = merge_dendrogram(pd, cd_list)
	dend = adjust_dend_by_x(dend)
	dend = dendextend::`order.dendrogram<-`(dend, value = structure(1:length(classes), names = names(classes))[labels(dend)])

	dend
}

# == title
# Max depth of the hierarchy
#
# == param
# -object A `HierarchicalPartition-class` object.
#
# == value
# A numeric value.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# max_depth(golub_cola_rh)
setMethod(f = "max_depth",
	signature = "HierarchicalPartition",
	definition = function(object) {

	if(has_hierarchy(object)) {
		max(nchar(object@hierarchy[, 2]))
	} else {
		1
	}
})


# == title
# Parameters to merge subgroup dendrogram.
#
# == param
# -depth Depth of the dendrogram.
# -min_n_signatures Minimal number of signatures for the partitioning on each node.
# -min_p_signatures Minimal fraction of sigatures compared to the total number of rows on each node.
# -node_height The height of the sub-dendrogram to cut
#
merge_node_param = function(depth = Inf, min_n_signatures = -Inf, 
	min_p_signatures = -Inf, node_height = -Inf) {
	
	list(depth = depth, min_n_signatures = min_n_signatures, 
		min_p_signatures = min_p_signatures, node_height = node_height)
}

# == title
# All nodes in the hierarchy
#
# == param
# -object A `HierarchicalPartition-class` object.
# -merge_node Parameters to merge sub-dendrograms, see `merge_node_param`.
#
# == value
# A vector of node ID.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# all_nodes(golub_cola_rh)
setMethod(f = "all_nodes",
	signature = "HierarchicalPartition",
	definition = function(object, merge_node = merge_node_param()) {

	if(has_hierarchy(object)) {
		hierarchy = get_hierarchy_table(object, merge_node)
		return(unique(as.vector(t(hierarchy))))
	} else {
		return(character(0))
	}
})

get_hierarchy_table = function(object, merge_node = merge_node_param()) {

	hierarchy = object@hierarchy

	n_signatures = object@node_level$n_signatures
	p_signatures = n_signatures/nrow(object)
	node_height = object@node_level$node_height

	if(is.null(node_height)) {
		hierarchy = hierarchy[ n_signatures[hierarchy[, 1]] >= merge_node$min_n_signatures &
		                       p_signatures[hierarchy[, 1]] >= merge_node$min_p_signatures, , drop = FALSE ]
	} else {
		hierarchy = hierarchy[ n_signatures[hierarchy[, 1]] >= merge_node$min_n_signatures &
		                       p_signatures[hierarchy[, 1]] >= merge_node$min_p_signatures &
		                       node_height[hierarchy[, 1]] >= merge_node$node_height, , drop = FALSE ]
	}
	
	hierarchy = hierarchy[nchar(hierarchy[, 2]) <= merge_node$depth, , drop = FALSE]
	hierarchy = hierarchy[ !(hierarchy[, 1] != "0" & !hierarchy[, 1] %in% hierarchy[, 2]), , drop = FALSE]
	
	hierarchy
}

# == title
# All leaves in the hierarchy
#
# == param
# -object A `HierarchicalPartition-class` object.
# -merge_node Parameters to merge sub-dendrograms, see `merge_node_param`.
#
# == value
# A vector of node ID.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# data(golub_cola_rh)
# all_leaves(golub_cola_rh)
setMethod(f = "all_leaves",
	signature = "HierarchicalPartition",
	definition = function(object, merge_node = merge_node_param()) {

	if(has_hierarchy(object)) {
		hierarchy = get_hierarchy_table(object, merge_node)
		tb = table(hierarchy)
		names(tb[tb <= 1])
	} else {
		"0"
	}
})

# == title
# Test whether a node is a leaf node
#
# == param
# -object A `HierarchicalPartition-class` object.
# -node A vector of node IDs.
# -merge_node Parameters to merge sub-dendrograms, see `merge_node_param`.
#
# == example
# data(golub_cola_rh)
# is_leaf_node(golub_cola_rh, all_leaves(golub_cola_rh))
setMethod(f = "is_leaf_node",
	signature = "HierarchicalPartition",
	definition = function(object, node, merge_node = merge_node_param()) {

	all_nodes = all_nodes(object, merge_node)
	all_leaves = all_leaves(object, merge_node)

	l = node %in% all_leaves
	l[!node %in% all_nodes] = NA
	l
})

get_children = function(object, node = "0") {
	hierarchy = unique(object@hierarchy)
	hierarchy[hierarchy[, 1] == node, 2]
}
