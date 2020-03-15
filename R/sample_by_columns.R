######################################
## only for column sampling
######################################

# for a single partition, all classes are assigned as 0 if the corresponding label is NA
# clue:::.cl_membership_from_class_ids
cl_membership2 = function(x, k) {
	m = matrix(0, nrow = n_of_objects(x), ncol = k)
	
	j = cl_class_ids(x)
	i = seq_along(j)
	
	nr = nrow(m)
    nc = ncol(m)

    ind = (j - 1) * nr + i

    m[ind[!is.na(j)]] = 1

	return(m)
}

.random_stochastic_matrix = getFromNamespace(".random_stochastic_matrix", ns = "clue")
.weighted_sum_of_matrices = getFromNamespace(".weighted_sum_of_matrices", ns = "clue")
.project_to_leading_columns = getFromNamespace(".project_to_leading_columns", ns = "clue")
.stochastify = getFromNamespace(".stochastify", ns = "clue")
.cl_membership_from_memberships = getFromNamespace(".cl_membership_from_memberships", ns = "clue")

cl_consensus2 = function(clusterings, k_max) {
    n <- n_of_objects(clusterings)

    k <- k_max

    maxiter <- 100

    reltol <- sqrt(.Machine$double.eps)

    nruns <- 1L

    M = .random_stochastic_matrix(n, k)
    w = rep(1, length(clusterings))
    k_all <- max(k, k_max)
    value <- function(M, memberships, w) {
        sum(w * sapply(memberships, function(u) sum((u - M)^2)))
    }
    match_memberships <- function(M, N) {
        M[, solve_LSAP(crossprod(N, M), maximum = TRUE), drop = FALSE]
    }
    fit_M <- function(memberships, w, k) {
        M <- .weighted_sum_of_matrices(memberships, w, nrow(M))
        if (k < ncol(M))
            M <- .project_to_leading_columns(M, k)
        M
    }
    memberships <- lapply(clusterings, cl_membership2, k_all)
    V_opt <- Inf
    M_opt <- NULL

    if (k < k_all) M <- cbind(M, matrix(0, nrow(M), k_all - k))
    memberships <- lapply(memberships, match_memberships, M)
    old_value <- value(M, memberships, w)
    iter <- 1L
    while (iter <= maxiter) {
        M <- fit_M(memberships, w, k)
        memberships <- lapply(memberships, match_memberships, M)
        new_value <- value(M, memberships, w)
        if (abs(old_value - new_value) < reltol * (abs(old_value) + reltol)) break
        old_value <- new_value
        iter <- iter + 1L
    }
    if (new_value < V_opt) {
        converged <- (iter <= maxiter)
        V_opt <- new_value
        M_opt <- M
    }
        
    M <- .stochastify(M_opt)
    rownames(M) <- rownames(memberships[[1L]])

    if(any(is.na(M))) {
        stop_wrap("There are at least one sample that have not been selected in all random resamplings. Either increase `p_sampling` or `partition_repeats`, or select another random seed.")
    }

    meta <- list(objval = value(M, memberships, w), converged = converged)
    M <- .cl_membership_from_memberships(M[, seq_len(k), drop = FALSE], k, meta)
    as.cl_partition(M)
}
