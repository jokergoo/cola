
set.seed(123)

partition_list = list(
	c(1, 1, 1, NA, 1, 2, 2, 2, 2, 2),
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1)
)
library(clue)
partition_list = lapply(partition_list, as.cl_hard_partition)

p = cl_ensemble(list = partition_list)
m1 = cl_membership(cola:::cl_consensus2(p, 2))
class(m1) = "matrix"
m1 = m1[, 1:2]


partition_list = list(
	c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1)
)
partition_list = lapply(partition_list, as.cl_hard_partition)

p = cl_ensemble(list = partition_list)
m2 = cl_membership(cl_consensus(p))
class(m2) = "matrix"
m2 = m2[, 1:2]

if(m1[1, 1] != m2[1, 1]) {
	m2 = m2[, 2:1]
}

test_that("test cl_consensus2", {
	expect_equal(m1, m2)
})

#############################################################

partition_list = list(
	c(2, 2, 2, NA,2, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, NA, 2, 1, 1, 1, 1, 1)
)
library(clue)
partition_list = lapply(partition_list, as.cl_hard_partition)

p = cl_ensemble(list = partition_list)
m1 = cl_membership(cola:::cl_consensus2(p, 2))
class(m1) = "matrix"
m1 = m1[, 1:2]


partition_list = list(
	c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 2, 1, 1, 1, 1, 1)
)
partition_list = lapply(partition_list, as.cl_hard_partition)

p = cl_ensemble(list = partition_list)
m2 = cl_membership(cl_consensus(p))
class(m2) = "matrix"
m2 = m2[, 1:2]

if(m1[1, 1] != m2[1, 1]) {
	m2 = m2[, 2:1]
}

test_that("test cl_consensus2", {
	expect_equal(m1, m2)
})


#############################################################

partition_list = list(
	c(2, 2, 2, NA,1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 1, 1, 1, 1, 1, 1, 1)
)
library(clue)
partition_list = lapply(partition_list, as.cl_hard_partition)

p = cl_ensemble(list = partition_list)
m1 = cl_membership(cola:::cl_consensus2(p, 2))
class(m1) = "matrix"
m1 = m1[, 1:2]


partition_list = list(
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1),
	c(2, 2, 2, 1, 1, 1, 1, 1, 1, 1)
)
partition_list = lapply(partition_list, as.cl_hard_partition)

p = cl_ensemble(list = partition_list)
m2 = cl_membership(cl_consensus(p))
class(m2) = "matrix"
m2 = m2[, 1:2]

if(m1[1, 1] != m2[1, 1]) {
	m2 = m2[, 2:1]
}

test_that("test cl_consensus2", {
	expect_equal(m1, m2)
})
