
#include <Rcpp.h>
using namespace Rcpp;

float logical_sum(LogicalVector v1, LogicalVector v2) {
	int s = 0;
	for(int i = 0; i < v1.size(); i ++) {
		s = s + (v1[i] & v2[i]);
	}
	return s;
}

// [[Rcpp::export]]
NumericMatrix gene_cooccurrence_in_geneset(LogicalMatrix match_mat) {

	int n_gene = match_mat.nrow();

	NumericMatrix m(n_gene, n_gene);

	for(int i = 0; i < n_gene; i ++) {
		for(int j = i + 1; j < n_gene - 1; j ++) {
			m(i, j) = logical_sum(match_mat(i, _), match_mat(j, _));
			m(j, i) = m(i, j);
		}
	}
	for(int i = 0; i < n_gene; i ++) {
		m(i, i) = 0;
	}
	return m;
}
