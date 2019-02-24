#include <Rcpp.h>
using namespace Rcpp;

float nequal(IntegerVector x1, IntegerVector x2) {
	int n = x1.size();
	int x = 0;
	for(int i = 0; i < n; i ++) {
		if(x1[i] == x2[i]) {
			x ++;
		}
	}
	return x;
}

// [[Rcpp::export]]
NumericMatrix get_consensus_matrix(IntegerMatrix membership_each) {
	int n = membership_each.nrow();
	int nc = membership_each.ncol();

	NumericMatrix consensus_mat(n, n);

	for(int i = 0; i < n - 1; i ++) {
		for(int j = i+1; j < n; j ++) {
			consensus_mat(i, j) = nequal(membership_each(i, _), membership_each(j, _))/(nc + 0.0);
			consensus_mat(j, i) = consensus_mat(i, j);
		}
	}
	for(int i = 0; i < n; i ++) {
		consensus_mat(i, i) = 1;
	}
	return consensus_mat;
}
