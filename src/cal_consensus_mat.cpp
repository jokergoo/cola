#include <Rcpp.h>
using namespace Rcpp;

float p_equal(IntegerVector x1, IntegerVector x2) {
	int n = x1.size();
	int x = 0;
	int s = 0;
	for(int i = 0; i < n; i ++) {
		if(x1[i] == 0 || x2[i] == 0) {
			continue;
		}
		if(x1[i] == x2[i]) {
			x ++;
		}
		s ++;
	}
	return x/(s + 0.0);
}

// [[Rcpp::export]]
NumericMatrix get_consensus_matrix(IntegerMatrix membership_each) {
	int n = membership_each.nrow();

	NumericMatrix consensus_mat(n, n);

	for(int i = 0; i < n - 1; i ++) {
		for(int j = i+1; j < n; j ++) {
			consensus_mat(i, j) = p_equal(membership_each(i, _), membership_each(j, _));
			consensus_mat(j, i) = consensus_mat(i, j);
		}
	}
	for(int i = 0; i < n; i ++) {
		consensus_mat(i, i) = 1;
	}
	return consensus_mat;
}

