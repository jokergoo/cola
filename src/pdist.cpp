#include <Rcpp.h>
using namespace Rcpp;


float euclidean(NumericVector x1, NumericVector x2) {
	float y = 0;
	int n = x1.size();
	
	for(int i = 0; i < n; i ++) {
		y += (x1[i] - x2[i])*(x1[i] - x2[i]);
	}

	y = std::sqrt(y);
	return y;
}


// [[Rcpp::export]]
NumericMatrix pdist(NumericMatrix m1, NumericMatrix m2) {
	// distance between every row in m1 to every row in m2

	int nr = m1.nrow();
	int nc = m2.nrow();

	NumericMatrix md(nr, nc);

	for(int i = 0; i < nr; i ++) {
		for(int j = 0; j < nc; j ++) {
			md(i, j) = euclidean(m1(i, _), m2(j, _));
		}
	}
	return md;
}
