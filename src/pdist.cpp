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

float cosine(NumericVector x1, NumericVector x2) {
	float y1 = 0;
	float y2 = 0;
	float y3 = 0;
	float y = 0;
	int n = x1.size();
	
	for(int i = 0; i < n; i ++) {
		y1 += x1[i]*x2[i];
		y2 += x1[i]*x1[i];
		y3 += x2[i]*x2[i];
	}

	y = 1 - y1/std::sqrt(y2)/std::sqrt(y3);
	return y;
}


// [[Rcpp::export]]
NumericMatrix pdist(NumericMatrix m1, NumericMatrix m2, int dm) {
	// distance between every row in m1 to every row in m2

	int nr = m1.nrow();
	int nc = m2.nrow();

	NumericMatrix md(nr, nc);

	if(dm == 1) {
		for(int i = 0; i < nr; i ++) {
			for(int j = 0; j < nc; j ++) {
				md(i, j) = euclidean(m1(i, _), m2(j, _));
			}
		}
	} else {
		for(int i = 0; i < nr; i ++) {
			for(int j = 0; j < nc; j ++) {
				md(i, j) = cosine(m1(i, _), m2(j, _));
			}
		}
	}
	return md;
}
