#include <Rcpp.h>

using namespace Rcpp;
	

float singleATC(NumericVector x, float min_cor, float power, int k_neighbours) {
	
	std::sort(x.begin(), x.end());

	int n = x.size();
		
	if(k_neighbours > 0) {
		for(int i = 0; i < n - k_neighbours; i ++) {
			x[i] = -1;
		}
	}

	int n2 = 0;
	for(int i = 0; i < n; i ++) {
		if(x[i] > -0.5) {
			n2 ++;
		}
	}

	float s = 0;
	int j = 0;
	for(int i = 0; i < n; i ++) {
		if(x[i] > -0.5) {
			if(i == n - 1) {
				if(x[i] >= min_cor) {
					s += (1 - pow(x[i], power))*1;
				} else {
					s += (1 - pow(min_cor, power))*1;
				}
			} else {
				if(x[i] >= min_cor) {
					s += (pow(x[i+1], power) - pow(x[i], power))*(j+1)/n2;
				} else if(x[i] < min_cor && x[i+1] >= min_cor) {
					s += (pow(x[i+1], power) - pow(min_cor, power))*(j+1)/n2;
				}
				j ++;
			}
		}
	}
	
	s = 1 - min_cor - s;
	return(s);
}

// [[Rcpp::export]]
NumericVector rowATC(NumericMatrix m, float min_cor, float power, int k_neighbours, IntegerVector self) {
	int n = m.nrow();
	NumericVector s(n);

	for(int i = 0; i < n; i ++) {
		NumericVector x = m(i, _);
		x[ self[i]-1 ] = -1;
		s[i] = singleATC(x, min_cor, power, k_neighbours);
	}
	return(s);
}

