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

// [[Rcpp::export]]
NumericMatrix cal_diff_ratio_r(NumericMatrix mat, NumericMatrix sig_mat, int n_perm, int dm) {
	// mat: a n_sample x n_sig matrix
	// sig_mat: a n_group x n_sig matrix
	// return diff_ratio_r: a n_sample x n_perm matrix
	int n_sample = mat.nrow();
	int n_sig = mat.ncol();
	int n_group = sig_mat.nrow();

	NumericMatrix diff_ratio_r(n_sample, n_perm);

	for(int i = 0; i < n_perm; i ++) {
		// a n_sample x n_group matrix
		NumericMatrix sig_mat_r(n_group, n_sig);
		for(int k = 0; k < n_group; k ++) {
			NumericVector foo = sig_mat(k, _);
			NumericVector v = clone(foo);
			std::random_shuffle(v.begin(), v.end());
			sig_mat_r(k, _) = v;
		}
		NumericMatrix dist_to_signatures_r = pdist(mat, sig_mat_r, dm);

		for(int j = 0; j < n_sample; j ++) {
			NumericVector x = dist_to_signatures_r(j, _);
			std::sort(x.begin(), x.end());
			diff_ratio_r(j, i) = std::abs(x[0] - x[1])/(std::accumulate(std::begin(x), std::end(x), 0.0) / x.size());
		}
	}

	return diff_ratio_r;
}
