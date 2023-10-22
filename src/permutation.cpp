#include <Rcpp.h>
using namespace Rcpp;
#include "utils.h"


// [[Rcpp::export]]
NumericMatrix cpp_random_aggregatioin(IntegerVector size, NumericVector value, int perm) {
	int n = size.size();

	NumericMatrix m(n, perm);
	NumericVector v2;

	for(int i = 0; i < perm; i ++) {

		if(i % 10 == 0) {
			message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
			message("going through " + std::to_string(i) + " / " + std::to_string(perm) + " permutations ...", false);
		}

		for(int j = 0; j < n; j ++) {
			v2 = sample(value, size[j]);
			m(j, i) = mean(v2);
		}
	}

	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
	message("going through " + std::to_string(perm) + " / " + std::to_string(perm) + " permutations ... Done.", true);


	return m;
}
