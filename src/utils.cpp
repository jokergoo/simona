
#include <Rcpp.h>
using namespace Rcpp;


void reset_logical_vector_to_false(LogicalVector& x) {
	int n = x.size();
	for(int i = 0; i < n; i ++) {
		x[i] = false;
	}
}

void reset_logical_vector_to_true(LogicalVector& x) {
	int n = x.size();
	for(int i = 0; i < n; i ++) {
		x[i] = true;
	}
}

void fill_logical_vector_to_true(LogicalVector&x, IntegerVector ind) {
	for(int i = 0; i < ind.size(); i ++) {
		x[ind[i]] = true;
	}
}

void reset_numeric_vector_to_zero(NumericVector& x) {
	for(int i = 0; i < x.size(); i ++) {
		x[i] = 0;
	}
}

void reset_integer_vector_to_zero(IntegerVector& x) {
	for(int i = 0; i < x.size(); i ++) {
		x[i] = 0;
	}
}

IntegerVector _dag_depth(S4 dag) {
	Environment term_env = dag.slot("term_env");
	IntegerVector depth = term_env["dag_depth"];
	return depth;
}


// return the indices where elements in `l` is true
IntegerVector _which(LogicalVector l) {
	int n = l.size();
	int n2 = sum(l);
	IntegerVector ind(n2);

	if(n2 == 0) {
		return(ind);
	}

	int i2 = 0;
	for(int i = 0; i < n; i ++) {
		if(l[i]) {
			ind[i2] = i;
			i2 ++;
		}
	}

	return ind;
}

LogicalVector integer_to_logical_vector(IntegerVector i, int n) {
	LogicalVector l(n);
	for(int k = 0; k < i.size(); k ++) {
		l[i[k]] = true;
	}
	return l;
}


// [[Rcpp::export]]
IntegerVector cpp_match_index(IntegerVector ind1, IntegerVector ind2) {

	// ind1 and ind2 are all sorted, for each i in ind1, assign the index in ind2 which is the nearest larger value to i
	int n1 = ind1.size();

	IntegerVector ind(n1);

	int i2 = 0;
	for(int i1 = 0; i1 < n1; i1 ++) {
		while( true ) {
			if(ind1[i1] < ind2[i2]) {
				ind[i1] = ind2[i2];
				break;
			} else {
				i2 ++;
			}
		}
	}

	return ind;
}


void message(String s, bool appendLF = true) {
	Function _message("message");
	_message(s, Named("appendLF", appendLF));
	return;
}
