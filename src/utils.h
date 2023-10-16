#ifndef __UTILS__
#define __UTILS__

void reset_logical_vector_to_false(LogicalVector& x);
void reset_logical_vector_to_true(LogicalVector& x);
void reset_numeric_vector_to_zero(NumericVector& x);
void reset_integer_vector_to_zero(IntegerVector& x);
IntegerVector _dag_depth(S4 dag);
IntegerVector _which(LogicalVector l);
LogicalVector integer_to_logical_vector(IntegerVector i, int n);
IntegerVector cpp_match_index(IntegerVector ind1, IntegerVector ind2);
void message(String s, bool appendLF = true);
LogicalVector merge_two_logical_vectors_by_or(LogicalVector l1, LogicalVector l2);
LogicalVector merge_two_logical_vectors_by_and(LogicalVector l1, LogicalVector l2);
IntegerVector _order(NumericVector x);
IntegerVector _order(IntegerVector x);

#endif

