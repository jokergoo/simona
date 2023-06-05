
#ifndef __DIST__
#define __DIST__

double cpp_find_path_length_single(S4 dag, int from, int to, NumericVector weight, int type);
IntegerVector cpp_find_path_single(S4 dag, int from, int to, NumericVector weight, int type);
int cpp_distance_single(S4 dag, int from, int to);
int cpp_longest_distance_single(S4 dag, int from, int to);;

#endif
