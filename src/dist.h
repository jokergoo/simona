#ifndef __DIST__
#define __DIST__

int cpp_tpl_shortest_path_length(S4 dag, int from, int to);
int cpp_tpl_longest_path_length(S4 dag, int from, int to);
IntegerVector cpp_tpl_shortest_path(S4 dag, int from, int to);
IntegerVector cpp_tpl_longest_path(S4 dag, int from, int to);

double cpp_tpl_shortest_path_sum_value(S4 dag, int from, int to, NumericVector value);
double cpp_tpl_longest_path_sum_value(S4 dag, int from, int to, NumericVector value);

#endif
