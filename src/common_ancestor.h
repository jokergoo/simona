#ifndef __ANCESTOR__
#define __ANCESTOR__


NumericMatrix cpp_max_ancestor_v(S4 dag, IntegerVector nodes, NumericVector v);
IntegerMatrix cpp_max_ancestor_id(S4 dag, IntegerVector nodes, NumericVector v);
IntegerMatrix cpp_distances(S4 dag, IntegerVector nodes, int type = 1);
IntegerMatrix cpp_shortest_distances_via_NCA(S4 dag, IntegerVector nodes);
IntegerMatrix cpp_distances_directed(S4 dag, IntegerVector nodes, int type = 1);
IntegerMatrix cpp_longest_distances_directed(S4 dag, IntegerVector nodes);
IntegerMatrix cpp_shortest_distances_directed(S4 dag, IntegerVector nodes);
IntegerMatrix cpp_nearest_common_ancestor(S4 dag, IntegerVector nodes);

#endif

