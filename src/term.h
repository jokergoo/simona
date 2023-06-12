#ifndef __TERMS__
#define __TERMS__


double _calc_wang_s(List lt_children, List lt_children_relations, NumericVector contribution, int i_node, int i_end, LogicalVector l_background);
NumericVector cpp_ic_wang(S4 dag, NumericVector contribution);
IntegerVector cpp_max_leaves_id(S4 dag, IntegerVector nodes, NumericVector v);

#endif

