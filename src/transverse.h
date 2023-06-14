#ifndef __TRANSVERSE__
#define __TRANSVERSE__

extern const int SET_UNION;
extern const int SET_INTERSECT;
extern const int SET_UNIQU_IN_1;
extern const int SET_UNIQU_IN_2;

void _add_parents(List lt_parents, int i_node, LogicalVector& l_ancestors);
void _add_parents_within_background(List lt_parents, int i_node, LogicalVector& l_ancestors, LogicalVector l_background);
void _find_ancestors(List lt_parents, int i_node, LogicalVector& l_ancestors, bool include_self = false);
void _find_ancestors_with_background(List lt_parents, int i_node, LogicalVector& l_ancestors, LogicalVector l_background, bool include_self = false);
IntegerVector cpp_ancestors(S4 dag, int node, bool include_self = false);
IntegerVector cpp_ancestors_within_background(S4 dag, int node, IntegerVector background, bool include_self = false);
void _add_children(List lt_children, int i_node, LogicalVector& l_offspring);
void _add_children_within_background(List lt_children, int i_node, LogicalVector& l_offspring, LogicalVector l_background);
void _find_offspring(List lt_children, int i_node, LogicalVector& l_offspring, bool include_self = false);
void _find_offspring_within_background(List lt_children, int i_node, LogicalVector& l_offspring, LogicalVector l_background, bool include_self = false);
IntegerVector cpp_offspring(S4 dag, int node, bool include_self = false);
IntegerVector cpp_offspring_within_background(S4 dag, int node, IntegerVector background, bool include_self = false);
void _add_leaves(List lt_children, int i_node, LogicalVector& l_offspring);
void _find_connected_leaves(List lt_children, int i_node, LogicalVector& l_offspring);

IntegerVector cpp_n_ancestors(S4 dag, bool include_self = false);
IntegerVector cpp_n_offspring(S4 dag, bool include_self = false);
IntegerVector cpp_n_leaves(S4 dag);
IntegerVector cpp_ancestors_of_a_group(S4 dag, IntegerVector nodes, int type = 1, bool include_self = false);
IntegerVector cpp_ancestors_of_a_group_within_background(S4 dag, IntegerVector nodes, IntegerVector background, int type = 1, bool include_self = false);
IntegerVector cpp_ancestors_of_two_groups(S4 dag, IntegerVector nodes1, IntegerVector nodes2, int type, bool include_self = false);
IntegerVector cpp_offspring_of_a_group(S4 dag, IntegerVector nodes, bool include_self = false);
LogicalMatrix cpp_is_reachable(S4 dag, IntegerVector nodes, bool directed = false);

IntegerVector cpp_dag_depth(S4 dag);
IntegerVector cpp_dag_longest_dist_to_offspring(S4 dag, IntegerVector from_node, LogicalVector l_background = LogicalVector(0));
IntegerVector cpp_dag_longest_dist_to_offspring(S4 dag, int from_node, LogicalVector l_background = LogicalVector(0));
IntegerVector cpp_dag_shortest_dist_to_offspring(S4 dag, IntegerVector from_node, LogicalVector l_background = LogicalVector(0));
IntegerVector cpp_dag_shortest_dist_to_offspring(S4 dag, int from_node, LogicalVector l_background = LogicalVector(0));
IntegerVector cpp_dag_height(S4 dag);
IntegerVector cpp_dag_longest_dist_to_ancestors(S4 dag, IntegerVector from_node, LogicalVector l_background = LogicalVector(0));
IntegerVector cpp_dag_longest_dist_to_ancestors(S4 dag, int from_node, LogicalVector l_background = LogicalVector(0));
IntegerVector cpp_dag_shortest_dist_to_ancestors(S4 dag, IntegerVector from_node, LogicalVector l_background = LogicalVector(0));
IntegerVector cpp_dag_shortest_dist_to_ancestors(S4 dag, int from_node, LogicalVector l_background = LogicalVector(0));

#endif
