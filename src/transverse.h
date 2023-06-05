#ifndef __TRANSVERSE__
#define __TRANSVERSE__


void _add_parents(List lt_parents, int i_node, LogicalVector& l_ancestor);
void _add_parents_within_background(List lt_parents, int i_node, LogicalVector& l_ancestor, LogicalVector l_background);
void _find_ancestor(List lt_parents, int i_node, LogicalVector& l_ancestor, bool include_self = false);
void _find_ancestor_with_background(List lt_parents, int i_node, LogicalVector& l_ancestor, LogicalVector l_background, bool include_self = false);
IntegerVector cpp_ancestor(S4 dag, int node, bool include_self = false);
IntegerVector cpp_ancestor_within_background(S4 dag, int node, IntegerVector background, bool include_self = false);
void _add_children(List lt_children, int i_node, LogicalVector& l_offspring);
void _add_children_within_background(List lt_children, int i_node, LogicalVector& l_offspring, LogicalVector l_background);
void _find_offspring(List lt_children, int i_node, LogicalVector& l_offspring, bool include_self = false);
void _find_offspring_within_background(List lt_children, int i_node, LogicalVector& l_offspring, LogicalVector l_background, bool include_self = false);
IntegerVector cpp_offspring(S4 dag, int node, bool include_self = false);
IntegerVector cpp_offspring_within_background(S4 dag, int node, IntegerVector background, bool include_self = false);
void _add_leaves(List lt_children, int i_node, LogicalVector& l_offspring);
void _find_connected_leaves(List lt_children, int i_node, LogicalVector& l_offspring);
IntegerVector cpp_ancestor_of_a_group(S4 dag, IntegerVector nodes, int type = 1, bool include_self = false);
IntegerVector cpp_ancestor_of_a_group_within_background(S4 dag, IntegerVector nodes, IntegerVector background, int type = 1, bool include_self = false);
IntegerVector cpp_ancestor_of_two_groups(S4 dag, IntegerVector nodes1, IntegerVector nodes2, int type, bool include_self = false);
IntegerVector cpp_offspring_of_a_group(S4 dag, IntegerVector nodes, bool include_self = false);
LogicalMatrix cpp_is_reachable(S4 dag, IntegerVector nodes, bool directed = false);

#endif
