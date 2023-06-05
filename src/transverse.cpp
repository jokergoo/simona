#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"

void _add_parents(List lt_parents, int i_node, LogicalVector& l_ancestor) {
	IntegerVector parents = lt_parents[i_node];
	if(parents.size() > 0) {
		for(int i = 0; i < parents.size(); i ++) {
			int i_parent = parents[i] - 1;
			l_ancestor[i_parent] = true;
			_add_parents(lt_parents, i_parent, l_ancestor);
		}
	}
}

void _add_parents_within_background(List lt_parents, int i_node, LogicalVector& l_ancestor, LogicalVector l_background) {
	if(l_background[i_node]) {
		IntegerVector parents = lt_parents[i_node];
		if(parents.size() > 0) {
			for(int i = 0; i < parents.size(); i ++) {
				int i_parent = parents[i] - 1;
				if(l_background[i_parent]) {
					l_ancestor[i_parent] = true;
					_add_parents_within_background(lt_parents, i_parent, l_ancestor, l_background);
				}
			}
		}
	}
}

// find all ancestor node for a given node
void _find_ancestor(List lt_parents, int i_node, LogicalVector& l_ancestor, bool include_self = false) {
	_add_parents(lt_parents, i_node, l_ancestor);
	if(include_self) {
		l_ancestor[i_node] = true;
	}
}

void _find_ancestor_with_background(List lt_parents, int i_node, LogicalVector& l_ancestor, LogicalVector l_background, bool include_self = false) {
	_add_parents_within_background(lt_parents, i_node, l_ancestor, l_background);
	if(include_self) {
		l_ancestor[i_node] = true;
	}
}


// [[Rcpp::export]]
IntegerVector cpp_ancestor(S4 dag, int node, bool include_self = false) {
	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestor(n);
	_find_ancestor(lt_parents, node - 1, l_ancestor, include_self);

	IntegerVector ancestor = _which(l_ancestor);
	if(ancestor.size() > 0) {
		ancestor = ancestor + 1;
	}
	return ancestor;
}

// [[Rcpp::export]]
IntegerVector cpp_ancestor_within_background(S4 dag, int node, IntegerVector background, bool include_self = false) {
	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestor(n);
	LogicalVector l_background(n);
	for(int i = 0; i < background.size(); i ++) {
		l_background[background[i] - 1] = true;
	}
	_find_ancestor_with_background(lt_parents, node - 1, l_ancestor, l_background, include_self);

	IntegerVector ancestor = _which(l_ancestor);
	if(ancestor.size() > 0) {
		ancestor = ancestor + 1;
	}
	return ancestor;
}

void _add_children(List lt_children, int i_node, LogicalVector& l_offspring) {
	IntegerVector children = lt_children[i_node];
	if(children.size() > 0) {
		for(int i = 0; i < children.size(); i ++) {
			int i_child = children[i] - 1;
			l_offspring[i_child] = true;
			_add_children(lt_children, i_child, l_offspring);
		}
	}
}

void _add_children_within_background(List lt_children, int i_node, LogicalVector& l_offspring, LogicalVector l_background) {
	if(l_background[i_node]) {
		IntegerVector children = lt_children[i_node];
		if(children.size() > 0) {
			for(int i = 0; i < children.size(); i ++) {
				int i_child = children[i] - 1;
				if(l_background[i_child]) {
					l_offspring[i_child] = true;
					_add_children_within_background(lt_children, i_child, l_offspring, l_background);
				}
			}
		}
	}
}

// find all offspring node for a given node
void _find_offspring(List lt_children, int i_node, LogicalVector& l_offspring, bool include_self = false) {
	_add_children(lt_children, i_node, l_offspring);
	if(include_self) {
		l_offspring[i_node] = true;
	}
}

void _find_offspring_within_background(List lt_children, int i_node, LogicalVector& l_offspring, LogicalVector l_background, bool include_self = false) {
	_add_children_within_background(lt_children, i_node, l_offspring, l_background);
	if(include_self) {
		l_offspring[i_node] = true;
	}
}

// [[Rcpp::export]]
IntegerVector cpp_offspring(S4 dag, int node, bool include_self = false) {
	List lt_children = dag.slot("lt_children");
	int n = lt_children.size();

	LogicalVector l_offspring(n);
	_find_offspring(lt_children, node - 1, l_offspring, include_self);

	IntegerVector offspring = _which(l_offspring);
	if(offspring.size() > 0) {
		offspring = offspring + 1;
	}

	return offspring;
}

// [[Rcpp::export]]
IntegerVector cpp_offspring_within_background(S4 dag, int node, IntegerVector background, bool include_self = false) {
	List lt_children = dag.slot("lt_children");
	int n = lt_children.size();

	LogicalVector l_offspring(n);
	LogicalVector l_background(n);
	for(int i = 0; i < background.size(); i ++) {
		l_background[background[i] - 1] = true;
	}
	_find_offspring_within_background(lt_children, node - 1, l_offspring, l_background, include_self);

	IntegerVector offspring = _which(l_offspring);
	if(offspring.size() > 0) {
		offspring = offspring + 1;
	}

	return offspring;
}


void _add_leaves(List lt_children, int i_node, LogicalVector& l_offspring) {
	IntegerVector children = lt_children[i_node];
	if(children.size() > 0) {
		for(int i = 0; i < children.size(); i ++) {
			int i_child = children[i] - 1;
			_add_leaves(lt_children, i_child, l_offspring);
		}
	} else {
		l_offspring[i_node] = true;
	}
}

void _find_connected_leaves(List lt_children, int i_node, LogicalVector& l_offspring) {
	_add_leaves(lt_children, i_node, l_offspring);
}

// [[Rcpp::export]]
IntegerVector cpp_connected_leaves(S4 dag, int node) {
	List lt_children = dag.slot("lt_children");
	int n = lt_children.size();

	LogicalVector l_leaf(n);
	_find_connected_leaves(lt_children, node - 1, l_leaf);

	IntegerVector leaves = _which(l_leaf);
	if(leaves.size() > 0) {
		leaves = leaves + 1;
	}

	return leaves;
}

// union of ancestors of a group of nodes
// [[Rcpp::export]]
IntegerVector cpp_ancestor_of_a_group(S4 dag, IntegerVector nodes, int type = 1, bool include_self = false) {  // type 1: union; 2: intersect
	int m = nodes.size();

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestor(n);
	if(type == 1) {
		for(int i = 0; i < m; i ++) {
			_find_ancestor(lt_parents, nodes[i] - 1, l_ancestor, include_self);
		}
	} else {
		LogicalVector l(n, true);
		LogicalVector l_an(n, false);
		for(int i = 0; i < m; i ++) {
			_find_ancestor(lt_parents, nodes[i] - 1, l_an, include_self);
			l = l & l_an;
			reset_logical_vector_to_false(l_an);
		}
		l_ancestor = l;
	}
	IntegerVector aid = _which(l_ancestor);
	if(aid.size() > 0) {
		aid = aid + 1;
	}
	return aid;
}

// [[Rcpp::export]]
IntegerVector cpp_ancestor_of_a_group_within_background(S4 dag, IntegerVector nodes, IntegerVector background, int type = 1, bool include_self = false) {  // type 1: union; 2: intersect
	int m = nodes.size();

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestor(n);
	LogicalVector l_background(n);
	for(int i = 0; i < background.size(); i ++) {
		l_background[background[i] - 1] = true;
	}
	if(type == 1) {
		for(int i = 0; i < m; i ++) {
			_find_ancestor_with_background(lt_parents, nodes[i] - 1, l_ancestor, l_background, include_self);
		}
	} else {
		LogicalVector l(n, true);
		LogicalVector l_an(n, false);
		for(int i = 0; i < m; i ++) {
			_find_ancestor_with_background(lt_parents, nodes[i] - 1, l_an, l_background, include_self);
			l = l & l_an;
			reset_logical_vector_to_false(l_an);
		}
		l_ancestor = l;
	}
	IntegerVector aid = _which(l_ancestor);
	if(aid.size() > 0) {
		aid = aid + 1;
	}
	return aid;
}

// [[Rcpp::export]]
IntegerVector cpp_ancestor_of_two_groups(S4 dag, IntegerVector nodes1, IntegerVector nodes2, int type, bool include_self = false) { // type 1: union, 2: intersect, 3: only in node1, 4: only in node2
	int m1 = nodes1.size();
	int m2 = nodes2.size();

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestor1(n);
	LogicalVector l_ancestor2(n);
	for(int i = 0; i < m1; i ++) {
		_find_ancestor(lt_parents, nodes1[i] - 1, l_ancestor1, include_self);
	}
	for(int i = 0; i < m2; i ++) {
		_find_ancestor(lt_parents, nodes2[i] - 1, l_ancestor2, include_self);
	}

	LogicalVector l_ancestor(n);
	if(type == 1) {
		l_ancestor = l_ancestor1 | l_ancestor2;
	} else if(type == 2) {
		l_ancestor = l_ancestor1 & l_ancestor2;
	} else if(type == 3) {
		l_ancestor = l_ancestor1 & !l_ancestor2;
	} else if(type == 4) {
		l_ancestor = !l_ancestor1 & l_ancestor2;
	}

	IntegerVector aid = _which(l_ancestor);
	if(aid.size() > 0) {
		aid = aid + 1;
	}
	return aid;
}


// union of ancestors of a group of nodes
// [[Rcpp::export]]
IntegerVector cpp_offspring_of_a_group(S4 dag, IntegerVector nodes, bool include_self = false) {
	int m = nodes.size();

	List lt_children = dag.slot("lt_children");
	int n = lt_children.size();

	LogicalVector l_offspring(n);
	for(int i = 0; i < m; i ++) {
		_find_offspring(lt_children, nodes[i] - 1, l_offspring, include_self);
	}

	IntegerVector aid = _which(l_offspring);
	if(aid.size() > 0) {
		aid = aid + 1;
	}
	return aid;
}


// whether node i and node j are ancestor/offspring
// [[Rcpp::export]]
LogicalMatrix cpp_is_reachable(S4 dag, IntegerVector nodes, bool directed = false) {
	List lt_parents = dag.slot("lt_parents");
	List lt_children = dag.slot("lt_children");
	int n_all = lt_parents.size();
	int root = dag.slot("root");

	IntegerVector depth = _dag_depth(dag);
	int max_depth = max(depth);

	int n = nodes.size();
	LogicalMatrix m(n, n);
	m.fill_diag(true);

	if(n <= 1) {
		return m;
	}

	LogicalVector l_ancestor(n_all);
	LogicalVector l_offspring(n_all);
	bool flag = false;

	for(int i = 0; i < n - 1; i ++) {
		if(nodes[i] == root) {
			for(int j = i+1; j < n; j ++) {
				m(i, j) = true;
				if(!directed) {
					m(j, i) = true;
				}
				flag = true;
			}
		}

		if(!flag) {
			
			// if i_node is closer to the root, first scan its ancestor nodes
			if(depth[i] < max_depth*0.5) {
				if(!directed) {
					_find_ancestor(lt_parents, nodes[i] - 1, l_ancestor);
					// whether j_node is an ancestor of i_node
					for(int j = i+1; j < n; j ++) {
						if(l_ancestor[ nodes[j] - 1 ]) {
							m(i, j) = true;
							if(!directed) {
								m(j, i) = true;
							}
							flag = true;
						}
					}
				}

				// whether j_node is an offspring of i_node
				if(!flag) {
					_find_offspring(lt_children, nodes[i] - 1, l_offspring);
					for(int j = i+1; j < n; j ++) {
						if(l_offspring[ nodes[j] - 1 ]) {
							m(i, j) = true;
							if(!directed) {
								m(j, i) = true;
							}
						}
					}
					reset_logical_vector_to_false(l_offspring);
				}
			} else {
				_find_offspring(lt_children, nodes[i] - 1, l_offspring);
				for(int j = i+1; j < n; j ++) {
					if(l_offspring[ nodes[j] - 1 ]) {
						m(i, j) = true;
						if(!directed) {
							m(j, i) = true;
						}
					}
				}

				// whether j_node is an offspring of i_node
				if(!flag & directed) {
					_find_ancestor(lt_parents, nodes[i] - 1, l_ancestor);
					// whether j_node is an ancestor of i_node
					for(int j = i+1; j < n; j ++) {
						if(l_ancestor[ nodes[j] - 1 ]) {
							m(i, j) = true;
							if(!directed) {
								m(j, i) = true;
							}
							flag = true;
						}
					}
					reset_logical_vector_to_false(l_ancestor);
				}
			}
		}

		flag = false;
		reset_logical_vector_to_false(l_ancestor);
		reset_logical_vector_to_false(l_offspring);
	}

	return m;
}
