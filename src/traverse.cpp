#include <Rcpp.h>
using namespace Rcpp;

#include "utils.h"


// used in _find_ancestors()
void _add_parents(List lt_parents, int i_node, LogicalVector& l_ancestors) {
	IntegerVector parents = lt_parents[i_node];
	if(parents.size() > 0) {
		for(int i = 0; i < parents.size(); i ++) {
			int i_parent = parents[i] - 1;
			if(!l_ancestors[i_parent]) {
				l_ancestors[i_parent] = true;
				_add_parents(lt_parents, i_parent, l_ancestors);
			}
		}
	}
}

// find all ancestor node for a given node
void _find_ancestors(List lt_parents, int i_node, LogicalVector& l_ancestors, bool include_self = false) {
	_add_parents(lt_parents, i_node, l_ancestors);
	if(include_self) {
		l_ancestors[i_node] = true;
	}
}


// used in _find_ancestors_with_background()
void _add_parents_within_background(List lt_parents, int i_node, LogicalVector& l_ancestors, LogicalVector l_background) {
	if(l_background[i_node]) {
		IntegerVector parents = lt_parents[i_node];
		if(parents.size() > 0) {
			for(int i = 0; i < parents.size(); i ++) {
				int i_parent = parents[i] - 1;
				if(l_background[i_parent] && !l_ancestors[i_parent]) {
					l_ancestors[i_parent] = true;
					_add_parents_within_background(lt_parents, i_parent, l_ancestors, l_background);
				}
			}
		}
	}
}

void _find_ancestors_with_background(List lt_parents, int i_node, LogicalVector& l_ancestors, LogicalVector l_background, bool include_self = false) {
	_add_parents_within_background(lt_parents, i_node, l_ancestors, l_background);
	if(include_self) {
		l_ancestors[i_node] = true;
	}
}


// [[Rcpp::export]]
IntegerVector cpp_ancestors(S4 dag, int node, bool include_self = false) {
	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestors(n);
	_find_ancestors(lt_parents, node - 1, l_ancestors, include_self);

	IntegerVector ancestors = _which(l_ancestors);
	if(ancestors.size() > 0) {
		ancestors = ancestors + 1;
	}
	return ancestors;
}

// [[Rcpp::export]]
IntegerVector cpp_ancestors_within_background(S4 dag, int node, IntegerVector background, bool include_self = false) {
	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestors(n);
	LogicalVector l_background = integer_to_logical_vector(background - 1, n);

	_find_ancestors_with_background(lt_parents, node - 1, l_ancestors, l_background, include_self);

	IntegerVector ancestors = _which(l_ancestors);
	if(ancestors.size() > 0) {
		ancestors = ancestors + 1;
	}
	return ancestors;
}

void _add_children(List lt_children, int i_node, LogicalVector& l_offspring) {
	IntegerVector children = lt_children[i_node];
	if(children.size() > 0) {
		for(int i = 0; i < children.size(); i ++) {
			int i_child = children[i] - 1;
			if(!l_offspring[i_child]) {
				l_offspring[i_child] = true;
				_add_children(lt_children, i_child, l_offspring);
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

void _add_children_within_background(List lt_children, int i_node, LogicalVector& l_offspring, LogicalVector l_background) {
	if(l_background[i_node]) {
		IntegerVector children = lt_children[i_node];
		if(children.size() > 0) {
			for(int i = 0; i < children.size(); i ++) {
				int i_child = children[i] - 1;
				if(l_background[i_child] && !l_offspring[i_child]) {
					l_offspring[i_child] = true;
					_add_children_within_background(lt_children, i_child, l_offspring, l_background);
				}
			}
		}
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
LogicalMatrix cpp_all_offspring(S4 dag, bool include_self = false) {

	// offsprings for all nodes,
	// it returns a binary matrix where m[i, j] = 1 means term j is a offspring of term i

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	// offspring on columns
	LogicalMatrix m_offspring(n, n);
	if(include_self) {
		m_offspring.fill_diag(true);
	}

	IntegerVector leaves = dag.slot("leaves");
	LogicalVector l_current_nodes(n);
	for(int i = 0; i < leaves.size(); i ++) {
		l_current_nodes[leaves[i]-1] = true;
	}
	LogicalVector l_current_nodes2(n);

	while(sum(l_current_nodes)) {
		for(int i = 0; i < n; i ++) {
			if(l_current_nodes[i]) {
				IntegerVector parents = lt_parents[i];

				for(int j = 0; j < parents.size(); j ++) {
					l_current_nodes2[parents[j]-1] = true;
					m_offspring(parents[j]-1, i) = true;

					for(int k = 0; k < n; k ++) {
						if(m_offspring(i, k)) {
							m_offspring(parents[j]-1, k) = true;
						}
					}
				}
			}
		}

		for(int i = 0 ; i < n; i ++) {
			l_current_nodes[i] = l_current_nodes2[i];
		}
		reset_logical_vector_to_false(l_current_nodes2);
	}

	return m_offspring;
}

// [[Rcpp::export]]
IntegerVector cpp_offspring_within_background(S4 dag, int node, IntegerVector background, bool include_self = false) {
	List lt_children = dag.slot("lt_children");
	int n = lt_children.size();

	LogicalVector l_offspring(n);
	LogicalVector l_background = integer_to_logical_vector(background - 1, n);

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


// -------------------------
// number of connected leaves
// number of ancestor nodes
// number of offspring nodes


// [[Rcpp::export]]
IntegerVector cpp_n_ancestors(S4 dag, bool include_self = false) {
	List lt_parents = dag.slot("lt_parents");

	int n = lt_parents.size();
	IntegerVector num(n, 0);

	LogicalVector l_ancestors(n, false);
	for(int i = 0; i < n; i ++) {
			
		_find_ancestors(lt_parents, i, l_ancestors, include_self);
		num[i] = sum(l_ancestors);

		reset_logical_vector_to_false(l_ancestors);

	}
	return num;
}


// [[Rcpp::export]]
IntegerVector cpp_n_ancestors_on_tree(S4 dag, bool include_self = false) {

	// faster than treating it as a DAG
	List lt_children = dag.slot("lt_children");

	int n = lt_children.size();
	IntegerVector num(n, 0);

	IntegerVector current = dag.slot("root");
	LogicalVector l_current = integer_to_logical_vector(current-1, n);
	bool has_child = true;
	while(has_child) {
		LogicalVector l_current2(n);
		has_child = false;
		for(int i = 0; i < n; i ++) {
			if(l_current[i]) {
				IntegerVector children = lt_children[i];
				if(children.size()) {

					for(int j = 0; j < children.size(); j ++) {
						num[ children[j] - 1 ] += num[ i ] + 1;
						l_current2[ children[j] -1] = true;
					}
					has_child = true;
				}
			}
		}
		l_current = l_current2;
	}

	if(include_self) {
		num = num + 1;
	}
	
	return num;
}


// [[Rcpp::export]]
IntegerVector cpp_n_offspring(S4 dag, bool include_self = false) {
	List lt_children = dag.slot("lt_children");

	int n = lt_children.size();
	IntegerVector num(n, 0);

	LogicalVector l_offspring(n, false);
	for(int i = 0; i < n; i ++) {
		_find_offspring(lt_children, i, l_offspring, include_self);
		num[i] = sum(l_offspring);

		reset_logical_vector_to_false(l_offspring);

	}
	return num;
}


// [[Rcpp::export]]
IntegerVector cpp_n_offspring_on_tree(S4 dag, bool include_self = false) {
	
	// faster than treating it as a DAG
	List lt_children = dag.slot("lt_children");
	IntegerVector depth = _dag_depth(dag);

	int max_depth = max(depth);

	int n = lt_children.size();
	IntegerVector num(n, 0);

	for(int i_depth = max_depth; i_depth >= 0; i_depth --) {

		for(int i = 0; i < n; i ++) {
			if(depth[i] == i_depth) {
				IntegerVector children = lt_children[i];
				if(children.size()) {

					for(int j = 0; j < children.size(); j ++) {
						num[ i ] += num[ children[j]-1 ] + 1;
					}
				}
			}
		}
	}

	if(include_self) {
		num = num + 1;
	}
	
	return num;
}

// [[Rcpp::export]]
IntegerVector cpp_n_offspring_with_intersect(S4 dag, IntegerVector nodes, bool include_self = false) {
	List lt_children = dag.slot("lt_children");

	int n = lt_children.size();
	IntegerVector num(n, 0);

	int m = nodes.size();

	if(m == 0) {
		return num;
	}

	LogicalVector l_offspring(n, false);
	for(int i = 0; i < n; i ++) {			
		_find_offspring(lt_children, i, l_offspring, include_self);

		// offspring overlap with `nodes`
		for(int j = 0; j < m; j ++) {
			if(l_offspring[nodes[j]-1]) {
				num[i] ++;
			}
		}

		reset_logical_vector_to_false(l_offspring);
		
	}
	return num;
}

// [[Rcpp::export]]
IntegerVector cpp_n_leaves(S4 dag) {
	List lt_children = dag.slot("lt_children");

	int n = lt_children.size();
	IntegerVector num(n, 0);

	LogicalVector l_leaf(n, false);
	for(int i = 0; i < n; i ++) {
		IntegerVector children = lt_children[i];
		if(children.size() > 0) {
			
			_find_connected_leaves(lt_children, i, l_leaf);
			num[i] = sum(l_leaf);

			reset_logical_vector_to_false(l_leaf);
		} else {
			num[i] = 1;
		}
	}
	return num;
}

// [[Rcpp::export]]
IntegerVector cpp_n_leaves_on_tree(S4 dag) {
	
	// faster than treating it as a DAG
	List lt_children = dag.slot("lt_children");
	IntegerVector depth = _dag_depth(dag);

	int max_depth = max(depth);

	int n = lt_children.size();
	IntegerVector num(n, 0);

	for(int i_depth = max_depth; i_depth >= 0; i_depth --) {

		for(int i = 0; i < n; i ++) {
			if(depth[i] == i_depth) {
				IntegerVector children = lt_children[i];
				if(children.size() == 0) {
					num[i] = 1;
				} else {
					for(int j = 0; j < children.size(); j ++) {
						num[ i ] += num[ children[j]-1 ];
					}
				}
			}
		}
	}
	
	return num;
}

const int SET_UNION = 1;
const int SET_INTERSECT = 2;
const int SET_UNIQU_IN_1 = 3;
const int SET_UNIQU_IN_2 = 4;


// [[Rcpp::export]]
IntegerVector cpp_ancestors_of_a_group(S4 dag, IntegerVector nodes, int type = 1, bool include_self = false) {  // type 1: union; 2: intersect
	// union/intersection of ancestors of a group of nodes

	int m = nodes.size();

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestors(n);
	if(type == SET_UNION) {
		for(int i = 0; i < m; i ++) {
			_find_ancestors(lt_parents, nodes[i] - 1, l_ancestors, include_self);
		}
	} else {
		LogicalVector l(n, true);
		LogicalVector l_an(n, false);
		for(int i = 0; i < m; i ++) {
			_find_ancestors(lt_parents, nodes[i] - 1, l_an, include_self);
			l = l & l_an;
			reset_logical_vector_to_false(l_an);
		}
		l_ancestors = l;
	}
	IntegerVector aid = _which(l_ancestors);
	if(aid.size() > 0) {
		aid = aid + 1;
	}
	return aid;
}

IntegerVector cpp_ancestors_of_a_group_within_background(S4 dag, IntegerVector nodes, IntegerVector background, int type = SET_UNION, bool include_self = false) {  // type 1: union; 2: intersect
	int m = nodes.size();

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestors(n);
	LogicalVector l_background(n);
	for(int i = 0; i < background.size(); i ++) {
		l_background[background[i] - 1] = true;
	}
	if(type == SET_UNION) {
		for(int i = 0; i < m; i ++) {
			_find_ancestors_with_background(lt_parents, nodes[i] - 1, l_ancestors, l_background, include_self);
		}
	} else {
		LogicalVector l(n, true);
		LogicalVector l_an(n, false);
		for(int i = 0; i < m; i ++) {
			_find_ancestors_with_background(lt_parents, nodes[i] - 1, l_an, l_background, include_self);
			l = l & l_an;
			reset_logical_vector_to_false(l_an);
		}
		l_ancestors = l;
	}
	IntegerVector aid = _which(l_ancestors);
	if(aid.size() > 0) {
		aid = aid + 1;
	}
	return aid;
}

// [[Rcpp::export]]
IntegerVector cpp_ancestors_of_two_groups(S4 dag, IntegerVector nodes1, IntegerVector nodes2, int type, bool include_self = false) { // type 1: union, 2: intersect, 3: only in node1, 4: only in node2
	int m1 = nodes1.size();
	int m2 = nodes2.size();

	List lt_parents = dag.slot("lt_parents");
	int n = lt_parents.size();

	LogicalVector l_ancestors1(n);
	LogicalVector l_ancestors2(n);
	for(int i = 0; i < m1; i ++) {
		_find_ancestors(lt_parents, nodes1[i] - 1, l_ancestors1, include_self);
	}
	for(int i = 0; i < m2; i ++) {
		_find_ancestors(lt_parents, nodes2[i] - 1, l_ancestors2, include_self);
	}

	LogicalVector l_ancestors(n);
	if(type == SET_UNION) {
		l_ancestors = l_ancestors1 | l_ancestors2;
	} else if(type == SET_INTERSECT) {
		l_ancestors = l_ancestors1 & l_ancestors2;
	} else if(type == SET_UNIQU_IN_1) {
		l_ancestors = l_ancestors1 & (!l_ancestors2);
	} else if(type == SET_UNIQU_IN_2) {
		l_ancestors = (!l_ancestors1) & l_ancestors2;
	}

	IntegerVector aid = _which(l_ancestors);
	if(aid.size() > 0) {
		aid = aid + 1;
	}
	return aid;
}


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

// [[Rcpp::export]]
NumericVector cpp_offspring_aggregate(S4 dag, NumericVector value, int method = 1) {
	int n = dag.slot("n_terms");
	List lt_children = dag.slot("lt_children");

	NumericVector s(n);
	LogicalVector l_offspring(n);
	NumericVector v2;
	for(int i = 0; i < n; i ++) {
		if(i % 1000 == 0) {
			message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
			message("going through " + std::to_string(i) + " / " + std::to_string(n) + " terms ...", false);
		}

		_find_offspring(lt_children, i, l_offspring, true);
		v2 = value[l_offspring];
		if(method == 1) {  // mean
			s[i] = sum(v2)/sum(l_offspring);
		} else {   // sum
			s[i] = sum(v2);
		}

		reset_logical_vector_to_false(l_offspring);
	}

	message("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", false);
	message("going through " + std::to_string(n) + " / " + std::to_string(n) + " terms ... Done.", true);


	return s;
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

	LogicalVector l_ancestors(n_all);
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
					_find_ancestors(lt_parents, nodes[i] - 1, l_ancestors);
					// whether j_node is an ancestor of i_node
					for(int j = i+1; j < n; j ++) {
						if(l_ancestors[ nodes[j] - 1 ]) {
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
					_find_ancestors(lt_parents, nodes[i] - 1, l_ancestors);
					// whether j_node is an ancestor of i_node
					for(int j = i+1; j < n; j ++) {
						if(l_ancestors[ nodes[j] - 1 ]) {
							m(i, j) = true;
							if(!directed) {
								m(j, i) = true;
							}
							flag = true;
						}
					}
					reset_logical_vector_to_false(l_ancestors);
				}
			}
		}

		flag = false;
		reset_logical_vector_to_false(l_ancestors);
		reset_logical_vector_to_false(l_offspring);
	}

	return m;
}


const int TRAVERSE_UPSTREAM = -1;
const int TRAVERSE_DOWNSTREAM = 1;

const bool USE_MAX_DIST = true;
const bool USE_MIN_DIST = false;

// assume l_background include from_node
IntegerVector cpp_dag_traverse_bfs(S4 dag, IntegerVector from_node = IntegerVector(0), bool use_max_dist = USE_MAX_DIST, 
	LogicalVector l_background = LogicalVector(0), int direction = TRAVERSE_DOWNSTREAM) {

	// it calcualtes distances for all nodes

	String slot_name;
	if(direction == TRAVERSE_DOWNSTREAM) {
		slot_name = "lt_children";
	} else {
		slot_name = "lt_parents";
	}

	List lt_nodes = dag.slot(slot_name);

	if(from_node.size() == 0) {
		if(direction == TRAVERSE_DOWNSTREAM) {
			from_node = dag.slot("root");
		} else {
			from_node = dag.slot("leaves");
		}
	}
	IntegerVector i_from = from_node - 1;

	int n = lt_nodes.size();
	IntegerVector d(n, -1);

	bool has_background = false;
	if(l_background.size() > 0) {
		has_background = true;
	}

	LogicalVector l_current_nodes(n, false);
	for(int i = 0; i < i_from.size(); i ++) {
		l_current_nodes[ i_from[i] ] = true;
		d[ i_from[i] ] = 0;
	}
	
	int n_current_nodes = sum(l_current_nodes);

	while(n_current_nodes) {
		for(int i = 0; i < n; i ++) {
			if(l_current_nodes[i]) {
				int cr = i;

				l_current_nodes[i] = false;
				
				IntegerVector nodes = lt_nodes[cr];

				if(nodes.size()) {
					for(int j = 0; j < nodes.size(); j ++) {
						int i_node = nodes[j]-1;

						if( (has_background && l_background[i_node]) || !has_background ) {
							if(d[i_node] == -1) { // if the node has not been visited, the depth/height is the previous + 1
								d[i_node] = d[cr] + 1;
							} else {
								if(use_max_dist) {  // if it is visited, compare to current value
									if(d[cr] + 1 > d[i_node]) {
										d[i_node] = d[cr] + 1;
									}
								} else {
									if(d[cr] + 1 < d[i_node]) {
										d[i_node] = d[cr] + 1;
									}
								}
							}
							l_current_nodes[i_node] = true;
						}
					}
				}
			}
		}
		n_current_nodes = sum(l_current_nodes);
	}

	return d;
}

// [[Rcpp::export]]
IntegerVector cpp_dag_depth(S4 dag) {
	IntegerVector from_node(1);
	from_node[0] = dag.slot("root");
	return cpp_dag_traverse_bfs(dag, from_node, USE_MAX_DIST, LogicalVector(0), TRAVERSE_DOWNSTREAM);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_dist_from_root(S4 dag) {
	IntegerVector from_node(1);
	from_node[0] = dag.slot("root");
	return cpp_dag_traverse_bfs(dag, from_node, USE_MIN_DIST, LogicalVector(0), TRAVERSE_DOWNSTREAM);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_longest_dist_to_offspring(S4 dag, IntegerVector from_node, LogicalVector l_background = LogicalVector(0)) {
	// for every node in the DAG, maximal dist from nodes in `from_nodes`
	return cpp_dag_traverse_bfs(dag, from_node, USE_MAX_DIST, l_background, TRAVERSE_DOWNSTREAM);
}

IntegerVector cpp_dag_longest_dist_to_offspring(S4 dag, int from_node, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector from_node2(1);
	from_node2[0] = from_node;
	return cpp_dag_traverse_bfs(dag, from_node2, USE_MAX_DIST, l_background, TRAVERSE_DOWNSTREAM);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_shortest_dist_to_offspring(S4 dag, IntegerVector from_node, LogicalVector l_background = LogicalVector(0)) {
	return cpp_dag_traverse_bfs(dag, from_node, USE_MIN_DIST, l_background, TRAVERSE_DOWNSTREAM);
}

IntegerVector cpp_dag_shortest_dist_to_offspring(S4 dag, int from_node, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector from_node2(1);
	from_node2[0] = from_node;
	return cpp_dag_traverse_bfs(dag, from_node2, USE_MIN_DIST, l_background, TRAVERSE_DOWNSTREAM);
}


// [[Rcpp::export]]
IntegerVector cpp_dag_height(S4 dag) {
	IntegerVector to_node = dag.slot("leaves");
	return cpp_dag_traverse_bfs(dag, to_node, USE_MAX_DIST, LogicalVector(0), TRAVERSE_UPSTREAM);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_dist_to_leaves(S4 dag) {
	IntegerVector to_node = dag.slot("leaves");
	return cpp_dag_traverse_bfs(dag, to_node, USE_MIN_DIST, LogicalVector(0), TRAVERSE_UPSTREAM);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_longest_dist_from_ancestors(S4 dag, IntegerVector to_node, LogicalVector l_background = LogicalVector(0)) {
	return cpp_dag_traverse_bfs(dag, to_node, USE_MAX_DIST, l_background, TRAVERSE_UPSTREAM);
}

IntegerVector cpp_dag_longest_dist_from_ancestors(S4 dag, int to_node, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector to_node2(1);
	to_node2[0] = to_node;
	return cpp_dag_traverse_bfs(dag, to_node2, USE_MAX_DIST, l_background, TRAVERSE_UPSTREAM);
}

// [[Rcpp::export]]
IntegerVector cpp_dag_shortest_dist_from_ancestors(S4 dag, IntegerVector to_node, LogicalVector l_background = LogicalVector(0)) {
	return cpp_dag_traverse_bfs(dag, to_node, USE_MIN_DIST, l_background, TRAVERSE_UPSTREAM);
}

IntegerVector cpp_dag_shortest_dist_from_ancestors(S4 dag, int to_node, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector to_node2(1);
	to_node2[0] = to_node;
	return cpp_dag_traverse_bfs(dag, to_node2, USE_MIN_DIST, l_background, TRAVERSE_UPSTREAM);
}

// -----------------------------
// cyclic node
void _go_child(List lt_children, int node, IntegerVector path, CharacterVector terms, List& cyclic_paths) {
	int i_node = node - 1;
	if(path.size() > 0) {
		for(int i = 0; i < path.size(); i ++) {
			if(path[i] == node) {
				IntegerVector path2;
				for(int j = i; j < path.size(); j ++) {
					path2.push_back(path[j]);
				}
				path2.push_back(node);
				cyclic_paths.push_back(path2);

				if(cyclic_paths.size() > 1000) {
					stop("Too many cyclic paths (> 1000).");
				}

				return;
			}
		}
	}

	IntegerVector children = lt_children[i_node];
	for(int i = 0; i < children.size(); i ++) {
		IntegerVector path2 = clone(path);
		path2.push_back(node);
		_go_child(lt_children, children[i], path2, terms, cyclic_paths);
	}
}

// [[Rcpp::export]]
List cpp_check_cyclic_node(S4 dag, int node = -1) {
	List lt_children = dag.slot("lt_children");
	if(node == -1) {
		node = dag.slot("root");
	}
	CharacterVector terms = dag.slot("terms");

	List cyclic_paths;

	IntegerVector path(0);
	_go_child(lt_children, node, path, terms, cyclic_paths);

	return cyclic_paths;
}


// ---------------------------
List cpp_dag_traverse_bfs_sum_value(S4 dag, IntegerVector from_node, NumericVector value, bool use_max_dist = USE_MAX_DIST, 
	LogicalVector l_background = LogicalVector(0), int direction = TRAVERSE_DOWNSTREAM) {

	// it calcualtes distances for all nodes

	String slot_name;
	if(direction == TRAVERSE_DOWNSTREAM) {
		slot_name = "lt_children";
	} else {
		slot_name = "lt_parents";
	}

	List lt_nodes = dag.slot(slot_name);

	if(from_node.size() == 0) {
		if(direction == TRAVERSE_DOWNSTREAM) {
			from_node = dag.slot("root");
		} else {
			from_node = dag.slot("leaves");
		}
	}
	IntegerVector i_from = from_node - 1;

	int n = lt_nodes.size();
	IntegerVector d(n, -1);
	NumericVector v(n, 0.0);

	bool has_background = false;
	if(l_background.size() > 0) {
		has_background = true;
	}

	LogicalVector l_current_nodes(n, false);
	for(int i = 0; i < i_from.size(); i ++) {
		l_current_nodes[ i_from[i] ] = true;
		d[ i_from[i] ] = 0;
	}

	int n_current_nodes = sum(l_current_nodes);

	while(n_current_nodes) {
		for(int i = 0; i < n; i ++) {
			if(l_current_nodes[i]) {
				int cr = i;

				l_current_nodes[i] = false;
				
				IntegerVector nodes = lt_nodes[cr];

				if(nodes.size()) {
					for(int j = 0; j < nodes.size(); j ++) {
						int i_node = nodes[j]-1;

						if( (has_background && l_background[i_node]) || !has_background ) {
							if(d[i_node] == -1) { // if the node has not been visited, the depth/height is the previous + 1
								d[i_node] = d[cr] + 1;
								v[i_node] = v[cr] + value[i_node];
							} else {
								if(use_max_dist) {  // if it is visited, compare to current value
									if(d[cr] + 1 > d[i_node]) {
										d[i_node] = d[cr] + 1;
										v[i_node] = v[cr] + value[i_node];
									}
								} else {
									if(d[cr] + 1 < d[i_node]) {
										d[i_node] = d[cr] + 1;
										v[i_node] = v[cr] + value[i_node];
									}
								}
							}
							l_current_nodes[i_node] = true;
						}
					}
				}
			}
		}
		n_current_nodes = sum(l_current_nodes);
	}

	List lt = List::create(_["d"] = d, _["v"] = v);
	return lt;
}


List cpp_dag_longest_path_to_offspring_sum_value(S4 dag, IntegerVector from_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	return cpp_dag_traverse_bfs_sum_value(dag, from_node, value, USE_MAX_DIST, l_background, TRAVERSE_DOWNSTREAM);
}

List cpp_dag_longest_path_to_offspring_sum_value(S4 dag, int from_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector from_node2(1);
	from_node2[0] = from_node;
	return cpp_dag_traverse_bfs_sum_value(dag, from_node2, value, USE_MAX_DIST, l_background, TRAVERSE_DOWNSTREAM);
}

List cpp_dag_shortest_path_to_offspring_sum_value(S4 dag, IntegerVector from_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	return cpp_dag_traverse_bfs_sum_value(dag, from_node, value, USE_MIN_DIST, l_background, TRAVERSE_DOWNSTREAM);
}

// [[Rcpp::export]]
List cpp_dag_shortest_path_to_offspring_sum_value(S4 dag, int from_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector from_node2(1);
	from_node2[0] = from_node;
	return cpp_dag_traverse_bfs_sum_value(dag, from_node2, value, USE_MIN_DIST, l_background, TRAVERSE_DOWNSTREAM);
}

List cpp_dag_longest_path_from_ancestors_sum_value(S4 dag, IntegerVector to_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	return cpp_dag_traverse_bfs_sum_value(dag, to_node, value, USE_MAX_DIST, l_background, TRAVERSE_UPSTREAM);
}

List cpp_dag_longest_path_from_ancestors_sum_value(S4 dag, int to_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector to_node2(1);
	to_node2[0] = to_node;
	return cpp_dag_traverse_bfs_sum_value(dag, to_node, value, USE_MAX_DIST, l_background, TRAVERSE_UPSTREAM);
}

List cpp_dag_shortest_path_from_ancestors_sum_value(S4 dag, IntegerVector to_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	return cpp_dag_traverse_bfs_sum_value(dag, to_node, value, USE_MIN_DIST, l_background, TRAVERSE_UPSTREAM);
}
List cpp_dag_shortest_path_from_ancestors_sum_value(S4 dag, int to_node, NumericVector value, LogicalVector l_background = LogicalVector(0)) {
	IntegerVector to_node2(1);
	to_node2[0] = to_node;
	return cpp_dag_traverse_bfs_sum_value(dag, to_node2, value, USE_MIN_DIST, l_background, TRAVERSE_UPSTREAM);
}
