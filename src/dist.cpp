#include <Rcpp.h>
using namespace Rcpp;

#include <limits.h>

const int USE_LONGEST_DISTANCE = 1;
const int USE_SHORTEST_DISTANCE = 2;


int cpp_tpl_path_length(S4 dag, int from, int to, int type = 2) {
	
	if(from == to) {
		return 0;
	}

	IntegerVector tpl_sorted = dag.slot("tpl_sorted");
	IntegerVector tpl_pos = dag.slot("tpl_pos"); 
	List lt_children = dag.slot("lt_children");
	int n = lt_children.size();

	if(from > n || from < 1) {
		stop("'from' node is not in the DAG.");
	}
	if(to > n || to < 1) {
		stop("'to' node is not in the DAG.");
	}

	int i_from = from - 1;
	int i_to = to - 1;

	int default_dist;
	if(type == USE_LONGEST_DISTANCE) {
		default_dist = INT_MIN;
	} else {
		default_dist = INT_MAX - 1;
	}

	if(tpl_pos[i_from] > tpl_pos[i_to]) {
		return -1;
	}

	int i_pos_from = tpl_pos[i_from] - 1;
	int i_pos_to = tpl_pos[i_to] - 1;

	int len = i_pos_to - i_pos_from + 1;
	
	IntegerVector dist(len, default_dist);
	dist[0] = 0; // from -> from

	for(int i_pos = i_pos_from; i_pos <= i_pos_to; i_pos ++) {

		int i_node = tpl_sorted[i_pos] - 1;  // the current node on the sorted list

		IntegerVector children = lt_children[i_node];
		
		for(int j = 0; j < children.size(); j ++) {
			int i_child = children[j] - 1;
			int i_pos_child = tpl_pos[i_child] - 1; // position of the child on the sorted list
			if(i_pos_child <= i_pos_to) {
				
				if(type == USE_LONGEST_DISTANCE) { // max
					if(dist[i_pos_child - i_pos_from] < dist[i_pos - i_pos_from] + 1) {
						dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
					}
				} else {
					if(dist[i_pos_child - i_pos_from] > dist[i_pos - i_pos_from] + 1) {
						dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
					}
				}
				
			}
		}
	}

	if(dist[len-1] == default_dist || dist[len-1] == default_dist + 1) {
		dist[len-1] = -1;
	}
	
	return dist[len-1];
}


IntegerVector cpp_tpl_find_path(S4 dag, int from, int to, int type = 2) {
	
	if(from == to) {
		IntegerVector path(1, from);
		return path;
	}

	IntegerVector tpl_sorted = dag.slot("tpl_sorted");
	IntegerVector tpl_pos = dag.slot("tpl_pos"); 
	List lt_children = dag.slot("lt_children");
	List lt_children_relations = dag.slot("lt_children_relations");
	int n = lt_children.size();

	if(from > n || from < 1) {
		stop("'from' node is not in the DAG.");
	}
	if(to > n || to < 1) {
		stop("'to' node is not in the DAG.");
	}

	int i_from = from - 1;
	int i_to = to - 1;

	if(tpl_pos[i_from] > tpl_pos[i_to]) {
		IntegerVector path(0);
		return path;
	}

	int i_pos_from = tpl_pos[i_from] - 1;
	int i_pos_to = tpl_pos[i_to] - 1;

	int len = i_pos_to - i_pos_from + 1; // include both from and to
	int default_dist;
	if(type == USE_LONGEST_DISTANCE) {
		default_dist = INT_MIN;
	} else {
		default_dist = INT_MAX - 1;
	}
	
	IntegerVector dist(len, default_dist);
	dist[0] = 0; // from -> from

	IntegerVector predecessor(len, -1);

	for(int i_pos = i_pos_from; i_pos <= i_pos_to; i_pos ++) {

		int i_node = tpl_sorted[i_pos] - 1;  // the current node on the sorted list

		IntegerVector children = lt_children[i_node];
		
		for(int j = 0; j < children.size(); j ++) {
			int i_child = children[j] - 1;
			int i_pos_child = tpl_pos[i_child] - 1; // position of the child on the sorted list
			if(i_pos_child <= i_pos_to) {
				
				if(type == USE_LONGEST_DISTANCE) {
					if(dist[i_pos_child - i_pos_from] < dist[i_pos - i_pos_from] + 1) {
						dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
						predecessor[i_pos_child - i_pos_from] = i_pos - i_pos_from;
					}
				} else {
					if(dist[i_pos_child - i_pos_from] > dist[i_pos - i_pos_from] + 1) {
						dist[i_pos_child - i_pos_from] = dist[i_pos - i_pos_from] + 1;
						predecessor[i_pos_child - i_pos_from] = i_pos - i_pos_from;
					}
				}
			}
		}
		
	}

	IntegerVector path(len, -1);
	path[len-1] = to;

	int i = len - 1;
	int i_previous;
	while(1) {
		i_previous = predecessor[i];
		path[i-1] = tpl_sorted[i_previous + i_pos_from];
		i = i_previous;

		if(i_previous < 0) {
			break;
		}
	}

	path = path[path > 0];
	
	return path;
}

// [[Rcpp::export]]
int cpp_tpl_shortest_path_length(S4 dag, int from, int to) {
	return cpp_tpl_path_length(dag, from, to, USE_SHORTEST_DISTANCE);
}

// [[Rcpp::export]]
int cpp_tpl_longest_path_length(S4 dag, int from, int to) {
	return cpp_tpl_path_length(dag, from, to, USE_LONGEST_DISTANCE);
}

// [[Rcpp::export]]
IntegerVector cpp_tpl_shortest_path(S4 dag, int from, int to) {
	return cpp_tpl_find_path(dag, from, to, USE_SHORTEST_DISTANCE);
}

// [[Rcpp::export]]
IntegerVector cpp_tpl_longest_path(S4 dag, int from, int to) {
	return cpp_tpl_find_path(dag, from, to, USE_LONGEST_DISTANCE);
}


