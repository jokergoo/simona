#ifndef __TERMS_
#define __TERMS__

IntegerVector cpp_n_offspring(S4 dag);
IntegerVector cpp_dag_depth_bfs(S4 dag, int from_node = 0, bool use_max = true, LogicalVector l_background = LogicalVector(0));

#endif

