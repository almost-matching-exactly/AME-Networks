// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
void tester2(IntegerMatrix A) {
    // Load relevant packages
  	Environment igraph("package:igraph");

    // Load relevant functions from those packages

  	Function graph_from_adjacency = igraph["graph_from_adjacency_matrix"];
  	Function get_neighbors = igraph["neighbors"];
  	Function is_iso = igraph["is_isomorphic_to"];
    Function induced_subgraph = igraph["induced_subgraph"];
    SEXP g = graph_from_adjacency(Named("adjmatrix", A), Named("mode", "undirected"));
    int n = A.nrow();
    for (int i = 1; i <= n; ++i) {
      IntegerVector neighbors = get_neighbors(Named("graph", g), Named("v", i));
      SEXP neighborhood_sg = induced_subgraph(Named("graph", g), Named("vids", neighbors));
      Rcpp::as<bool>(is_iso(Named("graph1", g),
                                Named("graph2", neighborhood_sg)));
    }
  	// RObject g1 = graph_from_adjacency(Named("adjmatrix", A1), Named("mode", "undirected"));
  	// RObject g2 = graph_from_adjacency(Named("adjmatrix", A2), Named("mode", "undirected"));
  	// bool iso = is_iso(A1, A2);
  return;
}

// threshold_list_subgraphs = function(V, k) {
//   if (k == 'max') {
//     k <- min(length(V), k)
//   }
//   else {###### I know this is stupid.
//     k <- min(length(V), k)
//   }
//   sgs = c()
//     for (j in 1:k) {
//       if (j == 1 && length(V) == 1) {
//         sgs = c(sgs, list(V)) # Can probably simplify
//       } else if (j <= length(V)) {
//         sgs = c(sgs, combn(V, j, simplify = FALSE))
//       }
//     }
//     sgs
// }
// 
// threshold_all_neighborhood_subgraphs = function(G, k = 'max') {
//   sgs = list()
//   for (i in V(G)) {
//     neighbs <- neighbors(G, i)
//     if (length(neighbs) == 0)
//       sgs[[i]] = numeric(0)
//       else
//         sgs[[i]] = threshold_list_subgraphs(neighbs, k)
//   }
//   sgs
// }
Function my_combn("my_combn");

// [[Rcpp::export]]
List get_neighb_subgraphs(IntegerMatrix A, IntegerVector ids, int threshold = 0) {
  int n = ids.length();
  if (threshold == 0) {
    threshold = A.nrow(); 
  }
  List all_neighb_subgraphs (n);
  for (int i = 0; i < n; ++i) { // for all units
    IntegerVector adj_row = A(_, ids(i) - 1);
    IntegerVector neighbors; 
    for (int k = 0; k < A.nrow(); ++k) {
      if (adj_row[k] == 1) {
        neighbors.push_back(k + 1);
      }
    }
    int n_neighbors = neighbors.length();
    if (n_neighbors == 0) {
      all_neighb_subgraphs[i] = neighbors;
    }
    else {
      List neighb_subgraphs (pow(2, n_neighbors) - 1);
      int counter = 0;
      for (int j = 1; j <= n_neighbors && j <= threshold; ++j) {
        List all_combs = my_combn(Named("x", neighbors), Named("m", j));
        int n_combs = all_combs.length();
        for (int l = 0; l < n_combs; ++l) {
          neighb_subgraphs[counter] = all_combs[l];
          counter += 1; 
        }
        // neighb_subgraphs.push_back(my_combn(Named("x", neighbors), Named("m", j)));
      }
      all_neighb_subgraphs[i] = neighb_subgraphs;
    }
  }
  return all_neighb_subgraphs;
}


// Load relevant packages
Environment igraph("package:igraph");
// Environment utils("package:utils");
// Environment combinat("package:combinat")
// Load relevant functions from those packages
// Function my_combn("my_combn");
// Function combn = utils["combn"];
Function graph_from_adjacency = igraph["graph_from_adjacency_matrix"];
Function get_neighbors = igraph["neighbors"];
Function is_iso = igraph["is_isomorphic_to"];
Function induced_subgraph = igraph["induced_subgraph"];
Function is_graph = igraph["is_igraph"];

// [[Rcpp::export]]
List get_node_subgraph_counts(IntegerMatrix A, IntegerVector Z) {
//   // Load relevant packages
// 	Environment igraph("package:igraph");
//   // Environment utils("package:utils");
//   // Environment combinat("package:combinat")
//   // Load relevant functions from those packages
//   Function my_combn("my_combn");
//   // Function combn = utils["combn"];
// 	Function graph_from_adjacency = igraph["graph_from_adjacency_matrix"];
// 	Function get_neighbors = igraph["neighbors"];
// 	Function is_iso = igraph["is_isomorphic_to"];
//   Function induced_subgraph = igraph["induced_subgraph"];
//   Function is_graph = igraph["is_igraph"];
	// Get number of treated individuals
	// int n_treated = sum(Z);
	// int n_treated = std::accumulate(Z.begin(), Z.end(), 0);
	int n = A.nrow();

	// Get the graph (will need this for determining isomorphisms)
	SEXP g = graph_from_adjacency(Named("adjmatrix", A), Named("mode", "undirected"));
	if (Rcpp::as<bool>(is_graph(g)) == false) {
	  // cout << "blygh";
	}
	// how to set vertex attributes?
	List all_subgraph_types;
	int n_seen_subgraphs; // total # non-isomorphic subgraphs seen (length of above)
	int sg_num = -1; // default value for not being isormophic to a previously seen graph

	// What we'll be returning: a list with a vector for each unit
  // corresponding to its (neighborhood) subgraph counts
	List all_nodes_subgraph_counts;

	for (int i = 1; i <= n; ++i) { // for each unit
	  // cout << "break here 1";
	  // Where we'll store the unit's counts; at least as many as the types of graphs we've already seen
	  IntegerVector subgraph_counts (all_subgraph_types.length(), 0);
	  // Get neighbors
	  IntegerVector neighbors = get_neighbors(Named("graph", g), Named("v", i)); // can in theory just go from the adj mat
	  // cout << "break here 2";
	  // for each subgraph size we're considering (< n_treated bc otherwise won't be entirely treated graph)
	  // in the case where color matters, should have j <= min(n_treated, neighbors.length())
  	for (int j = 1; j <= neighbors.length(); ++j) { // Start at 1 because we don't care about 0 node subgraphs
  	  // cout << "break here 9 ";
  	  // cout << i << " ";
  	  // cout << j << " ";
  	  // cout << neighbors.length();
  	  // cout << neighbors; 
  	  // would be faster (and harder) computing combinations ourselves and then doing what follows one at a time
  	  List all_combs = my_combn(Named("x", neighbors), Named("m", j));
  	
	    // cout << "Unit: " << i << " \n";
	    // cout << "Neighbors: " << neighbors << " \n";
	    // cout << "Subgraph Size: " << j << " \n";
  	  
  	  // List all_combs = combn(Named("x", neighbors), Named("m", j), Named("simplify", false));
  	  // cout << "break here 10";
  	  int n_combs = all_combs.length(); // number of subgraphs with j many nodes
  	  for (int k = 0; k < n_combs; ++k) { // for each of these j-node subgraphs
  	    sg_num = -1; // resets 'not seen' flag
  	    // For now assuming coloring and the fact that we only want fully treated graphs
  	    // Ideally would initialize this beforehand but I don't know its length...
  	    // cout << "break here 3";
  	    IntegerVector current_subgraph = all_combs[k]; // nodes in current subgraph
  	    // cout << "break here 4";
  	    // not sure if reinitializing each time is super slow; check that.
  	    SEXP neighborhood_sg = induced_subgraph(Named("graph", g), Named("vids", current_subgraph));
  	    // cout << "break here 5";
  	    // check that subgraph indices don't start at 0
  	    // ok temporarily don't do this and assume it's all treated
  	    // if (sum(Z[current_subgraph - 1]) != current_subgraph.length()) { // If not every node treated. -1 bc 0 based indexing.
  	    //   continue; // move on to the next j-node subgraph
  	    // }
  	    n_seen_subgraphs = all_subgraph_types.length();
  	    // if (i == 2 && j == 2) {
  	    //   cout << n_seen_subgraphs;
  	    // }
  	    if (n_seen_subgraphs == 0) { // Haven't seen any graphs at all
  	      subgraph_counts.push_back(1); // Start the count at 1 for first subgraph ever seen.
  	      all_subgraph_types.push_back(neighborhood_sg); // Add this first one to the list
  	    }
  	    else { // We've seen at least 1 graph before
  	      // cout << n_seen_subgraphs;
  	      for (int l = 0; l < n_seen_subgraphs; ++l) { // Loop through previously seen subgraphs
  	        // cout << "break here 6 ";
  	        if (Rcpp::as<bool>(is_iso(Named("graph1", all_subgraph_types[l]),
                       Named("graph2", neighborhood_sg)))) { // take coloring into account!!
              sg_num = l;
  	          break;
  	        }
  	        // cout << l;
  	      }
  	      // cout << "flag ";
  	      if (sg_num == -1) { // haven't seen this graph before
  	        subgraph_counts.push_back(1); // Start my count at 1
            all_subgraph_types.push_back(neighborhood_sg); // Add it to the list
  	      }
          else {
            subgraph_counts[sg_num] += 1;
          }
  	    }
  	  }
  	}
  	all_nodes_subgraph_counts.push_back(subgraph_counts);
	}
	return all_nodes_subgraph_counts;
}
/***R
my_combn <- function(x, m) {
  if (length(x) == 1) {
    return(list(x))
  }
  return(utils::combn(as.integer(x), m, simplify = FALSE))
}
A <- matrix(c(0, 0, 0, 0, 0,
              0, 0, 0, 1, 1,
              0, 0, 0, 1, 0,
              0, 1, 1, 0, 1,
              0, 1, 0, 1, 0),
            nrow = 5)
out <- get_neighb_subgraphs(A, c(1, 4))
# G <- erdos.renyi.game(50, 0.08)
# A <- get.adjacency(G, type = 'both', sparse = FALSE)
# for (i in 1:5) {
#   n <- dim(A)[1]
#   Z <- rep(1, n)
#   dta <- get_node_subgraph_counts(A, Z)
# }
*/
