#include <ANN/ANN.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/graph/incremental_components.hpp>
#include <boost/graph/edge_list.hpp>
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <sstream>
#include <ext/hash_map>
#include <algorithm>

using namespace boost;

double** pointgenerator(int d, int n,  unsigned int rseed);

typedef std::pair<int, int> Edge;

struct isLonger{
    bool operator() (std::pair<Edge, double> a, std::pair<Edge, double> b) {return (a.second > b.second);}
  };



int main(int argc, char* argv[]) {

  /* READING ARGS & ERROR CHECKING */
  if( argc != 5) {
    std::cerr << "Incorrect number of arguments. Usage " << argv[0] << " num_dimensions num_points num_neighbors seed" << std::endl; 
    return 1;
  }

  int d, n, k;
  unsigned int s;
  int i, j;

  std::istringstream(argv[1]) >> d;
  std::istringstream(argv[2]) >> n;
  std::istringstream(argv[3]) >> k;
  std::istringstream(argv[4]) >> s;
  k = k+1;

  if( n <= 0 || d <= 0 ) {
    std::cerr << "Incorrect arguments.\n";
    return 1;
  }


  /* POINT GENERATION */
  double** points = pointgenerator(d, n, s);
  
  /*
  for(i=0; i < n; i++) {
    for(j=0; j < d; j++) 
      printf("%lf\t", points[i][j]);
    printf("\n");
  }
  */

  ANNpointArray pointlist = (ANNpointArray) points;
  
  ANNkd_tree tree(pointlist, n, d);

  ANNidxArray indices = new ANNidx[k];
  ANNdistArray distances = new ANNdist[k];


  /* DEFINE EDGES AND EDGELISTS */
  typedef std::vector<std::pair<Edge, double> > EdgeList;



  /* CREATE EDGE LIST AND LENGTH VECTOR*/
  EdgeList AllowedEdges(n*k);


  
  int counter = 0;
  /* PREPARE EDGES */
  for(i=0; i < n; i++) {
    tree.annkSearch(pointlist[i], k, indices, distances);
    for(j=0; j < k; j++) {
      if (indices[j] != i) {
	Edge e= Edge(i, indices[j]);
	AllowedEdges[counter].first = e;
	AllowedEdges[counter].second = distances[j];
	counter++;
      }
    }				         
  }

  AllowedEdges.resize(counter);


  /* CREATE HEAP OF EDGES */
  
 

  std::make_heap(AllowedEdges.begin(), AllowedEdges.end(), isLonger());


  /* INITALIZE INCREMENTAL COMPONENTS */

  typedef adjacency_list <hash_setS, vecS, undirectedS> Graph;
  typedef graph_traits<Graph>::vertex_descriptor Vertex;
  typedef graph_traits<Graph>::vertices_size_type size_type;
  

  Graph G(n);

  std::vector<size_type> rank(num_vertices(G));
  std::vector<Vertex> parent(num_vertices(G));
  typedef size_type* Rank;
  typedef Vertex* Parent;
  disjoint_sets<Rank, Parent>  ds(&rank[0], &parent[0]);

  initialize_incremental_components(G, ds);
  incremental_components(G, ds);

  graph_traits<Graph>::edge_descriptor edge;
  bool flag;

  typedef graph_traits<Graph>::vertex_iterator vertex_iter;
  std::pair<vertex_iter, vertex_iter> vp;
  vp = vertices(G); 



  /* POP SMALLER EDGE, ADD TO GRAPH UNTIL CONNECTED */

  for(std::vector<std::pair<Edge, double> >::iterator iter = AllowedEdges.end(); iter > AllowedEdges.begin(); iter--) { 
    int num  = ds.count_sets(vp.first, vp.second);
    if(num < 2)
      break;
    

    std::pop_heap(AllowedEdges.begin(), iter, isLonger());
    Edge e = (iter-1)->first;
   

    boost::tie(edge,flag) = add_edge(e.first, e.second, G);
    if(flag == true) {
      ds.union_set(e.first, e.second);
    }

    /*
    printf("%d=(%lf, %lf)\t%d=(%lf, %lf)\t%lf\n", 
	   e.first, points[e.first][0], points[e.first][1], 
	   e.second, points[e.second][0], points[e.second][1],
	   (iter-1)->second);
    */

  }
  typedef property_map<Graph, vertex_index_t>::type IndexMap;
  IndexMap index = get(vertex_index, G);


  // print parameters first
  std::cout << n << "\t" << d << "\t" << num_edges(G) << std::endl;
  
  // print list of node coordinates
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < d; j++) {
      std::cout << points[i][j] << "\t";
    }
    std::cout << std::endl;
  }
  
  // print edge list
  graph_traits<Graph>::edge_iterator ei, ei_end;
  for (tie(ei, ei_end) = edges(G); ei != ei_end; ++ei) {
    if (index[source(*ei, G)] < index[target(*ei, G)])
      std::cout << index[source(*ei, G)] << "\t" << index[target(*ei, G)] << std::endl;
    else
      std::cout << index[target(*ei, G)] << "\t" << index[source(*ei, G)] << std::endl;
  }

  return 0;
}


// GENERATES n POINTS UNIFORMLY IN THE UNIT SPHERE IN d-dimensions
double** pointgenerator(int d, int n, unsigned int rseed) {
  double**  pointlist;

  pointlist = new double* [n];

  double pi = atan(1.0f) * 4.0f;

  srand(rseed);

  double random, temp;

  for(int i=0; i < n; i++) {    
    /*random = ((double) rand())/RAND_MAX;
    temp = pow(random, 1.0/((double) d));
    pointlist[i] = new double[d];

    for(int k=0; k < d-1; k++) {      
      random  =((double) rand())/RAND_MAX*2*pi;
      pointlist[i][k] = temp * sin(random);
      temp *= cos(random);
    }
    pointlist[i][d-1] = temp;*/

    pointlist[i] = new double[d];

    for(int k=0; k < d; k++) {      
      double magn = ((double) rand())/RAND_MAX;
      double sign = ((double) rand())/RAND_MAX > 0.5 ? 1.0 : -1.0;
      pointlist[i][k] = sign * magn;
    }
  }

  return pointlist;
}
