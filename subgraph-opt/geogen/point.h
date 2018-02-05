#include <stdlib.h>
#include <stdio.h>
using namespace std;

class Point {
 public:
  double* coordinates;
  int dim;

  Point(){};

  Point(int k = 2) {
    dim = k;
    coordinates = new double[k];
  }

  Point(const Point& t) {
    dim = t.dim;
    coordinates = new double[dim];
    for(int i=0; i < dim; i++) 
      coordinates[i] = t.coordinates[i];
  }

  void simpleprint () {
    for(int i=0; i < dim; i++) 
      printf("%f\t", coordinates[i]);
    printf("\n");
  }

  ~Point() {
    delete coordinates;
  }
};


