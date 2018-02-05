#include <stdlib.h>
#include <list>
#include<functional>
#include <cmath>
#include "point.h"

class Box {
 public:
  Point center;
  double sidelength;
  double estimate;
  int numpoints;
  int dim;

  list< Box* > neighbors;   // container neighbors
  list< Box* > attractors;  //container Attractors
  list< Box* > successors; //container successors
  list< Point >* contained; //points contained - CHANGE TO LIST OF POINTERS, MODIFY GENERATION

  Box(list< Point >* ptrcontained, int);
  ~Box();

  double volume();

  void split();
  void simpleprintvertices();
 
 private:
  void shrink();
}
;


struct less_volume : public binary_function<Box*, Box*, bool> {
  bool operator() (Box* a, Box* b) { return a->volume() < b->volume();}
};
