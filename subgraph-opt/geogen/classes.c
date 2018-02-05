#include "classes.h"
#include <cfloat>

// BOX CONSTRUCTOR
Box::Box(list< Point >* ptrcontained = new list< Point > , int d=2) {
  contained = ptrcontained; // CONTAINED INITIALIZED
  numpoints = contained->size(); // NUMPOINTS INITIALIZED
  dim = d; //DIM INITIALIZED
  estimate = DBL_MAX; // ESTIMATE INITIALIZED
  // NEIGHBORS, ATTRACTORS AND SUCCESSORs ARE EMPTY TO START WITH
  // SHOULD SOMEBODY CHECK ALL POINTS HAVE SAME DIM?
  

  shrink(); // WILL INITIALIZE 
            // SIDELENGTH AND CENTER
  
} 

//BOX COPY CONSTRUCTOR - NOT NECESSARY SO FAR

// BOX DESTRUCTOR
Box::~Box(){
  delete contained;
}


//SHRINK FUNCTION
void Box::shrink() {
  double max[dim];
  double min[dim];
  
  for(int i=0; i < dim; i++) {
    max[i] = DBL_MIN;
    min[i] = DBL_MAX;
  }
  
  
  for(list< Point >::iterator it = contained->begin();
      it != contained->end(); ++it) 
    for(int i=0; i < dim ; i++) {
      if ((it->coordinates)[i] > max[i])
	max[i] = (it->coordinates)[i];
      if ((it->coordinates)[i] <  min[i])
	min[i] = (it->coordinates)[i];
    }
  double maxlength = 0;
  for(int i=0; i < dim; i++) {
    if(max[i] - min[i] > maxlength)
      maxlength = max[i] - min[i];
  }
  
  center.coordinates = new double[dim];

  for(int i=0; i < dim; i++) {
    center.coordinates[i] = max[i] - maxlength/2;
  }
  
  sidelength = maxlength;
}
    
 

// DISTANCE FUNCTIONS
double squaredmax(Box& b) {
  return (b.sidelength*b.sidelength)*b.dim;
}
 
double squaredmax(Box& b, Box& c) {
  if(int dim = b.dim != c.dim)
    return -1;

  double value = 0;
  for(int i=0; i < b.dim; i++) 
    value += b.center.coordinates[i]*b.center.coordinates[i] + c.center.coordinates[i]*c.center.coordinates[i] + (b.sidelength + c.sidelength)*(b.sidelength + c.sidelength);
  
}

// PRINTFUNCTIONS
void Box::simpleprintvertices() {
  int k, j, temp;
  int end = 1;
  
  for(j=0; j < dim; j++) 
    end = 2 * end;


  for(j=0; j < end; j++) {
    temp = j;
    for(k=0; k < dim; k++) {
      if(temp%2)
	printf("%f\t", center.coordinates[k] + sidelength/2);
      else
	printf("%f\t", center.coordinates[k] - sidelength/2);
      
      temp=temp/2;
    }
    printf("\n");
  }
       
}


void Box::split() {
  


double Box::volume() {
  return sidelength;
  //return  pow(sidelength, dim);
}
