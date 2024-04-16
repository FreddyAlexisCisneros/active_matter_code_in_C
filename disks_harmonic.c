// need
// to
// edit
// this 
// code
// and
// import more parts of the original code.



#include <stdlib.h>
#include <stdio.h>
#include <time.h> 
#include <math.h>
#include <complex.h>
#include <string.h>
#include "mt19937-64.c"
#define PI (3.141592653589793)

typedef struct {
 double x,y;
 double theta;
 double new_theta;
 double Fx,Fy,E;
 double dx,dy;
} particle;
