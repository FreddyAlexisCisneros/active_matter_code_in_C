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
typedef struct elem elem;
struct elem
{
 int num;
 struct elem *next;
};  
typedef elem* site;
typedef struct elem_n elem_n;
struct elem_n
{
 int x;
 int y;
 struct elem_n *next;
};

typedef struct { 
 int n;
} Cells;
