#include <stdio.h>
#include <math.h>

// -------- set values --------
#define L 1.0
#define N_SIDE 100
#define N_NODE_INIT ((N_SIDE+1)*(N_SIDE+1)-(N_SIDE-1)*(N_SIDE-1))
#define N_EDGE_INIT 0
double boundary = 1.0;
int nums_node_init[N_NODE_INIT];//    = { 0, 1, 2, 3, 4, 7, 8, 11, 12, 13, 14, 15 }; 
double phis_node_init[N_NODE_INIT];// = { 5, 5, 5, 5, 0, 0, 0,  0,  0,  0,  0,  0 };

int nums_node_from[N_EDGE_INIT];
int nums_node_to[N_EDGE_INIT];
double qs[N_EDGE_INIT];

#define N_SERIES 100 // value for analytical solution plotting

// ---------------

#define l (L/N_SIDE)
#define N_node_side (N_SIDE+1)
#define N_node (N_node_side*N_node_side)
#define N_unknown_node (N_node-N_NODE_INIT)
double Q[N_node][N_node];


typedef struct {
    double x;
    double y;
} Node;

Node init_Node(double x, double y) {
    Node node;
    node.x = x;
    node.y = y;
    return node;
}

// ESM is the abbreviation of ElementStiffnessMatrix
typedef struct {
  double array[3][3];
} ESM;

ESM init_ESM(Node node1, Node node2, Node node3) {
  ESM esm;
  double x1,x2,x3,y1,y2,y3;
  double a1,a2,a3,b1,b2,b3,c1,c2,c3;
  double A; // area of a target element
  
  // initialize
  x1 = node1.x;
  y1 = node1.y;
  x2 = node2.x;
  y2 = node2.y;
  x3 = node3.x;
  y3 = node3.y;

  A = (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2;

  // calculate a,b,c
  a1 = (x2 * y3 - x3 * y2) / (2 * A);
  b1 = (y2 - y3) / (2 * A);
  c1 = (x3 - x2) / (2 * A);
  a2 = (x3 * y1 - x1 * y3) / (2 * A);
  b2 = (y3 - y1) / (2 * A);
  c2 = (x1 - x3) / (2 * A);
  a3 = (x1 * y2 - x2 * y1) / (2 * A);
  b3 = (y1 - y2) / (2 * A);
  c3 = (x2 - x1) / (2 * A);

  // initialize elemetns of array
  esm.array[0][0] = b1 * b1 + c1 * c1;
  esm.array[0][1] = esm.array[1][0] = b1 * b2 + c1 * c2;
  esm.array[0][2] = esm.array[2][0] = b1 * b3 + c1 * c3;
  esm.array[1][1] = b2 * b2 + c2 * c2;
  esm.array[1][2] = esm.array[2][1] = b2 * b3 + c2 * c3;
  esm.array[2][2] = b3 * b3 + c3 * c3;
  
  return esm;
}

// eesm is the abbreviation of ExtentedElemetsSNiffnessMatrix
typedef struct {
  double array[N_node][N_node];
} EESM;

EESM init_EESM(ESM esm, int nums_node[3]) {
  EESM eesm;
  
  // initialize eesm as zero elemetns
  for(int i=0; i<N_node; i++) {
    for(int j=0; j<N_node; j++) {
      eesm.array[i][j] = 0.0;
    }
  }

  // extend esm to eesm
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      eesm.array[nums_node[i]][nums_node[j]] = esm.array[i][j];
    }
  }

  return eesm;
}

// osm is the abbreviation of OverallSniffnessMatrix
typedef struct {
  double array[N_node][N_node];
} OSM;

OSM init_OSM(void) {
  OSM osm;
  for(int i=0; i<N_node; i++) {
    for(int j=0; j<N_node; j++) {
      osm.array[i][j] = 0.0;
    }
  }

  return osm;
}

// sosm is the abbreviation of SmallOverallSniffnessMatrix
typedef struct {
  double array[N_unknown_node][N_unknown_node];
} SOSM;

SOSM init_SOSM(OSM osm, int nums_unknown_node[N_unknown_node]) {
  SOSM sosm;

 for(int i=0; i<N_unknown_node; i++) {
    for(int j=0; j<N_unknown_node; j++) {
      sosm.array[i][j] = osm.array[nums_unknown_node[i]][nums_unknown_node[j]];
    }
  }   

  return sosm;
}

typedef struct {
  double array[N_unknown_node][N_unknown_node];
} INV_SOSM;

// Gaussian elimination
INV_SOSM init_INV_SOSM(SOSM sosm) {
  INV_SOSM inv_sosm;

  // initialize inv_sosm as identity matrix
  for(int i=0; i<N_unknown_node; i++) {
    for(int j=0; j<N_unknown_node; j++) {
      inv_sosm.array[i][j] = (i==j);
    }
  }

  for(int i=0; i<N_unknown_node; i++) {
    double tmp_normalize = 1.0 / sosm.array[i][i];
    // normalize elements along column based on diagonal element
    for(int j=0; j<N_unknown_node; j++) {
      sosm.array[i][j] *= tmp_normalize;
      inv_sosm.array[i][j] *= tmp_normalize;
    }
  // convert into zero except diagonal element
  for(int m=0; m<N_unknown_node; m++) {
    if(m!=i) {
      double tmp_toZero = sosm.array[m][i];
      for (int n=0; n<N_unknown_node; n++) {
        sosm.array[m][n] -= sosm.array[i][n] * tmp_toZero;
        inv_sosm.array[m][n] -= inv_sosm.array[i][n] * tmp_toZero;
      }
    }
  }
  }
  return inv_sosm;
}

// EF is the abbreviation of ExtendedF
typedef struct {
  double array[N_node];
} EF;

EF init_EF(int nums_node[3]) {
  EF ef;
  
  for(int i=0; i<N_node; i++) {
    ef.array[i] = 0.0;
  }

  double q1,q2,q3,L1,L2,L3;
  L1 = L2 = l;
  L3 = sqrt(2) * l;
  
  q1 = Q[nums_node[0]][nums_node[1]];
  q2 = Q[nums_node[1]][nums_node[2]];
  q3 = Q[nums_node[2]][nums_node[0]];

  ef.array[nums_node[0]] = (-1 / 2) * (q1 * L1 + q3 * L3);
  ef.array[nums_node[1]] = (-1 / 2) * (q1 * L1 + q2 * L2);
  ef.array[nums_node[2]] = (-1 / 2) * (q2 * L2 + q3 * L3);

  return ef;
}

// OF is the abbreviation of OverallF
typedef struct {
  double array[N_node];
} OF;

OF init_OF(void) {
  OF of;
  for(int i=0; i<N_node; i++) {
    of.array[i] = 0.0;
  }

  return of; 
}

OF update_OF(OF of, OSM osm) {
  for(int i=0; i<N_NODE_INIT; i++) {
    for(int j=0; j<N_node; j++) {
      of.array[j] -= phis_node_init[i] * osm.array[j][nums_node_init[i]];
    }
  }

  return of;
}

// SOF is the abbreviation of SmallOverallF
typedef struct {
  double array[N_unknown_node];
} SOF;

SOF init_SOF(OF of, int nums_unknown_node[N_unknown_node]) {
  SOF sof;

  for(int i=0; i<N_unknown_node; i++) {
    sof.array[i] = of.array[nums_unknown_node[i]];
  }  

  return sof; 
}


int main(void) {
  // initialize boundaries
  int p = 0;
  for(int i=0; i<N_node; i++) {
    if(i%N_node_side==(N_node_side-1)) {
      nums_node_init[p] = i;
      phis_node_init[p] = boundary;
      p++;
    }
    else if((i>=0 && i<N_node_side-1) || (i%N_node_side==0) || (i>(N_node-1-N_node_side) && i<N_node-1)) {
      nums_node_init[p] = i;
      phis_node_init[p] = 0;
      p++;
    }
  }

  // sorting
  for (int i=0; i<N_NODE_INIT; i++) {
    for (int j=i+1; j<N_NODE_INIT; j++) {
      if (nums_node_init[i] > nums_node_init[j]) {
        double tmp_nums =  nums_node_init[i];
        nums_node_init[i] = nums_node_init[j];
        nums_node_init[j] = tmp_nums;
        double tmp_phis =  phis_node_init[i];
        phis_node_init[i] = phis_node_init[j];
        phis_node_init[j] = tmp_phis;
      }
    }
  }

  //initialize Q
  for(int i=0; i<N_node; i++) {
    for(int j=0; j<N_node; j++) {
      Q[i][j] = 0.0;
    }
  }

  if(N_EDGE_INIT>0) {
    for(int i=0; i<N_EDGE_INIT; i++) {
      Q[nums_node_from[i]][nums_node_from[i]] = qs[i];
    }
  }
  // initialize nodes
  Node nodes[N_node_side][N_node_side];
  for(int i=0; i<N_node_side; i++) {
    for(int j=0; j<N_node_side; j++) {
      nodes[i][j] = init_Node(l*i, l*j);
    }
  }

  // start calculation
  OSM osm = init_OSM();
  OF of = init_OF();
  // uppser triangle
  for(int i=0; i<N_node_side-1; i++) {
    for(int j=0; j<N_node_side-1; j++) {
      ESM esm = init_ESM(nodes[i][j], nodes[i+1][j], nodes[i+1][j+1]);
      int num_node1 = N_node_side * i + j;
      int num_node2 = N_node_side * (i+1) + j;
      int num_node3 = N_node_side * (i+1) + (j+1);
      int nums_node[3] = { num_node1, num_node2, num_node3 };
      EESM eesm = init_EESM(esm, nums_node);
      EF ef = init_EF(nums_node);
      
      // update osm
      for(int m=0; m<N_node; m++) {
        for(int n=0; n<N_node; n++) {
          osm.array[m][n] += eesm.array[m][n];
        }
      }
      // update of
      for(int m=0; m<N_node; m++) {
        of.array[m] += ef.array[m];
      }
    }
  }

  // lower triangle
  for(int i=1; i<N_node_side; i++) {
    for(int j=1; j<N_node_side; j++) {
      ESM esm = init_ESM(nodes[i][j], nodes[i-1][j], nodes[i-1][j-1]);
      int num_node1 = N_node_side * i + j;
      int num_node2 = N_node_side * (i-1) + j;
      int num_node3 = N_node_side * (i-1) + (j-1);
      int nums_node[3] = { num_node1, num_node2, num_node3 };
      EESM eesm = init_EESM(esm, nums_node);
      EF ef = init_EF(nums_node);
      
      // update osm & of
      for(int m=0; m<N_node; m++) {
        for(int n=0; n<N_node; n++) {
          osm.array[m][n] += eesm.array[m][n];
        }
      }
      // update of
      for(int m=0; m<N_node; m++) {
        of.array[m] += ef.array[m];
      }
    }
  }

  // extract numbers of unknown nodes
  int nums_unknown_node[N_unknown_node];
  int m = 0;
  int n = 0;
  for(int i=0; i<N_node; i++) {
    if(i==nums_node_init[m]) {
      m++;
    } else {
      nums_unknown_node[n] = i;
      n++;
    }
  }
  
  // apply initial node conditions
  of = update_OF(of, osm);
  SOSM sosm = init_SOSM(osm, nums_unknown_node);
  INV_SOSM inv_sosm = init_INV_SOSM(sosm);
  SOF sof = init_SOF(of, nums_unknown_node);

  // output result
  double phis_unknown_node[N_unknown_node];
  double sum_element;
  for(int i=0; i<N_unknown_node; i++) {
    sum_element = 0;
    for(int j=0; j<N_unknown_node; j++) {
      sum_element += inv_sosm.array[i][j] * sof.array[j];
    }
    phis_unknown_node[i] = sum_element;
  }
  
  // write output on dat file
  double phis_node[N_node];
  m = 0;
  n = 0;
  for(int i=0; i<N_node; i++) {
    if(i==nums_node_init[m]) {
      phis_node[i] = phis_node_init[m];
      m++;
    }else {
      phis_node[i] = phis_unknown_node[n];
      n++;
    }
  }

  FILE *file;
  file = fopen("output/output.dat","w");
  m = 0;
  for(int i=0; i<N_node_side; i++) {
    for(int j=0; j<N_node_side; j++) {
      fprintf(file, "%f %f %f\n", nodes[i][j].x, nodes[i][j].y, phis_node[m]); 
      m++;
    }
    fprintf(file, "\n");
  }


  fclose(file);

  // plot by gnuplot
  FILE *gp;
  gp = popen("gnuplot -persist", "w");
  fprintf(gp, "u(x,y,n) = sum[k=1:n] u0*2.0*(1-(-1)**k)*sin(k*pi*x/L)*(exp(k*pi*y/L)-exp(-k*pi*y/L))/(k*pi*(exp(k*pi)-exp(-k*pi)))\n");
  fprintf(gp, "u0 = %f\n", boundary);
  fprintf(gp, "L = %f\n", L);
  fprintf(gp, "set xrange [0:%f]\n", L);
  fprintf(gp, "set yrange [0:%f]\n", L);
  //fprintf(gp, "splot u(x,y,%d), \"output/output.dat\" with lines\n", N_SERIES);
  fprintf(gp, "splot \"output/output.dat\" with lines\n", N_SERIES);
  pclose(gp);

 // for(int i=0; i<N_node; i++) {
 //   printf("%f\n", phis_node[i]);
 // }

  return 0;
}
