#include <stdio.h>

// set values
#define L 10.0
#define N_SIDE 2

#define l (L/N_SIDE)
#define N_node_side (N_SIDE+1)
#define N_node (N_node_side*N_node_side)

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

int main(void) {
  // initialize nodes
  Node nodes[N_node_side][N_node_side];
  for(int i=0; i<N_node_side; i++) {
    for(int j=0; j<N_node_side; j++) {
      nodes[i][j] = init_Node(l*i, l*j);
    }
  }

  OSM osm = init_OSM();
  // uppser triangle
  for(int i=0; i<N_node_side-1; i++) {
    for(int j=0; j<N_node_side-1; j++) {
      ESM esm = init_ESM(nodes[i][j], nodes[i+1][j], nodes[i+1][j+1]);
      int num_node1 = N_node_side * i + j;
      int num_node2 = N_node_side * (i+1) + j;
      int num_node3 = N_node_side * (i+1) + (j+1);
      int nums_node[3] = { num_node1, num_node2, num_node3 };
      EESM eesm = init_EESM(esm, nums_node);
      
      // update osm
      for(int m=0; m<N_node; m++) {
        for(int n=0; n<N_node; n++) {
          osm.array[i][j] += eesm.array[i][j];
        }
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
      
      // update osm
      for(int m=0; m<N_node; m++) {
        for(int n=0; n<N_node; n++) {
          osm.array[i][j] += eesm.array[i][j];
        }
      }
    }
  }

   
  for(int i=1; i<N_node; i++) {
    for(int j=1; j<N_node; j++) {
      printf("%f\n", osm.array[i][j]);
    }
  }
  
  return 0;
}
