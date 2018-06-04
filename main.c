#include <stdio.h>

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
  double x1,x2,x3,y1,y2,y3
  double a1,a2,a3,b1,b2,b3,c1,c2,c3;
  double A; // area of a target element
  
  // initialize
  x1 = node1.x;
  y1 = node1.y;
  x2 = node2.x;
  y2 = node2.y;
  x3 = node3.x;
  y3 = node3.y;

  A = (x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;

  // calculate a,b,c
  a1 = (x2*y3-x3*y2)/(2*A);
}

int main(void) {
  double L = 10.0;
  int N_side = 2;

  double l = L / N_side;
  int N_node_side = N_side + 1;

  // initialize
  Node nodes[N_node_side][N_node_side];
  for(int i=0; i<N_node_side; i++) {
    for(int j=0; j<N_node_side; j++) {
      nodes[i][j] = init_Node(l*i, l*j);
    }
  }


  return 0;
}
