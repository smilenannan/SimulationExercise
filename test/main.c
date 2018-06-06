#include <stdio.h>

int main(void) {
  int a[3][3];
  int e = 0;
  for(int i=0; i<3; i++) {
    for(int j=0; j<3; j++) {
      a[i][j] = e;
      e++;
    }
  }


  return 0;
}
