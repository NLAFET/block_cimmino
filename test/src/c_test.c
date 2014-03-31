#include "abcd_c.h"
#include <stdio.h>

int main()
{
  typedef struct abcd_solver abcd;
  int i;
  abcd *obj = new_solver();
  printf("%p\n", obj->icntl);
  for(i = 0; i < 20; i++)
    printf("%d   %d\n", i, obj->icntl[i]);
  return 0;
}
