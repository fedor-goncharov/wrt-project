#include <stdio.h>
#include <stdlib.h>


int init_parameters_alloc(int** parameters) {
  *parameters = (int*)malloc(4 * sizeof(int));
  
  return 0;
}

int init_parameters(char* filename, int* parameters) {
  
  FILE *f;
  char buffer[128];
  
  f = fopen(filename, "r");
  int i;
  
  for (i = 0; i < 4; ++i) {
    fscanf(f, "%d %[^\n]", &(parameters[i]), buffer);
    fgetc(f);
  }
  fclose(f);
  return 0;
  
}


int main(void) {
  
  FILE *f;
  char buffer[128];
  
  f = fopen("config.txt", "r");
  int i;
  
  int* parameters;
  init_parameters_alloc(&parameters);
  init_parameters("config.txt", parameters);
  
  for (i = 0; i < 4; ++i) {
     printf("value:%d\n", parameters[i]);
  }
  
  free(parameters);
  fclose(f);
  
}