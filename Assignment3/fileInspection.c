#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

  if (argc != 2) {
    printf("Usage: ./fileInspection fileName\n");
    return -1;
  }

  char *filename = argv[1];
  FILE *fp = fopen(filename, "rb");
  fseek(fp, 0L, SEEK_END);
  size_t sz = ftell(fp);
  fseek(fp, 0L, SEEK_SET);
  printf("Size of the file: %zu\n", sz);

  int i;
  double f;
  for (i=0; i<sz/48; i++) {
    fread(&f, sizeof(double), 1, fp);
    printf("%f ", f);
    fread(&f, sizeof(double), 1, fp);
    printf("%f ", f);
    fread(&f, sizeof(double), 1, fp);
    printf("%f ", f);
    fread(&f, sizeof(double), 1, fp);
    printf("%f ", f);
    fread(&f, sizeof(double), 1, fp);
    printf("%f\n", f);
  }

  fclose(fp);



  return 0;
}