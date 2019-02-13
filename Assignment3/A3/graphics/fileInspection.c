#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]){

    char *filename = argv[2];
    FILE *fp = fopen(filename, "rb");
    fseek(fp, 0L, SEEK_END);
    int sz = ftell(fp);
    printf("Size of the file: %d\n", sz);


    return 0;
}