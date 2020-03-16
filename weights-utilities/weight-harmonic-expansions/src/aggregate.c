#include <stdlib.h>
#include <stdio.h>

void copy(FILE* in, FILE* out) {
  
  char* buffer = malloc(sizeof(char)*10*1024*1024); // 10mb chunk
  int size;
  do {
    size = fread(buffer, sizeof(char), 10*1024*1024, in);
    if (size <= 0) break;
    fwrite(buffer, sizeof(char), size, out);
  } while (size == sizeof(buffer));
  //reached EOF, exit copying process
  free(buffer);
}

void aggregate_chunks(char*** output_filenames, char**** chunks_filenames, int nchunks, int ndegrees) {
  
  int i;
  for (i = 0; i < ndegrees; ++i) {
    FILE* fout_real = fopen(output_filenames[i][0], "w");
    FILE* fout_imag = fopen(output_filenames[i][1], "w");
     
    int j; 
    for (j = 0; j < nchunks; ++j) {
      FILE *fin_chunk_real = fopen(chunks_filenames[j][i][0], "r");
      FILE *fin_chunk_imag = fopen(chunks_filenames[j][i][1], "r");
      copy(fin_chunk_real, fout_real);
      copy(fin_chunk_imag, fout_imag);

      fclose(fin_chunk_real);
      fclose(fin_chunk_imag);   
    }
    fclose(fout_real);
    fclose(fout_imag);
  }
}
