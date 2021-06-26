/*
  The MIT License

   SCSit  (A high-efficiency preprocessing tool for single-cell sequencing data from SPLiT-seq),


   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

   Copyright (C) 2021 Intel Corporation, Bioinformatics Lab, Hainan University.
*/

/*
*  io_c: Basic I/O of FASTQ format
*/

#include "io.h"
#include <string.h>

char buffer[BUFFERLEN];

int getPrimerString(char* result, char* file_name){
    FILE* f = fopen(file_name, "r");
    if(f == NULL){
        printf("Primer file is not exist.\n");
        return 0;
    }
    unsigned long i = 0;
    for(char temp_c = fgetc(f); temp_c != EOF; temp_c = fgetc(f)){
        if( temp_c == '\n'){
            continue;
        }
        *(result+i) = temp_c;
        i++;
    }
    *(result+i) = '\0';
    fclose(f);
    return 1;
}

void write2file(struct FASTQ fastq, FILE* f){
    fprintf(f, "%s\n", fastq.info);
    fprintf(f, "%s\n", fastq.data);
    fprintf(f, "%s\n", fastq.comment);
    fprintf(f, "%s\n", fastq.sanger);
}


int getRead(struct FASTQ* fastq, FILE* f){
    strcpy(fastq->info, getNext(buffer, f));
    fastq->info[strlen(buffer)-1] = 0;
    strcpy(fastq->data, getNext(buffer, f));
    fastq->data[strlen(buffer)-1] = 0;
    strcpy(fastq->comment, getNext(buffer, f));
    fastq->comment[strlen(buffer)-1] = 0;
    strcpy(fastq->sanger, getNext(buffer, f));
    fastq->sanger[strlen(buffer)-1] = 0;
    if(fastq->info[0] && fastq->data[0] && fastq->comment[0] && fastq->sanger[0]){
        return 1;
    }
    return 0;
}

void printRead(struct FASTQ read){
    printf("%s\n", read.info);
    printf("%s\n", read.data);
    printf("%s\n", read.comment);
    printf("%s\n", read.sanger);
}


void* getNext(char* buffer_arg, FILE* f){
    return fgets(buffer_arg, BUFFERLEN, f);
}
