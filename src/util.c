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
    util_c: Basic basic transformation between base and number
 */

#include "util.h"



int base2int(char c){
    if (c == 'A'){
        return 0b00;
    }
    if (c == 'T'){
        return 0b01;
    }
    if (c == 'G'){
        return 0b10;
    }
    if (c == 'C'){
        return 0b11;
    }
     
    return 0;
}

//a sequence is converted to a decimal number

int decodeSqe(char* seq){
    int len = 0;
    len = strlen(seq);
    int num = 0;

    for(int i = 0; i < len; i++){
        num += base2int(seq[i]) * pow(4, len - i - 1);                               //4 ** (len-i-1) 
    }
    return num;
}

int bitCompare(char* seq1, char* seq2){
    int s1, s2;
    s1 = decodeSqe(seq1);
    s2 = decodeSqe(seq2);
    return s1^s2?0:1;                                                                //bitwise xor，if buffer1 is equal to buffer2，the value is 0, otherwise 1
}

void strCopy(char str[], char template[], int start, int end){
    if(template[0] == '\0'){
        return ;
    }
    if(end < start){
        return ;
    }
    char* ptr;
    ptr = template + start;
    strncpy(str, ptr, end - start);
    *(str + end - start) = '\0';
}

unsigned long long getFileSize(char* filename)
{
    struct stat statbuf;
    stat(filename,&statbuf);
    long long size=statbuf.st_size;
    return size;
}
