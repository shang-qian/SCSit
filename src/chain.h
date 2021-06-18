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
    chain_h: init chain list
 */


#ifndef chain_h
#define chain_h
#include <stdio.h>

//  typedef struct LNode 

struct Node{
    struct Node* next;
    struct Node* end;
    int value;
    int score;
};

void initChain(struct Node* *node, int value, int score);

void addNode(struct Node* head, int value, int score);

void printChain(struct Node* node);

struct Node* getNode(struct Node* node , int index);

struct Node* getMaxOfScore(struct Node* node);

int getMaxValueByScoreAndReset(struct Node* node);

void deleteChain(struct Node* node);

#endif /* chain_h */
