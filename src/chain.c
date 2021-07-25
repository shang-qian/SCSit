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
*   chain_c:   primer table stored in chain list
*/

#include "chain.h"
#include <stdlib.h>

void deleteChain(struct Node* node){
    while (node) {
        struct Node* next = node->next;
		free(node);
        node = next;
    }

}

int getMaxValueByScoreAndReset(struct Node* node){
    int max = 0;
    int value = -1;
    for(int i = 0; node; node = node->next){
        if(node->score > max){
            max = node->score;
            value = node->value;
        }
        i++;
        node->score = 0;
    }
    return value;
}

void initChain(struct Node* *node, int value, int score){
    *node = (struct Node *)malloc(1 * sizeof(struct Node));
    (*node)->value = value;
    (*node)->score = score;
}

void addNode(struct Node* head, int value, int score){
    if(head->next == NULL){
        head->value = value;
        head->score = score;
        return ;
    }
    while (head->next) {
        head = head->next;
    }
    head->next = (struct Node*)malloc(1 * sizeof(struct Node));
    head->next->value = value;
    head->next->score = score;
    return ;
}

void printChain(struct Node* node){
    for(int i = 0; node; node = node->next){
        
        i++;
    }
}

struct Node* getNode(struct Node* node , int index){
    for(int i = 0; node; node = node->next){
        if(i == index){
            return node;
        }
        i++;
    }
    return NULL;
}

struct Node* getMaxOfScore(struct Node* node){
    struct Node* max_node = node;
    for(int i = 0; node; node = node->next){
        if(node->score > max_node->score){
            max_node = node;
        }
        i++;
    }
    return max_node;
}

/*  chain_c    */
