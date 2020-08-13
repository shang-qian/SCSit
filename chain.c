//
//  chain.c
//  单细胞测序数据分类2
//
//  Created by SALNO3 on 2018/5/23.
//  Copyright © 2018年 SALNO3. All rights reserved.
//

#include "chain.h"
#include <stdlib.h>

void deleteChain(struct Node* node){
    while (node) {
        struct Node* next = node->next;
//        printf("v = %d, s = %d, ptr = %p\n", node->value, node->score, node);
        free(node);
        node = next;
    }
//    if(node != NULL){
//        deleteChain(node->next);
//        printf("v = %d, s = %d, ptr = %p\n", node->value, node->score, node);
//        free(node);
//    }
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
        printf("index: %d, value: %d, score: %d\n", i, node->value, node->score);
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
