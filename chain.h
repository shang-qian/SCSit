//
//  chain.h
//  单细胞测序数据分类2
//
//  Created by SALNO3 on 2018/5/23.
//  Copyright © 2018年 SALNO3. All rights reserved.
//

#ifndef chain_h
#define chain_h
#include <stdio.h>
// 链表数据结构
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
