//
//  util.h
//  单细胞测序数据分类2
//
//  Created by SALNO3 on 2018/5/23.
//  Copyright © 2018年 SALNO3. All rights reserved.
//

#ifndef util_h
#define util_h

#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <math.h>

int base2int(char c);
// 将碱基转为数字

int decodeSqe(char* seq);
// 编码一段DNA序列

void strCopy(char* str, char* template, int start, int end);
// 字符串复制, 将 template[start:end]的复制入 str

unsigned long long getFileSize(char* filename);
// 这个函数有个bug, 无法获取大型文件的大小

int bitCompare(char* seq1, char* seq2);
// 按位比较两个DNA序列
#endif /* util_h */
