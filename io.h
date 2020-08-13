//
//  io.h
//  单细胞测序数据分类2
//
//  Created by SALNO3 on 2018/5/23.
//  Copyright © 2018年 SALNO3. All rights reserved.
//

#ifndef io_h
#define io_h

#include <stdio.h>

#define BUFFERLEN 210



struct FASTQ {
    char info[200];
    char data[200];
    char comment[200];
    char sanger[200];
};

void write2file(struct FASTQ fastq, FILE* f);

int getPrimerString(char result[], char* file_name);

int getRead(struct FASTQ* fastq, FILE* f);

void printRead(struct FASTQ read);

void* getNext(char* buffer_arg, FILE* f);


#endif /* io_h */
