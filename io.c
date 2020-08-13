//
//  io.c
//  单细胞测序数据分类2
//
//  Created by SALNO3 on 2018/5/23.
//  Copyright © 2018年 SALNO3. All rights reserved.
//

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
    //printf("result:%s\n",result);
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
    //printf("string:%d\t",strlen(buffer));
    //printf("%d\t%s\n",strlen(fastq->info),fastq->info);
    strcpy(fastq->data, getNext(buffer, f));
    fastq->data[strlen(buffer)-1] = 0;
    //printf("string:%d\n",strlen(buffer));
    strcpy(fastq->comment, getNext(buffer, f));
    fastq->comment[strlen(buffer)-1] = 0;
    //printf("string:%d\n",strlen(buffer));
    strcpy(fastq->sanger, getNext(buffer, f));
    fastq->sanger[strlen(buffer)-1] = 0;
    //printf("string:%d\n",strlen(buffer));
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
