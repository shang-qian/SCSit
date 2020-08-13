//
//  util.c
//  单细胞测序数据分类2
//
//  Created by SALNO3 on 2018/5/23.
//  Copyright © 2018年 SALNO3. All rights reserved.
//

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

int decodeSqe(char* seq){
    int len = 0;
    len = strlen(seq);
    int num = 0;
	//一条序列转化成一个十进制的数
    for(int i = 0; i < len; i++){
        num += base2int(seq[i]) * pow(4, len - i - 1); //4 的 len-i-1 次方的值
    }
    return num;
}

int bitCompare(char* seq1, char* seq2){
    int s1, s2;
    s1 = decodeSqe(seq1);
    s2 = decodeSqe(seq2);
    return s1^s2?0:1;    //按位异或，如果buffer1和buffer2相同，则为0，否则为1
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
