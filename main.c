//
//  main.c
//  单细胞测序分类
//
//  Created by SALNO3 on 2018/5/19.
//  Copyright © 2018年 SALNO3. All rights reserved.
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>
#include "chain.h"
#include "util.h"
#include "io.h"

#define KMODEL 8
#define PRIMERLEN 65536 // 4 ** 8
#define BARCODE_TABLE_SIZE 65536
#define BARCODE_LEN 8
#define FEATURE_TABLE_SIZE 65536
#define BARC_LIST_LEN 97

//#define OUTR1FILE "result.R1.fastq"
//#define OUTR2FILE "result.R2.fastq"

//#define OUTR1FILE "/Users/salno3/Desktop/单细胞测序数据分类2/debug/debug.R1.fastq"
//#define OUTR2FILE "/Users/salno3/Desktop/单细胞测序数据分类2/debug/debug.R2.fastq"
int c = 0;  //记录行数

char R1FILE[128] = {};
char R2FILE[128] = {};
char OUTPREFIX[128] = {};
char BARC_LIST_FILE[128] = {};
char PRIMER_LIST_FILE[128] = {};

int barc_ro1_tab[BARCODE_TABLE_SIZE];
// barc_ro1_tab == barcode round1 table
int barc_ro2_tab[BARCODE_TABLE_SIZE];
int barc_ro3_tab[BARCODE_TABLE_SIZE];

char barc_ro1_list[BARC_LIST_LEN][BARCODE_LEN+2];
// barc_ro2_list == barcode round2 list
char barc_ro2_list[BARC_LIST_LEN][BARCODE_LEN+2];
char barc_ro3_list[BARC_LIST_LEN][BARCODE_LEN+2];

int feature32_tab[FEATURE_TABLE_SIZE];
int feature21_tab[FEATURE_TABLE_SIZE];

char feature32[] = "GTGGCCGATGTTTCGCATCGGCGTACGACT";
char feature21[] = "ATCCACGTGCTTGAGAGGCCAGAGCATTCG";
// feature32 和 feature21 取 KMODEL == 8 时, 没有重复的 K 段
// feature32 指的是round3的feature + round2的feature
// feature21 指的是round2的feature + round1的feature

struct Node* primerTable[PRIMERLEN];

char* primer;

FILE* r1_file;
FILE* r2_file;
FILE* out_r1_file;
FILE* out_r2_file;
FILE* out_err1_file;
FILE* out_err2_file;
FILE* out_der1_file;
FILE* out_der2_file;
FILE* id_list;

unsigned long long end;
// 文件的末尾
unsigned long long start;
// 文件当前读到的位置
pthread_mutex_t read_mutex;
pthread_mutex_t write_mutex;

int t_number = 2;
// 默认情况下的线程数, 为2


int distance(char a, char b){
    // 距离函数
    
    if(a == b){
        return 0;
    }
    else{
        return 1;
    }
    
}

double oDistance(char* str1, char* str2){
    //    求两个barcode的欧氏距离(稍微做了些改动)
    //    计算公式: ∑(√(xi^2 - yi^2)*(1+(n-i)/n))
    //    越往后越容易错, 用上面的式子来保证 碱基错误带来的距离 大于 错误位置带来的距离
    //    即: 保证oDistance("AACTGGAC", "AACTGGAG") < oDistance("AACTGGAC", "TTCTGGAC")
    //           1.875000   一位错,位置为末尾           2.125000   两位错, 位置在起始
    double len = strlen(str1);
    
    float dist = 0;
    
    for(int i = 0; i < len; i++){
        //dist = dist + sqrt(distance(str1[i], str2[i]))*(1+(len-i-1)/len);
        dist = dist + sqrt(distance(str1[i], str2[i]))*(1+i/len);
        //printf("%c=%c %d,%f\t",str1[i],str2[i],i,dist);
    }
    //printf("\n");
    return dist;
}

void barcodeTrans(char* barc1, char barc_list[BARC_LIST_LEN][BARCODE_LEN+2]){
    // 这个函数在这里实际上做的是对barcode的一次修正,
    // 遍历barc_list, 找到与barc1最相似的一个barcode, 且距离小于limit(如果大于limit则不做任何处理),
    // 作为修正后的barc1的值
    double limit = 2;
    //    d < 2 时容错1位
    // 这里的limit作为两个barcode的欧氏距离的阈值
    double min_dist = DBL_MAX;
    //printf("%f\t",min_dist); 
    int index = -1;
    char buffer[BARCODE_LEN+5];
    for(int i = 0; i < BARC_LIST_LEN; i++){
        //printf("%d:",i+1);
        double d = oDistance(barc1, barc_list[i]);
        //printf("%d\t%f\n",i,d);
        if(min_dist > d){
            min_dist = d;
            index = i;
            if(min_dist == 0){
                break;
            }
        }
    }
    if(min_dist < limit){
        strcpy(barc1, barc_list[index]);
    }
}
// 因为有3轮测序, 所以有3张barcode表
void barc1Trans(char* barc1){
    //printf("bac01:%s\n",barc1);
    barcodeTrans(barc1, barc_ro1_list);
}

void barc2Trans(char* barc1){
    barcodeTrans(barc1, barc_ro2_list);
}

void barc3Trans(char* barc1){
    barcodeTrans(barc1, barc_ro3_list);
}

int barcode_judge(char barcode[]){
    barc1Trans(barcode);
    int n_barc = decodeSqe(barcode);
    int i_barc = barc_ro1_tab[n_barc];
    return i_barc;
}

int * featureposi(int *point, int istart, int iend, int posi[], char feature[], struct FASTQ *fastq){
    //int point[2] = *point;
    int len = strlen(feature32);
    int flag = -1;
    int count = 0;
    int KUP = 1;
    char feature_string[KMODEL+KUP]; 
    char fastqbuffer[KMODEL+KUP]; 
    char feature_stringtab[2][KMODEL+KUP];
    char fastq_stringtab[2][KMODEL+KUP];
    char buffer_bar[KMODEL+KUP];
      
    //printf("i-s:%d\te:%d\n",istart,iend);
    for(int k = istart; k <= iend; k++){    //记录出现错配的次数
        if(posi[k] == -1){
            count++;
        }
    }
    //printf("count:%d\n",count);
    // 特征30个碱基分开三段，0-7,8-21，22-29
    //printf("f32posi-s:%d\te:%d\n",posi[istart],posi[iend]);
    if((posi[istart] == 0) && (posi[iend] == 22)){    // 前后两段没有错误
        if(((posi[iend] - posi[istart]) == (iend - istart)) && (count <= 8)){  //中间段允许一个错配
            point[0] = istart;     // fastq的起始位置
            point[1] = iend + KMODEL;         // fastq的终止位置
        }
        if(((posi[iend] - posi[istart]) != (iend - istart)) && (count <= 8)){  //中间段存在插入或缺失
            point[0] = istart;
            point[1] = iend + KMODEL;   
        }
    }
    else if((posi[istart] != 0) && (posi[iend] == 22) && (count == 0)){    // 前段有错误
            //printf("i-s:%d\te:%d\n",istart,iend);
            if(istart < 8){
            	point[0] = 0;
            }
            else{
                strCopy(feature_string, feature, 0, posi[istart] - 1);
                strCopy(fastqbuffer, fastq->data, istart - posi[istart], istart - 1);
                if(posi[istart] == 2){
                    strCopy(feature_stringtab[0], feature, 0, posi[istart] - 1);
                    strCopy(feature_stringtab[1], feature, posi[istart] - 1, posi[istart]);
                    strCopy(fastq_stringtab[0], fastq->data, istart - posi[istart], istart - posi[istart] + 1);
                    strCopy(fastq_stringtab[1], fastq->data, istart - posi[istart] + 1, istart - posi[istart] + 2);
                    
                    if(strcmp(feature_stringtab[1],fastq_stringtab[0]) == 0){
                        point[0] = istart - posi[istart] - 1;
                    }
                    else if(strcmp(feature_stringtab[0],fastq_stringtab[1]) == 0){
                        point[0] = istart - posi[istart] + 1;
                    }
                    //printf("fea0:%s\tfea1:%s\tfq0:%s\tfq1:%s\n",feature_stringtab[0],feature_stringtab[1],fastq_stringtab[0],fastq_stringtab[1]);
                }
                
                else if(posi[istart] == 1){
                    strCopy(feature_stringtab[0], feature, 0, posi[istart]);
                    strCopy(feature_stringtab[1], feature, posi[istart], posi[istart] + 1);
                    strCopy(fastq_stringtab[0], fastq->data, istart - posi[istart] - 1, istart - posi[istart]);
                    strCopy(fastq_stringtab[1], fastq->data, istart - posi[istart] + 1, istart - posi[istart] + 2);

                    if((strcmp(feature_stringtab[0],fastq_stringtab[0]) == 0)&&(strcmp(feature_stringtab[1],fastq_stringtab[1]) == 0)){
                        strCopy(buffer_bar, fastq->data, istart - posi[istart] - 7, istart - posi[istart] + 1);
                        flag = barcode_judge(buffer_bar);
                        if (flag == -1){
                            point[0] = istart - posi[istart] - 1;
                        }
                        else {
                            point[0] = istart - posi[istart];
                        }
                        
                    }
                    else if(strcmp(feature_stringtab[1],fastq_stringtab[1]) == 0){
                        strCopy(buffer_bar, fastq->data, istart - posi[istart] - 7, istart - posi[istart] + 1);
                        flag = barcode_judge(buffer_bar);
                        if (flag == -1){
                            point[0] = istart - posi[istart];
                        }
                        else {
                            point[0] = istart - posi[istart] + 1;
                        }
                    }
                    //printf("fea0:%s\tfea1:%s\tfq0:%s\tfq1:%s\n",feature_stringtab[0],feature_stringtab[1],fastq_stringtab[0],fastq_stringtab[1]);
                }
                else if(strcmp(feature_string,fastqbuffer) == 0){
                    point[0] = istart - posi[istart];
                }
                else{
                    strCopy(feature_stringtab[0], feature_string, 0, strlen(feature_string)-1);
                    strCopy(feature_stringtab[1], feature_string, 1, strlen(feature_string));
                    strCopy(fastq_stringtab[0], fastqbuffer, 0, strlen(fastqbuffer)-1);
                    strCopy(fastq_stringtab[1], fastqbuffer, 1, strlen(fastqbuffer));
                    if(strcmp(feature_stringtab[1],fastq_stringtab[0]) == 0){
                        point[0] = istart - posi[istart] - 1;
                    }
                    else if(strcmp(feature_stringtab[0],fastq_stringtab[1]) == 0){
                        point[0] = istart - posi[istart] + 1;
                    }
                    else{
                        point[0] = -1;
                    }
                }
            }
            point[1] = iend + KMODEL;
     }
     else if((posi[istart] == 0) && (posi[iend] != 22) && (count == 0)){    // 后段有错误
            //printf("i-s:%d\te:%d\n",istart,iend);
            point[0] = istart;
            if(iend > 85){
            	point[1] = 94;
            }
            else{
                strCopy(feature_string, feature, posi[iend] + KMODEL + 1, 30);
                strCopy(fastqbuffer, fastq->data, iend + KMODEL + 1, istart + 30);
                //printf("f32:%s\tfq:%s\n",feature_string,fastqbuffer);
                if((strcmp(feature_string,fastqbuffer) == 0)&&(strcmp(feature_string,"") != 0)){
                    point[1] = iend  - posi[iend] + strlen(feature);
                }
                else if(strcmp(feature_string,"") == 0){
                    strCopy(buffer_bar, fastq->data, istart - posi[istart] + len, istart - posi[istart] + len + 8);
                    flag = barcode_judge(buffer_bar);
                    if (flag > -1){
                        point[1] = istart - posi[istart] + len ;                 
                    }
                    else{
                        point[1] = istart - posi[istart] + len - 1;
                    }
                }
                else{
                    strCopy(feature_stringtab[0], feature_string, 0, strlen(feature_string)-2);
                    strCopy(feature_stringtab[1], feature_string, 1, strlen(feature_string)-1);
                    strCopy(fastq_stringtab[0], fastqbuffer, 0, strlen(fastqbuffer)-2);
                    strCopy(fastq_stringtab[1], fastqbuffer, 1, strlen(fastqbuffer)-1);
                    
                    if(strcmp(feature_stringtab[0],fastq_stringtab[1]) == 0){
                        
                        point[1] = iend  - posi[iend] + strlen(feature) + 1;
                    }
                    else if(strcmp(feature_stringtab[1],fastq_stringtab[0]) == 0){
                        //printf("iend:%d\tf32end:%d\tf32len:%d\n",iend,posi[iend],strlen(feature));
                        point[1] = iend  - posi[iend] + strlen(feature) - 1;
                    }
                    else{
                        point[1] = -1;
                    }
                }
            }

      }
        
      //printf("point-s:%d\te:%d\n",point[0],point[1]);
      //printf("line:%d\n",c); 
      return point;
}


int transR1(struct FASTQ *fastq){
    // 处理R1的每条read
    unsigned long len = strlen(fastq->data);
    char buffer[KMODEL+1];
    int cut_posi = -1;
    // 切点位置
    //printf("%d\n",len);
    struct Node* record = (struct Node*) malloc(sizeof(struct Node));
    // 用来记录位置和得分
    record->value = -1;
    record->score = -1;
    
    
    for(int i = 0; i != len ; i += KMODEL){
        if(i > len - KMODEL){
            i = len - KMODEL;
        }
        strCopy(buffer, fastq->data, i, i+KMODEL);
        //printf("%d\t%s\n",i,buffer);
        struct Node* temp_node = primerTable[decodeSqe(buffer)];
        int primer_posi = temp_node->value;
        //printf("primer_posi:%d\n",primer_posi);
        while (primer_posi != -1) {
            int mistake = 0;
            int threshold = 1;
            
            for(int j = 0; j < len - i - KMODEL; j++){
                if(primer[primer_posi+KMODEL+j] != fastq->data[i+KMODEL+j]){  //序列与primer发生错配的碱基
                    mistake++;
                    if(mistake > threshold){
                        addNode(record, i, i+KMODEL+j-1);
                        //printf("%d\t%d\n",i, i+KMODEL+j);
                        break;
                    }
                }
                
            }
            temp_node = temp_node->next;
            if(temp_node == NULL){
                break;
            }
            primer_posi = temp_node->value;
        }
        
    }
    //printf("recode:%d\n",record->score);
    int max_step = 0; // 0
    struct Node* temp_node = record;
    while (temp_node != NULL) {
        int step = record->score - record->value; //> max_step;
        //printf("score:%d\tvalue:%d\tmax:%d\tstep:%d\n",record->score,record->value,max_step,step);
        if(step > max_step){
            max_step = step;
            cut_posi = record->score;
            //printf("%d\n",cut_posi);
        }
        temp_node = temp_node->next;
    }
    
    if(cut_posi != -1){
        char t_data[210];
        strCopy(t_data, fastq->data, cut_posi, strlen(fastq->data));
        strcpy(fastq->data, t_data);
        //printf("%s\n",fastq->data);
        strCopy(t_data, fastq->sanger, cut_posi, strlen(fastq->sanger));
        strcpy(fastq->sanger, t_data);
        deleteChain(record);
        record = NULL;
        return 1;
    }
    
    
    return 0;
    
}


int transR2(struct FASTQ *fastq){
    if(fastq == NULL){
        return 0;
    }
    //c++;
    int KUP = 1;
    //char *null_ptr = "NNNNNNNN";
//    kmodel 扩展值, 这里这个值是用来扩展buffer的,
//    比如 buffer[8] = "AAAATTTT";
//    在调用 strCopy时会在最后加上 '\0', 这时候就会产生溢出,
//    KUP在这里的作用是增加buffer的大小(例如: buffer[8+KUP]), 防止数组溢出
    char umi_buffer[15];
    char kbuffer[KMODEL+KUP];
    strCopy(umi_buffer, fastq->data, 0, 10);
    //printf("%s\n",umi_buffer);
    int len = strlen(fastq->data);
    //printf("%d\t%s\n", len,fastq->data);
    int f32posi[len+KUP] ;
    int f21posi[len+KUP] ;
    //  f32posi == feature32 position
    //  f21posi == feature21 position
    //  这两个量表示在 feature 上的位点

    int f32point[2] = {-1, -1};    //  这两个量表示 feature 在 fastq->data 上的位点
    int f21point[2] = {-1, -1};
    //  这两个量表示在 fastq->data 上的位点
    
    int f32_len = strlen(feature32);
    int f21_len = strlen(feature21);
    
    int i = 0;    //可以从0开始，barcode3为缺失

    int c_start = -1;
    int c_end = -1;
    int temp_posi[2] = {-1,-1};
    int i_start = -1;
    int i_end = -1;   
    // 记录i（feature比对在fastq-data的断点位置）的起始和终止位置     
   
    //  接下来的这个循环是找到fastq->data上有关feature32的断点
   
    for(; i + KMODEL < 66; i += 1){               //可以到 65，barcode1,2为缺失  94-29=65（0-64）， 如果i=0,i+KMODEL=8（0-7）;i=57,i+KMODEL=65（57-64）。
        strCopy(kbuffer, fastq->data, i, i + KMODEL);
        int temp = decodeSqe(kbuffer);
        f32posi[i] = feature32_tab[temp];
        //printf("f32posi:%d\ti:%d\n",f32posi[i],i);
        if(f32posi[i] != -1){
            if(temp_posi[0] == -1){
                temp_posi[0] = i;
                temp_posi[1] = i; 
                  
                //printf("i-s:%d\te:%d\n",i_start,i_end);
            }
            else {
                 if((i - temp_posi[0]) > (f32posi[i] - f32posi[temp_posi[0]] + 2)){
                     i_start = temp_posi[0] ;
                     i_end = temp_posi[1]; 
                     temp_posi[0] = i;
                 }
                 //i_start = (f32posi[i] <= f32posi[i_start]) ? i:i_start;
                 if((i - temp_posi[0] <= 23)){                // 超出特征范围，直接忽略
                     temp_posi[1] = (i > temp_posi[1]) ? i:temp_posi[1];    // 取出特征比对最后的位置
                         
                     //printf("i-s:%d\te:%d\n",i_start,i_end);
                       
                 }
                 else{
                      break;
                 }
            }
                
        }
    }
    //printf("is:%d\tie:%d\nts:%d\tte:%d\n",i_start,i_end,temp_posi[0],temp_posi[1]);
    if(i_end - i_start >= temp_posi[1] - temp_posi[0]){
        c_start = i_start;
        c_end = i_end;   
    }
    else{
        c_start = temp_posi[0] ;
        c_end = temp_posi[1];
    }
    //printf("start:%d\tend:%d\n",c_start,c_end);
    if (c_start == -1 || c_end == -1 ){
        //  如果两个断点中有任何一个没有找到, 则视为丢弃当前这条r2 read, 函数返回
        return 0;
    }
    featureposi(f32point, c_start, c_end, f32posi, feature32, fastq);
    //printf("f32point-s:%d\te:%d\n",f32point[0],f32point[1]);
    
    i = 29;     //从 29开始，barcode3,2为缺失  94-29=65（29-93）， 如果i=29,i+KMODEL=37（29-36）;i=86,i+KMODEL=94（86-93）。
    i_start = -1;  
    i_end = -1; 
    c_start = -1;
    c_end = -1;
    temp_posi[0] = -1;
    temp_posi[1] = -1;
    for(; i + KMODEL <= len; i += 1){
        strCopy(kbuffer, fastq->data, i, i + KMODEL);
        int temp = decodeSqe(kbuffer);
        f21posi[i] = feature21_tab[temp];
        //printf("f21posi:%d\ti:%d\n",f21posi[i],i);
        if(f21posi[i] != -1){
            if(temp_posi[0] == -1){
                temp_posi[0] = i;
                temp_posi[1] = i; 
                   
                //printf("i-s:%d\te:%d\n",i_start,i_end);
            }
            else {
                 if((i - temp_posi[0]) > (f21posi[i] - f21posi[temp_posi[0]] + 2)){
                     i_start = temp_posi[0] ;
                     i_end = temp_posi[1]; 
                     temp_posi[0] = i;
                 }
                 //i_start = (f32posi[i] <= f32posi[i_start]) ? i:i_start;
                 if((i - temp_posi[0] <= 23)){                // 超出特征范围，直接忽略
                     temp_posi[1] = (i > temp_posi[1]) ? i:temp_posi[1];    // 取出特征比对最后的位置
                         
                     //printf("i-s:%d\te:%d\n",i_start,i_end);
                       
                 }
                 else{
                     break;
                 }
             }
                
        }
    }
    if(i_end - i_start >= temp_posi[1] - temp_posi[0]){
        c_start = i_start;
        c_end = i_end;   
    }
    else{
        c_start = temp_posi[0] ;
        c_end = temp_posi[1];
    }
    //printf("start:%d\tend:%d\n",c_start,c_end);
    if (c_start == -1 || c_end == -1 ){
        //  如果两个断点中有任何一个没有找到, 则视为丢弃当前这条r2 read, 函数返回
        return 0;
    }
    //printf("f21point:%d\n",f21point); 
    featureposi(f21point, c_start, c_end, f21posi, feature21, fastq);
    //printf("f21point-s:%d\te:%d\n",f21point[0],f21point[1]);
    
    if (f32point[0] == -1 || f32point[1] == -1 || f21point[0] == -1 || f21point[1] == -1){
        //  如果两个断点中有任何一个没有找到, 则视为丢弃当前这条r2 read, 函数返回
        return 0;
    }

    char barc01[BARCODE_LEN+KUP];
    char barc02[BARCODE_LEN+KUP];
    char barc03[BARCODE_LEN+KUP];
    char cDNA[len/2];
    

    
    if(f32point[0]-BARCODE_LEN == -1){
        strCopy(barc03, fastq->data, 0, f32point[0]);
    }
    else{
        strCopy(barc03, fastq->data, f32point[0]-BARCODE_LEN, f32point[0]);
    }
    strCopy(barc02, fastq->data, f32point[1], f21point[0]);
    strCopy(barc01, fastq->data, f21point[1], f21point[1]+BARCODE_LEN);
    strCopy(cDNA, fastq->data, f21point[1]+BARCODE_LEN, len);
    //printf ("%s\n",fastq->info);
    //printf("start:%d\tend:%d\n",f21point+f21_len, f21point+f21_len+BARCODE_LEN);
    //printf("barc01:%s\tlen:%d\n",barc01,strlen(barc01));
    //printf("barc02:%s\tlen:%d\n",barc02,strlen(barc02));
    //printf("barc03:%s\tlen:%d\n",barc03,strlen(barc03));
    
    
    if(strlen(barc01) != 8){
        //printf("barc01:%s\tlen:%d\n",barc01,strlen(barc01));
        if(strlen(barc01) == 7) {
           char *temp=(char *) malloc(strlen(barc01) + 1);
           strcpy(temp, barc01); 
           strcat(temp, "N");
	   strcpy(barc01, temp);
           //printf("barc01:%s\tlen:%d\n",barc01,strlen(barc01)); 
        }
        else{
           strcpy(barc01, "NNNNNNNN");
        }
        cDNA[0] = '\0';
    }
    //printf("barc01:%s\tlen:%d\n",barc01,strlen(barc01));
    //printf("cDNA:%s\n",cDNA);
    /*if(strlen(barc02) != 8){
        strcpy(barc02, "NNNNNNNN");
    }*/
    int i_barc2 = -1;
    int n_barc2 = -1;
    if(strlen(barc02) != 8){
        //printf("barc01:%s\tlen:%d\n",barc01,strlen(barc01));
        if(strlen(barc02) == 7) {
           char *temp = (char *) malloc(strlen(barc02) + 1);
           char *barc02_temp = (char *) malloc(strlen(barc02) + 1);
           strcpy(barc02_temp, barc02);
           strcpy(temp, barc02_temp); 
           strcat(temp, "N");
	   strcpy(barc02_temp, temp);
           barc2Trans(barc02_temp);
           n_barc2 = decodeSqe(barc02_temp);
           i_barc2 = barc_ro2_tab[n_barc2]; 
           //printf("barcode2: %s\n",barc02);
           //printf("barcode2: %s\n",barc02_temp);
           //printf("barcode2position: %d\n",i_barc2);
           if(i_barc2 == -1){
               //temp=(char *) malloc(strlen(barc02) + 1);
               strcpy(temp,"N"); 
               strcat(temp, barc02);
	       strcpy(barc02, temp);
               barc2Trans(barc02);
               n_barc2 = decodeSqe(barc02);
               i_barc2 = barc_ro2_tab[n_barc2]; 
               //printf("barcode2: %s\n",barc02);
               
               //printf("barcode2position: %d\n",i_barc2);
           }
        }
        else if(strlen(barc02) >= 9){
           for(int l = 0; l < strlen(barc02) - BARCODE_LEN + 1; l++){
               char *temp=(char *) malloc(BARCODE_LEN + 1); 
               strCopy(temp, barc02, l, l + BARCODE_LEN);
               barc2Trans(temp);
               n_barc2 = decodeSqe(temp);
               i_barc2 = barc_ro2_tab[n_barc2];
               //printf("barcode2position: %d\n",i_barc2);
               if(i_barc2 != -1){break;}
           }
        }
        else{
            strcpy(barc02, "NNNNNNNN");
            barc2Trans(barc02);
            n_barc2 = decodeSqe(barc02);
            i_barc2 = barc_ro2_tab[n_barc2];
        }

    }
    else{
        barc2Trans(barc02);
        n_barc2 = decodeSqe(barc02);
        i_barc2 = barc_ro2_tab[n_barc2];
    }
    //printf("barc02:%s\tlen:%d\n",barc02,strlen(barc02));
    //printf("barcode2position: %d\n",i_barc2);
    if(strlen(barc03) != 8){
        //printf("barc03:%s\tlen:%d\n",barc03,strlen(barc03));
        if(strlen(barc03) == 7) {
           char *temp=(char *) malloc(strlen(barc03) + 2);
           strcpy(temp,"N"); 
           strcat(temp, barc03);
	   strcpy(barc03, temp);
           //printf("barc03:%s\tlen:%d\n",barc03,strlen(barc03)); 
        }
        else{
           strcpy(barc03, "NNNNNNNN");
        }
    }
    //printf("barc03:%s\tlen:%d\n",barc03,strlen(barc03));

    barc1Trans(barc01);
    
    barc3Trans(barc03);
    //    barcode 容错
    //printf("barc01:%s\tlen:%d\n",barc01,strlen(barc01));
    //printf("barc03:%s\tlen:%d\n",barc03,strlen(barc03));
    char info[210];
    
    int n_barc1 = decodeSqe(barc01);
    
    int n_barc3 = decodeSqe(barc03);
    
    int i_barc1 = barc_ro1_tab[n_barc1];
    
    int i_barc3 = barc_ro3_tab[n_barc3];
    //printf("barcode 1 2 3 position: %d\t%d\t%d\n",i_barc1,i_barc2,i_barc3);
    if(i_barc1 == -1 || i_barc2 == -1 || i_barc3 == -1){
        // 如果容错后的barcode不再barcodeList中,
        // 即错误的barcode
        // 则丢弃当前r2 read
        return 2;
    }
    
    //  这时候查到的index是从0开始的序列中的, 所以index整个向右偏移1
    i_barc1++;
    i_barc2++;
    i_barc3++;
    
    //printf("barcode 1 2 3 position: %d\t%d\t%d\n",i_barc1,i_barc2,i_barc3);
    fastq->info[0] = ' ';
    //    这行的是将文件原本的第一行中的 @ 删掉
    sprintf(info, "@%d_%d_%d:%s%s",
            i_barc1,
            i_barc2,
            i_barc3,
            umi_buffer,
            fastq->info);
    
    //    例如: @21_54_85:32141 表示:
    //    UMI 为 ATCTCGGCAT
    //    round1 中的 barcode 序为 21
    //    round2 中的 barcode 序为 54
    //    round3 中的 barcode 序为 85
    strcpy(fastq->info, info);
    strcpy(fastq->data, cDNA);
    strCopy(info, fastq->sanger, len - strlen(cDNA), len);
    //printf("%s\n",info);
    //    这里用 info 只是临时把它当做一个缓冲区来使用
    strcpy(fastq->sanger, info);
    
    
    return 1;
}

void generateFeatureTable(){
    for(int i = 0; i < FEATURE_TABLE_SIZE; i++){
        feature32_tab[i] = -1;
        feature21_tab[i] = -1;
    }
    //    unsigned long len32 = strlen(feature32);
    //    unsigned long len21 = strlen(feature21);
    //    这应该是为了更好地扩展性的常规写法, 但是在本例中:
    //    feature32 和 feature21 一样长, 所以这里就暂时用一个 len 长度
    unsigned long len = strlen(feature32);
    char buffer[KMODEL + 1];
    for(int i = 0; i < len - KMODEL + 1; i++){
        strCopy(buffer, feature32, i, i+KMODEL);
        //printf("%d->f32:%s\n",i,buffer);
        feature32_tab[decodeSqe(buffer)] = i;
        strCopy(buffer, feature21, i, i+KMODEL);
        feature21_tab[decodeSqe(buffer)] = i;
    }
}

void generatePrimerTable(){
    for(unsigned int i = 0; i < PRIMERLEN; i++){
        initChain(&primerTable[i], -1, 0);
    }
    char primer[20000];
    // 这里有一个巨大无比的bug, 这里暂时先写死,
    // 即无法获取所有primer加起来的总和的长度, 这是一个外部文件的输入
    // 自己写的获取文件大小的函数始终有问题, 所以这里就偷了个懒
    // 将primer的空间开的大了一点, 并且写死,
    // 当primer文件过大时, 这里就会产生一个bug
    getPrimerString(primer, PRIMER_LIST_FILE);
    
    long len = strlen(primer);
    //printf("%ld\n",len);
    char buffer[KMODEL];
    for(unsigned int i = 0; i < len - KMODEL + 1 ; i++){
        strCopy(buffer, primer, i, i+KMODEL);
        addNode(primerTable[decodeSqe(buffer)], i, 0);
        //printf("%d\n",i);
    }

    return ;
    
}


void getIndexInfo(char* buffer, char* info){
    //  获取index信息, 例: @23_74_12
    int len = strlen(info);
    for(int i = 0; i < len; i++){
        if(info[i] == ' '){
            strCopy(buffer, info, 0, i);
            break;
        }
    }
}




int generateBarcTable(char* filename){
    // 生成 Barcode Table
    // Barcode Table 有3张:
    // barc_ro1_tab
    // barc_ro2_tab
    // barc_ro3_tab
    FILE* f;
    f = fopen(filename, "r");
    
    if(f == NULL){
        printf("Barcode list file is not exist.\n");
        return 0;
    }
    
    for(int i = 0; i < BARCODE_TABLE_SIZE; i++){
        barc_ro1_tab[i] = -1;
        barc_ro2_tab[i] = -1;
        barc_ro3_tab[i] = -1;
    }
    char buffer[BARCODE_LEN+2];
    
    
    //    int i = 0;
    int flag = 0;
    
    for(int i = 0; fgets(buffer, BARCODE_LEN+2, f); ){
        if(buffer[0] == '@' && flag == 0){
            flag = 1;
            continue;
        }else if(buffer[0] == '@' && flag == 1){
            break;
        }else{
            buffer[BARCODE_LEN] = '\0';
            barc_ro1_tab[decodeSqe(buffer)] = i;
            strcpy(barc_ro1_list[i], buffer);
            i++;
        }
    }
    
    for(int i = 0; fgets(buffer, BARCODE_LEN+2, f); ){
        if(buffer[0] == '@'){
            break;
        }else{
            buffer[BARCODE_LEN] = '\0';
            barc_ro2_tab[decodeSqe(buffer)] = i;
            strcpy(barc_ro2_list[i], buffer);
            i++;
        }
    }
    
    for(int i = 0; fgets(buffer, BARCODE_LEN+2, f); ){
        if(buffer[0] == '@' || buffer[0] == '\n'){
            break;
        }else{
            buffer[BARCODE_LEN] = '\0';
            barc_ro3_tab[decodeSqe(buffer)] = i;
            strcpy(barc_ro3_list[i], buffer);
            i++;
        }
    }
    fclose(f);
    return 1;
}





void sub_main(){
    struct FASTQ readr1;
    struct FASTQ readr2;
    //printf("%d\t%d\n",start,end);
    
    while (start < end) {
        pthread_mutex_lock(&read_mutex);
        getRead(&readr2, r2_file);
        getRead(&readr1, r1_file);
        start = ftell(r2_file);
        pthread_mutex_unlock(&read_mutex);
        //printf("start:%d\n",start);
        int flag = transR2(&readr2);
        //printf("flag:%d\n",flag);
        if(flag == 1){
            //printRead(readr2);
            //            write2file(readr2, out_r2_file);
            
            transR1(&readr1);
            char buffer_1[20];
            char buffer_2[210];
            getIndexInfo(buffer_1, readr2.info);
            strcpy(buffer_2, readr1.info);
            fprintf(id_list, "%s\n", buffer_2);
            //printf("%s\t%s\n", buffer_1, buffer_2);
            buffer_2[0] = ' ';
            sprintf(readr1.info, "%s%s", buffer_1, buffer_2);
            pthread_mutex_lock(&write_mutex);
            write2file(readr1, out_r1_file);
            write2file(readr2, out_r2_file);
            pthread_mutex_unlock(&write_mutex);
        }
        else if(flag == 0){        // feature断点没有找到
            pthread_mutex_lock(&write_mutex);
            write2file(readr1, out_err1_file);
            write2file(readr2, out_err2_file);
            pthread_mutex_unlock(&write_mutex);
        }
        else if(flag == 2){       // 无对应的barcode
            pthread_mutex_lock(&write_mutex); 
            write2file(readr1, out_der1_file);
            write2file(readr2, out_der2_file);
            pthread_mutex_unlock(&write_mutex);
        }
        
        printf("start: %lld, end: %lld\n", start, end);
        printf("%.2lf%% finishing...\n", ((double)start)/end*100);
    }
}

void printHelp(){
    printf("-r1 R1File\n");
    printf("-r2 R2File\n");
    printf("-p primerList\n");
    printf("-b barcodeList\n");
    printf("-t number of threads\n");
    printf("-o output prefix\n");
}


void analyzeArgs(int args, const char * argv[]){
    int i = 1;

    char* _r1 = "-r1";
    char* _r2 = "-r2";
    char* _b = "-b";
    char* _p = "-p";
    char* _t = "-t";
    char* _o = "-o";
    char* _h = "-h";
    
    if(args == 1){
        printHelp();
    }
    
    for(; i < args; ){
        if(!strcmp(argv[i], _r1)){
            i++;
            strcpy(R1FILE, argv[i]);
            i++;
        }else if(!strcmp(argv[i], _r2)){
            i++;
            strcpy(R2FILE, argv[i]);
            i++;
        }else if(!strcmp(argv[i], _b)){
            i++;
            strcpy(BARC_LIST_FILE, argv[i]);
            i++;
        }else if(!strcmp(argv[i], _p)){
            i++;
            strcpy(PRIMER_LIST_FILE, argv[i]);
            i++;
        }else if(!strcmp(argv[i], _t)){
            i++;
            t_number = atoi(argv[i]);
            i++;
        }else if(!strcmp(argv[i], _o)){
            i++;
            strcpy(OUTPREFIX, argv[i]);
            i++;
        }else if(!strcmp(argv[i], _h)){
            printHelp();
            exit(5);
        }
        else{
            i++;
            printf("Illegal parameters.\n");
            i++;
        }
    }
}


int main(int argc, const char * argv[]) {
    
    analyzeArgs(argc, argv);
    
    r1_file = fopen(R1FILE, "r");
    if(r1_file == NULL){
        printf("R1 file does not exist.\n");
        return 1;
    }
    r2_file = fopen(R2FILE, "r");
    if(r2_file == NULL){
        printf("R2 file does not exist.\n");
        return 2;
    }
    
    char OUTR1FILE[128]={} ;
    strcpy(OUTR1FILE,OUTPREFIX);
    //printf("%s\n",OUTR1FILE); 
    strcat(OUTR1FILE,".R1.fastq");
    
    char OUTR2FILE[128]={} ;
    strcpy(OUTR2FILE,OUTPREFIX);
    //printf("%s\n",OUTR2FILE); 
    strcat(OUTR2FILE,".R2.fastq");
    
    char OUTERR1FILE[128]={} ;
    strcpy(OUTERR1FILE,OUTPREFIX);
    //printf("%s\n",OUTR1FILE); 
    strcat(OUTERR1FILE,".ERR1.fastq");
    
    char OUTERR2FILE[128]={} ;
    strcpy(OUTERR2FILE,OUTPREFIX);
    //printf("%s\n",OUTR2FILE); 
    strcat(OUTERR2FILE,".ERR2.fastq");
    
    char OUTDER1FILE[128]={} ;
    strcpy(OUTDER1FILE,OUTPREFIX);
    //printf("%s\n",OUTR1FILE); 
    strcat(OUTDER1FILE,".DER1.fastq");
    
    char OUTDER2FILE[128]={} ;
    strcpy(OUTDER2FILE,OUTPREFIX);
    //printf("%s\n",OUTR2FILE); 
    strcat(OUTDER2FILE,".DER2.fastq");
    
    char ID_LIST[128]={} ;
    strcpy(ID_LIST,OUTPREFIX);
    //printf("%s\n",OUTR1FILE); 
    strcat(ID_LIST,".id.list");
    
    //char* OUTR2FILE = OUTPREFIX;
    //strcat(OUTR2FILE,".R2.fastq");
    //printf("%s\n",OUTR1FILE);  
    //printf("%s\n",OUTR2FILE);  
  
    out_r1_file = fopen(OUTR1FILE, "w");
    out_r2_file = fopen(OUTR2FILE, "w");
    out_err1_file = fopen(OUTERR1FILE, "w");
    out_err2_file = fopen(OUTERR2FILE, "w");
    out_der1_file = fopen(OUTDER1FILE, "w");
    out_der2_file = fopen(OUTDER2FILE, "w");
    id_list = fopen(ID_LIST, "w");

    primer = malloc(getFileSize(PRIMER_LIST_FILE));
    int ret = getPrimerString(primer, PRIMER_LIST_FILE);
    if(ret == 0){
        exit(3);
    }
    generatePrimerTable();
    ret = generateBarcTable(BARC_LIST_FILE);
    if(ret == 0){
        exit(4);
    }
    generateFeatureTable();
    
    fseek(r2_file, 0L, SEEK_END);
    end = ftell(r2_file);
    fseek(r2_file, 0L, SEEK_SET);
    start = ftell(r2_file);
    
    pthread_mutex_init(&read_mutex, NULL);
    pthread_mutex_init(&write_mutex, NULL);
    
    pthread_t thread_pool[t_number];
    int tid_pool[t_number];
    
    for(int i = 0; i < t_number; i++){
        tid_pool[i] = pthread_create(&thread_pool[i], NULL, (void*) &sub_main , NULL);
        sleep(5);
    }
    
    for(int i = 0; i < t_number; i++){
        pthread_join(thread_pool[i], NULL);
    }

    
    fclose(r1_file);
    fclose(r2_file);
    fclose(out_r1_file);
    fclose(out_r2_file);
    fclose(out_err1_file);
    fclose(out_err2_file);
    fclose(out_der1_file);
    fclose(out_der2_file);
    return 0;
}
