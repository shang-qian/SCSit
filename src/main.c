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
    main_c : mainly function of SCSit
    label three round barcodes and UMI
    filter out primers
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>
#include <ctype.h>
#include "chain.h"
#include "util.h"
#include "io.h"

#define KMODEL 8
#define PRIMERLEN 65536                                                       // 4 ** 8
#define BARCODE_TABLE_SIZE 65536
#define BARCODE_LEN 8
#define FEATURE_TABLE_SIZE 65536
#define BARC_LIST_LEN 97

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "1.0.0"
#endif


char R1FILE[128] = {};
char R2FILE[128] = {};
char OUTPREFIX[128] = "output";
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

/*

feature32 refers to feature of round3 and round2
feature21 refers to feature of round2 and round1

*/

struct Node* primerTable[PRIMERLEN];

char* primer;

FILE* r1_file;
FILE* r2_file;
FILE* out_r1_file;
FILE* out_r2_file;


unsigned long long end;                                                       // end of file
unsigned long long start;                                                     // start of file
pthread_mutex_t read_mutex;
pthread_mutex_t write_mutex;

int t_number = 2;                                                             // default: 2


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
    /*    
	     Euclidean distance between two barcodes
         distance formula: ∑(√(xi^2 - yi^2)*(1+(n-i)/n)) 
		 Errors tend to occur later in the sequence, so it makes sure that the distance from base error is greater than the distance from positional fault using above formula.
         e.g.: oDistance("AACTGGAC", "AACTGGAG") < oDistance("AACTGGAC", "TTCTGGAC") 
	*/
    double len = strlen(str1);
    
    float dist = 0;
    
    for(int i = 0; i < len; i++){
        //dist = dist + sqrt(distance(str1[i], str2[i]))*(1+(len-i-1)/len);
        dist = dist + sqrt(distance(str1[i], str2[i]))*(1+i/len);
        
    }
    
    return dist;
}

void barcodeTrans(char* barc1, char barc_list[BARC_LIST_LEN][BARCODE_LEN+2]){
	/*
	    barcode fault-tolerant and correction
        limit is the threshhold of euclidean distance comparing two barcodes
		Traversing from barcode list, if barcode is most similar to barc1/barc2/barc3 and distance between barcodes is less than limit that barcode used as corrected barc1/barc2/barc3.
		If distance between barcodes is larger than the limit, then no processing will be done.
	*/
   
    double limit = 2;
    
    double min_dist = DBL_MAX;
   
    int index = -1;
    char buffer[BARCODE_LEN+5];
    for(int i = 0; i < BARC_LIST_LEN; i++){
        
        double d = oDistance(barc1, barc_list[i]);                         // d < 2, fault-tolerant one base
        
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

// Three rounds barcodes
void barc1Trans(char* barc1){
    
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
    
    int len = strlen(feature32);
    int flag = -1;
    int count = 0;
    int KUP = 1;
    char feature_string[KMODEL+KUP]; 
    char fastqbuffer[KMODEL+KUP]; 
    char feature_stringtab[2][KMODEL+KUP];
    char fastq_stringtab[2][KMODEL+KUP];
    char buffer_bar[KMODEL+KUP];

    //Record the number of mismatches  
    for(int k = istart; k <= iend; k++){                                           
        if(posi[k] == -1){
            count++;
        }
    }
    
    // Feature of 30 bases divided into three segments，0-7,8-21，22-29


    // The first and third is no mistakes.
    if((posi[istart] == 0) && (posi[iend] == 22)){                                  
        if(((posi[iend] - posi[istart]) == (iend - istart)) && (count <= 8)){      // The second segment allows a mismatch
            point[0] = istart;                                                     // the start of reads
            point[1] = iend + KMODEL;                                              // the end of reads
        }
        if(((posi[iend] - posi[istart]) != (iend - istart)) && (count <= 8)){      // The second segment allows a insertion and deletion
            point[0] = istart;
            point[1] = iend + KMODEL;   
        }
    }

	// The first is a mistakes.

    else if((posi[istart] != 0) && (posi[iend] == 22) && (count == 0)){    
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

	// The third is a mistakes.

     else if((posi[istart] == 0) && (posi[iend] != 22) && (count == 0)){    
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
    /*  process Read1 and identify primers  */

    unsigned long len = strlen(fastq->data);
    char buffer[KMODEL+1];
    int cut_posi = -1;                                                            // cut position 
    
    struct Node* record = (struct Node*) malloc(sizeof(struct Node));             // record position and score
    record->value = -1;
    record->score = -1;
    
    
    for(int i = 0; i != len ; i += KMODEL){
        if(i > len - KMODEL){
            i = len - KMODEL;
        }
        strCopy(buffer, fastq->data, i, i+KMODEL);
        
        struct Node* temp_node = primerTable[decodeSqe(buffer)];
        int primer_posi = temp_node->value;
        
        while (primer_posi != -1) {
            int mistake = 0;
            int threshold = 1;
            
            for(int j = 0; j < len - i - KMODEL; j++){
                if(primer[primer_posi+KMODEL+j] != fastq->data[i+KMODEL+j]){       // The base is mismatched comparing read to a primer
                    mistake++;
                    if(mistake > threshold){
                        addNode(record, i, i+KMODEL+j-1);
                        
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
    int max_step = 0;                                                            
    struct Node* temp_node = record;
    while (temp_node != NULL) {
        int step = record->score - record->value;                                 
        //printf("score:%d\tvalue:%d\tmax:%d\tstep:%d\n",record->score,record->value,max_step,step);
        if(step > max_step){
            max_step = step;
            cut_posi = record->score;
            
        }
        temp_node = temp_node->next;
    }
    
    if(cut_posi != -1){
        char t_data[210];
        strCopy(t_data, fastq->data, cut_posi, strlen(fastq->data));
        strcpy(fastq->data, t_data);
        
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
    
    int KUP = 1;                                                                        // extended value of KMODEL 
    char umi_buffer[15];
    char kbuffer[KMODEL+KUP];
    strCopy(umi_buffer, fastq->data, 0, 10);
    
    int len = strlen(fastq->data);
   
    //  f32posi and f21posi refer to position of feature: f32posi == feature32 position, f21posi == feature21 position
    int f32posi[len+KUP] ;
    int f21posi[len+KUP] ;
    

    //  f32point[2] and f21point[2] refer to position of feature existing on fastq->data
    int f32point[2] = {-1, -1};    
    int f21point[2] = {-1, -1};
    
	
    
    int f32_len = strlen(feature32);
    int f21_len = strlen(feature21);
    
    int i = 0;                                                                           // Begin with 0，barcode3 is missing

    int c_start = -1;
    int c_end = -1;
    int temp_posi[2] = {-1,-1};
    int i_start = -1;
    int i_end = -1;   
    //  find breakpoints of feature32 in fastq->data
   
    for(; i + KMODEL < 66; i += 1){                                                      //  Until 65，barcode2 and barcode1 are missing， if i=0,i+KMODEL=8(0-7);if i=57,i+KMODEL=65(57-64)
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
                 if((i - temp_posi[0] <= 23)){                
                     temp_posi[1] = (i > temp_posi[1]) ? i:temp_posi[1];    
                         
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
        return 0;
    }
    featureposi(f32point, c_start, c_end, f32posi, feature32, fastq);
    //printf("f32point-s:%d\te:%d\n",f32point[0],f32point[1]);
    
    i = 29;                                                                              //Begin with 29，barcode3 and barcode2 are missing， if i=29,i+KMODEL=37(29-36);if i=86,i+KMODEL=94(86-93)。
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
                 if((i - temp_posi[0] <= 23)){                                         // Out of the range of features, directly ignored
                     temp_posi[1] = (i > temp_posi[1]) ? i:temp_posi[1];               // Take out the last position of feature alignment
                         
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
        return 0;
    }
    
    featureposi(f21point, c_start, c_end, f21posi, feature21, fastq);
    //printf("f21point-s:%d\te:%d\n",f21point[0],f21point[1]);

    //  If either of the two breakpoints is not found, Read2 is given up.
    if (f32point[0] == -1 || f32point[1] == -1 || f21point[0] == -1 || f21point[1] == -1){
        
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

    //    fault-tolerant barcode 
    
    char info[210];
    
    int n_barc1 = decodeSqe(barc01);
    
    int n_barc3 = decodeSqe(barc03);
    
    int i_barc1 = barc_ro1_tab[n_barc1];
    
    int i_barc3 = barc_ro3_tab[n_barc3];
    
   // If fault-tolerant barcode is not exist on barcode list, Read2 is given up
	if(i_barc1 == -1 || i_barc2 == -1 || i_barc3 == -1){
        
        return 0;
    }
    
    //  The index is in the sequence starting at 0, so the whole index is offset to the right by 1
    i_barc1++;
    i_barc2++;
    i_barc3++;
    
    //printf("barcode 1 2 3 position: %d\t%d\t%d\n",i_barc1,i_barc2,i_barc3);
    fastq->info[0] = ' ';
    
    sprintf(info, "@%d_%d_%d:%s%s",
            i_barc1,
            i_barc2,
            i_barc3,
            umi_buffer,
            fastq->info);
    
    /*    
	    e.g.: @21_54_85:32141 :
        UMI is ATCTCGGCAT
        round1 barcode is 21
        round2 barcode is 54
        round3 barcode is 85
	*/
    strcpy(fastq->info, info);
    strcpy(fastq->data, cDNA);
    strCopy(info, fastq->sanger, len - strlen(cDNA), len);
    strcpy(fastq->sanger, info);
    
    
    return 1;
}

void generateFeatureTable(){

	/* generate feature table */

    for(int i = 0; i < FEATURE_TABLE_SIZE; i++){
        feature32_tab[i] = -1;
        feature21_tab[i] = -1;
    }
    unsigned long len = strlen(feature32);                                             // length of feature32 and feature21 is equal
    char buffer[KMODEL + 1];
    for(int i = 0; i < len - KMODEL + 1; i++){
        strCopy(buffer, feature32, i, i+KMODEL);
        feature32_tab[decodeSqe(buffer)] = i;
        strCopy(buffer, feature21, i, i+KMODEL);
        feature21_tab[decodeSqe(buffer)] = i;
    }
}

void generatePrimerTable(){

	/* genetate primer table */

    for(unsigned int i = 0; i < PRIMERLEN; i++){
        initChain(&primerTable[i], -1, 0);
    }
    char primer[20000];
    
    getPrimerString(primer, PRIMER_LIST_FILE);
    
    long len = strlen(primer);
    
    char buffer[KMODEL];
    for(unsigned int i = 0; i < len - KMODEL + 1 ; i++){
        strCopy(buffer, primer, i, i+KMODEL);
        addNode(primerTable[decodeSqe(buffer)], i, 0);
        
    }

    return ;
    
}


void getIndexInfo(char* buffer, char* info){
    /*  get index information, e.g.: @23_74_12  */
    int len = strlen(info);
    for(int i = 0; i < len; i++){
        if(info[i] == ' '){
            strCopy(buffer, info, 0, i);
            break;
        }
    }
}




int generateBarcTable(char* filename){
    /*  
	    generated Barcode Table
     
        barc_ro1_tab
        barc_ro2_tab
        barc_ro3_tab
	*/
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
    
    
    while (start < end) {
        pthread_mutex_lock(&read_mutex);
        getRead(&readr2, r2_file);
        getRead(&readr1, r1_file);
        start = ftell(r2_file);
        pthread_mutex_unlock(&read_mutex);
        
        int flag = transR2(&readr2);
        
        if(flag == 1){
            
            transR1(&readr1);
            char buffer_1[20];
            char buffer_2[210];
            getIndexInfo(buffer_1, readr2.info);
            strcpy(buffer_2, readr1.info);
            
            
            buffer_2[0] = ' ';
            sprintf(readr1.info, "%s%s", buffer_1, buffer_2);
            pthread_mutex_lock(&write_mutex);
            write2file(readr1, out_r1_file);
            write2file(readr2, out_r2_file);
            pthread_mutex_unlock(&write_mutex);
        }
        
        
        printf("start: %lld, end: %lld\n", start, end);
        printf("%.2lf%% finishing...\n", ((double)start)/end*100);
    }
}

void printHelp(){
    printf("Program: SCSit (A high-efficiency preprocessing tool for single-cell sequencing data from SPLiT-seq)\n");
    printf("Version: %s\n", PACKAGE_VERSION);
    printf("Contact: Shangqian Xie <sqianxie@foxmail.com>\n\n");
    printf("Usage: SCSit [options] -r1 input_r1.fastq -r2 input_r2.fastq -p primer.list -b barcode.list -o output \n\n");
    printf("\nInput/Output :\n\n");
    printf("       -r1 R1File        Paths to files that contain input read1 of SPLiT-seq pair-end files <str>\n");
    printf("       -r2 R2File        Paths to files that contain input read2 of SPLiT-seq pair-end files <str>\n");
    printf("       -p primerList     Primer list of all oligonucleotide sequences used <str>\n");
    printf("       -b barcodeList    The 96 well plate oligonucleotides used for each round of barcodes <str>\n");
    printf("Options:\n\n");
    printf("       -t int            Maximum number of threads to use [default: 2]\n");
    printf("       -o output_prefix  Paths to files that contain output file prefix [default: output]\n");
    printf("       -h                show this help message and exit\n");
    
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
        printf("[Error] R1 file does not exist.\n");
        return 1;
    }
    r2_file = fopen(R2FILE, "r");
    if(r2_file == NULL){
        printf("[Error] R2 file does not exist.\n");
        return 2;
    }
    
    if(BARC_LIST_FILE == NULL){
        printf("[Error] Barcode list is missing.\n");
    	return 1;
    }

    if(PRIMER_LIST_FILE == NULL){
    	printf("[Error] Primer list is missing.\n");
    	return 1;
    }
    
    for(int i=0; i<sizeof(OUTPREFIX); i++ ){
        if(OUTPREFIX[i] == '\0'){
            if (OUTPREFIX[i-1] == '/') {
                strcat(OUTPREFIX,"output");
            }
            break;
        }
        
    }
    //printf("%s\n",OUTPREFIX);
        

    char OUTR1FILE[128]={} ;
    strcpy(OUTR1FILE,OUTPREFIX);
    //printf("%s\n",OUTR1FILE); 
    strcat(OUTR1FILE,".R1.fastq");
    
    char OUTR2FILE[128]={} ;
    strcpy(OUTR2FILE,OUTPREFIX);
    //printf("%s\n",OUTR2FILE); 
    strcat(OUTR2FILE,".R2.fastq");
        
    
    //char* OUTR2FILE = OUTPREFIX;
    //strcat(OUTR2FILE,".R2.fastq");
    //printf("%s\n",OUTR1FILE);  
    //printf("%s\n",OUTR2FILE);  
  
    out_r1_file = fopen(OUTR1FILE, "w");
    out_r2_file = fopen(OUTR2FILE, "w");
    
   
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
    
    return 0;
}
