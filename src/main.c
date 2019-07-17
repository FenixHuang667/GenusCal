//
//  main.c
//  ucmap
//
//  Created by Fenix Huang on 5/7/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "uc_map_generator.h"
#include "ucmap.h"
#include "uc_utility.h"

#include "permutation.h"

char Inputfilename[500];
char Outputfilename[500];
char **Infile;
char *Outfile;
FILE *Input;
FILE *Output;


void usage()
{
    printf("usage:\n"
           "Gcal [-OPTIONS}\n"
           "[-i inputfilename] Specify the input file (default input.in)\n"
           "[-o outputfilename] Specify the output file (default output.out)\n"
           "[-n] Use dot/bracket input form (default enable)\n"
           );
    exit(0);
}


int main(int argc, const char * argv[]) {
    
    int *result, weight[1000], label[1000], i, j, StringInput = 1, length, *ptable, NumInput = 0;
    uc_map *map;
    char structure[1000];
    
    srand48(time(NULL)); //random seed
    
    strcpy(Inputfilename, "./input.in");
    strcpy(Outputfilename, "./output.out");
    
    // Uniform_sampling (16,2);
    
    
    for (i=1; i<argc; i++) {
        if (argv[i][0]=='-') {
            switch (argv[i][1]) {
                case 'i':
                    if (i==argc-1) usage();
                    Infile = argv[++i];
                    strcpy(Inputfilename, Infile);
                    break;
                case 'o':
                    if (i==argc-1) usage();
                    Outfile = argv[++i];
                    strcpy(Outputfilename, Outfile);
                    break;
                case 'n':
                    StringInput = 0;
                    break;
                default: usage();
            }
        }
    }
    
    Input=fopen(Inputfilename,"r");
    if (Input==NULL) {
        printf("Input file error!\n");
        exit(0);
    }
    
    Output=fopen(Outputfilename,"w");
    if (Output==NULL) {
        printf("Output file error!\n");
        exit(0);
    }
    
    while(!feof(Input)) {
        if (StringInput) {
            fscanf(Input, "%s", structure);
            
            //printf("%s\n", structure);
            //debug input
            
            ptable = structure2pair(structure);
            
            //for (i=0; i<=ptable[0]; i++)
            //    printf("%d ", ptable[i]);
            //printf("\n");
             
            //debug input
            map = structure_to_map(ptable, weight);
            // printf("Genus: %d\n", map->genus);
            fprintf(Output, "%s Genus: %d\n", structure, map->genus);
            
            freemap(map);
            free(ptable);
            
        } else {
            fscanf(Input, "%d", &length);
            ptable = (int *) calloc (length + 5, sizeof(int));
            ptable[0] = length;
            for (i=1; i<length; i++) {
                fscanf(Input, "%d ", &ptable[i]);
            }
            
            for (i=0; i<=ptable[0]; i++)
                printf("%d ", ptable[i]);
            printf("\n");
            //debug input
        }
        NumInput ++;
        
    }
    
    fclose(Input);
    fclose(Output);
    
    printf("done!\n");
    
    return 0;
}
