//
//  uc_utility.h
//  ucmap
//
//  Created by Fenix Huang on 5/14/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

//Description:
//This is a pacakge for functions which deal with unicellular map
//
// To be added..




//Unicellular map data structure

//Converting a structure (in base pair style) to a unicellular map

//Converting a unicellular map to a structure (in bp style (p[i] pair p[j] then p[i]=j and p[j]=i. p[0] is the length)

//Find trisection

//Glue process


#ifndef ucmap_uc_utility_h
#define ucmap_uc_utility_h


#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "ucmap.h"


extern uc_map* newmap (int l);
// initial a new uc_map structure
extern void freemap(uc_map *s);
//free an uc_map structure
extern void circ(int *a, int *b, int *result);
// apply b \circ a
// Note!! be aware of the order
extern uc_map *copy_uc(uc_map *s);
//copy an uc_map
extern void printmap(uc_map *s);
extern void printstructure(int *s);
extern void printcycle(int *a);



extern int countbc(uc_map *s); 
//This function is to identify boundary components in an uc_map.
//bc is the number of boundary component
//By computing bc, the genus can also be computed in this function
extern slist *compare_new(int *p, slist **L);
//Check whether the input uc_map is already in the list
extern uc_map *relabel(uc_map *s, int *weight);
// relabel the uc-map as canonical backbone, need also consider the weight array
extern uc_map *I_glue(uc_map *s, int *label); 


extern uc_map *structure_to_map (int *struc, int *weight);
// map a structure to a uc-map with weights
extern int *map_to_structure(uc_map *s, int *weight);
//from uc-map to structure



extern void n_slice_donw (uc_map *input, int * label);



extern char *pair2structure(int *pair);
extern int *structure2pair(char *s);
