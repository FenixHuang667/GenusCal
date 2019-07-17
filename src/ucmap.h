//
//  ucmap.h
//  ucmap
//
//  Created by Fenix Huang on 5/18/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#ifndef ucmap_ucmap_h
#define ucmap_ucmap_h

#define MAXG 10
#define MAXN 1000

//Data structure for uc_map.

typedef struct uc_map uc_map;
struct uc_map
{
    int *tour; // boundary component
    int *sigma; // vertices
    int *alpha; // edges (not
    //These are three permutation presenting a unicellular map
    int *trisection; // number of trisection
    int *min; // rememeber the minimum sector/half-edge in each vertex
    int genus; // store the genus of the uc_map. g = (e+1-v)/2 for orienated and g = (e+1-v) for non-orienated
    int vertex; // number of vertices
};

//Structure of a uc_map (oriented or non-oriented)




typedef struct slist slist;
struct slist
{
    int *pair;
    int cnt;
    slist *Next;
};
// List of structure in bp style

#endif