//
//  permutation.c
//  ucmap
//
//  Created by Fenix Huang on 5/18/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#include "permutation.h"
#include "uc_utility.h"
#include "ucmap.h"


int tour[100];
int sigma[100];
int alpha[100];

int label[100][3];
int label_trisection;
int depth=0;
int UD[100];
int cnt=0, over_cnt=0,cnt_origin=0;
int tri;


uc_map *permuttaion_to_ucmap (int *signed_p)
{
    uc_map *s;
    int i, n=signed_p[0];
    
    s=newmap(2*n);
    
    
    for (i=2; i<2*n; i+=2) {
        s->alpha[i]=i+1;
        s->alpha[i+1]=i;
    }
    s->alpha[1]=2*n;
    s->alpha[2*n]=1;
    
    for (i=1;i<=n;i++) {
        if (signed_p[i]>0) {
            s->tour[2*signed_p[i]-1] = 2*signed_p[i];
            if (signed_p[i+1]>0) {
                s->tour[2*signed_p[i]] = 2*signed_p[i+1]-1;
            } else {
                s->tour[2*signed_p[i]] = -2*signed_p[i+1];
            }
        } else {
            s->tour[-2*signed_p[i]] = -2*signed_p[i]-1;
            if (signed_p[i+1]>0) {
                s->tour[-2*signed_p[i]-1] = 2*signed_p[i+1]-1;
            } else {
                s->tour[-2*signed_p[i]-1] = -2*signed_p[i+1];
            }
        }
    }
    
    s->tour[2*signed_p[n]] = 2*signed_p[1]-1;
    
    circ(s->tour, s->alpha, s->sigma);
    
    countbc(s);
    return s;
    
}

//convert a permutation, signed or unsiged to a ucmap

uc_map *relabel_from_root(uc_map *s, int root)
{
    uc_map *res;
    int *map, l, i, k;
    
    l=s->alpha[0];
    
    res=newmap(l);
    
    res->tour[0]=l;
    for (i=1;i<l;i++) res->tour[i]=i+1;
    res->tour[l]=1;
    
    res->alpha[0]=l;
    res->alpha[1]=l; res->alpha[l]=1;
    
    map=(int *) calloc (l+10, sizeof(int));
    
    k=1;
    i=1;
    do {
        map[k]=i;
        i++;
        k=s->tour[k];
    } while (k!=1);
    
    
    for (i=2; i<l; i++) {
        res->alpha[map[i]] = map[s->alpha[i]];
    }
    
    circ(res->tour, res->alpha, res->sigma);
    countbc(res);
    
    free(map);
    return res;
}

int odd_length_of_bc (uc_map *s)
{
    int i, k, n=s->alpha[0], bc_length;
    
    for (i=2; i<=n; i+=2) {
        bc_length = 0;
        k=i;
        do {
            k=s->sigma[k];
            bc_length++;
        } while (k!=i);
        
        if (bc_length % 2 == 0) return 0;
    }
    return 1;
}

int odd_length_tri_of_bc (uc_map *s)
{
    int i, k, n=s->alpha[0], bc_length;
    
    for (i=2; i<=n; i+=2) {
        bc_length = 0;
        k=i;
        do {
            k=s->sigma[k];
            bc_length++;
        } while (k!=i);
        
        if (i!=n && bc_length != 3) return 0;
    }
    return 1;
}



void unsigned_permutation(int n, int step)
{
    int i;
    int ex_tri_cnt;
    uc_map *s1 ,*s2;
    
    if (step>n+1) {
        
        
        alpha[0]=n+2;
        alpha[1]=1;
        alpha[n+2]=n+2;
        
        ex_tri_cnt=0;
        
        s2=permuttaion_to_ucmap(alpha);
        s1=relabel_from_root(s2, s2->alpha[0]);
        
        for (i=1; i<=s1->trisection[0]; i++) {
            if (s1->trisection[i] % 2==0) ex_tri_cnt++;
        }
        
        if (odd_length_tri_of_bc (s1) ) {
            //          if (shape(s1)) {
            
            
            //  	    printf(">genome%d\n", cnt+1);
            //  	    for (i=1;i<=n+2;i++)
            //  	      printf("%d ", alpha[i]);
            //  	    printf("\n");
            //  	    printf("%d ", s1->genus);
            //   	    printcycle(s1->tour);
            
            
            for (i=1; i<=s1->trisection[0]; i++) {
                if (s1->trisection[i] % 2==0)
                    // 		printf("%d ", s1->trisection[i]);
                    ex_tri_cnt++;
            }
            
            if (ex_tri_cnt==0) {
                for (i=2; i<=n+1; i++) printf("%d ", alpha[i]-1);
                printf("\n");
                printf("%d ", s1->genus);
                printcycle(s1->sigma);
                printcycle(s1->alpha);
                printf("\n");
                cnt++;
            }
            
            freemap(s1);
            freemap(s2);
        }
        
        return;
    }
    for (i=2;i<=n+1;i++) { 
        if (UD[i]==0) { 
            UD[i]=1;
            
            alpha[step]=i;
            unsigned_permutation(n, step+1);
            alpha[step]=0;
            
            UD[i]=0;
        }
    }
    return;
}
