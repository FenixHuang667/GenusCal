//
//  uc_map_generator.c
//  ucmap
//
//  Created by Fenix Huang on 5/15/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#include "uc_map_generator.h"
#include "ucmap.h"
#include "uc_utility.h"


double *E[MAXG]; // the number of uc_map for the current genus of size n.
double bc[MAXN][MAXN]; // binomial coefficients, (n \choose k)
double Secondary[MAXN]; // number of secondary structure of length n (including unpaired bases)
float Trans[MAXG][MAXG]; // Transition probability
int frequency[10000];
int fre_loop_standard[10000];
int fre_loop_pk[10000];
int fre;

slist *List;

// initial start
// ------------------------------------------------------------
void initial_matrix (int n)
{
    int i;
    for (i=0;i<MAXG;i++) {
        E[i]=(double *) calloc (n+10,sizeof(double));
    }
}
//initial the E matrix

double compute_binomial (int n)
{
    int i, j;
    for (i=0;i<=n;i++) {
        bc[i][0]=1;
        bc[i][i]=1;
    }
    
    for (i=1;i<=n;i++) {
        for (j=1;j<i;j++) {
            bc[i][j] = bc[i-1][j-1] + bc[i-1][j];
        }
    }
    return 0;
}
//compute the binomial coefficients using a 2-dim table


void Catalan(int n)
{
    int i,j;
    E[0][0]=1;
    for (i=1;i<=n;i++) {
        E[0][i]=0;
        for (j=0;j<i;j++) {
            E[0][i]+=E[0][j]*E[0][i-j-1];
        }
    }
    return;
}
// compute the number of secondary structure/tree/uc_map with genus 0.

void intitial_compute(int n) //initial the number of uc-map
{
    int i,j;
    Catalan(n);
    for (i=1;i<MAXG;i++) {
        E[i][n]=0;
        if (n<2*i) continue;
        for (j=0; j<i; j++) {
            E[i][n]+=bc[n+1-2*j][2*(i-j)+1]*E[j][n];
        }
        E[i][n]=E[i][n]/(2*i);
    }
    return;
}
//Fill E. Compute the number of uc_maps with genus > 0

void compute_transition_matrix (int n, int target_g)
{
    int i,j;
    float scale;
    
    for (i=0;i<MAXG;i++) {
        for (j=0;j<MAXG;j++) {
            Trans[i][j]=0;
        }
    }
    
    for (i=target_g; i>0; i--) {
        scale=0;
        if (i==target_g) scale=1;
        else {
            for (j=target_g; j>i; j--) scale+=Trans[j][i];
        }
        for (j=i-1; j>=0; j--) {
            Trans[i][j]+=scale*bc[n+1-2*j][2*(i-j)+1]*E[j][n]/(2*i*E[i][n]);
        }
    }
    
    /*
     for (i=0; i<target_g; i++) {
     for (j=i+1; j<=target_g; j++) {
     printf("%f  ", Trans[j][i]);
     }
     printf("\n");
     }
     */ //debug
    
    return;
}
// Compute the transition probability for selecting the next genus. For detail please read the paper.


void initial_secondary_structure (int l) // l is the length
{
    int i, j;
    
    for (i=0;i<=l;i++) Secondary[i]=0;
    Secondary[0]=1;
    Secondary[1]=0; //disallow unpaired vertices
    
    for (i=2;i<=l;i++) {
        //     Secondary[i]+=Secondary[i-1]; // allow unpaired vertices
        for (j=0; j<=i-2; j++) {
            Secondary[i]+=Secondary[j]*Secondary[i-2-j];
        }
    }
    //   for (i=0;i<=l;i++) printf("%lf\n", Secondary[i]);
    return;
}
//Compute the number of secondary structures (with unpaired bases) recursively.


// ------------------------------------------------------------
// initial end


void precompute(int l)
{
    int n=l/2, i;
    
    initial_matrix(n+1);
    compute_binomial(l+3);
    initial_secondary_structure(l);
    
    for (i=1;i<=n;i++) {
        intitial_compute(i);
    } 
}


int determin_genus (int g, int target_g, int n)
{
    int new_g=0, i;
    float I_array[MAXG];
    float seed;
    
    I_array[g]=0;
    for (i=g+1; i<=target_g; i++) {
        I_array[i]=I_array[i-1]+Trans[i][g];
        //     printf("%f\n", I_array[i]); //debug
    }
    
    seed=(float)((1-drand48())*I_array[target_g]);
    
    for (i=g+1; i<=target_g; i++) {
        if (seed>I_array[i-1] && seed<=I_array[i]) new_g=i;
    }
    
    return new_g;
}
//Compute the next genus via traisition probability


int select_vertex(int n) // select a labeled vertex from n with uniform probability
{
    float seed;
    int i, s=0;
    
    seed=(float)((1-drand48())*n);
    
    //   printf("%f ", seed);
    
    for (i=0;i<n;i++) {
        if (seed > (float) i && seed <= (float) i+1 ) {
            s=i;
        }
    }
    return s+1;
}

int *label_vertex(uc_map *s, int number) //pcik up n=2k+1 vertices from uc-map s with uniform probability
{
    int *label, *res, *choose_set;
    int n, i, n_V, k, j, l, choice;
    
    n=s->vertex;
    l=s->alpha[0];
    
    label=(int *) calloc (l+10, sizeof (int));
    res=(int *) calloc (n+10, sizeof (int));
    choose_set=(int *) calloc (n+10, sizeof (int));
    
    for (i=0;i<=n;i++) label[i]=0;
    
    n_V=n-1;
    for (i=0; i<number; i++) {
        k=1;
        j=1;
        while (j<l) {
            if (s->min[j]==1 && label[j]==0) {
                choose_set[k]=j;
                k++;
            }
            j=s->tour[j];
        }
        
        choice=select_vertex(n_V);
        label[choose_set[choice]]=1;
        n_V--;
    }
    
    res[0]=number;
    k=1;
    for (j=1;j<=l;j++) {
        if (label[j]==1) {
            //       printf("%d ", j); //debug
            res[k]=j;
            k++;
        }
    }
    //   printf("\n"); //debug
    
    return res;
} 

int loop_statistic(uc_map *s)
{
    int *used,i,l=s->sigma[0],bc=0, start, bc_length, pk;
    used=(int *) calloc (l+2,sizeof(int));
    for (i=1;i<=l;i++) {
        used[i]=1;
    }
    start=1;
    do {
        i=s->sigma[start];
        bc_length=1;
        pk=0;
        
        while (i!=start) {
            used[i]=0;
            bc_length++;
            if (s->sigma[i]<i && s->sigma[i]!=start) pk=1;
            i=s->sigma[i];
        }
        
        if (pk) {
            fre_loop_pk[bc_length]++;
        } else {
            fre_loop_standard[bc_length]++;
        }
        
        used[start]=0;
        
        start=s->tour[start];
        while (used[start]!=1 && start!=1) start=s->tour[start];
    } while (start!=1);
    
    free(used);
    return bc;
    
}

// Note!! Remeber to initiate the array fre...


int *assemble (int *initial, int target_g) // satrt from genus 0 to target_g
{
    int n, l, i, genus, next_genus; //n: #edges l:length
    int *weight, *result, *label;
    uc_map *s1, *s2, *t, *s3;
    
    
    l=initial[0];
    n=0;
    for (i=1; i<=l; i++) if (initial[i]!=0) n++;
    if (n<4*target_g) {
        printf("Not enough edges to form genus %d structure\n", target_g);
        return NULL;
    }
    
    compute_transition_matrix(n/2, target_g);
    //initialization
    
    weight=(int *) calloc (l+10, sizeof(int));
    //apply array
    
    
    s1=structure_to_map(initial, weight);
    genus=s1->genus;
    while (genus<target_g) {
        next_genus=determin_genus(genus, target_g, n/2);
        //     printf("Next Genus: %d\n", next_genus); //debug
        
        label=label_vertex(s1, 2*(next_genus-genus)+1);
        s2 = I_glue(s1, label);
        s3=relabel(s2, weight);
        
        //     printmap(s3); //debug
        
        t=s1;
        s1=s3;
        freemap(t);
        freemap(s2);
        free(label);
        genus=next_genus;
    }
    
    //   printmap(s1); //debug
    
    loop_statistic(s1);  //statisc of the result
    
    result = map_to_structure(s1, weight);
    
    return result;
    
}


// external : secondary structure sampler

void determin_arc(int i, int length, int *array)
{
    int k;
    double seed, a1, a2;
    
    if (array[i]!=0 || length<1) return;
    seed=(double)((1-drand48()))*Secondary[length];
    //   if (seed<=Secondary[length-1]) {
    //     determin_arc(i+1, length-1, array); //unpaired vertices
    //     return;
    //   }
    //   a1=Secondary[length-1]; //allow unpaired vertices
    a1=0; // disallow unpaired vertices
    for (k=0; k<=length-2; k++) {
        a2=a1+Secondary[k]*Secondary[length-2-k];
        if (seed>a1 && seed<=a2) {
            array[i]=i+k+1;
            array[i+k+1]=i;
            determin_arc(i+1, k, array);
            determin_arc(i+k+2, length-2-k, array);
            return;
        }
        a1=a2;
    }
}


int *tree_sampler (int length) //without unpaired vertices
{
    int i, *res;
    
    res=(int *) calloc (length+10, sizeof(int));
    
    res[0]=length;
    for (i=1;i<=length;i++) res[i]=0;
    
    determin_arc(1, length, res);
    
    //   for (i=1;i<=length;i++) printf("%d ", res[i]); printf("\n"); //debug
    
    //   printtree(res); // deubug
    return res;
    
}

int determin_unpaired (int length, int target_g)
{
    int i, res;
    double seed, delta=0, a1, a2;
    
    for (i=2*target_g; i<=length/2; i++) {
        delta+=bc[length][length-2*i]*E[target_g][i];
    }
    //   printf("%lf\n", delta); //debug
    
    seed = (double)((1-drand48()))*delta;
    
    a1=0;
    for (res=2*target_g; res<=length/2; res++) {
        a2=a1+bc[length][length-2*res]*E[target_g][res];
        if (a1 < seed && seed <= a2) return (length-2*res);
        a1=a2;
    }
    
    return 0;
}


int *pick_unpaired (int *matching, int up)
{
    int *res, i, l, p, j ,t, *position;
    
    l=matching[0]+up;
    
    res=(int *) calloc (l+10, sizeof(int));
    position=(int *) calloc (l+10, sizeof(int));
    res[0]=l;
    
    for (i=1;i<=l;i++) position[i]=i;
    for (i=1;i<=l;i++) res[i]=1;
    
    
    for (i=1;i<=up;i++) {
        p=select_vertex(l-i+1);
        res[position[p]]=0;
        for (j=p;j<l-i+1;j++)
            position[j]=position[j+1];
    }
    
    i=1;
    t=0;
    j=1;
    while (i<=matching[0]) {
        while (res[j]==0) {
            j++;
            t++;
        }
        position[i]=t;
        i++;
        j++;
    }
    
    i=1;
    j=1;
    while (i<=matching[0]) {
        while (res[j]==0) j++;
        res[j]=matching[i]+position[matching[i]];
        i++;
        j++;
    }
    
    //   for (i=1;i<=l;i++) printf("%d ", res[i]); printf("\n");
    //   for (i=1;i<=matching[0];i++) printf("%d ", position[i]); printf("\n");
    
    free(matching);
    return res;
}


int *G_sampler (int length, int target_g)
{
    int *struc, *struc_g, up;
    
    up=determin_unpaired(length, target_g);
    
    initial_secondary_structure(length-up);
    
    struc = tree_sampler(length-up);
    
    struc_g = assemble(struc, target_g);
    
    struc_g = pick_unpaired(struc_g, up);
    
    //   for (i=1; i<=struc_g[0]; i++) printf("%d ", struc_g[i]); printf("\n"); //debug
    
    free(struc);
    return struc_g;
}


int Uniform_sampling (int n, int target_g)
{
    int i, j, cnt;
    int *result;
    slist *index;
    
    srand48(time(NULL)); //random seed
    
    //   printf("%lf\n", E[4][n]);
    cnt=0;
    List=NULL;
    
    for (i=0;i<10000;i++) {
        frequency[i]=0;
        fre_loop_pk[i]=0;
        fre_loop_standard[i]=0;
    }
    
    precompute(n);
    
    
    
    /* debug
     initial_matrix(n); // intiial E matrix
     compute_binomial(n+1); //compute the binomial coeeficients
     intitial_compute(n); //compute the E matrix
     
     s1=structure_to_map(struc, weight);
     printmap(s1);
     
     for (i=0;i<=weight[0];i++) printf("%d ", weight[i]); printf("\n"); //debug
     
     label=label_vertex(s1, 5);
     for (i=1;i<=label[0];i++) printf("%d ", label[i]);
     printf("\n");
     
     s2=I_glue(s1, label);
     
     printmap(s2);
     
     s3=relabel(s2, weight);
     printmap(s3);
     
     for (i=0;i<=weight[0];i++) printf("%d ", weight[i]); printf("\n"); //debug
     
     map_to_structure(s3, weight);
     */
    
    
    for (i=1; i<=100000; i++) {
        result=G_sampler(n,target_g);
        //     result=tree_sampler(n);
        //     shape = structure_to_shape(result);
        compare_new(result, &List);
        free(result);
        
    }
    
    for (i=2*target_g; i<=6*target_g-1; i++) {
        cnt=0;
        fre=0;
        index=List;
        for (j=0; j<1000; j++) frequency[j]=0;
        while (index!=NULL) {
            if (index->cnt>=1 && index->pair[0] == 2*i) {
                frequency[index->cnt]++;
                fre+=index->cnt;
                //      for (i=1;i<=index->pair[0];i++) printf("%d ", index->pair[i]);
                //         printf("%d,", index->cnt);
                cnt++;
            }
            index=index->Next;
        }
        printf("\n");
        for (j=1; j<=200;j++)printf("%d,", frequency[j]);
        printf("\nCnt: %d \n", cnt);
        printf("Fre: %d \n", fre);
    }
    
    /*
     for (i=0; i<100; i++) printf("%d,", fre_loop_standard[i]);
     printf("\n");
     for (i=0; i<100; i++) printf("%d,", fre_loop_pk[i]);
     printf("\n");
     printf("%d \n", cnt);
     
     */
    //   result=assemble(struc, 3);
    
    //   for (i=0;i<1000;i++) {
    //     genus=determin_genus(2, 5, n);
    //     genus=select_vertex(10);
    //     printf("%d ", genus);
    //     if (genus == 2) cnt++;
    //   }
    
    
    
    return 0;
}


