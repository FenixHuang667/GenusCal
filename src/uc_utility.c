//
//  uc_utility.c
//  ucmap
//
//  Created by Fenix Huang on 5/14/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#include "uc_utility.h"
#include "ucmap.h"

uc_map* newmap (int l)
{
    uc_map *result;
    int i;
    
    result=(uc_map *) calloc (1,sizeof(uc_map));
    result->tour=(int *) calloc (l+10,sizeof(int));
    result->tour[0]=l;
    result->sigma=(int *) calloc (l+10,sizeof(int));
    result->sigma[0]=l;
    result->alpha=(int *) calloc (l+10,sizeof(int));
    result->alpha[0]=l;
    result->trisection=(int *) calloc (l+10,sizeof(int));
    result->trisection[0]=0;
    result->min=(int *) calloc (l+10,sizeof(int));
    result->genus = 0;
    result->vertex = 0;
    for (i=0;i<=l;i++) result->min[i]=0;
    return result;
}

// initial a new uc_map structure

void freemap(uc_map *s)
{
    free(s->tour);
    free(s->sigma);
    free(s->alpha);
    free(s->trisection);
    free(s->min);
    free(s);
}

//free an uc_map structure



void circ(int *a, int *b, int *result)
{
    int l, i;
    l=a[0];
    result[0]=l;
    for (i=1;i<=l;i++) {
        result[i]=b[a[i]];
    }
}
// apply b \circ a
// Note!! be aware of the order


void printcycle(int *a)
{
    int l, i, n, *used;
    
    n=a[0];
    
    used=(int *) calloc (n+2,sizeof(int));
    for (i=0;i<=n;i++) used[i]=0;
    l=1;
    while (l<=n)
    {
        printf("(%d", l);
        used[l]=1;
        i=a[l];
        used[i]=1;
        while (i!=l) {
            printf(",");
            printf("%d",i);
            used[i]=1;
            i=a[i];
        }
        printf(")");
        while (used[l]!=0) l++;
    }
    
    printf("\n");
    free(used);
}
// for debug purpose, display the cycle structure of a permutation.



void printmap(uc_map *s)
{
    int i;
    printf("Length: %d\n", s->alpha[0]);
    printf("Tour: "); printcycle(s->tour);
    printf("Sigma: "); printcycle(s->sigma);
    printf("Alpha: "); printcycle(s->alpha);
    printf("Trisections: ");
    for (i=1;i<=s->trisection[0];i++)
        printf("%d ", s->trisection[i]);
    printf("\n");
    printf("Genus: %d\n", s->genus);
    printf("#Vertex: %d\n", s->vertex);
}
// for debug purpose, display the structure of a uc_map
// Note!! This function may be changed frequently for experimental purpose





void copy(int *a, int *b) // copy an array
{
    int i;
    for (i=0;i<=a[0];i++) {
        b[i]=a[i];
    }
}

uc_map *copy_uc(uc_map *s)
{
    uc_map *res;
    int i, l;
    
    l=s->alpha[0];
    
    res=newmap(l);
    
    copy(s->alpha, res->alpha);
    copy(s->tour, res->tour);
    copy(s->sigma, res->sigma);
    copy(s->trisection, res->trisection);
    
    for (i=0;i<=l;i++) {
        res->min[i]=s->min[i];
    }
    
    res->genus=s->genus;
    res->vertex=s->vertex;
    return res;
}
//copy an uc_map


int find_trisection(uc_map *s, int a)
{
    int b,i,l, start;
    b=s->sigma[a];
    l=s->sigma[0];
    if (a==b) return 0;
    
    start=1;
    i=s->tour[start];
    while (i!=start) {
        if (i==a) return 0;
        if (i==b && s->min[b]!=1) return 1;
        i=s->tour[i];
    }
    return 0;
}
//Identify trisections in an uc_map


int countbc(uc_map *s)
{
    int *used,i,l=s->sigma[0],bc=0, start;
    used=(int *) calloc (l+2,sizeof(int));
    for (i=1;i<=l;i++) {
        used[i]=1;
    }
    s->trisection[0]=0;
    start=1;
    s->min[1]=1;
    do {
        i=s->sigma[start];
        while (i!=start) {
            if (find_trisection(s, i)) {
                //  	printf("tri: %d\n", i);
                s->trisection[0]+=1;
                s->trisection[s->trisection[0]]=i;
            }
            used[i]=0;
            i=s->sigma[i];
        }
        used[start]=0;
        bc+=1;
        
        start=s->tour[start];
        while (used[start]!=1 && start!=1) start=s->tour[start];
        s->min[start]=1;
        s->vertex++;
    } while (start!=1);
    
    s->genus=(1-bc+l/2)/2;
    
    free(used);
    return bc;
}

//This function is to identify boundary components in an uc_map.
//bc is the number of boundary component
//By computing bc, the genus can also be computed in this function




int same(int *a, int *b)
{
    int i;
    if (a[0]!=b[0]) return 0;
    for (i=1;i<=a[0];i++)
        if (a[i]!=b[i]) return 0;
    return 1;
}
//If two array is the same return 1 otherwise return 0


slist *compare_new(int *p, slist **L)
{
    slist *index;
    index=*L;
    // 	if (index->pair==NULL) index=NULL;
    while (index!=NULL) {
        if (same(index->pair,p)) {
            index->cnt++;
            return NULL;
        }
        index=index->Next;
    }
    index=(slist *) calloc (1,sizeof(slist));
    index->pair=(int *) calloc (p[0]+10,sizeof(int));
    copy(p, index->pair);
    index->cnt=1;
    index->Next=*L;
    *L=index;
    return index;
}
//Check whether the input uc_map is already in the list


uc_map *glue (uc_map *s, int a1, int a2, int a3) // gluing operation
{
    
    int l=s->sigma[0],a,b,c,t, i;
    uc_map *res;
    
    i=1;
    while (i!=a1 && i!=a2 && i!=a3) i=s->tour[i];
    a=i;
    i=s->tour[i];
    while (i!=a1 && i!=a2 && i!=a3) i=s->tour[i];
    b=i;
    i=s->tour[i];
    while (i!=a1 && i!=a2 && i!=a3) i=s->tour[i];
    c=i;
    
    res=newmap(l);
    
    copy(s->alpha, res->alpha);
    copy(s->sigma, res->sigma);
    
    t=res->sigma[a];
    res->sigma[a]=res->sigma[b];
    res->sigma[b]=res->sigma[c];
    res->sigma[c]=t;
    
    circ(res->sigma, res->alpha, res->tour);
    countbc(res);
    
    return res;
    
}
// glue three vertex in an uc_map. The three vertices are labeld by their min half-edge.
// a3 can be a trisection


uc_map *I_glue(uc_map *s, int *label)
{
    int i, number, k;
    uc_map *s1, *s2, *t;
    
    number=label[0];
    k=(number-1)/2;
    
    s1=copy_uc(s);
    
    for (i=k; i>0; i--) {
        s2=glue(s1, label[2*i-1], label[2*i], label[2*k+1]);
        t=s1;
        s1=s2;
        freemap(t);
    }
    
    return s1;
}
// interativly glue vertices in label array

uc_map *relabel(uc_map *s, int *weight)
// relabel the uc-map as canonical backbone, need also consider the weight array
{
    uc_map *res;
    int *map, l, i, k, *new_weight;
    
    l=s->alpha[0];
    
    res=newmap(l);
    
    res->tour[0]=l;
    for (i=1;i<l;i++) res->tour[i]=i+1;
    res->tour[l]=1;
    
    res->alpha[0]=l;
    res->alpha[1]=l; res->alpha[l]=1;
    
    map=(int *) calloc (l+10, sizeof(int));
    new_weight=(int *) calloc (l+10, sizeof(int));
    
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
    for (i=1; i<=l; i++){
        new_weight[map[i]] = weight[i];
    }
    for (i=1; i<=l; i++) {
        weight[i]=new_weight[i];
    }
    
    circ(res->tour, res->alpha, res->sigma);
    countbc(res);
    
    free(map);
    free(new_weight);
    return res;
}


void printstructure (int *s)
{
    int i;
    for (i=1; i<=s[0]; i++) {
        printf("%d, ", s[i]);
    }
    printf("\n");
}

// map a structure to a uc-map with weights
uc_map *structure_to_map (int *struc, int *weight)
{
    int l=struc[0], i, n=0, *map, unpair;
    uc_map *res;
    
    map = (int *) calloc (l+10, sizeof(int));
    
    n=1;
    unpair=0;
    
    for (i=1;i<=l;i++) {
        if (struc[i]==0) unpair++;
        else if (struc[i]!=0) {
            weight[n]=unpair;
            n++;
            unpair=0;
            map[i]=n;
        }
    }
    weight[n]=unpair;
    n++;
    weight[n]=0;
    
    if (n==2) return NULL;
    weight[0]=n;
    
    res=newmap(n);
    
    res->alpha[0]=n;
    res->alpha[1]=n; res->alpha[n]=1;
    for (i=1;i<=l;i++) {
        if (struc[i]!=0 && i<struc[i]) {
            res->alpha[map[i]] = map[struc[i]];
            res->alpha[map[struc[i]]] = map[i];
        }
    }
    
    res->tour[0]=n;
    
    for (i=1;i<n;i++) {
        res->tour[i]=i+1;
    }
    res->tour[n]=1;
    
    circ(res->tour, res->alpha, res->sigma);
    countbc(res);
    
    //   printmap (res);  //debug
    //   for (i=1;i<=n;i++) printf("%d ", weight[i]); printf("\n"); //debug
    free(map);
    return res;
    
}

int *map_to_structure(uc_map *s, int *weight) //from uc-map to structure
{
    int *res, l, i, *map, k, n;
    
    l=n=s->alpha[0];
    for (i=1;i<=weight[0];i++) l+=weight[i];
    l=l-2;
    res=(int *) calloc (l+10, sizeof(int));
    map=(int *) calloc (n+10, sizeof(int));
    res[0]=l;
    
    k=weight[1];
    for (i=2;i<n;i++) {
        map[i]=i+k-1;
        k+=weight[i];
    }
    
    //   for (i=2; i<n; i++) printf("%d ", map[i]); printf("\n"); //debug
    for (i=1;i<=l;i++) res[i]=0;
    for (i=2;i<n;i++) {
        res[map[i]]=map[s->alpha[i]];
    }
    
    //   for (i=0; i<=l; i++) printf("%d ", res[i]);printf("\n"); //debug
    
    
    free(map);
    return res;
}



int small_tour(uc_map *s, int a, int b)
{
    int i;
    i=1;
    while (i!=s->sigma[0]) {
        if (s->tour[i] == a) return 1;
        if (s->tour[i] == b) return 0;
        i=s->tour[i];
    }
    return 0;
}
//compqre two half-edge a and b respec to gamma in s a<b return 1, otherwise 0

uc_map * n_slice (uc_map *s, int tri, int * label)
{
    uc_map *res;
    int l=s->alpha[0], a, b, c, i;
    
    res=newmap(l);
    c = s->sigma[tri];
    
    a = tri;
    while (s->min[a]!=1) {
        a=s->sigma[a];
    }
    // find the minimum half-edge on the vertex that tri belongs to
    
    i = a; b=tri;
    while (s->sigma[i]!=c) {
        if (small_tour(s, c, s->sigma[i]) && small_tour(s, s->sigma[i] ,b)) {
            b = s->sigma[i];
        }
        i = s->sigma[i];
    }
    //find b
    
    label[a]=1; label[b]=1; label[c]=1;
    
    copy(s->alpha, res->alpha);
    copy(s->sigma, res->sigma);
    
    res->sigma[a] = s->sigma[c];
    res->sigma[b] = s->sigma[a];
    res->sigma[c] = s->sigma[b];
    
    circ(res->sigma, res->alpha, res->tour);
    countbc(res);
    
    //make res
    
    return res;
}


uc_map * n_I_slice (uc_map *s, int tri, int * label)
{
    uc_map *s1, *res;
    
    s1 = copy_uc(s);
    
    do{
        res = n_slice(s1, tri, label);
        freemap(s1);
        s1 = res;
    } while (find_trisection(res, tri) == 1);
    
    return res;
}


void n_slice_donw(uc_map *input, int *label)
{
    int i,j;
    uc_map *s1;
    
    
    if (0) {
        printmap(input);
        for (i=1;i<=input->alpha[0]; i++) {
            if (label[i]==1) printf("%d,", i);
        }
        printf("\n");
    } else {
        j=1;
        for (i=2; i<=input->trisection[0]; i++)
            if (small_tour(input, input->trisection[i], input->trisection[j]))
                j=i;

        s1 = n_I_slice(input, input->trisection[j], label);
        printf("Tri_select: %d\n", input->trisection[j]);
//        n_slice_donw(s1, label);
        printmap(s1);
        for (i=1;i<=s1->alpha[0]; i++) {
            if (label[i]==1) printf("%d,", i);
        }
        printf("\n");

    }
}

char *pair2structure(int *pair)
{
    int l,i,top=0;
    int *stack1, *stack2, *stack3, mov1=0, mov2=0, mov3=0;
    char *res;
    
    l=pair[0];
    
    stack1=(int *) calloc((l+2), sizeof(int));
    stack2=(int *) calloc((l+2), sizeof(int));
    stack3=(int *) calloc((l+2), sizeof(int));
    res=(char *) calloc((l+2), sizeof(int));
    
    for (i=1;i<=l;i++) {
        if (pair[i]!=0 && pair[i]>i) {
            if (mov1==0) {
                stack1[mov1]=pair[i];
                mov1++;
            } else if (mov1>0 && mov2==0) {
                if (pair[i]>stack1[mov1-1]) {
                    stack2[mov2]=pair[i];
                    mov2++;
                } else {
                    stack1[mov1]=pair[i];
                    mov1++;
                }
            } else if (mov1>0 && mov2>0 && mov3==0) {
                if (stack1[mov1-1]<stack2[mov2-1] && pair[i]>stack2[mov2-1]) {
                    stack3[mov3]=pair[i];
                    mov3++;
                } else if (pair[i]>stack1[mov1-1] && pair[i]<stack2[mov2-1]) {
                    stack2[mov2]=pair[i];
                    mov2++;
                } else if (pair[i]<stack1[mov1-1]) {
                    stack1[mov1]=pair[i];
                    mov1++;
                }
            } else {
                stack3[mov3]=pair[i];
                mov3++;
            }
        } else if (pair[i]!=0 && pair[i]<i) {
            if (i==stack1[mov1-1]) {
                res[i-1]=')';
                res[pair[i]-1]='(';
                mov1--;
            } else if (i==stack2[mov2-1]) {
                res[i-1]=']';
                res[pair[i]-1]='[';
                mov2--;
            } else if (i==stack3[mov3-1]) {
                res[i-1]='}';
                res[pair[i]-1]='{';
                mov3--;
            }
        } else if (pair[i]==0) {
            res[i-1]='.';
        }
    }
    res[i-1]=0;
    free(stack1);
    free(stack2);
    free(stack3);
    return res;
}

int *structure2pair(char *s)
{
    int *res, length, i;
    int *stack1, *stack2, *stack3, mov1=0, mov2=0, mov3=0;
    
    length=strlen(s);
    stack1=(int *) calloc (length+2, sizeof(int));
    stack2=(int *) calloc(length+2, sizeof(int));
    stack3=(int *) calloc (length+2, sizeof(int));
    res=(int *) calloc (length+2, sizeof(int));
    res[0]=length;
    
    for (i=1;i<=length;i++) res[i]=0;
    for (i=0;i<length;i++) {
        if (s[i]==':' || s[i]=='.') continue;
        else if (s[i]=='(') {
            stack1[mov1]=i+1;
            mov1++;
        } else if (s[i]==')') {
            mov1--;
            res[i+1]=stack1[mov1];
            res[stack1[mov1]]=i+1;
        } else if (s[i]=='[') {
            stack2[mov2]=i+1;
            mov2++;
        } else if (s[i]==']') {
            mov2--;
            res[i+1]=stack2[mov2];
            res[stack2[mov2]]=i+1;
        } else if (s[i]=='{') {
            stack3[mov3]=i+1;
            mov3++;
        } else if (s[i]=='}') {
            mov3--;
            res[i+1]=stack3[mov3];
            res[stack3[mov3]]=i+1;
        }
    }
    free(stack1);
    free(stack2);
    free(stack3);
    return res;
}

