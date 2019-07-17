//
//  uc_map_generator.h
//  ucmap
//
//  Created by Fenix Huang on 5/15/15.
//  Copyright (c) 2015 Fenix Huang. All rights reserved.
//

#ifndef __ucmap__uc_map_generator__
#define __ucmap__uc_map_generator__

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>


#endif /* defined(__ucmap__uc_map_generator__) */



extern int *G_sampler (int length, int target_g);
//uniform sample a structue (including unpaired bases with genus g)
extern int Uniform_sampling (int n, int target_g);
//main function
extern void precompute(int l);
