//include "../include/read_fof_snapshot.h"
#include "../include/adp2d.h"
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*
#define STB_IMAGE_IMPLEMENTATION
#include "../include/stb_image.h"
*/

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../include/stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "../include/stb_image_resize2.h"

#define MAX_SERIAL_MERGING 40000
#define MAX_N_NGBH 1000
#define PREALLOC_BORDERS 10

#define MAX(x,y) x > y ? x : y
#define MIN(x,y) x < y ? x : y

unsigned int data_dims;
idx_t Npart;
const border_t border_null = {.density = -1.0, .error = 0, .idx = NOBORDER};
const SparseBorder_t SparseBorder_null = {.density = -1.0, .error = 0, .idx = NOBORDER, .i = NOBORDER, .j = NOBORDER};

/*****************************
 * Clusters object functions *
 *****************************/

void Clusters_allocate(Clusters * c, int s)
{

    /*************************************
     * allocate additional resources and *
     * pointers for Clusters object      *
     *************************************/
    if(c -> centers.count == 0)
    {
        printf("Provide a valid cluster centers list\n");
        return;
    }

    idx_t nclus = c -> centers.count;
    
    if(s)
    {
	    //printf("Using sparse implementation\n");
	    c -> UseSparseBorders = 1;
	    c -> SparseBorders = (AdjList_t*)malloc(nclus*sizeof(AdjList_t));
	    for(idx_t i = 0; i < nclus; ++i)
	    {
		    c -> SparseBorders[i].count = 0;
		    c -> SparseBorders[i].size  = PREALLOC_BORDERS;
		    c -> SparseBorders[i].data  = (SparseBorder_t*)malloc(PREALLOC_BORDERS*sizeof(SparseBorder_t));
	    }

    }
    else
    {
	    //printf("Using dense implementation\n");
	    c -> UseSparseBorders = 0;
	    c -> __borders_data         = (border_t*)malloc(nclus*nclus*sizeof(border_t)); 
	    c -> borders                = (border_t**)malloc(nclus*sizeof(border_t*));

	    #pragma omp parallel for
	    for(idx_t i = 0; i < nclus; ++i)
	    {
		c -> borders[i]         = c -> __borders_data + i*nclus;
		for(idx_t j = 0; j < nclus; ++j)
		{
		    c -> borders[i][j] = border_null;
		}
	    }
    }
}

void AdjList_Insert(AdjList_t* l, SparseBorder_t b)
{
	if(l -> count < l -> size)
	{
		l -> data[l -> count] = b;
		l -> count++;
	}
	else
	{
		l -> size += PREALLOC_BORDERS; 
		l -> data = realloc( l -> data, sizeof(SparseBorder_t) * ( l -> size));
		l -> data[l -> count] = b;
		l -> count++;
	}
}

void AdjList_reset(AdjList_t* l)
{
	free(l -> data);
	l -> count = 0;
	l -> size  = 0;
	l -> data  = NULL;
}

void Clusters_Reset(Clusters * c)
{
	if(c -> UseSparseBorders)
	{
		for(idx_t i = 0; i < c -> centers.count; ++i)
		{
			AdjList_reset((c -> SparseBorders) + i);
		
		}
		free(c -> SparseBorders);
		c -> SparseBorders = NULL;
	}
	else
	{
		free(c -> __borders_data);
		free(c -> borders);
	}
    free(c -> centers.data);
}

void Clusters_free(Clusters * c)
{

    Clusters_Reset(c);
}


void SparseBorder_Insert(Clusters *c, SparseBorder_t b)
{
	idx_t i = b.i;
	AdjList_t l = c -> SparseBorders[i];
	int check = 1;
	for(idx_t k = 0; k < l.count; ++k)
	{
		SparseBorder_t p = l.data[k];
		if(p.i == b.i && p.j == b.j)
		{
			if( b.density > p.density)
			{
				l.data[k] = b;
			}
			check = 0;
		}
	}
	if(check) AdjList_Insert(c -> SparseBorders + i, b);
	return;
}

SparseBorder_t SparseBorder_get(Clusters* c, idx_t i, idx_t j)
{
	SparseBorder_t b = SparseBorder_null;
	AdjList_t l = c -> SparseBorders[i];
	for(idx_t el = 0; el < l.count; ++el)
	{
		SparseBorder_t candidate = l.data[el];
		if(candidate.i == i && candidate.j == j)
		{
			b = candidate;
		}
	}

	return b;
}

/*****************
 * Dyanmic Array *
 *****************/

void DynamicArray_allocate(lu_dynamicArray * a)
{
    a -> data = (idx_t*)malloc(ARRAY_INCREMENT*sizeof(idx_t));
    a -> count = 0;
    a -> size = ARRAY_INCREMENT;
}

void DynamicArray_pushBack(lu_dynamicArray * a, idx_t p)
{
    if(a -> count < a -> size)
    {
        a -> data[a -> count] =  p;
        a -> count += 1;
    }
    else{
        a -> size += ARRAY_INCREMENT;
        a -> data = realloc(a -> data, a -> size * sizeof(idx_t));
        a -> data[a -> count] =  p;
        a -> count += 1;
    }
}

void DynamicArray_Reset(lu_dynamicArray * a){
    a -> count = 0;
}

void DynamicArray_Reserve(lu_dynamicArray * a, idx_t n)
{
    a -> data = realloc(a -> data, n*sizeof(idx_t));
    a -> size = n;
}

void DynamicArray_Init(lu_dynamicArray * a)
{
    a -> data = NULL;
    a -> count = 0;
    a -> size = 0;
}


/*******************
 * Clustering part *
 *******************/


int cmp(const void * a, const void * b){
    FLOAT_TYPE aa = *((FLOAT_TYPE*)a);
    FLOAT_TYPE bb = *((FLOAT_TYPE*)b);
    return 2*(aa > bb ) - 1; 
}



FLOAT_TYPE avg(const FLOAT_TYPE * x, const idx_t n)
{
    FLOAT_TYPE f = 0;
    for(idx_t i = 0; i < n; ++i)
    {
        f += x[i];
    }
    return f/(FLOAT_TYPE)n;
}



FLOAT_TYPE mEst2(FLOAT_TYPE * x, FLOAT_TYPE *y, idx_t n)
{

    /********************************************
     * Estimate the m coefficient of a straight *
     * line passing through the origin          *
     * params:                                  *
     * - x: x values of the points              *
     * - y: y values of the points              *
     * - n: size of the arrays                  *
     ********************************************/
     

    //FLOAT_TYPE x_avg, y_avg;
    FLOAT_TYPE num = 0;
    FLOAT_TYPE den = 0;
    FLOAT_TYPE dd;
    for(idx_t i = 0; i < n; ++i)
    {
        FLOAT_TYPE xx = x[i];
        FLOAT_TYPE yy = y[i];

        dd = xx;
        num += dd*yy;
        den += dd*dd;

    }
  
    return num/den;
}
FLOAT_TYPE mEst(FLOAT_TYPE * x, FLOAT_TYPE *y, idx_t n)
{
    FLOAT_TYPE x_avg, y_avg;
    x_avg = avg(x,n);
    y_avg = avg(y,n);
    FLOAT_TYPE num = 0;
    FLOAT_TYPE den = 0;
    FLOAT_TYPE dd;
    for(idx_t i = 0; i < n - 1; ++i)
    {
        FLOAT_TYPE xx = x[i];
        FLOAT_TYPE yy = y[i];

        dd = (xx - x_avg);
        num += dd*(yy - y_avg);
        den += dd*dd;

    }
  
    return num/den;
}

int cmpPP(const void* p1, const void *p2)
{
    /***********************************************
     * Utility function to perform quicksort then  *
     * when clustering assignment is performed     *
     ***********************************************/
    Datapoint_info* pp1 = *(Datapoint_info**)p1;
    Datapoint_info* pp2 = *(Datapoint_info**)p2;
	float_t g1 = pp1 -> g;
	float_t g2 = pp2 -> g;
  	//return - ( DensA > DensB) + (DensA < DensB);
    return - (g1 > g2) + (g1 < g2);
    //return 2*(pp1 -> g < pp2 -> g) - 1;
}

void computeCorrection(Datapoint_info* dpInfo, int* mask, idx_t n, FLOAT_TYPE Z)
{
    /*****************************************************************************
     * Utility function, find the minimum value of the density of the datapoints *
     * and shift them up in order to further work with values greater than 0     *
     *****************************************************************************/
    FLOAT_TYPE min_log_rho = 999999.9;
    

    #pragma omp parallel
    {
        FLOAT_TYPE thread_min_log_rho = 9999999.;
        #pragma omp for
        for(idx_t i = 0; i < n; ++i)
        {
            FLOAT_TYPE tmp = dpInfo[i].log_rho - Z*dpInfo[i].log_rho_err;
            if(tmp < thread_min_log_rho && mask[i]){
                thread_min_log_rho = tmp;
            }
        }
        #pragma omp critical
        if(thread_min_log_rho < min_log_rho) min_log_rho = thread_min_log_rho;

        #pragma omp barrier 
        #pragma omp for
        for(idx_t i = 0; i < n; ++i)
        {
            dpInfo[i].log_rho_c = dpInfo[i].log_rho - min_log_rho + 1;
            dpInfo[i].g = dpInfo[i].log_rho_c - dpInfo[i].log_rho_err;
        }
    }
    //printf("%lf\n",min_log_rho);
}

//Clusters Heuristic1(Datapoint_info* dpInfo, int* mask, int nrows, int ncols)
Clusters Heuristic1(Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols)
{
    /**************************************************************
     * Heurisitc 1, from paper of Errico, Facco, Laio & Rodriguez *
     * ( https://doi.org/10.1016/j.ins.2021.01.010 )              *
     *                                                            *
     * args:                                                      *
     * - dpInfo: array of Datapoint structures                 *
     * - data: pointer to the dataset                             *
     * - n: number of Datapoints                                  *
     **************************************************************/

    struct timespec start_tot, finish_tot;
    double elapsed_tot;

    printf("H1: Preliminary cluster assignment\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);

    //idx_t ncenters = 0;
    //idx_t putativeCenters = n;
    lu_dynamicArray allCenters, removedCenters, actualCenters, max_rho;
    DynamicArray_allocate(&allCenters);
    DynamicArray_allocate(&removedCenters);
    DynamicArray_allocate(&actualCenters);
    DynamicArray_allocate(&max_rho);

    Datapoint_info** dpInfo_ptrs = (Datapoint_info**)malloc(nrows*ncols*sizeof(Datapoint_info*));

    struct timespec start, finish;
    double elapsed;


    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

    for(int i = 0; i < (int)nrows; ++i)
    for(int j = 0; j < (int)ncols; ++j)
    {   
        /*

        Find the centers of the clusters as the points of higher density in their neighborhoods
        A point is tagged as a putative center if it is the point of higer density of its neighborhood 
        
        */

        dpInfo_ptrs[i*ncols + j] = dpInfo + i*ncols + j;
        int r = (int)dpInfo[i*ncols + j].kstar;
        //int r = 50; 
        FLOAT_TYPE gi = dpInfo[i*ncols + j].g;
        dpInfo[i*ncols + j].is_center = mask[i*ncols + j] ? 1 : 0;
        dpInfo[i*ncols + j].cluster_idx = -1;
        //printf("%lf\n",p -> g);
		int jjmin = j - r > 0 			? j - r : 0;  
		int jjmax = j + r + 1 < (int)ncols 	? j + r + 1 : (int)ncols;  

		int iimin = i - r > 0 	 		? i - r : 0;  
		int iimax = i + r + 1 < (int)nrows 	? i + r + 1 : (int)nrows;  
		
		if(mask[i*ncols + j])
		{
			for(int ii = iimin; ii < iimax; ++ii)
			for(int jj = jjmin; jj < jjmax; ++jj)
			{
				idx_t ngbh_index = (idx_t)ii*ncols + jj; 
				FLOAT_TYPE gj = dpInfo[ngbh_index].g;
				if(gj > gi && mask[ngbh_index] && ((int)ngbh_index != (int)(i*ncols + j) )){
					dpInfo[i*ncols + j].is_center = 0;
					break;
				}
			}
		}
        if(dpInfo[i*ncols + j].is_center && mask[i*ncols + j]){
                DynamicArray_pushBack(&allCenters, i*ncols + j);
        }


    }

    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinding putative centers: %.3lfs\n",elapsed);
        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

	qsort(dpInfo_ptrs, nrows*ncols, sizeof(Datapoint_info*), cmpPP);

    idx_t * to_remove = (idx_t*)malloc(allCenters.count*sizeof(idx_t));
    for(idx_t c = 0; c < allCenters.count; ++c) {to_remove[c] = MY_SIZE_MAX;}

	idx_t* to_remove_mask = (idx_t*)malloc(nrows*ncols*sizeof(idx_t));
    for(idx_t p = 0; p < nrows*ncols; ++p) {to_remove_mask[p] = MY_SIZE_MAX;}

	
    #pragma omp parallel shared(to_remove_mask)
    {
        #pragma omp for
        for(idx_t p = 0; p < nrows*ncols; ++p)
        {
        	Datapoint_info pp = *(dpInfo_ptrs[p]);
			int i = (int)pp.array_idx / (int)ncols;
			int j = (int)pp.array_idx % (int)ncols;
			int r = (int)pp.kstar; //ATTENTION

			int jjmin = j - r > 0 				? j - r : 0;  
			int jjmax = j + r + 1 < (int)ncols 	? j + r + 1 : (int)ncols;  

			int iimin = i - r > 0 	 			? i - r : 0;  
			int iimax = i + r + 1 < (int)nrows 	? i + r + 1 : (int)nrows;  
			int flag = 0;
			idx_t ppp = 0;
			
			if(mask[i*ncols + j])
			{
				for(int ii = iimin; ii < iimax; ++ii)
				for(int jj = jjmin; jj < jjmax; ++jj)
				{
					idx_t jidx = ii*ncols + jj;
					if(dpInfo[jidx].is_center && pp.g > dpInfo[jidx].g && mask[jidx])
					{
						
						#pragma omp critical 
						{
							ppp = to_remove_mask[jidx];
							flag = ppp != MY_SIZE_MAX;							
							to_remove_mask[jidx] = flag ? (pp.g > dpInfo[ppp].g ? pp.array_idx : ppp) : pp.array_idx; 
						}
						
						//#pragma omp atomic read 
						//ppp = to_remove_mask[jidx];

						//flag = ppp != MY_SIZE_MAX;							
						//
						//#pragma omp atomic write
						//to_remove_mask[jidx] = flag ? (pp.g > dpInfo[ppp].g ? pp.array_idx : ppp) : pp.array_idx; 
					}
				}
			}
		}
	}
    
    

    for(idx_t p = 0; p < allCenters.count; ++p)
    {
        idx_t i = allCenters.data[p];
        int e = 0;
        //FLOAT_TYPE gi = dpInfo[i].g;
        idx_t mr = to_remove_mask[i];
        if(mr != MY_SIZE_MAX)
        {
            //if(dpInfo[mr].g > gi) e = 1;
			e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    DynamicArray_pushBack(&removedCenters,i);
                    dpInfo[i].is_center = 0;
                    for(idx_t c = 0; c < removedCenters.count - 1; ++c)
                    {
                        if(mr == removedCenters.data[c])
                        {
                            mr = max_rho.data[c];
                        }
                    }
                    DynamicArray_pushBack(&max_rho,mr);
                    
                }
                break;
            case 0:
                {
                    DynamicArray_pushBack(&actualCenters,i);
                    dpInfo[i].cluster_idx = actualCenters.count - 1;
                }
                break;
            default:
                break;
        }
    }


	
/*	
    #pragma omp parallel
    {
            
        idx_t * to_remove_private = (idx_t*)malloc(allCenters.count*sizeof(idx_t));
    	for(idx_t c = 0; c < allCenters.count; ++c) {to_remove_private[c] = MY_SIZE_MAX;}

        #pragma omp for
        for(idx_t p = 0; p < nrows*ncols; ++p)
        {
        	Datapoint_info pp = *(dpInfo_ptrs[p]);
			int i = (int)pp.array_idx / (int)ncols;
			int j = (int)pp.array_idx % (int)ncols;
			int r = (int)pp.kstar; //ATTENTION

			int jjmin = j - r > 0 			? j - r : 0;  
			int jjmax = j + r + 1 < ncols 	? j + r + 1 : ncols;  

			int iimin = i - r > 0 	 		? i - r : 0;  
			int iimax = i + r + 1 < nrows 	? i + r + 1 : nrows;  
			
			if(mask[i*ncols + j])
			{
				for(int ii = iimin; ii < iimax; ++ii)
				for(int jj = jjmin; jj < jjmax; ++jj)
				{
					idx_t jidx = ii*ncols + jj;
					if(dpInfo[jidx].is_center && pp.g > dpInfo[jidx].g && mask[jidx])
					{
						//dpInfo[jidx].is_center = 0;
						for(idx_t c = 0; c < allCenters.count; ++c)
						{
							if(allCenters.data[c] == jidx)
							{

								if(to_remove_private[c] != MY_SIZE_MAX)
								{
									to_remove_private[c] = pp.g > 	dpInfo[to_remove_private[c]].g  ? pp.array_idx : to_remove_private[c];
								}
								else
								{
									to_remove_private[c] = pp.array_idx;
								}
							}
						}
					}
				}
			}
        }
        #pragma omp critical
        {
        	for(idx_t c = 0; c < allCenters.count; ++c)
        	{
        		if(to_remove_private[c] != MY_SIZE_MAX)
			{
				if(to_remove[c] != MY_SIZE_MAX)
				{
					to_remove[c] = dpInfo[to_remove_private[c]].g > dpInfo[to_remove[c]].g ?
						       to_remove_private[c] : to_remove[c];
				}
				else
				{
					to_remove[c] = to_remove_private[c];
				}
			}
        	}
        }

            free(to_remove_private);
    }
	

	
	
    for(idx_t p = 0; p < allCenters.count; ++p)
    {
        idx_t i = allCenters.data[p];
        int e = 0;
        //FLOAT_TYPE gi = dpInfo[i].g;
        idx_t mr = to_remove[p];
        if(mr != MY_SIZE_MAX)
        {
            //if(dpInfo[mr].g > gi) e = 1;
	    e = 1;
        }
        switch (e)
        {
            case 1:
                {
                    DynamicArray_pushBack(&removedCenters,i);
                    dpInfo[i].is_center = 0;
                    for(idx_t c = 0; c < removedCenters.count - 1; ++c)
                    {
                        if(mr == removedCenters.data[c])
                        {
                            mr = max_rho.data[c];
                        }
                    }
                    DynamicArray_pushBack(&max_rho,mr);
                    
                }
                break;
            case 0:
                {
                    DynamicArray_pushBack(&actualCenters,i);
                    dpInfo[i].cluster_idx = actualCenters.count - 1;
                }
                break;
            default:
                break;
        }
    }
	
	*/
	free(to_remove);
	free(to_remove_mask);


    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinding actual centers:   %.3lfs\n",elapsed);

        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif


    //idx_t nclusters = 0;


    /*****************************************************************************
     * Sort all the dpInfo based on g and then perform the cluster assignment *
     * in asceding order                                                         *
     * UPDATE: dpInfo already sorted                                          *
     *****************************************************************************/
                                                                                

    //qsort(dpInfo_ptrs, n, sizeof(Datapoint_info*), cmpPP);
	
	
	idx_t* fromWho = (idx_t*)malloc(nrows*ncols*sizeof(idx_t));
    for(idx_t pidx = 0; pidx < nrows*ncols; ++pidx) fromWho[pidx] = SIZE_MAX;

	
	
	
	#pragma omp parallel for schedule(dynamic)
    for(idx_t pidx = 0; pidx < nrows*ncols; ++pidx)
    {   
        Datapoint_info* p = dpInfo_ptrs[pidx];
		int i = (int)(p -> array_idx) / (int)ncols;
		int j = (int)(p -> array_idx) % (int)ncols;
		int r = p -> kstar ; //ATTENTION
		//int r = 5; //ATTENTION
		int iimin, iimax, jjmin, jjmax;
        //idx_t ele = p -> array_idx;
        //fprintf(f,"%lu\n",ele);
        if(!(p -> is_center) && mask[i*ncols + j])
        {
            //int cluster = -1;
            //idx_t k = 0;
            //assign each particle at the same cluster as the nearest particle of higher density
			jjmin = j - r > 0 			? j - r : 0;  
			jjmax = j + r + 1 < (int)ncols 	? j + r + 1 : (int)ncols;  

			iimin = i - r > 0 	 		? i - r : 0;  
			iimax = i + r + 1 < (int)nrows 	? i + r + 1 : (int)nrows;  
			
			//int ii_toTakeFrom, jj_toTakeFrom;
			long int minNgbhDist = nrows*nrows*ncols*ncols;
			float_t g_current = dpInfo[i*ncols + j].g;
			int foundFlag = 0;
			for(int ii = iimin; ii < iimax; ++ii)
			for(int jj = jjmin; jj < jjmax; ++jj)
			//Take the same cluster as the nearest neighbor with higher density
			{
				//take the ngbh
				long int currentDist = (ii-i)*(ii-i) + (jj-j)*(jj-j);	
				int notMySelf = (ii != i) || (jj != j);
				idx_t ngbhIdx = ii*ncols + jj;
				float_t g_ngbh = dpInfo[ngbhIdx].g; 
				//if(ii*ncols + j == 302000)
				//{
				//	printf("Nopeh\n");
				//}
				if(g_ngbh > g_current && notMySelf && currentDist < minNgbhDist)
				{
					minNgbhDist = currentDist;
					//ii_toTakeFrom = ii;
					//jj_toTakeFrom = jj;
					//cluster = dpInfo[p_idx].cluster_idx; 
					fromWho[i*ncols + j] = (idx_t)(ii*ncols+jj);
					foundFlag = 1;
				}
			}


            //
            if(!foundFlag)
            {
                FLOAT_TYPE gmax = -99999.;               
                idx_t gm_index = 0;

				for(int ii = iimin; ii < iimax; ++ii)
				for(int jj = jjmin; jj < jjmax; ++jj)
                {
                    idx_t ngbh_index = ii*ncols + jj;
                    for(idx_t m = 0; m < removedCenters.count; ++m)
                    {
                        FLOAT_TYPE gcand = dpInfo[max_rho.data[m]].g;
                        if(ngbh_index == removedCenters.data[m] && gcand > gmax)
                        {   
                            //printf("%lu -- %lu\n", ele, m);
                            gmax = gcand;
                            gm_index = max_rho.data[m];
							fromWho[i*ncols + j] = gm_index;
							foundFlag = 1;
                        }
                    }
                }

                //cluster = dpInfo[gm_index].cluster_idx;

            }
            //p -> cluster_idx = cluster;
			if(!foundFlag) mask[i*ncols + j] = 0;

		}
	}


	//printf("aa\n");
	//#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < (int)nrows; ++i)
	for(int j = 0; j < (int)ncols; ++j)
	{
		idx_t pidx = dpInfo_ptrs[i*ncols + j] -> array_idx;
		if(mask[pidx] && !(dpInfo[pidx].is_center))
		{
			idx_t idxToTakeFrom = fromWho[pidx];
			//int cluster = dpInfo[idxToTakeFrom].cluster_idx;				
			int cluster = -1;
			while(cluster == -1)
			{
				cluster = dpInfo[idxToTakeFrom].cluster_idx;				
				idxToTakeFrom = fromWho[idxToTakeFrom];
			}
			dpInfo[pidx].cluster_idx = cluster;
		}
	}
	
	
/*	
	
    for(idx_t pidx = 0; pidx < nrows*ncols; ++pidx)
    {   
        Datapoint_info* p = dpInfo_ptrs[pidx];
		int i = (int)(p -> array_idx) / ncols;
		int j = (int)(p -> array_idx) % ncols;
		//int r = p -> kstar + 1; //ATTENTION
		int r = p -> kstar; //ATTENTION
		//int r = 5; //ATTENTION
		idx_t ggg;
		int iimin, iimax, jjmin, jjmax;
        //idx_t ele = p -> array_idx;
        //fprintf(f,"%lu\n",ele);
        if(!(p -> is_center) && mask[i*ncols + j])
        {
            int cluster = -1;
            idx_t k = 0;
            idx_t p_idx;
            //assign each particle at the same cluster as the nearest particle of higher density
			jjmin = j - r > 0 				? j - r : 0;  
			jjmax = j + r + 1 < (int)ncols 	? j + r + 1 : (int)ncols;  

			iimin = i - r > 0 	 			? i - r : 0;  
			iimax = i + r + 1 < (int)nrows 	? i + r + 1 : (int)nrows;  
			
			int ii_toTakeFrom, jj_toTakeFrom;
			long int minNgbhDist = nrows*nrows*ncols*ncols;
			for(int ii = iimin; ii < iimax; ++ii)
			for(int jj = jjmin; jj < jjmax; ++jj)
			//Take the same cluster as the nearest neighbor with higher density
			{
				//take the ngbh
				long int currentDist = (ii-i)*(ii-i) + (jj-j)*(jj-j);	
				int notMySelf = (ii != i) || (jj != j);
				//if(ii*ncols + j == 302000)
				//{
				//	printf("Nopeh\n");
				//}
				if(dpInfo[ii*ncols + jj].cluster_idx != -1 && notMySelf && currentDist < minNgbhDist)
				{
					minNgbhDist = currentDist;
					p_idx = ii*ncols + jj;
					ii_toTakeFrom = ii;
					jj_toTakeFrom = jj;
					cluster = dpInfo[p_idx].cluster_idx; 
					ggg = p_idx;
				}
			}


            //
            if(cluster == -1)
            {
                FLOAT_TYPE gmax = -99999.;               
                idx_t gm_index = 0;

				for(int ii = iimin; ii < iimax; ++ii)
				for(int jj = jjmin; jj < jjmax; ++jj)
                {
                    idx_t ngbh_index = ii*ncols + jj;
                    for(idx_t m = 0; m < removedCenters.count; ++m)
                    {
                        FLOAT_TYPE gcand = dpInfo[max_rho.data[m]].g;
                        if(ngbh_index == removedCenters.data[m] && gcand > gmax)
                        {   
                            //printf("%lu -- %lu\n", ele, m);
                            gmax = gcand;
                            gm_index = max_rho.data[m];
							ggg = p_idx;
                        }
                    }
                }

                cluster = dpInfo[gm_index].cluster_idx;

            }
            p -> cluster_idx = cluster;
//			if(fromWho[i*ncols + j] != ggg) printf("Nope in %lu got m1 %lu mog %lu\n", 
					i*ncols + j, 
					fromWho[i*ncols + j], ggg); 
			if(cluster == -1) mask[i*ncols + j] = 0;


        }
	}


*/
	

    

    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tTentative clustering:     %.3lfs\n",elapsed);

        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

    free(dpInfo_ptrs);
    free(max_rho.data);
    free(removedCenters.data);
    free(allCenters.data);


    Clusters c_all;
    c_all.centers = actualCenters;


    #ifdef VERBOSE
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinalizing clustering:    %.3lfs\n",elapsed);
        printf("\n");
    #endif

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;


    printf("\tFound %lu clusters\n",(uint64_t)actualCenters.count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    c_all.n = nrows*ncols;
    return c_all;
}

void Heuristic2(Clusters* cluster, Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols)
{

    #define borders cluster->borders

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    //idx_t n = cluster -> n;

    printf("H2: Finding border points\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);


    idx_t nclus = cluster->centers.count; 
    //idx_t max_k = dpInfo[0].ngbh.N;


    for(int i = 0; i < (int)nrows; ++i)
    for(int j = 0; j < (int)ncols; ++j)
    {
            idx_t pp = NOBORDER;
            /*loop over n neighbors*/
            int c = dpInfo[i*ncols + j].cluster_idx;
            if((!(dpInfo[i*ncols + j].is_center)) && mask[i*ncols + j])
            {
				int r = (int)dpInfo[i*ncols + j].kstar;				
				int jjmin = j - r > 0 				? j - r : 0;  
				int jjmax = j + r + 1 < (int)ncols 	? j + r + 1 : (int)ncols;  

				int iimin = i - r > 0 	 			? i - r : 0;  
				int iimax = i + r + 1 < (int)nrows 	? i + r + 1 : (int)nrows;  
				
				long int minNgbhDist = nrows*nrows*ncols*ncols;
				for(int ii = iimin; ii < iimax; ++ii)
				for(int jj = jjmin; jj < jjmax; ++jj)
                {
                    /*index of the kth ngbh of n*/
                    idx_t jidx = ii*ncols + jj;
					long int currentNgbhDist = (ii-i)*(ii-i) + (jj - j)*(jj - j);
                    /*Loop over kn neigbhours to find if n is the nearest*/
                    /*if cluster of the particle in nbhg is c then check is neighborhood*/                                                
                    if(dpInfo[jidx].cluster_idx != -1 
							&& dpInfo[jidx].cluster_idx != c 
							&& !dpInfo[jidx].is_center  
							&& mask[jidx] 
							&& currentNgbhDist < minNgbhDist)
                    {
						minNgbhDist = currentNgbhDist;
                        pp = jidx;
                    }

                }
            }

            if(pp != NOBORDER)
            {
				int r = (int)dpInfo[pp].kstar;				
				int ngbh_i = (int)pp / (int)ncols; 
				int ngbh_j = (int)pp % (int)ncols; 
				int jjmin = ngbh_j - r > 0 				? ngbh_j - r : 0;  
				int jjmax = ngbh_j + r + 1 < (int)ncols ? ngbh_j + r + 1 : (int)ncols;  

				int iimin = ngbh_i - r > 0 	 			? ngbh_i - r : 0;  
				int iimax = ngbh_i + r + 1 < (int)nrows ? ngbh_i + r + 1 : (int)nrows;  
				
				long int minNgbhDist = nrows*nrows*ncols*ncols;
				idx_t nearestBelongingToC = NOBORDER;
				for(int ii = iimin; ii < iimax; ++ii)
				for(int jj = jjmin; jj < jjmax; ++jj)
                {
					idx_t pp_ngbh_idx = ii*ncols + jj;
					long int currentNgbhDist = (ii-i)*(ii-i) + (jj - j)*(jj - j);
					//find if the nearest is the starting point
                    if(dpInfo[pp_ngbh_idx].cluster_idx == c && currentNgbhDist < minNgbhDist )
                    {
						minNgbhDist = currentNgbhDist;
						nearestBelongingToC = pp_ngbh_idx;
                    }
                }
				if(nearestBelongingToC != i*ncols + j)
				{
					pp = NOBORDER;
				}
            }
                            /*if it is the maximum one add it to the cluster*/
            if(pp != NOBORDER)
            {
				int ppc = dpInfo[pp].cluster_idx;
				if(cluster -> UseSparseBorders)
				{
					//insert one and symmetric one
					SparseBorder_t b = {.i = c, .j = ppc, .idx = i*ncols + j, .density = dpInfo[i*ncols + j].g, .error = dpInfo[i*ncols + j].log_rho_err}; 
					SparseBorder_Insert(cluster, b);
					//get symmetric border
					SparseBorder_t bsym = {.i = ppc, .j = c, .idx = i*ncols + j, .density = dpInfo[i*ncols + j].g, .error = dpInfo[i*ncols + j].log_rho_err}; 
					SparseBorder_Insert(cluster, bsym);

				}
				else
				{
					if(dpInfo[i*ncols + j].g > borders[c][ppc].density)
					{
						borders[c][ppc].density = dpInfo[i*ncols + j].g;
						borders[ppc][c].density = dpInfo[i*ncols + j].g;
						borders[c][ppc].idx = i*ncols + j;
						borders[ppc][c].idx = i*ncols + j;
					}
				}
			}

}


	if(cluster -> UseSparseBorders)
	{
		for(idx_t c = 0; c < nclus; ++c)
		{
			for(idx_t el = 0; el < cluster -> SparseBorders[c].count; ++el)
			{
				//fix border density, write log rho c
				idx_t idx = cluster -> SparseBorders[c].data[el].idx; 
				cluster -> SparseBorders[c].data[el].density = dpInfo[idx].log_rho_c;
			}
		}

	}
	else
	{
		for(idx_t bi = 0; bi < nclus - 1; ++bi)
		{
		for(idx_t bj = bi + 1; bj < nclus; ++bj)
		{
			idx_t p = borders[bi][bj].idx;
			if(p != NOBORDER)
			{   

			borders[bi][bj].density = dpInfo[p].log_rho_c;
			borders[bj][bi].density = dpInfo[p].log_rho_c;

			borders[bi][bj].error = dpInfo[p].log_rho_err;
			borders[bj][bi].error = dpInfo[p].log_rho_err;
			}
		}
		}

		for(idx_t dd = 0; dd < nclus; ++dd)
		{
		borders[dd][dd].density = -1.0;
		borders[dd][dd].error = 0.0;
		}
    }

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    return;
    #undef borders
   }



void Merge_A_into_B(idx_t* who_amI, idx_t cluster_A, idx_t cluster_B, idx_t n)
{
    #pragma omp parallel if(n > MAX_SERIAL_MERGING)
    {
	    idx_t tmp;
	    #pragma omp for
	    for(idx_t i = 0; i < n; ++i)
	    {   
		//substitute occurencies of b with a 
		tmp = who_amI[i] == cluster_A ? cluster_B : who_amI[i];
		who_amI[i] = tmp;
	    }
    }
    return;
}


int compare_merging_density( const void *A, const void *B)
{
  FLOAT_TYPE DensA = ((merge_t*)A)->density;
  FLOAT_TYPE DensB = ((merge_t*)B)->density;

  return - ( DensA > DensB) + (DensA < DensB);
}


inline int is_a_merging( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
			 FLOAT_TYPE dens2, FLOAT_TYPE dens2_err,
			 FLOAT_TYPE dens_border, FLOAT_TYPE dens_border_err,
			 FLOAT_TYPE Z)
/*
 * dens1 : the density of the particle that is the center of the first cluster
 * dens2 : the density of the particle that is the center of the second cluster
 * dens_border : the density of the border btw the cluster 1 and the cluster 2
 * *_err : the errors on the densities
 * Z     : the desired accuracy
 */
{
  /* in the original code it was:
   *
  FLOAT_TYPE a1 = dpInfo[cluster->centers.data[i]].log_rho_c - border_density[i][j];
  FLOAT_TYPE a2 = dpInfo[cluster->centers.data[j]].log_rho_c - border_density[i][j];
  
  FLOAT_TYPE e1 = Z*(dpInfo[cluster->centers.data[i]].log_rho_err + border_err[i][j]);
  FLOAT_TYPE e2 = Z*(dpInfo[cluster->centers.data[j]].log_rho_err + border_err[i][j]);
  */

  FLOAT_TYPE a1 = dens1 - dens_border;
  FLOAT_TYPE a2 = dens2 - dens_border;

  FLOAT_TYPE e1 = Z*(dens1_err + dens_border_err);
  FLOAT_TYPE e2 = Z*(dens2_err + dens_border_err);

  return (a1 < e1 || a2 < e2);
}


int merging_roles( FLOAT_TYPE dens1, FLOAT_TYPE dens1_err,
			  FLOAT_TYPE dens2, FLOAT_TYPE dens2_err,
			  FLOAT_TYPE dens_border, FLOAT_TYPE dens_border_err )
{
      
  FLOAT_TYPE c1 = (dens1 - dens_border) / (dens1_err + dens_border_err); 
  FLOAT_TYPE c2 = (dens2 - dens_border) / (dens2_err + dens_border_err);
  //printf("%.10lf %.10lf %d\n",c1,c2, c1 > c2);
  
  return ( c1 < c2 );     // if 1, this signal to swap 1 and 2
}

void fix_borders_A_into_B(idx_t A, idx_t B, border_t** borders, idx_t n)
{
   #pragma omp parallel for if(n > MAX_SERIAL_MERGING)
   for(idx_t i = 0; i < n; ++i) 
   {
        if(borders[A][i].idx != NOBORDER )
        {
            if(borders[B][i].idx != NOBORDER)
            {
                int mb = (borders[A][i].density > borders[B][i].density); 

                borders[B][i] = mb ? borders[A][i] : borders[B][i];
                borders[i][B] = borders[B][i];
            }
            else
            {
                borders[B][i] = borders[A][i];
                borders[i][B] = borders[B][i];
            }
        } 
        borders[A][i] = border_null;
        borders[i][A] = border_null;
   }
}

void Delete_adjlist_element(Clusters * c, const idx_t list_idx, const idx_t el)
{
	//swap last element with 
	idx_t count = c -> SparseBorders[list_idx].count;
	c -> SparseBorders[list_idx].data[el] = c -> SparseBorders[list_idx].data[count-1];
	c -> SparseBorders[list_idx].data[count-1] = SparseBorder_null;
	c -> SparseBorders[list_idx].count -= 1;
}

void fix_SparseBorders_A_into_B(idx_t s,idx_t t,Clusters* c)
{
	//delete border trg -> src
	
	//idx_t nclus = c -> centers.count;
	
	{
		{
			for(idx_t el = 0; el < c -> SparseBorders[t].count; ++el)
			{
				SparseBorder_t b = c -> SparseBorders[t].data[el];
				if(b.i == t && b.j == s)
				{
					//delete the border src trg
					Delete_adjlist_element(c, t, el);
				}
			}
		}
		//find the border and delete it, other insert them in correct place
		for(idx_t el = 0; el < c -> SparseBorders[s].count; ++el)
		{
			SparseBorder_t b = c -> SparseBorders[s].data[el];
		//	idx_t ii = b.i;
			if(b.j != t)
			{
				//insert these borders as trg -> j and j -> trg
				b.i = t;
				SparseBorder_Insert(c, b);
				SparseBorder_t bsym = b;
				bsym.i = b.j;
				bsym.j = b.i;
				SparseBorder_Insert(c, bsym);
				for(idx_t dl = 0; dl < c -> SparseBorders[b.j].count; ++dl)
				{
					SparseBorder_t b_del = c -> SparseBorders[b.j].data[dl];
					if(b_del.j == s)
					{
						//delete the border src trg
						Delete_adjlist_element(c, b.j, dl);
					}
				}
						
			}
		}
		//clean up all borders
		//delete the src list
		{
			AdjList_reset((c->SparseBorders) + s);
		}
		//delete all borders containing src
	//	for(idx_t i = 0; i < nclus; ++i)
	//	{
	//		for(idx_t el = 0; el < c -> SparseBorders[i].count; ++el)
	//		{
	//			SparseBorder_t b = c -> SparseBorders[i].data[el];
	//			if(b.j == s)
	//			{
	//				//delete the border src trg
	//				Delete_adjlist_element(c, i, el);
	//			}
	//		}
	//			
	//	}
	}


}

void Heuristic3_sparse(Clusters* cluster, Datapoint_info* dpInfo, FLOAT_TYPE Z, int halo)
{
  printf("Using sparse implementation\n");
  #define borders cluster->borders

  struct timespec start_tot, finish_tot;
  double elapsed_tot;

  struct timespec start, finish;
  double elapsed;

  printf("H3: Merging clusters\n");
  clock_gettime(CLOCK_MONOTONIC, &start_tot);
  #ifdef VERBOSE
 	 clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif

  idx_t nclus                 = cluster -> centers.count;  
  idx_t *  surviving_clusters = (idx_t*)malloc(nclus*sizeof(idx_t));
  for(idx_t i = 0; i < nclus; ++i)
    { 
        surviving_clusters[i] = i; 
    }

  idx_t   merge_count        = 0;
  idx_t   merging_table_size = 1000;
  merge_t *merging_table      = (merge_t*)malloc(sizeof(merge_t)*merging_table_size);
  
  /*Find clusters to be merged*/
  for(idx_t i = 0; i < nclus - 1; ++i)   
  {
    idx_t count = cluster -> SparseBorders[i].count;
    for(idx_t el = 0; el < count; ++el)   
    {
	      SparseBorder_t b = cluster -> SparseBorders[i].data[el];
	      if( b.j > b.i)
	      {
		      FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[b.i]].log_rho_c;
		      FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[b.i]].log_rho_err;
		      FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[b.j]].log_rho_c;
		      FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[b.j]].log_rho_err;
		      FLOAT_TYPE dens_border     = b.density;
		      FLOAT_TYPE dens_border_err = b.error;
	      
		      if ( is_a_merging( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err, Z ) )
			{
			  
			  if ( merge_count == merging_table_size ) {
			    merging_table_size *= 1.1;
			    merging_table = (merge_t*)realloc( merging_table, sizeof(merge_t) * merging_table_size ); }

			  idx_t src = b.j;
			  idx_t trg = b.i;

			  merging_table[merge_count].source = src;
			  merging_table[merge_count].target = trg;
			  merging_table[merge_count].density = b.density;
			  ++merge_count;
			}
	      }

	}
            
  }

  qsort( (void*)merging_table, merge_count, sizeof(merge_t), compare_merging_density);

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinding merges:   %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
  
    for( idx_t m = 0; m < merge_count; m++ )
    {
      
        #define src surviving_clusters[merging_table[m].source]
        #define trg surviving_clusters[merging_table[m].target]
        //printf("Found: %lu, %lu which now is %lu, %lu\n",merging_table[m].source, merging_table[m].target, src,trg);

        //int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
	//if(re_check)
	{
		idx_t new_src = (src < trg) ? src : trg;
		idx_t new_trg = (src < trg) ? trg : src;

                //pick who am I

                FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[new_src]].log_rho_c;
                FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[new_src]].log_rho_err;
                FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[new_trg]].log_rho_c;
                FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[new_trg]].log_rho_err;

		//borders get
		SparseBorder_t b 	   = SparseBorder_get(cluster, new_src, new_trg);
                FLOAT_TYPE dens_border     = b.density;
                FLOAT_TYPE dens_border_err = b.error;

                int i_have_to_merge = is_a_merging(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err,Z);            
                switch (i_have_to_merge && src != trg)
                {
                case 1:
                    {
                        int side = merging_roles(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err);
                        if(!side)
                        {
                            idx_t tmp;
                            tmp = new_src;
                            new_src = new_trg;
                            new_trg = tmp;
                        }

                        //borders[new_src][new_trg] = border_null;
                        //borders[new_trg][new_src] = border_null;
                        //printf("Merging %lu into %lu\n",new_src,new_trg);
                        fix_SparseBorders_A_into_B(new_src,new_trg,cluster);
                        Merge_A_into_B ( surviving_clusters, new_src, new_trg, nclus );	  
                    }
                    break;
                
                default:
                    break;
                }
	}
        
        #undef src
        #undef trg
    }

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tCluster merging:  %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
    /*Finalize clustering*/
    /*Acutally copying */
    lu_dynamicArray tmp_centers;
    lu_dynamicArray tmp_cluster_idx;


    DynamicArray_Init(&tmp_centers);
    DynamicArray_Init(&tmp_cluster_idx);

    DynamicArray_Reserve(&tmp_centers, nclus);
    DynamicArray_Reserve(&tmp_cluster_idx, nclus);

    idx_t final_cluster_count = 0;

    idx_t* old_to_new = (idx_t*)malloc(nclus*sizeof(idx_t));
    idx_t incremental_k = 0;
    for(idx_t i = 0; i < nclus; ++i)
    {
        
        if(surviving_clusters[i] == i){
            DynamicArray_pushBack(&tmp_centers, cluster->centers.data[i]);
            DynamicArray_pushBack(&tmp_cluster_idx, i);
            old_to_new[i] = incremental_k;
            ++incremental_k;
            ++final_cluster_count;
        }
    }

    //fill the rest of old_to_new
    for(idx_t i = 0; i < nclus; ++i)
    {
		idx_t cidx_to_copy_from = surviving_clusters[i];
		old_to_new[i] = old_to_new[cidx_to_copy_from];
    }

    /*allocate auxiliary pointers to store results of the finalization of the procedure*/

    AdjList_t* tmp_borders      = (AdjList_t*)malloc(final_cluster_count*sizeof(AdjList_t));

    //initialize temporary borders
    for(idx_t i = 0; i < final_cluster_count; ++i)
    {
	    tmp_borders[i].count = 0;
	    tmp_borders[i].size  = PREALLOC_BORDERS;
	    tmp_borders[i].data  = (SparseBorder_t*)malloc(PREALLOC_BORDERS*sizeof(SparseBorder_t));
    }

    /*initialize all pointers*/

    /*Fix cluster assignment*/
    #pragma omp parallel for
    for(idx_t i = 0; i < cluster -> n; ++i)
    {
        dpInfo[i].is_center = 0;
        int old_cidx = dpInfo[i].cluster_idx;
		if(old_cidx != -1)
		{
			dpInfo[i].cluster_idx = old_to_new[old_cidx];
		}

    }

    
    #pragma omp parallel for
    for(idx_t c = 0; c < final_cluster_count; ++c)
    {
        idx_t c_idx = tmp_cluster_idx.data[c];
		for(idx_t el = 0; el < cluster -> SparseBorders[c_idx].count; ++el)
		{
			//retrieve border
			SparseBorder_t b = cluster -> SparseBorders[c_idx].data[el];
			//change idexes of clusters
			b.i = old_to_new[b.i];
			b.j = old_to_new[b.j];

			AdjList_Insert(tmp_borders + c, b);
		}
    }

    Clusters_Reset(cluster);
    /*pay attention to the defined borders*/
    /*copy into members*/
    cluster -> SparseBorders = tmp_borders;


    cluster -> centers = tmp_centers;
    /**
     * Fix center assignment
    */
    for(idx_t i = 0; i < cluster -> centers.count; ++i)
    {
        int idx = cluster -> centers.data[i];
        dpInfo[idx].is_center = 1;
    }
    /*Halo*/
    switch (halo)
    {
    case 1:
	{
		FLOAT_TYPE* max_border_den_array = (FLOAT_TYPE*)malloc(final_cluster_count*sizeof(FLOAT_TYPE));
		#pragma omp parallel
		{
		    #pragma omp for
		    for(idx_t c = 0; c < final_cluster_count; ++c)
		    {
				FLOAT_TYPE max_border_den = -2.;
				for(idx_t el = 0; el < cluster -> SparseBorders[c].count; ++el)
				{
					SparseBorder_t b = cluster -> SparseBorders[c].data[el];
					if(b.density > max_border_den)
					{
						max_border_den = b.density;
					}
				}
				max_border_den_array[c] = max_border_den;
		    }

		    #pragma omp barrier

		    #pragma omp for
		    for(idx_t i = 0; i < cluster -> n; ++i)
		    {
				int cidx = dpInfo[i].cluster_idx;
				//int halo_flag;
				if(cidx != -1)
				{
					int halo_flag = dpInfo[i].log_rho_c < max_border_den_array[cidx] && !dpInfo[i].is_center; 
					dpInfo[i].cluster_idx = halo_flag ? -1 : cidx;
				}
		    }
		}
		free(max_border_den_array);
	}
        break;
    
    default:
        break;
    }    

    /*free memory and put the correct arrays into place*/
    free(tmp_cluster_idx.data);
    free(merging_table);
    //free(ipos.data);
    //free(jpos.data);
    free(surviving_clusters);
    free(old_to_new);

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinal operations: %.3lfs\n\n", elapsed);
  #endif

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tFound %lu possible merges\n",(uint64_t)merge_count);
    printf("\tSurviving clusters %lu\n",(uint64_t)final_cluster_count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

  #undef  borders  
}


void Heuristic3_dense(Clusters* cluster, Datapoint_info* dpInfo, FLOAT_TYPE Z, int halo)
{
  printf("Using dense implementation\n");
  #define borders cluster->borders

  struct timespec start_tot, finish_tot;
  double elapsed_tot;

  struct timespec start, finish;
  double elapsed;

  printf("H3: Merging clusters\n");
  clock_gettime(CLOCK_MONOTONIC, &start_tot);
  #ifdef VERBOSE
 	 clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif

  idx_t nclus              = cluster -> centers.count;  
  idx_t *  surviving_clusters = (idx_t*)malloc(nclus*sizeof(idx_t));
  for(idx_t i = 0; i < nclus; ++i)
    { 
        surviving_clusters[i] = i; 
    }

  idx_t   merge_count        = 0;
  idx_t   merging_table_size = 1000;
  merge_t *merging_table      = (merge_t*)malloc(sizeof(merge_t)*merging_table_size);
  
  /*Find clusters to be merged*/
  for(idx_t i = 0; i < nclus - 1; ++i)   
    for(idx_t j = i + 1; j < nclus; ++j)   
    {
	switch(borders[i][j].idx != NOBORDER)
	{
                    
	  case 1:		
	    {
	      FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[i]].log_rho_c;
	      FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[i]].log_rho_err;
	      FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[j]].log_rho_c;
	      FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[j]].log_rho_err;
	      FLOAT_TYPE dens_border     = borders[i][j].density;
	      FLOAT_TYPE dens_border_err = borders[i][j].error;
	      
	    if ( is_a_merging( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err, Z ) )
		{
		  
		  if ( merge_count == merging_table_size ) {
		    merging_table_size *= 1.1;
		    merging_table = (merge_t*)realloc( merging_table, sizeof(merge_t) * merging_table_size ); }

		  //int swap = merging_roles( dens1, dens1_err, dens2, dens2_err, dens_border, dens_border_err);
		  idx_t src = j;
		  idx_t trg = i;
		  //switch ( swap )
		  //  {
		  //  case 0: { src = j; trg = i;} break;
		  //  case 1: { src = i; trg = j;} break;
		  //  }

		  merging_table[merge_count].source = src;
		  merging_table[merge_count].target = trg;
		  merging_table[merge_count].density = borders[src][trg].density;
          ++merge_count;
		}
	      break;
	    }
        default:
	    {
	      break;
	    }
            
	  }
      }

  qsort( (void*)merging_table, merge_count, sizeof(merge_t), compare_merging_density);
  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinding merges:   %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
  
    for( idx_t m = 0; m < merge_count; m++ )
    {
      
        #define src surviving_clusters[merging_table[m].source]
        #define trg surviving_clusters[merging_table[m].target]
        //printf("Found: %lu, %lu which now is %lu, %lu\n",merging_table[m].source, merging_table[m].target, src,trg);

        //int re_check = ( (src != merging_table[m].source) || (trg != merging_table[m].target) );
	//if(re_check)
	{
		idx_t new_src = (src < trg) ? src : trg;
		idx_t new_trg = (src < trg) ? trg : src;

                //pick who am I

                FLOAT_TYPE dens1           = dpInfo[cluster->centers.data[new_src]].log_rho_c;
                FLOAT_TYPE dens1_err       = dpInfo[cluster->centers.data[new_src]].log_rho_err;
                FLOAT_TYPE dens2           = dpInfo[cluster->centers.data[new_trg]].log_rho_c;
                FLOAT_TYPE dens2_err       = dpInfo[cluster->centers.data[new_trg]].log_rho_err;

                FLOAT_TYPE dens_border     = borders[new_src][new_trg].density;
                FLOAT_TYPE dens_border_err = borders[new_src][new_trg].error;

                int i_have_to_merge = is_a_merging(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err,Z);            
                switch (i_have_to_merge && src != trg)
                {
                case 1:
                    {
                        int side = merging_roles(dens1,dens1_err,dens2,dens2_err,dens_border,dens_border_err);
                        if(!side)
                        {
                            idx_t tmp;
                            tmp = new_src;
                            new_src = new_trg;
                            new_trg = tmp;
                        }

                        borders[new_src][new_trg] = border_null;
                        borders[new_trg][new_src] = border_null;
                        //printf("Merging %lu into %lu\n",new_src,new_trg);
                        fix_borders_A_into_B(new_src,new_trg,borders,nclus);
                        Merge_A_into_B ( surviving_clusters, new_src, new_trg, nclus );	  
                    }
                    break;
                
                default:
                    break;
                }
	}
        
        #undef src
        #undef trg
    }

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tCluster merging:  %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 
  #endif
  
    /*Finalize clustering*/
    /*Acutally copying */
    lu_dynamicArray tmp_centers;
    lu_dynamicArray tmp_cluster_idx;


    DynamicArray_Init(&tmp_centers);
    DynamicArray_Init(&tmp_cluster_idx);

    DynamicArray_Reserve(&tmp_centers, nclus);
    DynamicArray_Reserve(&tmp_cluster_idx, nclus);

    idx_t final_cluster_count = 0;

    idx_t* old_to_new = (idx_t*)malloc(nclus*sizeof(idx_t));
    idx_t incremental_k = 0;
    for(idx_t i = 0; i < nclus; ++i)
    {
        
        if(surviving_clusters[i] == i){
            DynamicArray_pushBack(&tmp_centers, cluster->centers.data[i]);
            DynamicArray_pushBack(&tmp_cluster_idx, i);
            old_to_new[i] = incremental_k;
            ++incremental_k;
            ++final_cluster_count;
        }
    }

    //fill the rest of old_to_new
    for(idx_t i = 0; i < nclus; ++i)
    {
        if(surviving_clusters[i] != i){
            idx_t cidx_to_copy_from = surviving_clusters[i];
            old_to_new[i] = old_to_new[cidx_to_copy_from];
        }
    }

    /*allocate auxiliary pointers to store results of the finalization of the procedure*/

    border_t** tmp_borders      = (border_t**)malloc(final_cluster_count*sizeof(border_t*));
    border_t*  tmp_borders_data = (border_t*)malloc(final_cluster_count*final_cluster_count*sizeof(border_t));

    /*initialize all pointers*/
    for(idx_t i = 0; i < final_cluster_count; ++i)
    {
        tmp_borders[i] = tmp_borders_data + i*final_cluster_count;
    }

    /*Fix cluster assignment*/
    #pragma omp parallel for
    for(idx_t i = 0; i < cluster -> n; ++i)
    {
        dpInfo[i].is_center = 0;
        int old_cidx = dpInfo[i].cluster_idx;
		if(old_cidx != -1 )
		{
			dpInfo[i].cluster_idx = old_to_new[old_cidx];
		}
    }

    
    #pragma omp parallel for
    for(idx_t c = 0; c < final_cluster_count; ++c)
    {
        idx_t c_idx = tmp_cluster_idx.data[c];
        for(idx_t d = c; d < final_cluster_count; ++d)
        {
            idx_t c_jdx = tmp_cluster_idx.data[d];
            tmp_borders[c][d].density = borders[c_idx][c_jdx].density;
            tmp_borders[d][c].density = borders[c_idx][c_jdx].density;

            tmp_borders[c][d].idx = borders[c_idx][c_jdx].idx;
            tmp_borders[d][c].idx = borders[c_idx][c_jdx].idx;


            tmp_borders[c][d].error = borders[c_idx][c_jdx].error;
            tmp_borders[d][c].error = borders[c_idx][c_jdx].error;
        } 
    }

    Clusters_Reset(cluster);
    /*pay attention to the defined borders*/
    /*copy into members*/
    borders = tmp_borders;

    cluster -> __borders_data = tmp_borders_data;

    cluster -> centers = tmp_centers;
    /**
     * Fix center assignment
    */
    for(idx_t i = 0; i < cluster -> centers.count; ++i)
    {
        int idx = cluster -> centers.data[i];
        dpInfo[idx].is_center = 1;
    }
    /*Halo*/
    switch (halo)
    {
    case 1:
	{
		FLOAT_TYPE* max_border_den_array = (FLOAT_TYPE*)malloc(final_cluster_count*sizeof(FLOAT_TYPE));
		#pragma omp parallel
		{
		    #pragma omp for
		    for(idx_t c = 0; c < final_cluster_count; ++c)
		    {
			FLOAT_TYPE max_border_den = -2.;
			for(idx_t d = 0; d < final_cluster_count; ++d)
			{
			    if(tmp_borders[c][d].density > max_border_den)
			    {
				max_border_den = tmp_borders[c][d].density;
			    }
			}
			max_border_den_array[c] = max_border_den;
		    }

		    #pragma omp barrier

		    #pragma omp for
		    for(idx_t i = 0; i < cluster -> n; ++i)
		    {
			int cidx = dpInfo[i].cluster_idx;
			if(cidx != -1)
			{
				int halo_flag = dpInfo[i].log_rho_c < max_border_den_array[cidx]; 
				dpInfo[i].cluster_idx = halo_flag ? -1 : cidx;
			}
		    }
		}
		free(max_border_den_array);
	}
        break;
    
    default:
        break;
    }    

    /*free memory and put the correct arrays into place*/
    free(tmp_cluster_idx.data);
    free(merging_table);
    //free(ipos.data);
    //free(jpos.data);
    free(surviving_clusters);
    free(old_to_new);

  #ifdef VERBOSE
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinal operations: %.3lfs\n\n", elapsed);
  #endif

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tFound %lu possible merges\n", (uint64_t)merge_count);
    printf("\tSurviving clusters %lu\n", (uint64_t)final_cluster_count);
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

  #undef  borders  
}


void Heuristic3(Clusters* cluster, Datapoint_info* dpInfo, FLOAT_TYPE Z, int halo)
{
	if(cluster -> UseSparseBorders)
	{
		Heuristic3_sparse(cluster, dpInfo,  Z,  halo);
	}
	else
	{
		Heuristic3_dense(cluster, dpInfo,  Z,  halo);
	}
}



void freeDatapointArray(Datapoint_info* d, size_t n)
{
    free(d);
}

int FloatAndUintSize()
{
	int v = 0;
	int vf = sizeof(FLOAT_TYPE) == 8 ? 1 : 0; 
	int vi = sizeof(idx_t) == 8 ? 1 : 0; 
	v = vf + vi*2;
	return v;
}


void setRhoErrK(Datapoint_info* points, FLOAT_TYPE* rho, FLOAT_TYPE* rhoErr, idx_t* k, size_t n)
{
	for(size_t i = 0; i < n; ++i)
	{
		points[i].log_rho = rho[i];
		points[i].log_rho_err = rhoErr[i];
		points[i].g = points[i].log_rho - points[i].log_rho_err;
		points[i].kstar = k[i];
	}
	return;
}


Datapoint_info* computeDensityFromImg(FLOAT_TYPE* vals, int* mask, int nrows, int ncols, int rmax) { //use it to prune isolated pixels
    struct timespec start_tot, finish_tot;
    double elapsed_tot;

    printf("Density estimation from image\n");
    //printf("Got: nrows %lu ncols %lu radius %lu\n", nrows, ncols, rmax);
    printf("Got: nrows %d ncols %d radius %d\n", nrows, ncols, rmax);
    clock_gettime(CLOCK_MONOTONIC, &start_tot);
	
    Datapoint_info* p = (Datapoint_info*)malloc(nrows*ncols*sizeof(Datapoint_info));
	int* tmp_mask = (int*)malloc(nrows*ncols*sizeof(int));
	for(int idx = 0; idx < nrows*ncols; ++idx) tmp_mask[idx] = 1;

	#pragma omp parallel for schedule(dynamic)
	for(int i = 0; i < nrows; ++i)			
	for(int j = 0; j < ncols; ++j)			
	{
		int n = 0;
		FLOAT_TYPE avg = 0;
		FLOAT_TYPE var = 0;
		int r = 1;
		if(mask[i*ncols + j])
		{
			for(r = 1; r < rmax; ++r)
			{
				FLOAT_TYPE tmp_avg = avg;
				FLOAT_TYPE tmp_var = var;
				int 	   tmp_n   = n;

				n = 0;
				avg = 0;
				var = 0;
				int jjmin = j - r > 0 			? j - r : 0;  
				int jjmax = j + r + 1 < ncols 	? j + r + 1 : ncols;  

				int iimin = i - r > 0 	 		? i - r : 0;  
				int iimax = i + r + 1 < nrows 	? i + r + 1 : nrows;  

				for(int ii = iimin; ii < iimax; ++ii)
					for(int jj = jjmin; jj < jjmax; ++jj)
					{
						int index = ii*ncols + jj;
						n 	+= (mask[index] ? 1 : 0);	
						avg += (mask[index] ? vals[index] : 0.);	
						var += (mask[index] ? vals[index]*vals[index] : 0.);	
					}
				if(n > 1)
				{
					avg = avg/(float_t)n;
					var = var/(float_t)(n-1) - avg*avg*(float_t)n/(float_t)(n-1); 	
					var = var/(float_t)(n);
				}

				if(tmp_n > 2)
				{
					float_t sigma_comp = sqrt(var + tmp_var);
					//float_t sigma_comp = sqrt(var);
					int compatibilityCondition = (avg - tmp_avg < sigma_comp) && (tmp_avg - avg < sigma_comp);
					//int compatibilityCondition = (var < tmp_var);

					if(!compatibilityCondition)
					{
						break;
						var = tmp_var;
						avg = tmp_avg;
					}


				}

		}
		}
		if(n > 1 && mask[i*ncols + j])
		{
			p[i*ncols + j].log_rho = log(avg);
			p[i*ncols + j].log_rho_err = sqrt(var)/avg;
			p[i*ncols + j].g = p[i*ncols + j].log_rho - p[i*ncols + j].log_rho_err;
			p[i*ncols + j].kstar = (idx_t)r;
			p[i*ncols + j].array_idx = i*ncols + j;
			p[i*ncols + j].cluster_idx = -1;
		}
		else
		{
			tmp_mask[i*ncols + j] = 0;
			p[i*ncols + j].log_rho = -99999.;
			p[i*ncols + j].g = -99999.; 
		}

	}
	for(int idx = 0; idx < nrows*ncols; ++idx) mask[idx] = mask[idx] && tmp_mask[idx];
	free(tmp_mask);

    clock_gettime(CLOCK_MONOTONIC, &finish_tot);
    elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
    elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
    printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

	return p;
}

#define RED(x)      (3 * x) 
#define GREEN(x)    (3 * x + 1) 
#define BLUE(x)     (3 * x + 2)



void tiny_colorize(
        const char* fname, 
        Datapoint_info* dp, 
        FLOAT_TYPE* data, 
        uint32_t n_clusters, 
        uint32_t og_width, 
        uint32_t og_height, 
        uint32_t target_width,
        uint32_t target_height)
{
    FLOAT_TYPE a = 1;
    FLOAT_TYPE c = 0.99;

    unsigned char* img_buffer = (unsigned char*)malloc(3 * (og_width * 2) * og_height);

    uint32_t stride = 2 * og_width;
    uint32_t offset = og_width;
   
    unsigned char* palette = (unsigned char*)malloc(3 * (n_clusters + 1));
    /* generate palette */
    palette[RED(0)]     = 255;
    palette[GREEN(0)]   = 255;
    palette[BLUE(0)]    = 255;
    for(uint32_t i = 1; i < n_clusters + 1; ++i)
    {
        palette[RED(i)]     = (unsigned char)(rand() % 256);
        palette[GREEN(i)]   = (unsigned char)(rand() % 256);
        palette[BLUE(i)]    = (unsigned char)(rand() % 256);
    }

    for(uint32_t i = 0; i < og_height; ++i)
        for(uint32_t j = 0; j < og_width; ++j)
        {
            uint32_t idx = i * stride + j;
            int cluster_idx = dp[i * og_width + j].cluster_idx + 1;

            img_buffer[RED(idx)] = palette[RED(cluster_idx)]; 
            img_buffer[GREEN(idx)] = palette[GREEN(cluster_idx)]; 
            img_buffer[BLUE(idx)] = palette[BLUE(cluster_idx)]; 
        }

    FLOAT_TYPE data_max = -9999999.f;
    FLOAT_TYPE data_min =  9999999.f;
    for(uint32_t i = 0; i < og_height; ++i)
        for(uint32_t j = 0; j < og_width; ++j)
        {
            data_max = MAX(data_max, data[i * og_width + j] );
            data_min = MIN(data_min, data[i * og_width + j] );
        }

    FLOAT_TYPE delta = 1./(data_max - data_min);
    for(uint32_t i = 0; i < og_height; ++i)
        for(uint32_t j = 0; j < og_width; ++j)
        {
            uint32_t idx = i * stride + j + offset;
            FLOAT_TYPE val = data[i * og_width + j];
            FLOAT_TYPE vnorm = (val - data_min)*delta;

            //unsigned char v = (unsigned char)(a * vnorm/(c * vnorm + (a - c))*255.);
            //unsigned char v = (unsigned char)((val - data_min)*delta*255.);
            //unsigned char v = (unsigned char)(0.5*(tanh(100 * vnorm - 2.5) + 1)*255.);
            FLOAT_TYPE v = (a * vnorm/(c * vnorm + (a - c)));

            img_buffer[RED(idx)]    = (unsigned char)(v * 255); 
            img_buffer[GREEN(idx)]  = (unsigned char)(v * 255); 
            img_buffer[BLUE(idx)]   = (unsigned char)(v * 255); 

            
        }
    
     unsigned char* out_pixels = stbir_resize_uint8_srgb( img_buffer,  2 * og_width,  og_height,  0,
                                                  NULL, 2 * target_width, target_height, 0,
                                                  STBIR_RGB);
     stbi_write_png(fname, 2 * target_width, target_height, 3, out_pixels, 0);
   
     free(out_pixels);
     free(img_buffer);
     free(palette);
    
}

#undef BLUE
#undef GREEN
#undef RED

void export_cluster_assignment(Datapoint_info* points, int* labels, idx_t n)
{
	for(idx_t i = 0; i < n; ++i) labels[i] = points[i].cluster_idx;
}

