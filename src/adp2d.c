//include "../include/read_fof_snapshot.h"
#include "../include/adp2d.h"
#include <float.h>
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

#define MAX(x,y) (x > y ? x : y)
#define MIN(x,y) (x < y ? x : y)

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
    float_t aa = *((float_t*)a);
    float_t bb = *((float_t*)b);
    return (aa > bb ) - (aa < bb); 
}



float_t avg(const float_t * x, const idx_t n)
{
    float_t f = 0;
    for(idx_t i = 0; i < n; ++i)
    {
        f += x[i];
    }
    return f/(float_t)n;
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
    return (g1 < g2) - (g1 > g2);
}

void computeCorrection(Datapoint_info* dpInfo, int* mask, idx_t n, float_t Z)
{
    /*****************************************************************************
     * Utility function, find the minimum value of the density of the datapoints *
     * and shift them up in order to further work with values greater than 0     *
     *****************************************************************************/
    float_t min_log_rho = FLT_MAX;
    

    #pragma omp parallel
    {
        float_t thread_min_log_rho = FLT_MAX;
        #pragma omp for
        for(idx_t i = 0; i < n; ++i)
        {
            float_t tmp = dpInfo[i].log_rho - Z*dpInfo[i].log_rho_err;
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

typedef struct {
    int lb_row;
    int lb_col;
    int ub_row;
    int ub_col;
    int label;
} bounding_box_t;

Clusters adpWrapper(Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols, int min_size, float Z, bool halo, bool split_per_thread) 
{
    #ifdef OPENMP
        int machine_threads = omp_get_num_threads(); 
    #else
        int machine_threads = 1; 
    #endif
    Clusters c = {0};
    // this has to be called 
    printf("Computing correction\n");
    computeCorrection(dpInfo, mask, nrows * ncols, Z);

    // then see if I can split the compute per thread
    //
    int aa = 0;
    for(int i = 0; i < nrows * ncols; ++i)
    {
        if(aa < mask[i]) aa = mask[i];
    }
    if(split_per_thread) 
    {
        idx_t n_labels_mask = 0;

        // compute number of clusters
        printf("Computing n detection\n");
        #pragma omp parallel 
        {
            idx_t pvt_max = 0;
            #pragma omp for 
            for(int i = 0; i < nrows; ++i)
                for(int j = 0; j < ncols; ++j)
                {
                    pvt_max = MAX(mask[i*ncols + j], pvt_max);
                }

            #pragma omp critical (max_n_labels)
            {
                n_labels_mask = MAX(n_labels_mask, pvt_max);
            }
        }

        // since they go from 0 to N
        n_labels_mask++;

        printf("Computing bouding boxes\n");
        bounding_box_t* bounding_boxes = (bounding_box_t* )malloc(n_labels_mask * sizeof(bounding_box_t));

        for(int i = 0; i < n_labels_mask; ++i)
        {
            bounding_boxes[i].label = i;
            bounding_boxes[i].lb_row = nrows;
            bounding_boxes[i].lb_col = ncols;
            bounding_boxes[i].ub_row = 0;
            bounding_boxes[i].ub_col = 0;
        }

        #pragma omp parallel 
        {
            bounding_box_t* pvt_bounding_boxes = (bounding_box_t* )malloc(n_labels_mask * sizeof(bounding_box_t));

            for(int lab = 0; lab < n_labels_mask; ++lab)
            {
                pvt_bounding_boxes[lab].label = lab;
                pvt_bounding_boxes[lab].lb_row = nrows;
                pvt_bounding_boxes[lab].lb_col = ncols;
                pvt_bounding_boxes[lab].ub_row = 0;
                pvt_bounding_boxes[lab].ub_col = 0;
            }

            #pragma omp for 
            for(int i = 0; i < nrows; ++i)
                for(int j = 0; j < ncols; ++j)
                {
                    int lab = mask[i * ncols + j];
                    if(lab)
                    {
                        pvt_bounding_boxes[lab].lb_row = MIN(i, pvt_bounding_boxes[lab].lb_row);
                        pvt_bounding_boxes[lab].lb_col = MIN(j, pvt_bounding_boxes[lab].lb_col);
                        pvt_bounding_boxes[lab].ub_row = MAX(i, pvt_bounding_boxes[lab].ub_row);
                        pvt_bounding_boxes[lab].ub_col = MAX(j, pvt_bounding_boxes[lab].ub_col);
                    }
                }

            #pragma omp critical 
            {
                for(int lab = 0; lab < n_labels_mask; ++lab)
                {
                    bounding_boxes[lab].lb_row = MIN(bounding_boxes[lab].lb_row, pvt_bounding_boxes[lab].lb_row);
                    bounding_boxes[lab].lb_col = MIN(bounding_boxes[lab].lb_col, pvt_bounding_boxes[lab].lb_col);
                    bounding_boxes[lab].ub_row = MAX(bounding_boxes[lab].ub_row, pvt_bounding_boxes[lab].ub_row);
                    bounding_boxes[lab].ub_col = MAX(bounding_boxes[lab].ub_col, pvt_bounding_boxes[lab].ub_col);
                }
            }
            free(pvt_bounding_boxes);

        }

        // avoid 0

        printf("Computing clustering\n");
        int* clusters_per_box = (int*)calloc(n_labels_mask, sizeof(int));

        for(int lab = 0; lab < 15; ++lab)
        {
            bounding_box_t box = bounding_boxes[lab];

        }
        // compute bounding boxes 
        #pragma omp parallel for schedule(dynamic)
        for(int lab = 1; lab < n_labels_mask; ++lab)
        {
            bounding_box_t box = bounding_boxes[lab];

            // check if the box have been setted
            bool box_is_valid = (box.lb_row < nrows) && 
                                (box.lb_col < ncols) &&
                                (box.ub_row > 0) &&
                                (box.ub_col > 0);


            if(box_is_valid)
            {
                // printf("Processing box row [%d %d] col [%d %d]\n",  box.lb_row, box.ub_row, 
                //                                                     box.lb_col, box.ub_col);
                int ncols_box = (box.ub_col - box.lb_col) + 1;
                int nrows_box = (box.ub_row - box.lb_row) + 1;
                int n_pixels_in_box =  nrows_box * ncols_box;

                int* tmp_mask = (int*)malloc(ncols_box * nrows_box * sizeof(int));
                Datapoint_info* tmp_dp = (Datapoint_info*)malloc(ncols_box * nrows_box * sizeof(Datapoint_info));


                // copy datapoints and everything into a temporary array
                // on these adp should run
                for(int i = 0; i < nrows_box; ++i)
                    for(int j = 0; j < ncols_box; ++j)
                    {
                        int ii = i + box.lb_row;
                        int jj = j + box.lb_col;

                        memcpy(tmp_dp + i*ncols_box + j, dpInfo + ii*ncols + jj, sizeof(Datapoint_info));
                        tmp_dp[i*ncols_box + j].array_idx = i*ncols_box + j;

                        bool mask_is_lab = mask[ii*ncols + jj] == lab;
                        tmp_mask[i*ncols_box + j] = mask_is_lab ? lab : 0;
                    }

                Clusters c_tmp = Heuristic1(tmp_dp, tmp_mask, nrows_box, ncols_box, 1, false);
                Clusters_allocate(&c_tmp, false);
                Heuristic2(&c_tmp, tmp_dp, tmp_mask, nrows_box, ncols_box, 1, false);
                Heuristic3(&c_tmp, tmp_dp, Z, halo, 1, false);

                // if(c_tmp.n > 1000) printf("c_tmp.n %lu\n", c_tmp.n);
                clusters_per_box[lab] = c_tmp.centers.count;

                //copy cluster assignment
                for(int i = 0; i < nrows_box; ++i)
                    for(int j = 0; j < ncols_box; ++j)
                    {
                        int ii = i + box.lb_row;
                        int jj = j + box.lb_col;
                        if(mask[ii * ncols + jj] == lab) dpInfo[ii * ncols + jj].cluster_idx = tmp_dp[i * ncols_box + j].cluster_idx;
                    }

                Clusters_free(&c_tmp);
                free(tmp_mask);
                free(tmp_dp);
            }
        }

        // compute correct cluster indices

        int n_clusters = 0;
        // exclusive prefix sum
        //
        for(int i = 0; i < n_labels_mask; ++i)
        {
            int tmp_n = clusters_per_box[i];
            clusters_per_box[i] = n_clusters;
            n_clusters += tmp_n; 
        }

        #pragma omp parallel for
        for(idx_t i = 0; i < nrows * ncols; ++i)
        {
             int lab = mask[i];
             dpInfo[i].cluster_idx += clusters_per_box[lab];
        }
        
        printf("Final n clusters %d\n", n_clusters);
        c.centers.count = n_clusters;
        free(bounding_boxes);
        free(clusters_per_box);

    }
    else 
    {
        c = Heuristic1(dpInfo, mask, nrows, ncols, machine_threads, true);
        Clusters_allocate(&c, true);
        Heuristic2(&c, dpInfo, mask, nrows, ncols, machine_threads, true);
        Heuristic3(&c, dpInfo, Z, halo, machine_threads, true);
    }

    // filter out clusters less than min size

    int* pixel_count_per_cluster = (int*)calloc(c.centers.count, sizeof(int));
    int* new_labels              = (int*)calloc(c.centers.count, sizeof(int));

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < nrows; ++i)
        for(int j = 0; j < ncols; ++j)
        {
            int cidx = dpInfo[i*ncols + j].cluster_idx;
            if(cidx != -1)
            {
                #pragma omp atomic update
                pixel_count_per_cluster[cidx]++;
            }
        }
    
    printf("Pruning clusters smaller than %d pixels of area\n", min_size);
    int label_count = 0;
    for(int i = 0; i < c.centers.count; ++i)
    {
        if(pixel_count_per_cluster[i] > min_size)
        {
            new_labels[i] = label_count;
            ++label_count;
        }
        else
        {
            new_labels[i] = -1;
        }
    }

    printf("Final cluster count after pruning %d\n", label_count);

    #pragma omp parallel for 
    for(int i = 0; i < nrows; ++i)
        for(int j = 0; j < ncols; ++j)
        {
            int cidx = dpInfo[i*ncols + j].cluster_idx;
            dpInfo[i*ncols + j].cluster_idx = cidx != -1 ? new_labels[cidx] : cidx;
        }

    c.centers.count = label_count;

    free(pixel_count_per_cluster);
    free(new_labels);
    return c;
}


//Clusters Heuristic1(Datapoint_info* dpInfo, int* mask, int nrows, int ncols)- 
Clusters Heuristic1(Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols, int num_threads, bool verbose)
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

    if(verbose) printf("H1: Preliminary cluster assignment\n");
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


    if(verbose) clock_gettime(CLOCK_MONOTONIC, &start);

    #pragma omp parallel for num_threads(num_threads)
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
        float_t gi = dpInfo[i*ncols + j].g;
        dpInfo[i*ncols + j].is_center = mask[i*ncols + j] ? 1 : 0;
        dpInfo[i*ncols + j].cluster_idx = -1;
        //printf("%lf\n",p -> g);
		int jjmin = j - r > 0 			    ? j - r : 0;  
		int jjmax = j + r + 1 < (int)ncols 	? j + r + 1 : (int)ncols;  

		int iimin = i - r > 0 	 		    ? i - r : 0;  
		int iimax = i + r + 1 < (int)nrows 	? i + r + 1 : (int)nrows;  
		
		if(mask[i*ncols + j])
		{
			for(int ii = iimin; ii < iimax; ++ii)
			for(int jj = jjmin; jj < jjmax; ++jj)
			{
				idx_t ngbh_index = (idx_t)ii*ncols + jj; 
				float_t gj = dpInfo[ngbh_index].g;
				if(gj > gi && mask[ngbh_index] && ((int)ngbh_index != (int)(i*ncols + j) )){
					dpInfo[i*ncols + j].is_center = 0;
					break;
				}
			}
		}
    }

    for(int i = 0; i < (int)nrows; ++i)
        for(int j = 0; j < (int)ncols; ++j)
        {
            if(dpInfo[i*ncols + j].is_center && mask[i*ncols + j]){
                DynamicArray_pushBack(&allCenters, i*ncols + j);
            }
        }


    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinding putative centers: %.3lfs\n",elapsed);
        clock_gettime(CLOCK_MONOTONIC, &start);

    }

	qsort(dpInfo_ptrs, nrows*ncols, sizeof(Datapoint_info*), cmpPP);

    idx_t * to_remove = (idx_t*)malloc(allCenters.count*sizeof(idx_t));
    for(idx_t c = 0; c < allCenters.count; ++c) {to_remove[c] = MY_SIZE_MAX;}

	idx_t* to_remove_mask = (idx_t*)malloc(nrows*ncols*sizeof(idx_t));
    for(idx_t p = 0; p < nrows*ncols; ++p) {to_remove_mask[p] = MY_SIZE_MAX;}

	
    #pragma omp parallel shared(to_remove_mask) num_threads(num_threads)
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
						
                        // this critical makes jupyter crash ostia
						#pragma omp critical 
						{
							ppp = to_remove_mask[jidx];
							if(ppp != MY_SIZE_MAX)					
                            {
                                to_remove_mask[jidx] = pp.g > dpInfo[ppp].g ? pp.array_idx : ppp;
                            }
                            else
                            {
                                to_remove_mask[jidx] = pp.array_idx; 
                            }
						}
						
					}
				}
			}
		}
	}
    
    

    for(idx_t p = 0; p < allCenters.count; ++p)
    {
        idx_t i = allCenters.data[p];
        int e = 0;
        //float_t gi = dpInfo[i].g;
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


	free(to_remove);
	free(to_remove_mask);


    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinding actual centers:   %.3lfs\n",elapsed);
        clock_gettime(CLOCK_MONOTONIC, &start);
    }


    //idx_t nclusters = 0;


    /*****************************************************************************
     * Sort all the dpInfo based on g and then perform the cluster assignment *
     * in asceding order                                                         *
     * UPDATE: dpInfo already sorted                                          *
     *****************************************************************************/
                                                                                

    //qsort(dpInfo_ptrs, n, sizeof(Datapoint_info*), cmpPP);
	
	
	idx_t* fromWho = (idx_t*)malloc(nrows*ncols*sizeof(idx_t));
    for(idx_t pidx = 0; pidx < nrows*ncols; ++pidx) fromWho[pidx] = SIZE_MAX;

	
	
	
	#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
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
                float_t gmax = -99999.;               
                idx_t gm_index = 0;

				for(int ii = iimin; ii < iimax; ++ii)
				for(int jj = jjmin; jj < jjmax; ++jj)
                {
                    idx_t ngbh_index = ii*ncols + jj;
                    for(idx_t m = 0; m < removedCenters.count; ++m)
                    {
                        float_t gcand = dpInfo[max_rho.data[m]].g;
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
	
    

    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tTentative clustering:     %.3lfs\n",elapsed);
        clock_gettime(CLOCK_MONOTONIC, &start);

    }

    free(dpInfo_ptrs);
    free(max_rho.data);
    free(removedCenters.data);
    free(allCenters.data);
    free(fromWho);

    Clusters c_all;
    c_all.centers = actualCenters;


    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish);
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinalizing clustering:    %.3lfs\n",elapsed);
        printf("\n");

        clock_gettime(CLOCK_MONOTONIC, &finish_tot);
        elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
        elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;


        printf("\tFound %lu clusters\n",(uint64_t)actualCenters.count);
        printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    }


    c_all.n = nrows*ncols;
    return c_all;
}

void Heuristic2(Clusters* cluster, Datapoint_info* dpInfo, int* mask, size_t nrows, size_t ncols, int num_threads, bool verbose)
{

    #define borders cluster->borders

    struct timespec start_tot, finish_tot;
    double elapsed_tot;
    //idx_t n = cluster -> n;

    if(verbose) printf("H2: Finding border points\n");
    clock_gettime(CLOCK_MONOTONIC, &start_tot);


    idx_t nclus = cluster->centers.count; 
    //idx_t max_k = dpInfo[0].ngbh.N;


    for(int i = 0; i < (int)nrows; i++)
    for(int j = 0; j < (int)ncols; j++)
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
				for(int ii = iimin; ii < iimax; ii++)
				for(int jj = jjmin; jj < jjmax; jj++)
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

    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish_tot);
        elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
        elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
        printf("\tTotal time: %.3lfs\n\n", elapsed_tot);

    }

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
  float_t DensA = ((merge_t*)A)->density;
  float_t DensB = ((merge_t*)B)->density;

  return - ( DensA > DensB) + (DensA < DensB);
}


static inline int is_a_merging( 
                float_t dens1, float_t dens1_err,
                float_t dens2, float_t dens2_err,
                float_t dens_border, float_t dens_border_err,
                float_t Z)
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
  float_t a1 = dpInfo[cluster->centers.data[i]].log_rho_c - border_density[i][j];
  float_t a2 = dpInfo[cluster->centers.data[j]].log_rho_c - border_density[i][j];
  
  float_t e1 = Z*(dpInfo[cluster->centers.data[i]].log_rho_err + border_err[i][j]);
  float_t e2 = Z*(dpInfo[cluster->centers.data[j]].log_rho_err + border_err[i][j]);
  */

  float_t a1 = dens1 - dens_border;
  float_t a2 = dens2 - dens_border;

  float_t e1 = Z*(dens1_err + dens_border_err);
  float_t e2 = Z*(dens2_err + dens_border_err);

  return (a1 < e1 || a2 < e2);
}


int merging_roles( float_t dens1, float_t dens1_err,
			  float_t dens2, float_t dens2_err,
			  float_t dens_border, float_t dens_border_err )
{
      
  float_t c1 = (dens1 - dens_border) / (dens1_err + dens_border_err); 
  float_t c2 = (dens2 - dens_border) / (dens2_err + dens_border_err);
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

void Heuristic3_sparse(Clusters* cluster, Datapoint_info* dpInfo, float_t Z, int halo, int num_threads, bool verbose)
{
  if(verbose) printf("Using sparse implementation\n");
  #define borders cluster->borders

  struct timespec start_tot, finish_tot;
  double elapsed_tot;

  struct timespec start, finish;
  double elapsed;

  if(verbose) printf("H3: Merging clusters\n");
  clock_gettime(CLOCK_MONOTONIC, &start_tot);
  if(verbose) clock_gettime(CLOCK_MONOTONIC, &start); 

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
		      float_t dens1           = dpInfo[cluster->centers.data[b.i]].log_rho_c;
		      float_t dens1_err       = dpInfo[cluster->centers.data[b.i]].log_rho_err;
		      float_t dens2           = dpInfo[cluster->centers.data[b.j]].log_rho_c;
		      float_t dens2_err       = dpInfo[cluster->centers.data[b.j]].log_rho_err;
		      float_t dens_border     = b.density;
		      float_t dens_border_err = b.error;
	      
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

  if(verbose)
  {
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinding merges:   %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 

  }
  
  
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

                float_t dens1           = dpInfo[cluster->centers.data[new_src]].log_rho_c;
                float_t dens1_err       = dpInfo[cluster->centers.data[new_src]].log_rho_err;
                float_t dens2           = dpInfo[cluster->centers.data[new_trg]].log_rho_c;
                float_t dens2_err       = dpInfo[cluster->centers.data[new_trg]].log_rho_err;

		//borders get
		SparseBorder_t b 	   = SparseBorder_get(cluster, new_src, new_trg);
                float_t dens_border     = b.density;
                float_t dens_border_err = b.error;

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

  if(verbose)
  {
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tCluster merging:  %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 

  }
  
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

    
    #pragma omp parallel for num_threads(num_threads)
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
		float_t* max_border_den_array = (float_t*)malloc(final_cluster_count*sizeof(float_t));
		#pragma omp parallel
		{
		    #pragma omp for
		    for(idx_t c = 0; c < final_cluster_count; ++c)
		    {
				float_t max_border_den = -2.;
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

    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish); 
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinal operations: %.3lfs\n\n", elapsed);

        clock_gettime(CLOCK_MONOTONIC, &finish_tot);
        elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
        elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
        printf("\tFound %lu possible merges\n",(uint64_t)merge_count);
        printf("\tSurviving clusters %lu\n",(uint64_t)final_cluster_count);
        printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    }


  #undef  borders  
}


void Heuristic3_dense(Clusters* cluster, Datapoint_info* dpInfo, float_t Z, int halo, int num_threads, bool verbose)
{
  if(verbose) printf("Using dense implementation\n");
  #define borders cluster->borders

  struct timespec start_tot, finish_tot;
  double elapsed_tot;

  struct timespec start, finish;
  double elapsed;

  if(verbose) printf("H3: Merging clusters\n");
  if(verbose) clock_gettime(CLOCK_MONOTONIC, &start_tot);
  if(verbose) clock_gettime(CLOCK_MONOTONIC, &start); 

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
	      float_t dens1           = dpInfo[cluster->centers.data[i]].log_rho_c;
	      float_t dens1_err       = dpInfo[cluster->centers.data[i]].log_rho_err;
	      float_t dens2           = dpInfo[cluster->centers.data[j]].log_rho_c;
	      float_t dens2_err       = dpInfo[cluster->centers.data[j]].log_rho_err;
	      float_t dens_border     = borders[i][j].density;
	      float_t dens_border_err = borders[i][j].error;
	      
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
  if(verbose)
  {
	clock_gettime(CLOCK_MONOTONIC, &finish); 
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	printf("\tFinding merges:   %.3lfs\n", elapsed);
	clock_gettime(CLOCK_MONOTONIC, &start); 

  }
  
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

                float_t dens1           = dpInfo[cluster->centers.data[new_src]].log_rho_c;
                float_t dens1_err       = dpInfo[cluster->centers.data[new_src]].log_rho_err;
                float_t dens2           = dpInfo[cluster->centers.data[new_trg]].log_rho_c;
                float_t dens2_err       = dpInfo[cluster->centers.data[new_trg]].log_rho_err;

                float_t dens_border     = borders[new_src][new_trg].density;
                float_t dens_border_err = borders[new_src][new_trg].error;

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

    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish); 
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tCluster merging:  %.3lfs\n", elapsed);
        clock_gettime(CLOCK_MONOTONIC, &start); 

    }
  
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
    #pragma omp parallel for num_threads(num_threads)
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
		float_t* max_border_den_array = (float_t*)malloc(final_cluster_count*sizeof(float_t));
		#pragma omp parallel
		{
		    #pragma omp for
		    for(idx_t c = 0; c < final_cluster_count; ++c)
		    {
			float_t max_border_den = -2.;
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

    if(verbose)
    {
        clock_gettime(CLOCK_MONOTONIC, &finish); 
        elapsed = (finish.tv_sec - start.tv_sec);
        elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("\tFinal operations: %.3lfs\n\n", elapsed);

        clock_gettime(CLOCK_MONOTONIC, &finish_tot);
        elapsed_tot = (finish_tot.tv_sec - start_tot.tv_sec);
        elapsed_tot += (finish_tot.tv_nsec - start_tot.tv_nsec) / 1000000000.0;
        printf("\tFound %lu possible merges\n", (uint64_t)merge_count);
        printf("\tSurviving clusters %lu\n", (uint64_t)final_cluster_count);
        printf("\tTotal time: %.3lfs\n\n", elapsed_tot);
    }

  #undef  borders  
}


void Heuristic3(Clusters* cluster, Datapoint_info* dpInfo, float_t Z, int halo, int num_threads, bool verbose)
{
	if(cluster -> UseSparseBorders)
	{
		Heuristic3_sparse(cluster, dpInfo,  Z,  halo, num_threads, verbose);
	}
	else
	{
		Heuristic3_dense(cluster, dpInfo,  Z,  halo, num_threads, verbose);
	}
}



void freeDatapointArray(Datapoint_info* d, size_t n)
{
    free(d);
}

int FloatAndUintSize()
{
	int v = 0;
	int vf = sizeof(float_t) == 8 ? 1 : 0; 
	int vi = sizeof(idx_t) == 8 ? 1 : 0; 
	v = vf + vi*2;
	return v;
}


void setRhoErrK(Datapoint_info* points, float_t* rho, float_t* rhoErr, idx_t* k, size_t n)
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

void convolve(
    Datapoint_info* p,
    int*            tmp_mask,
    const float_t*  vals,
    const int*      mask,
    const float_t*  kernel,
    int             nrows,
    int             ncols,
    int             rmax,
    bool            use_log)
{
    idx_t row_len = 2 * rmax + 1;

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < nrows; ++i)
        for(int j = 0; j < ncols; ++j)
        {
            float_t avg    = 0;
            float_t var    = 0;
            float_t w_sum  = 0;
            float_t w2_sum = 0;

            int idx = i*ncols + j;

            if(mask[idx])
            {
                int jjmin = j - rmax > 0          ? j - rmax     : 0;
                int jjmax = j + rmax + 1 < ncols   ? j + rmax + 1 : ncols;
                int iimin = i - rmax > 0           ? i - rmax     : 0;
                int iimax = i + rmax + 1 < nrows   ? i + rmax + 1 : nrows;

                for(int ii = iimin; ii < iimax; ++ii)
                    for(int jj = jjmin; jj < jjmax; ++jj)
                    {
                        int index = ii*ncols + jj;
                        if(!mask[index]) continue;

                        float_t w  = kernel[(ii - i + rmax) * row_len + (jj - j + rmax)];
                        avg       += w * vals[index];
                        w_sum     += w;
                    }

                if(w_sum > 0) avg /= w_sum;

                for(int ii = iimin; ii < iimax; ++ii)
                    for(int jj = jjmin; jj < jjmax; ++jj)
                    {
                        int index = ii*ncols + jj;
                        if(!mask[index]) continue;

                        float_t w  = kernel[(ii - i + rmax) * row_len + (jj - j + rmax)];
                        float_t d  = vals[index] - avg;
                        var       += w * d * d;
                        w2_sum    += w * w;
                    }

                float_t denom = w_sum - w2_sum / w_sum;
                if(denom > 0) var /= denom;
                else          var  = 0;
            }

            if(w_sum > 0 && mask[idx])
            {
                p[idx].log_rho     = use_log ? log(avg) : avg;
                p[idx].log_rho_err = use_log ? sqrt(var) / avg : sqrt(avg);
                p[idx].g           = p[idx].log_rho - p[idx].log_rho_err;
                p[idx].kstar       = (idx_t)rmax;
                p[idx].array_idx   = idx;
                p[idx].cluster_idx = -1;
            }
            else
            {
                tmp_mask[idx]      = 0;
                p[idx].log_rho     = -FLT_MAX;
                p[idx].g           = -FLT_MAX;
                p[idx].array_idx   = idx;
                p[idx].cluster_idx = -1;
            }
        }
}

Datapoint_info* computeDensityFromImg(float_t* vals, int* mask, int nrows, int ncols, int rmax,
                                      density_alg_t algorithm, bool use_log, bool use_adaptive_radius, int param) 
{ 
    struct timespec start_tot, finish_tot;
    double elapsed_tot;

    printf("Density estimation from image\n");
    //printf("Got: nrows %lu ncols %lu radius %lu\n", nrows, ncols, rmax);
    printf("Got: nrows %d ncols %d radius %d\n", nrows, ncols, rmax);
    clock_gettime(CLOCK_MONOTONIC, &start_tot);
	
    Datapoint_info* p = (Datapoint_info*)malloc(nrows*ncols*sizeof(Datapoint_info));
	int* tmp_mask = (int*)malloc(nrows*ncols*sizeof(int));
	for(int idx = 0; idx < nrows*ncols; ++idx) tmp_mask[idx] = 1;

    switch(algorithm)
    {
        case MEAN:
            #pragma omp parallel for schedule(dynamic)
            for(int i = 0; i < nrows; ++i)			
                for(int j = 0; j < ncols; ++j)			
                {
                    int n = 0;
                    float_t avg = 0;
                    float_t var = 0;
                    int r = use_adaptive_radius ? 1 : rmax - 1;
                    if(mask[i*ncols + j])
                    {
                        for(r = 1; r < rmax; ++r)
                        {
                            float_t tmp_avg = avg;
                            float_t tmp_var = var;
                            int 	tmp_n   = n;

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

                            if(tmp_n > 1)
                            {
                                float_t sigma_comp = sqrt(var + tmp_var);
                                int compatibilityCondition = (avg - tmp_avg < sigma_comp) && (tmp_avg - avg < sigma_comp);
                                if(!compatibilityCondition)
                                {
                                    var = tmp_var;
                                    avg = tmp_avg;
                                    n   = tmp_n;
                                    break;
                                }
                            }

                        }
                    }
                    if(n > 1 && mask[i*ncols + j])
                    {

                        // Local density contrast test
                        // int idx     = i*ncols + j;
                        // float_t x   = vals[idx];
                        // float_t eps = 1e-8;
                        // float_t rho = (x - avg)/(avg + eps);

                        // p[idx].log_rho     = rho;
                        // p[idx].log_rho_err = sqrt(var);
                        // p[idx].g           = p[idx].log_rho - p[idx].log_rho_err;
                        // p[idx].kstar       = (idx_t)r;
                        // p[idx].array_idx   = idx;
                        // p[idx].cluster_idx = -1;
                        
                        //

                        p[i*ncols + j].log_rho = use_log ? log(avg) : avg;
                        p[i*ncols + j].log_rho_err = use_log ? sqrt(var)/avg : sqrt(avg);
                        p[i*ncols + j].g = p[i*ncols + j].log_rho - p[i*ncols + j].log_rho_err;
                        p[i*ncols + j].kstar = (idx_t)r;
                        p[i*ncols + j].array_idx = i*ncols + j;
                        p[i*ncols + j].cluster_idx = -1;
                    }
                    else
                    {
                        tmp_mask[i*ncols + j] = 0;
                        p[i*ncols + j].log_rho = -FLT_MAX;
                        p[i*ncols + j].g = -FLT_MAX; 
                        p[i*ncols + j].array_idx = i*ncols + j;
                        p[i*ncols + j].cluster_idx = -1;
                    }
                }
            break;
        case MEDIAN:
            #pragma omp parallel 
            {
                float_t* vals_for_median = (float_t*)calloc((2*rmax+1)*(2*rmax+1), sizeof(float_t)); 
                #pragma omp for schedule(dynamic)
                for(int i = 0; i < nrows; ++i)			
                    for(int j = 0; j < ncols; ++j)			
                    {
                        int n = 0;
                        float_t avg = 0;
                        float_t var = 0;
                        if(mask[i*ncols + j])
                        {
                            n = 0;
                            avg = 0;
                            var = 0;
                            int jjmin = j - rmax > 0 			? j - rmax : 0;  
                            int jjmax = j + rmax + 1 < ncols 	? j + rmax + 1 : ncols;  

                            int iimin = i - rmax > 0 	 		? i - rmax : 0;  
                            int iimax = i + rmax + 1 < nrows 	? i + rmax + 1 : nrows;  

                            for(int ii = iimin; ii < iimax; ++ii)
                                for(int jj = jjmin; jj < jjmax; ++jj)
                                {
                                    int index = ii*ncols + jj;
                                    vals_for_median[n] = vals[index];
                                    n 	+= (mask[index] ? 1 : 0);	
                                    avg += (mask[index] ? vals[index] : 0.);	
                                    var += (mask[index] ? vals[index]*vals[index] : 0.);	
                                }
                        }
                        if(n > 1 && mask[i*ncols + j])
                        {
                            avg = avg/(float_t)n;
                            var = var/(float_t)(n-1) - avg*avg*(float_t)n/(float_t)(n-1); 	
                            var = var/(float_t)(n);

                            qsort(vals_for_median, n, sizeof(float_t), cmp);
                            float_t median = 0.;
                            if(n % 2 == 0)
                            {   
                                median = vals_for_median[n/2 - 1] + vals_for_median[n / 2]; 
                                median = median/2.;
                            }
                            else 
                            {
                                median = vals_for_median[n/2]; 
                            }

                            p[i*ncols + j].log_rho = use_log ? log(median) : median;
                            p[i*ncols + j].log_rho_err = use_log ? sqrt(var)/avg : sqrt(var);
                            p[i*ncols + j].g = p[i*ncols + j].log_rho - p[i*ncols + j].log_rho_err;
                            p[i*ncols + j].kstar = (idx_t)rmax;
                            p[i*ncols + j].array_idx = i*ncols + j;
                            p[i*ncols + j].cluster_idx = -1;
                        }
                        else
                        {
                            tmp_mask[i*ncols + j] = 0;
                            p[i*ncols + j].log_rho = -FLT_MAX;
                            p[i*ncols + j].g = -FLT_MAX; 
                            p[i*ncols + j].array_idx = i*ncols + j;
                            p[i*ncols + j].cluster_idx = -1;
                        }

                    }
                free(vals_for_median);
            }
            break;
        case GAUSSIAN:
            {
                float_t* gaussian_weights = (float_t*)calloc((2 * rmax + 1) * (2 * rmax + 1), sizeof(float_t));
                // produce gaussian_weights
                idx_t row_len = 2 * rmax + 1;
                idx_t center = rmax;
                float_t sigma_sq = (float_t)param * (float_t)param;

                float_t tot_weights = 0.;
                for(idx_t i = 0; i < 2 * rmax + 1; ++i)
                    for(idx_t j = 0; j < 2 * rmax + 1; ++j)
                    {
                        float_t dist_sq = (float_t)(i - center) * (float_t)(i - center) + (float_t)(j - center) * (float_t)(j - center);
                        gaussian_weights[i * row_len + j] = exp(- dist_sq/sigma_sq);
                        tot_weights += gaussian_weights[i * row_len + j];
                    }

                // normalize kernel
                for(idx_t i = 0; i < (2*rmax + 1) * (2*rmax + 1); ++i) gaussian_weights[i] = gaussian_weights[i]/tot_weights;

                convolve(p, tmp_mask, vals, mask, gaussian_weights, nrows, ncols, rmax, use_log);

                free(gaussian_weights);
            }
            break;
        case SPLINE:
            {
                // Precompute SPH cubic spline weights 
                float_t* sph_weights = (float_t*)calloc((2 * rmax + 1) * (2 * rmax + 1), sizeof(float_t));
                idx_t row_len = 2 * rmax + 1;
                idx_t center  = rmax;
                float_t tot_weights = 0.;

                for(idx_t i = 0; i < 2 * rmax + 1; ++i)
                    for(idx_t j = 0; j < 2 * rmax + 1; ++j)
                    {
                        float_t dist = sqrt((float_t)(i - center)*(i - center)
                                        + (float_t)(j - center)*(j - center));
                        float_t q = dist / (float_t)param;   // normalised distance

                        float_t w = 0;
                        if     (q < 1.f) w = 1.f - 1.5f*q*q + 0.75f*q*q*q;
                        else if(q < 2.f) w = 0.25f * (2.f - q)*(2.f - q)*(2.f - q);
                        // q >= 2: w = 0 (compact support)

                        sph_weights[i * row_len + j] = w;
                        tot_weights += w;
                    }

                // Normalize
                for(idx_t i = 0; i < (2*rmax+1)*(2*rmax+1); ++i) sph_weights[i] /= tot_weights;

                convolve(p, tmp_mask, vals, mask, sph_weights, nrows, ncols, rmax, use_log);

                free(sph_weights);
            }
            break;
        default:
            printf("Select a valid algorithm `MEAN` or `MEDIAN` for the density computation\n");
            break;
    }
	for(int idx = 0; idx < nrows*ncols; ++idx) mask[idx] = mask[idx] * tmp_mask[idx];
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
        float_t* data, 
        uint32_t n_clusters, 
        uint32_t og_width, 
        uint32_t og_height, 
        uint32_t target_width,
        uint32_t target_height)
{
    float_t a = 1;
    float_t c = 0.99;

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

    float_t data_max = -9999999.f;
    float_t data_min =  9999999.f;
    for(uint32_t i = 0; i < og_height; ++i)
        for(uint32_t j = 0; j < og_width; ++j)
        {
            data_max = MAX(data_max, data[i * og_width + j] );
            data_min = MIN(data_min, data[i * og_width + j] );
        }

    float_t delta = 1./(data_max - data_min);
    for(uint32_t i = 0; i < og_height; ++i)
        for(uint32_t j = 0; j < og_width; ++j)
        {
            uint32_t idx = i * stride + j + offset;
            float_t val = data[i * og_width + j];
            float_t vnorm = (val - data_min)*delta;

            //unsigned char v = (unsigned char)(a * vnorm/(c * vnorm + (a - c))*255.);
            //unsigned char v = (unsigned char)((val - data_min)*delta*255.);
            //unsigned char v = (unsigned char)(0.5*(tanh(100 * vnorm - 2.5) + 1)*255.);
            float_t v = (a * vnorm/(c * vnorm + (a - c)));

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

void compute_covs(float_t* image, int* segmentation_map, int* mask, 
                  int nrows, int ncols, int nclusters, 
                  float_t* centers_of_mass, 
                  float_t* cov_matrices, 
                  float_t* flux,
                  int* areas,
                  float_t* rmax,
                  int* parent_id,
                  int* x_limits,
                  int* y_limits)
{
    #define LOWER_BOUND(x) (2*x) 
    #define UPPER_BOUND(x) (2*x + 1) 
    // initialization
    #pragma omp parallel for
    for(int i = 0; i < nclusters; ++i)
    {
        centers_of_mass[2*i]     = 0;
        centers_of_mass[2*i + 1] = 0;

        cov_matrices[4*i]     = 0.;
        cov_matrices[4*i + 1] = 0.;
        cov_matrices[4*i + 2] = 0.;
        cov_matrices[4*i + 3] = 0.;

        parent_id[i]    = -1;

        x_limits[LOWER_BOUND(i)] = ncols; 
        x_limits[UPPER_BOUND(i)] = 0;

        y_limits[LOWER_BOUND(i)] = nrows; 
        y_limits[UPPER_BOUND(i)] = 0;

        flux[i]    = 0.;
        areas[i]   = 0.;
        rmax[i]    = 0.;
    }

    #pragma omp parallel
    {
        float_t* pvt_centers_of_mass = (float_t*)calloc(2 * nclusters, sizeof(float_t));
        float_t* pvt_cov_matrices    = (float_t*)calloc(4 * nclusters, sizeof(float_t));
        float_t* pvt_flux            = (float_t*)calloc(nclusters, sizeof(float_t));
        float_t* pvt_r_max           = (float_t*)calloc(nclusters, sizeof(float_t));
        
        int* pvt_areas                  = (int*)calloc(nclusters, sizeof(int));
        int* pvt_x_limits               = (int*)calloc(2 * nclusters, sizeof(int));
        int* pvt_y_limits               = (int*)calloc(2 * nclusters, sizeof(int));
        int* pvt_parent_id              = (int*)malloc(nclusters * sizeof(int));

        for(int i = 0; i < nclusters; ++i)
        {
            pvt_x_limits[LOWER_BOUND(i)] = ncols; 
            pvt_x_limits[UPPER_BOUND(i)] = 0;
            pvt_y_limits[LOWER_BOUND(i)] = nrows; 
            pvt_y_limits[UPPER_BOUND(i)] = 0;
            pvt_parent_id[i] = -1;
        }

        #pragma omp for
        for(int yy = 0; yy < nrows; ++yy)
            for(int xx = 0; xx < ncols; ++xx)
            {
                int lab = segmentation_map[yy*ncols + xx]; 
                if (lab != -1)
                {   
                    float_t pix_flux = image[yy*ncols + xx];
                    if (pix_flux < 0) pix_flux = fabs(pix_flux);
                    {
                        pvt_centers_of_mass[2*lab]     += (float_t)xx * pix_flux;
                        pvt_centers_of_mass[2*lab + 1] += (float_t)yy * pix_flux;
                        pvt_flux[lab]                  += pix_flux;
                    }
                    
                    pvt_areas[lab] += 1;

                    // if(parent_id[lab] == -1) parent_id[lab] = mask[yy*ncols + xx];
                    
                    int det_id = mask[yy*ncols + xx];
                    if (det_id != -1 && pvt_parent_id[lab] == -1) {
                        pvt_parent_id[lab] = det_id;
                    }

                    // pvt_x_limits[LOWER_BOUND(lab)] = MIN(yy, pvt_x_limits[LOWER_BOUND(lab)]);
                    // pvt_x_limits[UPPER_BOUND(lab)] = MAX(yy, pvt_x_limits[UPPER_BOUND(lab)]);
                    pvt_x_limits[LOWER_BOUND(lab)] = MIN(xx, pvt_x_limits[LOWER_BOUND(lab)]);
                    pvt_x_limits[UPPER_BOUND(lab)] = MAX(xx, pvt_x_limits[UPPER_BOUND(lab)]);

                    // pvt_y_limits[LOWER_BOUND(lab)] = MIN(xx, pvt_y_limits[LOWER_BOUND(lab)]);
                    // pvt_y_limits[UPPER_BOUND(lab)] = MAX(xx, pvt_y_limits[UPPER_BOUND(lab)]);
                    pvt_y_limits[LOWER_BOUND(lab)] = MIN(yy, pvt_y_limits[LOWER_BOUND(lab)]);
                    pvt_y_limits[UPPER_BOUND(lab)] = MAX(yy, pvt_y_limits[UPPER_BOUND(lab)]);
                }
            }

        // reduction
        #pragma omp critical (merging_coms)
        {
            for(int i = 0; i < nclusters; ++i)
            {
                centers_of_mass[2*i]     += pvt_centers_of_mass[2*i];
                centers_of_mass[2*i + 1] += pvt_centers_of_mass[2*i + 1];

                areas[i]   += pvt_areas[i];
                flux[i]    += pvt_flux[i];

                x_limits[LOWER_BOUND(i)] = MIN(x_limits[LOWER_BOUND(i)], pvt_x_limits[LOWER_BOUND(i)]);
                x_limits[UPPER_BOUND(i)] = MAX(x_limits[UPPER_BOUND(i)], pvt_x_limits[UPPER_BOUND(i)]);

                y_limits[LOWER_BOUND(i)] = MIN(y_limits[LOWER_BOUND(i)], pvt_y_limits[LOWER_BOUND(i)]);
                y_limits[UPPER_BOUND(i)] = MAX(y_limits[UPPER_BOUND(i)], pvt_y_limits[UPPER_BOUND(i)]);

                if (parent_id[i] == -1 && pvt_parent_id[i] != -1) {
                    parent_id[i] = pvt_parent_id[i];
                }
            }
        }

        #pragma omp barrier

        #pragma omp single
        {
            int max_det_id = -1;
            for (int lab = 0; lab < nclusters; ++lab) {
                if (parent_id[lab] > max_det_id) max_det_id = parent_id[lab];
            }

            int* det_child_count = (int*)calloc(max_det_id + 1, sizeof(int));

            for (int lab = 0; lab < nclusters; ++lab) {
                int det_id = parent_id[lab];
                if (det_id != -1) {
                    det_child_count[det_id] += 1;
                }
            }

            for (int lab = 0; lab < nclusters; ++lab) {
                int det_id = parent_id[lab];
                if (det_id == -1) continue;

                if (det_child_count[det_id] <= 1)
                    parent_id[lab] = -1;      // not deblended
                else
                    parent_id[lab] = det_id;  // deblended child
            }

            free(det_child_count);
        }

        // #pragma omp barrier

        #pragma omp for
        for(int i = 0; i < nclusters; ++i)
        {
            // centers_of_mass[2*i]     = centers_of_mass[2*i]    /areas[i];
            // centers_of_mass[2*i + 1] = centers_of_mass[2*i + 1]/areas[i];
            centers_of_mass[2*i]     = centers_of_mass[2*i]    /flux[i];
            centers_of_mass[2*i + 1] = centers_of_mass[2*i + 1]/flux[i];
        }


        #pragma omp for
        for(int yy = 0; yy < nrows; ++yy)
            for(int xx = 0; xx < ncols; ++xx)
            {
                int lab = segmentation_map[yy*ncols + xx]; 
                if(lab != -1)
                {
                    float x_n = (float_t)xx - centers_of_mass[2*lab];
                    float y_n = (float_t)yy - centers_of_mass[2*lab + 1];

                    float_t pix_flux = image[yy*ncols + xx];
                    float_t r        = sqrt(x_n * x_n + y_n * y_n);
                    if (r > pvt_r_max[lab]) pvt_r_max[lab] = r;
                    
                    // pvt_cov_matrices[4*lab    ] += x_n * x_n;
                    // pvt_cov_matrices[4*lab + 1] += x_n * y_n;
                    // pvt_cov_matrices[4*lab + 2] += x_n * y_n;
                    // pvt_cov_matrices[4*lab + 3] += y_n * y_n;
                    if (pix_flux < 0) pix_flux = fabs(pix_flux);
                    {
                        pvt_cov_matrices[4*lab    ] += x_n * x_n * pix_flux;
                        pvt_cov_matrices[4*lab + 1] += x_n * y_n * pix_flux;
                        pvt_cov_matrices[4*lab + 2] += x_n * y_n * pix_flux;
                        pvt_cov_matrices[4*lab + 3] += y_n * y_n * pix_flux;
                    }
                }
            }

        #pragma omp critical (merging_cov)
        {
            for(int lab = 0; lab < nclusters; ++lab)
            {
                    cov_matrices[4*lab]     += pvt_cov_matrices[4*lab];
                    cov_matrices[4*lab + 1] += pvt_cov_matrices[4*lab + 1];
                    cov_matrices[4*lab + 2] += pvt_cov_matrices[4*lab + 2];
                    cov_matrices[4*lab + 3] += pvt_cov_matrices[4*lab + 3];

                    x_limits[2*lab    ] = MIN(pvt_x_limits[2*lab], x_limits[2*lab]);
                    x_limits[2*lab + 1] = MAX(pvt_x_limits[2*lab + 1], x_limits[2*lab + 1]);

                    y_limits[2*lab    ] = MIN(pvt_y_limits[2*lab], y_limits[2*lab]);
                    y_limits[2*lab + 1] = MAX(pvt_y_limits[2*lab + 1], y_limits[2*lab + 1]);
                    
                    rmax[lab]           = MAX(rmax[lab], pvt_r_max[lab]);
            }
        }

        #pragma omp barrier

        #pragma omp for
        for(int lab = 0; lab < nclusters; ++lab)
        {
                // cov_matrices[4*lab]     = cov_matrices[4*lab]    /(areas[lab]-1);
                // cov_matrices[4*lab + 1] = cov_matrices[4*lab + 1]/(areas[lab]-1);
                // cov_matrices[4*lab + 2] = cov_matrices[4*lab + 2]/(areas[lab]-1);
                // cov_matrices[4*lab + 3] = cov_matrices[4*lab + 3]/(areas[lab]-1);

                cov_matrices[4*lab]     = cov_matrices[4*lab]    /flux[lab];
                cov_matrices[4*lab + 1] = cov_matrices[4*lab + 1]/flux[lab];
                cov_matrices[4*lab + 2] = cov_matrices[4*lab + 2]/flux[lab];
                cov_matrices[4*lab + 3] = cov_matrices[4*lab + 3]/flux[lab];
        }
 
        free(pvt_centers_of_mass);
        free(pvt_cov_matrices);
        free(pvt_flux);
        free(pvt_x_limits);
        free(pvt_y_limits);
        free(pvt_areas);
        free(pvt_parent_id);
        free(pvt_r_max);
    }

    #undef LOWER_BOUND
    #undef UPPER_BOUND
}


// FIX: to check
void compute_eigensystem_2x2(const double *A, double *lambda, double *V) {
    // A is symmetric: [[a, b], [b, d]]
    double a = A[0];
    double b = A[1];
    double d = A[3]; // Note: A[2] is also 'b'

    // Use a small tolerance for comparison with zero
    const double EPS = 1e-9;

    // 1. Compute Eigenvalues (lambda)
    // Formula: lambda = ( (a+d) +/- sqrt((a-d)^2 + 4*b^2) ) / 2
    double sum = a + d;
    double diff_sq = (a - d) * (a - d);
    double four_b_sq = 4.0 * b * b;

    // Discriminant Delta = (a-d)^2 + 4*b^2. Always non-negative.
    double delta = diff_sq + four_b_sq;
    double sqrt_delta = sqrt(delta);

    // Store eigenvalues (lambda[0] = lambda1, lambda[1] = lambda2)
    lambda[0] = (sum + sqrt_delta) / 2.0; // Larger eigenvalue
    lambda[1] = (sum - sqrt_delta) / 2.0; // Smaller eigenvalue

    // 2. Compute Eigenvectors (V)
    double lambda1 = lambda[0];
    double lambda2 = lambda[1];

    // --- Eigenvector 1 (for lambda1) ---
    // General case: v1 = [ b, lambda1 - a ]^T
    if (fabs(b) < EPS) {
        // Case: Diagonal matrix (b=0)
        // Eigenvectors are [1, 0] and [0, 1].
        V[0] = 1.0;
        V[1] = 0.0;
    } else {
        // Standard case: [ b, lambda1 - a ]^T
        V[0] = b;
        V[1] = lambda1 - a;
    }
    
    // --- Eigenvector 2 (for lambda2) ---
    // V[2], V[3] store v2
    // General case: v2 = [ b, lambda2 - a ]^T
    if (fabs(b) < EPS) {
        // Case: Diagonal matrix (b=0)
        // Eigenvectors are [1, 0] and [0, 1].
        // If lambda1 == lambda2, use the orthogonal basis [0, 1] for v2
        if (fabs(lambda1 - lambda2) < EPS) {
            V[2] = 0.0;
            V[3] = 1.0;
        } else {
             // Distinct eigenvalues: The eigenvectors are [1, 0] and [0, 1].
            V[2] = 0.0;
            V[3] = 1.0;
        }
    } else {
        // Standard case: [ b, lambda2 - a ]^T
        V[2] = b;
        V[3] = lambda2 - a;
    }

    // normalize eigenvector 1
    double n1 = sqrt(V[0]*V[0] + V[1]*V[1]);
    if (n1 > EPS) {
        V[0] /= n1;
        V[1] /= n1;
    }

    // normalize eigenvector 2
    double n2 = sqrt(V[2]*V[2] + V[3]*V[3]);
    if (n2 > EPS) {
        V[2] /= n2;
        V[3] /= n2;
    }
}

void compute_eigensystems(float_t* cov_matrices, float_t* lambdas, float_t* vs, int nclusters)
{
    #pragma omp parallel for
    for(int lab = 0; lab < nclusters; ++lab)
    {
        compute_eigensystem_2x2(cov_matrices + lab*4, lambdas + 2*lab, vs + 4*lab);
    }
}

void export_cluster_assignment(Datapoint_info* points, int* labels, idx_t n)
{
	for(idx_t i = 0; i < n; ++i) labels[i] = points[i].cluster_idx;
}

