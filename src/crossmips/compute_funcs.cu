//------------------------------------------------------------------------------------------------
// Copyright (c) 2012  Alessandro Bria and Giulio Iannello (University Campus Bio-Medico of Rome).
// Copyright (c) 2017  Massimo Bernaschi (IAC-CNR Roma).  
// All rights reserved.
//------------------------------------------------------------------------------------------------

/*******************************************************************************************************************************************************************************************
*    LICENSE NOTICE
********************************************************************************************************************************************************************************************
*    By downloading/using/running/editing/changing any portion of codes in this package you agree to this license. If you do not agree to this license, do not download/use/run/edit/change
*    this code.
********************************************************************************************************************************************************************************************
*    1. This material is free for non-profit research, but needs a special license for any commercial purpose. Please contact Alessandro Bria at a.bria@unicas.it or Giulio Iannello at
*       g.iannello@unicampus.it for further details.
*    2. You agree to appropriately cite this work in your related studies and publications.
*
*       Bria, A., Iannello, G., "TeraStitcher - A Tool for Fast 3D Automatic Stitching of Teravoxel-sized Microscopy Images", (2012) BMC Bioinformatics, 13 (1), art. no. 316.
*
*    3. This material is provided by  the copyright holders (Alessandro Bria  and  Giulio Iannello),  University Campus Bio-Medico and contributors "as is" and any express or implied war-
*       ranties, including, but  not limited to,  any implied warranties  of merchantability,  non-infringement, or fitness for a particular purpose are  disclaimed. In no event shall the
*       copyright owners, University Campus Bio-Medico, or contributors be liable for any direct, indirect, incidental, special, exemplary, or  consequential  damages  (including, but not
*       limited to, procurement of substitute goods or services; loss of use, data, or profits;reasonable royalties; or business interruption) however caused  and on any theory of liabil-
*       ity, whether in contract, strict liability, or tort  (including negligence or otherwise) arising in any way out of the use of this software,  even if advised of the possibility of
*       such damage.
*    4. Neither the name of University  Campus Bio-Medico of Rome, nor Alessandro Bria and Giulio Iannello, may be used to endorse or  promote products  derived from this software without
*       specific prior written permission.
********************************************************************************************************************************************************************************************/

/******************
*    CHANGELOG    *
*******************
* 2019-04-25. Giulio.     @ADDED using subsampling to compute the first NCC map
* 2019-04-08. Giulio.     @FIXED introduced optimization in loops in function compute_NCC
* 2018-08-29. Giulio.     @CHECKED checked the code implementing the optimization
* 2018.08-27. Giulio.     @FIXED it could happen that auxiliary matrices ps1 and ps2 were not allocated because images are too small and 'seq_cpu_compute_partial_sums' raises an error 
* 2018-08-26. Giulio.     @ADDED code and comments concerning the optimization in computing multiple NCCs on the same images with different overlap
* 2018-06-05. Giulio.     @FIXED bug in allocation of d_missrv which was not reallocated when dimu and/or dimv change
* 2018-05-05. Giulio.     @ADDED conditional code for excluding timing instructions under Windows
* 2018-04-15. Massimo.    @FIXED bugs in CUDA code
* 2018-02-17. Giulio.     @FIXED a bug in 'compute_Neighborhood' when NCCs to be reused are moved internally to the NCC map
* 2017-08-21. Massimo.    @ADDED   CUDA implementation of NCC
* 2017-04-01. Giulio.     @CHANGED the algorithm that computes the peak width
* 2015-04-06. Giulio.     @CHANGED corrected compute_NCC_alignment to deal with the case widthX = 1 which likely to be an anomaly
* 2015-03-20. Giulio.     @CHANGED newu and newv have been moved as parameters in compute_Neighborhood
* 2014-10-31. Giulio.     @CHANGED computations in compute_NCC are performed in double precision (and not in single precision) to avoit roundoff errors
*/

/*
 * compute_funcs.cpp
 *
 *  Created on: September 2010
 *      Author: iannello
 *
 *  Last revision: May, 31 2011
 */

# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include <assert.h>

#ifdef USE_OPENMP
        #include <omp.h>
#endif

#ifndef TIMING_H
#define TIMING_H

#ifdef _WIN32

#define TIMER_DEF     

#define TIMER_START   

#define TIMER_STOP    

#define TIMER_ELAPSED 0    // this must be an expression

#else

#include <sys/time.h>

#define TIMER_DEF     struct timeval temp_1, temp_2

#define TIMER_START   gettimeofday(&temp_1, (struct timezone*)0)

#define TIMER_STOP    gettimeofday(&temp_2, (struct timezone*)0)

#define TIMER_ELAPSED ((temp_2.tv_sec-temp_1.tv_sec)*1.e6+(temp_2.tv_usec-temp_1 .tv_usec))

#endif

#endif // TIMING_H

//extern double simpletimer;
double simpletimer;

/********* CONTROL SYMBOLS AND VARIABLES FOR CUDA CODE ****************************************************/

/* the USECUDA symbol controls conditional compilation of CUDA code 
 * it is defined by CMake when the WITH_CUDA option is checked
 * when WITH_CUDA option is unchecked, this file is included by compute_funcs.cpp to force C++ compilation
 */

// NCC is computed on NVIDIA card if the environment variable USECUDA_X_NCC exists 
#define CUDAENVVAR "USECUDA_X_NCC"

/* this variable control actual execution on NVIDIA card It can assume the following values: 
 * -1: non initialized
 * 0: do not use NVIDIA card
 * 1: use NVIDIA card
 */
static int usecuda=-1; 
/*********************************************************************************************************/


# include "compute_funcs.h"

# define LOG2(V)   (log((double)V)/log(2.0))


/**************************** INTERNAL AUXILIARY FUNCTION *******************************/

/* search x in the sorted list (a,n) and returns if x has been found; if x is not found,
 * returns in pos the index at which x should be inserted (i.e. a[pos] is the first element
 * of the list that is greater than x)
 * specialized for use by enhance: returns the index of the linear transformation to be
 * applied to value x
 */
static
void binary_search ( iom::real_t *a, int n, iom::real_t x, bool &found, int &pos ) {
        int f = 0;
        int l = n-1;
        int m;

        found = false;
        while ( f <= l ) {
                m = (f+l)/2;
                if ( x == a[m] ) {
                        found = true;
                        pos = m+1; // if found, the index to be returned is the next one
                        return;
                }
                else if ( x < a[m] )
                        l = m - 1;
                else   // V > a[m]
                        f = m + 1;
        }
        pos = l+1;
        return;
}

/* given general parameters, the NCC map, its second dimension dimj, the index of NCC maximum
 * and the extension of the map along the vertical (wRangeThr1) and horizontal (wRangeThr2) directions
 * returns a measure of the vertical (width1) and horizontal (width2) half widths of the NCC peak
 * (i.e. the largest distance of the pixel equal to a given fraction of the maximum and the
 * maximum itself)
 */
static
void compute_NCC_width ( NCC_parms_t *NCC_params, iom::real_t *NCC, int dimj, int ind, int wRangeThr1, int wRangeThr2, bool failed, int &width1, int &width2 ) {
        bool found;
        iom::real_t prec_val;
        int dist;

        iom::real_t thr = NCC_params->widthThr * NCC[ind];

        if(failed)
        {
                width1 = width2 = NCC_params->INF_W;
        }
        else
        {
                // evaluates first maximum width parallel to second dimension (horizontal)
                if ( wRangeThr2 < NCC_params->minDim_NCCmap ) {
                        // the map is too narrow in horizontal direction
                        width2 = NCC_params->INF_W;
                }
                else {
                        // try to find a very clear peak
                        found = false;
                        width2 = 1;
                        while ( width2<=wRangeThr2 && !found )
                                if ( NCC[ind-width2] <= thr )
                                        found = true;
                                else
                                        width2++;
                        found = false;
                        while ( width2<=wRangeThr2 && !found )
                                if ( NCC[ind+width2] <= thr )
                                        found = true;
                                else
                                        width2++;
                        if ( !found ) { // try to find if there is a peak anyway
                                // skip NCC_params->minPoints points
                                prec_val = NCC[ind-NCC_params->minPoints];
                                dist = NCC_params->minPoints + 1;
                                while ( dist<=wRangeThr2 && !found )
                                        if ( NCC[ind-dist] >= prec_val )
                                                found = true;
                                        else {
                                                prec_val = NCC[ind-dist];
                                                dist++;
                                        }
                                if ( dist < (2*NCC_params->minPoints) ) // not enough points
                                        width2 = NCC_params->INF_W;
                                else
                                        // project the profile to compute the equivalent width at thr
                                        width2 = (int)floor((dist-1) * (NCC[ind] - thr) / (NCC[ind] - prec_val));

                                found = false;
                                prec_val = NCC[ind+NCC_params->minPoints];
                                dist = NCC_params->minPoints + 1;
                                while ( dist<=wRangeThr2 && !found ) {
                                        if ( NCC[ind+dist] >= prec_val ) // NCC increases
                                                found = true;
                                        else {
                                                prec_val = NCC[ind+dist];
                                                dist++;
                                        }
                                }
                                if ( dist < (2*NCC_params->minPoints) ) // not enough points
                                        width2 = NCC_params->INF_W;
                                else
                                        // project the profile to compute the equivalent width at thr and compare with the largest width that is not 'infinite'
                                        width2 = MIN(MAX(width2,(int)floor((dist-1) * (NCC[ind] - thr) / (NCC[ind] - prec_val))),NCC_params->INF_W-1);
                        }
                }

                // evaluates maximum width parallel to first dimension (vertical)
                if ( wRangeThr1 < NCC_params->minDim_NCCmap ) {
                        // the map is too narrow in vertical direction
                        width1 = NCC_params->INF_W;
                }
                else {
                        // try to find a very clear peak
                        found = false;
                        width1 = 1;
                        while ( width1<=wRangeThr1 && !found )
                                if ( NCC[ind-width1*dimj] <= thr )
                                        found = true;
                                else
                                        width1++;
                        found = false;
                        while ( width1<=wRangeThr1 && !found )
                                if ( NCC[ind+width1*dimj] <= thr )
                                        found = true;
                                else
                                        width1++;
                        if ( !found ) { // try to find if there is a peak anyway
                                prec_val = NCC[ind-NCC_params->minPoints*dimj];
                                dist = NCC_params->minPoints + 1;
                                while ( dist<=wRangeThr2 && !found )
                                        if ( NCC[ind-dist*dimj] >= prec_val )
                                                found = true;
                                        else {
                                                prec_val = NCC[ind-dist*dimj];
                                                dist++;
                                        }
                                if ( dist < (2*NCC_params->minPoints) ) // not enough points
                                        width1 = NCC_params->INF_W;
                                else
                                        width1 = (int)floor((dist-1) * (NCC[ind] - thr) / (NCC[ind] - prec_val));

                                found = false;
                                prec_val = NCC[ind+NCC_params->minPoints*dimj];
                                dist = NCC_params->minPoints + 1;
                                while ( dist<=wRangeThr2 && !found )
                                        if ( NCC[ind+dist*dimj] >= prec_val )
                                                found = true;
                                        else {
                                                prec_val = NCC[ind+dist*dimj];
                                                dist++;
                                        }
                                if ( dist < (2*NCC_params->minPoints) ) // not enough points
                                        width1 = NCC_params->INF_W;
                                else
                                        // project the profile to compute the equivalent width at thr and compare with the largest width that is not 'infinite'
                                        width1 = MIN(MAX(width1,(int)floor((dist-1) * (NCC[ind] - thr) / (NCC[ind] - prec_val))),NCC_params->INF_W-1);
                        }
                }
        }
}

/*************************** COMPUTE FINAL ALIGNMENT *****************************************/

/* given the general parameters, the index i of a dimension (0 for vertical, 1 for horizontal
 * and 2 for depth), and two alignments d1 and d2 with the corresponding NCC maxima (peak_val1
 * and peak_val2) and half widths (width1 and width2), returns in result the aligment for that
 * with a measure of its reliability and potential error
 *
 * 2015-04-06. Giulio. The algorithm has been vcorrected as follows:
 * when widthX is 1 the alignment is unreliable since it is very likely that it is due to
 * a spike or a too little NCC map (e.g. because the stack is very thin); for this reason
 * parametes width1 and width2 are first checked and changes if equal to 1
 */
static
void compute_NCC_alignment ( NCC_parms_t *NCC_params, NCC_descr_t *result, int i,
                                                        int d1, iom::real_t peak_val1, int width1, int d2, iom::real_t peak_val2, int width2 ) {

        // check width1 and widthw
        if ( width1 == 1 ) // alignment 1 is unreliable
                width1 = NCC_params->INF_W;
        if ( width2 == 1 ) // alignment 1 is unreliable
                width2 = NCC_params->INF_W;

        // check how many values contribute to final alignment
        if ( peak_val1 >= NCC_params->maxThr && width1 < NCC_params->INF_W ) // first value may be considered
                if ( peak_val2 >= NCC_params->maxThr && width2 < NCC_params->INF_W ) // second value may be considered too
                        if ( abs(d1 - d2) < MIN(width1,width2) ) { // both values must be considered
                                result->coord[i] = (int) floor((peak_val1 * d1 + peak_val2 * d2) / (peak_val1 + peak_val2) + 0.5); // weighted mean of alignments
                                result->NCC_maxs[i] = (peak_val1 * peak_val1 + peak_val2 * peak_val2) / (peak_val1 + peak_val2); // weighted mean of reliabilities
                                result->NCC_widths[i] = MAX(width1,width2); // maximum width
                        }
                        else { // only one value should be considered: take into account both peak value and peak width
                                if ( peak_val1/width1 > peak_val2/width2 ) { // first value should be considered
                                        result->coord[i] = d1;
                                        result->NCC_maxs[i] = peak_val1;
                                        result->NCC_widths[i] = width1;
                                }
                                else { // second value should be considered
                                        result->coord[i] = d2;
                                        result->NCC_maxs[i] = peak_val2;
                                        result->NCC_widths[i] = width2;
                                }
                        }
                else { // only first value should be considered
                        result->coord[i] = d1;
                        result->NCC_maxs[i] = peak_val1;
                        result->NCC_widths[i] = width1;
                }
        else
                if ( peak_val2 >= NCC_params->maxThr && width2 < NCC_params->INF_W ) { // only second value should be considered
                        result->coord[i] = d2;
                        result->NCC_maxs[i] = peak_val2;
                        result->NCC_widths[i] = width2;
                }
                else { // none value is reliable
                        result->coord[i] = NCC_params->INV_COORD;  // invalid coordinate
                        result->NCC_maxs[i] = NCC_params->UNR_NCC; // unreliable NCC
                        result->NCC_widths[i] = NCC_params->INF_W; // peak of infinite width
                }
}
/****************************************************************************************/



/*********** THREADS PAREMETERS AND CODE ************/

# if defined(_PAR_VERSION) && 0  // WARNING: the code enabled by symbol _PAR_VERSION is obsolete and not fully checked

# include <windows.h>

// parallel configuration ------------------------------
int n_procs = 0;
int par_degree = 0;

void init_configuration ( ) {
        fprintf(stdout,"--- PARALLEL VERSION RUNNING ---\n");
        n_procs = atoi(getenv("NUMBER_OF_PROCESSORS"));
        par_degree = n_procs;
}

// compute_3_MIPs --------------------------------------

typedef struct{
        // input parametres
        iom::real_t *A;
        iom::real_t *B;
        int dimi_v;
        int dimj_v;
        int dimk_v;
        int MIP_stridek;
        int stridei;
        int stridek;
        // input/ouput parameters
        // output parametres
        iom::real_t *MIP_xy1;
        iom::real_t *MIP_xz1;
        iom::real_t *MIP_yz1;
        iom::real_t *MIP_xy2;
        iom::real_t *MIP_xz2;
        iom::real_t *MIP_yz2;
} compute_3_MIPs_params_t;

DWORD WINAPI compute_3_MIPs_worker ( LPVOID lpParam ) {
        iom::real_t *A       = ((compute_3_MIPs_params_t *) lpParam)->A;
        iom::real_t *B       = ((compute_3_MIPs_params_t *) lpParam)->B;
        int dimi_v      = ((compute_3_MIPs_params_t *) lpParam)->dimi_v;
        int dimj_v      = ((compute_3_MIPs_params_t *) lpParam)->dimj_v;
        int dimk_v      = ((compute_3_MIPs_params_t *) lpParam)->dimk_v;
        int stridei     = ((compute_3_MIPs_params_t *) lpParam)->stridei;
        int stridek     = ((compute_3_MIPs_params_t *) lpParam)->stridek;
        int MIP_stridek = ((compute_3_MIPs_params_t *) lpParam)->MIP_stridek;
        iom::real_t *MIP_xy1 = ((compute_3_MIPs_params_t *) lpParam)->MIP_xy1;
        iom::real_t *MIP_xz1 = ((compute_3_MIPs_params_t *) lpParam)->MIP_xz1;
        iom::real_t *MIP_yz1 = ((compute_3_MIPs_params_t *) lpParam)->MIP_yz1;
        iom::real_t *MIP_xy2 = ((compute_3_MIPs_params_t *) lpParam)->MIP_xy2;
        iom::real_t *MIP_xz2 = ((compute_3_MIPs_params_t *) lpParam)->MIP_xz2;
        iom::real_t *MIP_yz2 = ((compute_3_MIPs_params_t *) lpParam)->MIP_yz2;

        iom::real_t *vol1, *vol2;
        int i, j, k;

        // calcola MIP su xy, xz, yz scandendo una sola volta i due volumi
        for ( k=0, vol1=A, vol2=B; k<dimk_v; k++, vol1+=stridek, vol2+=stridek )
                for ( i=0; i<dimi_v; i++, vol1+=stridei, vol2+=stridei )
                        for ( j=0; j<dimj_v; j++, vol1++, vol2++ ) {
                                MIP_xy1[i*dimj_v+j] = MAX(MIP_xy1[i*dimj_v+j],*vol1);
                                MIP_xz1[i*MIP_stridek+k] = MAX(MIP_xz1[i*MIP_stridek+k],*vol1); // MIP stride along k dimension is the original MIP k dimension
                                MIP_yz1[j*MIP_stridek+k] = MAX(MIP_yz1[j*MIP_stridek+k],*vol1); // MIP stride along k dimension is the original MIP k dimension
                                MIP_xy2[i*dimj_v+j] = MAX(MIP_xy2[i*dimj_v+j],*vol2);
                                MIP_xz2[i*MIP_stridek+k] = MAX(MIP_xz2[i*MIP_stridek+k],*vol2); // MIP stride along k dimension is the original MIP k dimension
                                MIP_yz2[j*MIP_stridek+k] = MAX(MIP_yz2[j*MIP_stridek+k],*vol2); // MIP stride along k dimension is the original MIP k dimension
                        }

        return 0;
}

// compute_NCC_map --------------------------------------

typedef struct{
        // input parametres
        iom::real_t *MIP_1;
        iom::real_t *MIP_2;
        int dimu;
        int dimv;
        int delayu;
        int delayv;
        int u_start;
        int v_start;
        int u_end;
        int v_end;
        // input/ouput parameters
        // output parametres
        iom::real_t *NCC_map;
} compute_NCC_map_params_t;

DWORD WINAPI compute_NCC_map_worker ( LPVOID lpParam ) {
        iom::real_t *MIP_1    = ((compute_NCC_map_params_t *) lpParam)->MIP_1;
        iom::real_t *MIP_2    = ((compute_NCC_map_params_t *) lpParam)->MIP_2;
        int dimu         = ((compute_NCC_map_params_t *) lpParam)->dimu;
        int dimv         = ((compute_NCC_map_params_t *) lpParam)->dimv;
        int delayu       = ((compute_NCC_map_params_t *) lpParam)->delayu;
        int delayv       = ((compute_NCC_map_params_t *) lpParam)->delayv;
        int u_start      = ((compute_NCC_map_params_t *) lpParam)->u_start;
        int v_start      = ((compute_NCC_map_params_t *) lpParam)->v_start;
        int u_end        = ((compute_NCC_map_params_t *) lpParam)->u_end;
        int v_end        = ((compute_NCC_map_params_t *) lpParam)->v_end;
        iom::real_t *NCC_map  = ((compute_NCC_map_params_t *) lpParam)->NCC_map;

        iom::real_t *im1, *im2;
        int u, v;

        // nel seguito u=0 rappresenta il massimo scostamento negativo del secondo MIP rispetto al primo
        // con riferimento alla prima coordinata; v=0 ha il medesimo significato con riferimento alla seconda
        // coordinata

        for ( u=u_start; u<=u_end; u++ )
                for ( v=v_start; v<=v_end; v++ ) {
                        im1 = MIP_1 + START_IND(u*dimv) + START_IND(v);
                        im2 = MIP_2 + START_IND(-u*dimv) + START_IND(-v);
                        NCC_map[(u+delayu)*(2*delayv+1)+(v+delayv)] = compute_NCC(im1,im2,dimu-abs(u),dimv-abs(v),abs(v));

                }

        return 0;
}

# endif // defined(_PAR_VERSION) && 0


/************ OPERATIONS IMPLEMENTATION *************/

void seq_cpu_compute_partial_sums(iom::real_t *image1, iom::real_t *image2, int height, int width, 
								  iom::real_t *ps1, iom::real_t *ps2){
	
	// 2018. Giulio. @ADDED if auxiliary matrices have not been allocated do nothing
	if ( !ps1 || !ps2 )
		return;
		
	// assumes that ps1 and ps2 are allocated and have size floor(height/TILE_SIDE) x floor(width/TILE_SIDE)
    int nh = height - (height % TILE_SIDE);
    int nw = width  - (width % TILE_SIDE);

    // int ph = nh / TILE_SIDE; // Giulio: not used?
    int pw = nw / TILE_SIDE;

    for(int i = 0; i < nh; i += TILE_SIDE){
        for(int j = 0; j < nw; j+= TILE_SIDE){
            ps1[(i / TILE_SIDE)*pw + (j/TILE_SIDE)] = 0;
            ps2[(i / TILE_SIDE)*pw + (j/TILE_SIDE)] = 0;
            for(int l = 0; l < TILE_SIDE; l++){
                for (int k = 0; k < TILE_SIDE; k++){                    
                    ps1[(i / TILE_SIDE)*pw + (j/TILE_SIDE)] += image1[(i + l) * width + (j + k)];
                    ps2[(i / TILE_SIDE)*pw + (j/TILE_SIDE)] += image2[(i + l) * width + (j + k)];
                }
            }
        }
    }
}

void compute_3_MIPs ( iom::real_t *A, iom::real_t *B,
                                          iom::real_t *MIP_xy1, iom::real_t *MIP_xz1, iom::real_t *MIP_yz1,
                                          iom::real_t *MIP_xy2, iom::real_t *MIP_xz2, iom::real_t *MIP_yz2,
                                          int dimi_v, int dimj_v, int dimk_v, int stridei, int stridek ) {
# ifndef _PAR_VERSION 

        iom::real_t *vol1, *vol2;
        int i, j, k;

        // calcola MIP su xy, xz, yz scandendo una sola volta i due volumi
        for ( k=0, vol1=A, vol2=B; k<dimk_v; k++, vol1+=stridek, vol2+=stridek )
                for ( i=0; i<dimi_v; i++, vol1+=stridei, vol2+=stridei )
                        for ( j=0; j<dimj_v; j++, vol1++, vol2++ ) {
                                MIP_xy1[i*dimj_v+j] = MAX(MIP_xy1[i*dimj_v+j],*vol1);
                                MIP_xz1[i*dimk_v+k] = MAX(MIP_xz1[i*dimk_v+k],*vol1);
                                MIP_yz1[j*dimk_v+k] = MAX(MIP_yz1[j*dimk_v+k],*vol1);
                                MIP_xy2[i*dimj_v+j] = MAX(MIP_xy2[i*dimj_v+j],*vol2);
                                MIP_xz2[i*dimk_v+k] = MAX(MIP_xz2[i*dimk_v+k],*vol2);
                                MIP_yz2[j*dimk_v+k] = MAX(MIP_yz2[j*dimk_v+k],*vol2);
                        }

# else // ndef _PAR_VERSION     WARNING: the code enabled by symbol _PAR_VERSION is obsolete and not fully checked

        HANDLE *workerHandles = new HANDLE[par_degree];
        compute_3_MIPs_params_t *compute_3_MIPs_params = new compute_3_MIPs_params_t[par_degree];
        int t, i, j;

        /*
         *  work decomposition is performed by partitioning the volum along the k (i.e. z) direction
         *  each thread compute a portion of MIPS in xz and yz planes and a partial MIP in xy plane
         *  partial MIPS are then merged
         */

        // partition dimk_v
        int n1 = dimk_v / par_degree;
        int n2 = dimk_v % par_degree;
        for ( t=0; t<n2; t++ )
                compute_3_MIPs_params[t].dimk_v = n1 + 1;
        for ( ; t<par_degree; t++ )
                compute_3_MIPs_params[t].dimk_v = n1;

        // allocate and initialize memory for partial MIP computation
        iom::real_t **MIP_xy1_lst = new iom::real_t *[par_degree];
        iom::real_t **MIP_xy2_lst = new iom::real_t *[par_degree];
        // first partial MIPs are stored in MIP_xy1 and MIP_xy2
        MIP_xy1_lst[0] = MIP_xy1;
        MIP_xy2_lst[0] = MIP_xy2;
        for ( t=1; t<par_degree; t++ ) {
                MIP_xy1_lst[t] = new iom::real_t[dimi_v*dimj_v];
                MIP_xy2_lst[t] = new iom::real_t[dimi_v*dimj_v];
        }

        int slice_dim = dimi_v*(dimj_v+stridei)+stridek;
        /*
         * number of pixels of one slice of volumes A and B
         *
         * case NORTH_SOUTH:
         *    dimi_v  = dimi - ni
         *    dimj_v  = dimj
         *    stridei = 0
         *    stridek = ni*dimj    ====>
         *      slice_dim = (dimi - ni) * dimj + ni * dimj = dimi * dimj
         *
         * case WEST_EAST:
         *    dimi_v  = dimi
         *    dimj_v  = dimj - nj
         *    stridei = nj
         *    stridek = 0          ====>
         *      slice_dim = dimi * ( dimj - nj + nj) + 0 = dimi * dimj
         */
        int depth = 0;
        for ( t=0; t<par_degree; t++ ) {
                compute_3_MIPs_params[t].A = A + depth*slice_dim;
                compute_3_MIPs_params[t].B = B + depth*slice_dim;
                compute_3_MIPs_params[t].dimi_v = dimi_v;
                compute_3_MIPs_params[t].dimj_v = dimj_v;
                compute_3_MIPs_params[t].stridei = stridei;
                compute_3_MIPs_params[t].stridek = stridek;
                compute_3_MIPs_params[t].MIP_stridek = dimk_v;
                compute_3_MIPs_params[t].MIP_xy1 = MIP_xy1_lst[t];
                compute_3_MIPs_params[t].MIP_xz1 = MIP_xz1 + depth;
                compute_3_MIPs_params[t].MIP_yz1 = MIP_yz1 + depth;
                compute_3_MIPs_params[t].MIP_xy2 = MIP_xy2_lst[t];
                compute_3_MIPs_params[t].MIP_xz2 = MIP_xz2 + depth;
                compute_3_MIPs_params[t].MIP_yz2 = MIP_yz2 + depth;

                workerHandles[t] = CreateThread( NULL, 0, compute_3_MIPs_worker, (compute_3_MIPs_params+t), 0, NULL);

                depth += compute_3_MIPs_params[t].dimk_v;
        }

        WaitForMultipleObjects(par_degree,workerHandles,TRUE,INFINITE);

        for ( t=0; t<par_degree; t++ )
                CloseHandle(workerHandles[t]);

        // compute global MIP_xy
        for ( t=1; t<par_degree; t++ )
                for ( i=0; i<dimi_v; i++ )
                        for ( j=0; j<dimj_v; j++ ) {
                                MIP_xy1[i*dimj_v + j] = MAX(MIP_xy1[i*dimj_v + j],MIP_xy1_lst[t][i*dimj_v + j]);
                                MIP_xy2[i*dimj_v + j] = MAX(MIP_xy2[i*dimj_v + j],MIP_xy2_lst[t][i*dimj_v + j]);
                        }

        // deallocation
        for ( t=1; t<par_degree; t++ ) {
                delete MIP_xy1_lst[t];
                delete MIP_xy2_lst[t];
        }
        delete MIP_xy1_lst;
        delete MIP_xy2_lst;

        delete compute_3_MIPs_params;
        delete workerHandles;
        
# endif // ndef _PAR_VERSION
}

#if defined(USECUDA)
static iom::real_t *d_im1=NULL, *d_im2=NULL, *d_rv=NULL, *d_missrv=NULL;
static iom::real_t *dev_ps1 = NULL; // partial sums of image 1 (sums optimization)
static iom::real_t *dev_ps2 = NULL; // partial sums of image 2 (sums optimization)
static int sizeps = -1;             // size auxiliary matrices (sums optimization)
static int *d_missu=NULL, *d_missv=NULL;
static iom::real_t *s_im1=NULL, *s_im2=NULL;
static unsigned int sizeimg=0, sizeout=0;
static int sizemiss = -1;
static int sizerv = -1; // 2018-06-05. Giulio. @ADDED variable d_missrv_size to correctly initialize d_missrv

#include "warp_reduce.h"

#define MY_CUDA_CHECK( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#define MY_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }

#define NSTREAMS 2
#define NTHREADS 1024
#define MAXBLOCKS 1
#define REAL float
//#define USE_LDG
#ifdef USE_LDG
#define LDG(x) (__ldg(&(x)))
#else
#define LDG(x) (x)
#endif

#define CALCSUM(v) (v)=warpReduceSumD((v)); \
        if(lane==0) shared[wid]=(v); \
        __syncthreads(); \
        (v) = (threadIdx.x<blockDim.x/warpSize) ? shared[lane] : 0; \
        if(wid==0) (v)=warpReduceSumD((v)); \
        if(tid==0) shared[0]=(v); \
        __syncthreads(); \
        (v)=shared[0]; \
        __syncthreads();

#define CALCAVE(v) CALCSUM(v) \
        (v) /= (dimi*dimj);

__global__ void compute_partial_sums(iom::real_t *image1, iom::real_t *image2, int height, int width, iom::real_t *ps1, iom::real_t *ps2, int tileSide){

        int tid_x = blockDim.x * blockIdx.x + threadIdx.x;
        int tid_y = blockDim.y * blockIdx.y + threadIdx.y;
    
        int th = blockDim.x * gridDim.x;
        int tw = blockDim.y * gridDim.y; 
        
        __shared__ double partials1[NTHREADS];
        __shared__ double partials2[NTHREADS];
        
        int nh = height - (height % tileSide); // assumption: height, width >= tileSide
        int nw = width  - (width % tileSide);
        int pw = nw / tileSide;
        int threadSpan = tileSide / blockDim.x;
    
        
        
    
        for (int x = tid_x * threadSpan; x < nh; x += th*threadSpan){
            for(int y = tid_y * threadSpan; y < nw; y += tw*threadSpan){
                // Each loop of this while computes the sum of a tile
                double v1 = 0;
                double v2 = 0;
                for (int i = x;  i < x + threadSpan; i += 1){
                    for(int j = y; j < y + threadSpan; j += 1){
                        v1 += image1[i * width + j];
                        v2 += image2[i * width + j];
                    }
                }
                            
                // local coordinates with respect to block
                int tid = threadIdx.x * blockDim.y + threadIdx.y;
                partials1[tid] = v1;
                partials2[tid] = v2;
                __syncthreads();
    
                // reduce sum insde a warp
                for(int i = NTHREADS / 2; i > 0; i /= 2){
                    if(tid < i){
                        partials1[tid] += partials1[tid + i];
                        partials2[tid] += partials2[tid + i];           
                    }
                    __syncthreads();
                }
            
                if(tid == 0){
                    ps1[(x / tileSide)*pw + (y/tileSide)] = (float) partials1[0];
                    ps2[(x / tileSide)*pw + (y/tileSide)] = (float) partials2[0];
                }
            } 
        }
    }
    
            
            
__global__ void
__launch_bounds__(1024, 1)
gpu_NCC_map ( iom::real_t *MIP_1, iom::real_t *MIP_2, iom::real_t *rv, int dimu, int dimv, int delayu, int delayv,
        iom::real_t *ps1, iom::real_t *ps2, int tileSide ) {

        static __shared__ double shared[32];
     
        double f_mean, t_mean, f_prime, t_prime, numerator, factor1, factor2;

        int u=-delayu+blockIdx.x;
        int v=-delayv+blockIdx.y;
        int dimi=dimu-abs(u);
        int dimj=dimv-abs(v);
        int stride=abs(v);

        int a_u = START_IND(u);
        int a_v = START_IND(v);
        int b_u = START_IND(-u);
        int b_v = START_IND(-v);
        
        iom::real_t *pxl1 = MIP_1 + dimv * a_u + a_v;
        iom::real_t *pxl2 = MIP_2 + b_u*dimv + b_v;

        /*
        iom::real_t *pxl1 = MIP_1 + START_IND(u*dimv) + START_IND(v);
        iom::real_t *pxl2 = MIP_2 + START_IND(-u*dimv) + START_IND(-v);
        */
        const unsigned int tid = threadIdx.x;

        int lane=tid%warpSize;
        int wid=tid/warpSize;

        unsigned int ij;

        f_mean = t_mean = 0;
       
        // --- Cristian's addition: START
        
        // see comments in function 'NNC_compute' for further documentation
        // WARNING: in this implementation it is assumed that ps1 and ps2 have been allocated and initialized
        //          actually ps1 and ps2 are initialized with global CUDA variables dev_ps1 and dev_ps2
        
        int pw = (dimv / tileSide);
        int blockSide = 32;
        int i, j;
        int tid_x = (threadIdx.x / blockSide);
        int tid_y = (threadIdx.x % blockSide);
        
         // start indices of the overlapping region of MIP_1 and MIP_2 assigned to this thread
        int start_a_i = tid_x + a_u;
        int start_a_j = tid_y + a_v;
        int start_b_i = tid_x + b_u;
        int start_b_j = tid_y + b_v;

        // Indexes defining the tiled region inside the overlapping region of MIP_1
        int start_tiled_block_a_u = a_u - (a_u % tileSide);
        int start_tiled_block_a_v = a_v - (a_v % tileSide);
        start_tiled_block_a_u = start_tiled_block_a_u == a_u ? start_tiled_block_a_u : start_tiled_block_a_u + tileSide;
        start_tiled_block_a_v = start_tiled_block_a_v == a_v ? start_tiled_block_a_v : start_tiled_block_a_v + tileSide;
        int end_tiled_block_a_u = a_u + dimi - ((a_u + dimi) % tileSide);
        int end_tiled_block_a_v = a_v + dimj - ((a_v + dimj) % tileSide);

        // Indexes defining the tiled region inside the overlapping region of MIP_2
        int start_tiled_block_b_u = b_u - (b_u % tileSide);
        int start_tiled_block_b_v = b_v - (b_v % tileSide);
        start_tiled_block_b_u = start_tiled_block_b_u == b_u ? start_tiled_block_b_u : start_tiled_block_b_u + tileSide;
        start_tiled_block_b_v = start_tiled_block_b_v == b_v ? start_tiled_block_b_v : start_tiled_block_b_v + tileSide;
        int end_tiled_block_b_u = b_u + dimi - ((b_u + dimi) % tileSide);
        int end_tiled_block_b_v = b_v + dimj - ((b_v + dimj) % tileSide);

        // Since continuous threads doesnt work on contiguous pixels in this part, we advance with a larger offset
        int newBlockSide = blockSide * tileSide;
        // First, sum all the values in the tiled region. We will move the threads only on special pixels of 
        // the overlapping region, that is, the ones with indexes that are multiples of tileSide and are inside the
        // tiled region of the overlapping region. (In a local coordinates of a tile, they are at position (0, 0))
        
        // first for image MIP_1
        for(i = start_tiled_block_a_u + tileSide*tid_x; i < end_tiled_block_a_u; i+= newBlockSide){
                int row = (i/tileSide)*pw;
                for(j = start_tiled_block_a_v + tileSide*tid_y; j < end_tiled_block_a_v; j+= newBlockSide){
                f_mean += ps1[ row + (j/tileSide)];
                }
        }
        // then for image MIP_2
        for(i = start_tiled_block_b_u + tileSide*tid_x; i < end_tiled_block_b_u; i+= newBlockSide){
                int row = (i/tileSide)*pw;
                for(j = start_tiled_block_b_v + tileSide*tid_y; j < end_tiled_block_b_v; j+= newBlockSide){
                t_mean += ps2[ row + (j/tileSide)];
                }
        }

        // now sum the pixels at the border, that is, all the ones that are not inside the tiled region
        for(i = start_a_i; i < dimi + a_u; i+= blockSide){
                j = start_a_j;
                while(j < dimj + a_v){
                    if(j < start_tiled_block_a_v || j >= end_tiled_block_a_v || i < start_tiled_block_a_u || i >= end_tiled_block_a_u){
                        f_mean += MIP_1[i*dimv + j];
                        j+= blockSide;
                    }else{
                        int pixels_to_jump = end_tiled_block_a_v - j;
                        j += ((pixels_to_jump / blockSide) + ((pixels_to_jump % blockSide) ? 1 : 0))*blockSide;
                    }
                }
            }
        
            for(i = start_b_i; i < dimi + b_u; i+= blockSide){
                j = start_b_j;
                while(j < dimj + b_v) {
                    if(j < start_tiled_block_b_v || j >= end_tiled_block_b_v || i < start_tiled_block_b_u || i >= end_tiled_block_b_u){
                        t_mean += MIP_2[i*dimv + j];
                        j+= blockSide;
                    }else{
                        int pixels_to_jump = end_tiled_block_b_v - j;
                        j += ((pixels_to_jump / blockSide) + ((pixels_to_jump % blockSide) ? 1 : 0))*blockSide;
                    }
                }
            }
        // --- Cristian's addition: END
        /*
        for ( ij=tid; ij<(dimi*dimj); ij+=blockDim.x) {
          f_mean+=LDG(pxl1[(ij%dimj)+(stride+dimj)*(ij/dimj)]);
          t_mean+=LDG(pxl2[(ij%dimj)+(stride+dimj)*(ij/dimj)]);
        }*/
        CALCAVE(f_mean);
        CALCAVE(t_mean);
        
      //if(tid==0)  printf("FMEAN: %f, TMEAN: %f\n",f_mean, t_mean);
        // applies the optimization at the beginning of section 5 of Lewis article (t_prime has zero mean)
        numerator = 0;
        factor1 = 0;
        factor2 = 0;
        for ( ij=tid; ij<(dimi*dimj); ij+=blockDim.x) {
          f_prime = LDG(pxl1[(ij%dimj)+(stride+dimj)*(ij/dimj)]) - f_mean;
          t_prime = LDG(pxl2[(ij%dimj)+(stride+dimj)*(ij/dimj)]) - t_mean;
          numerator += LDG(pxl1[(ij%dimj)+(stride+dimj)*(ij/dimj)]) * t_prime;
          factor1 += f_prime * f_prime;
          factor2 += t_prime * t_prime;
        }
        CALCSUM(numerator);
        CALCSUM(factor1);
        CALCSUM(factor2);

        if(tid==0) {
//        printf("NUMERATOR: %f, FACTOR1: %f, FACTOR2: %f\n",numerator, factor1, factor2);
          rv[(u+delayu)*(2*delayv+1)+(v+delayv)]=((float) (numerator / sqrt(factor1*factor2))); // the result is converted to single precision
        }

}

__global__ void
__launch_bounds__(1024, 1)
gpu_NCC_miss ( iom::real_t *MIP_1, iom::real_t *MIP_2, iom::real_t *missrv,
              int dimu, int dimv, int du, int dv,
              int newu, int newv, int *missu, int *missv) {

        static __shared__ double shared[32];

        double f_mean, t_mean, f_prime, t_prime, numerator, factor1, factor2;

        int u = missu[blockIdx.x] - newu + du;
        int v = missv[blockIdx.x] - newv + dv;
        int dimi=dimu-abs(u);
        int dimj=dimv-abs(v);

        if(START_IND(dimi) == 0 || START_IND(dimj) == 0)
                return;

        int stride=abs(v);
        iom::real_t *pxl1 = MIP_1 + START_IND(u*dimv) + START_IND(v);
        iom::real_t *pxl2 = MIP_2 + START_IND(-u*dimv) + START_IND(-v);

        const unsigned int tid = threadIdx.x;

        int lane=tid%warpSize;
        int wid=tid/warpSize;

        unsigned int ij;
        f_mean = t_mean = 0;
        for ( ij=tid; ij<(dimi*dimj); ij+=blockDim.x) {
          f_mean+=pxl1[(ij%dimj)+(stride+dimj)*(ij/dimj)];
          t_mean+=pxl2[(ij%dimj)+(stride+dimj)*(ij/dimj)];
        }
        CALCAVE(f_mean);
        CALCAVE(t_mean);
//      if(tid==0)  printf("FMEAN: %f, TMEAN: %f\n",f_mean, t_mean);
        // applies the optimization at the beginning of section 5 of Lewis article (t_prime has zero mean)
        numerator = 0;
        factor1 = 0;
        factor2 = 0;
        for ( ij=tid; ij<(dimi*dimj); ij+=blockDim.x) {
          f_prime = pxl1[(ij%dimj)+(stride+dimj)*(ij/dimj)] - f_mean;
          t_prime = pxl2[(ij%dimj)+(stride+dimj)*(ij/dimj)] - t_mean;
          numerator += pxl1[(ij%dimj)+(stride+dimj)*(ij/dimj)] * t_prime;
          factor1 += f_prime * f_prime;
          factor2 += t_prime * t_prime;
        }
        CALCSUM(numerator);
        CALCSUM(factor1);
        CALCSUM(factor2);

        if(tid==0) {
//        printf("NUMERATOR: %f, FACTOR1: %f, FACTOR2: %f\n",numerator, factor1, factor2);
          missrv[missu[blockIdx.x]*(2*newv+1)+missv[blockIdx.x]]=((float) (numerator / sqrt(factor1*factor2))); // the result is converted to single precision
        }

}

#endif

void compute_NCC_map ( iom::real_t *NCC_map, iom::real_t *MIP_1, iom::real_t *MIP_2,
                                               int dimu, int dimv, int delayu, int delayv ) {
	//# ifndef _PAR_VERSION

# if !defined(_PAR_VERSION)

        if(usecuda<0) { // check if NVIDIA card should be used
           usecuda=getenv(CUDAENVVAR)?1:0;
        }

#if defined(USECUDA) // 2019-04-25. Giulio. Check 'usecuda' variable only if CUDA code is enabled

        if(usecuda>0) { // use the NVIDIA card
        
			unsigned int nthreads=NTHREADS;
			dim3 dimGrid((2*delayu)+1,(2*delayv+1));
			// computes means
			cudaEvent_t start, stop;
			cudaEventCreate(&start);
			cudaEventCreate(&stop);
			cudaEventRecord( start, 0 );
		
			// ------ cristian's addition START
			dim3 dimThread(32, 32);
			dim3 dimGridPS(8, 8);
			int tileSide = 32;
		
			int ph = dimu / tileSide;
			int pw = dimv / tileSide;
			int newSizePS =  pw * ph;  // total size of the matrix containing the partial sums (i.e. sums of each tile)
			if(sizeps < newSizePS){
					if(sizeps > 0){
							cudaFree(dev_ps1);
							cudaFree(dev_ps2);
					}
					sizeps = newSizePS;
					MY_CUDA_CHECK ( cudaMalloc((void **)&dev_ps1, sizeof(iom::real_t) * sizeps));
					MY_CUDA_CHECK ( cudaMalloc((void **)&dev_ps2, sizeof(iom::real_t) * sizeps));
			}               
			// ------ cristian's addition END

			if(sizeimg<(dimu*dimv)) {
			  if(sizeimg>0) {
				MY_CUDA_CHECK(cudaFree(d_im1));
				MY_CUDA_CHECK(cudaFree(d_im2));
			  }
			  sizeimg=dimu*dimv;
			  MY_CUDA_CHECK ( cudaMalloc (( void **) &d_im1, sizeimg*sizeof ( iom::real_t ) ));
			  MY_CUDA_CHECK ( cudaMalloc (( void **) &d_im2, sizeimg*sizeof ( iom::real_t ) ));
			}
			if(sizeout<((2*delayu+1)*(2*delayv+1))) {
			  if(sizeout>0) {
				MY_CUDA_CHECK(cudaFree(d_rv));
			  }
			  sizeout=(2*delayu+1)*(2*delayv+1);
			  MY_CUDA_CHECK ( cudaMalloc (( void **) &d_rv, sizeout*sizeof ( iom::real_t ) ));

			}
			if(s_im1!=MIP_1) {
					MY_CUDA_CHECK( cudaMemcpy(d_im1, MIP_1,
								   sizeof(iom::real_t)*(dimu*dimv), cudaMemcpyHostToDevice) );
					s_im1=MIP_1;
			}
			if(s_im2!=MIP_2) {
					MY_CUDA_CHECK( cudaMemcpy(d_im2, MIP_2,
								   sizeof(iom::real_t)*(dimu*dimv), cudaMemcpyHostToDevice) );
					s_im2=MIP_2;
			}
			MY_CUDA_CHECK( cudaMemcpy(d_rv, NCC_map,
					 sizeof(iom::real_t)*(2*delayu+1)*(2*delayv+1), cudaMemcpyHostToDevice) );

			// --- Cristian's addition START 
			compute_partial_sums<<<dimGridPS, dimThread>>>(d_im1, d_im2, dimu, dimv, dev_ps1, dev_ps2, tileSide);
			gpu_NCC_map<<<dimGrid,nthreads>>>(d_im1, d_im2, d_rv, dimu, dimv, delayu, delayv, dev_ps1, dev_ps2, tileSide);
			// ---- Cristian's addition END

			MY_CUDA_CHECK( cudaMemcpy(NCC_map, d_rv,
					 sizeof(iom::real_t)*(2*delayu+1)*(2*delayv+1), cudaMemcpyDeviceToHost) );
			cudaEventRecord( stop, 0 );
			cudaEventSynchronize( stop );
			float elapsedTime;
			MY_CUDA_CHECK( cudaEventElapsedTime( &elapsedTime, start, stop ) );
			simpletimer += elapsedTime*1000;        
        } 
        else 
        
#endif //defined(USECUDA)
        { // 2019-04-25. Giulio. If CUDA code is not enabled this code should be executed in any case
        	TIMER_DEF;   
        	TIMER_START;
        	int ph = dimu / TILE_SIDE;
        	int pw = dimv / TILE_SIDE;
			int psmatrix_size =  ph * pw; // total size of the matrix containing the partial sums (i.e. sums of each tile)
			iom::real_t *ps1 = (iom::real_t *) 0;
			iom::real_t *ps2 = (iom::real_t *) 0;
        	if(psmatrix_size > 0){
				ps1 = (iom::real_t *) malloc(sizeof(iom::real_t) * psmatrix_size);
				ps2 = (iom::real_t *) malloc(sizeof(iom::real_t) * psmatrix_size);
				if(!ps1 || !ps2){
					// TODO: handle malloc failures according to policy

				}
			}
	
			seq_cpu_compute_partial_sums(MIP_1, MIP_2, dimu, dimv, ps1, ps2);

#ifdef USE_OPENMP // WARNING: OpenMP not sufficiently tested and not maintaine
        	
        	// assume processors are a power of 2
			int nProcessors = omp_get_max_threads();
    	    int a = 2;
        	while(a*a < nProcessors) a*=2;
        	int b = nProcessors / a;
        	omp_set_num_threads(a*b);
        	#pragma omp parallel
        	{
        		int u, v;
        		int start_u = omp_get_thread_num() / b;
        		int start_v = omp_get_thread_num() % b;

        		for ( u=-delayu + start_u; u<=delayu; u += a)
                	for ( v=-delayv + start_v; v<=delayv; v += b) {
                        NCC_map[(u+delayu)*(2*delayv+1)+(v+delayv)] = compute_NCC(MIP_1,MIP_2,dimu,dimv, u, v, ps1, ps2);
                	}
        	}

#else // USE_OPENMP
        	
	        int u, v;
        	for ( u=-delayu; u<=delayu; u++ ) 
            	for ( v=-delayv; v<=delayv; v++ ) {
// 2019-04-25. Giulio. @CHANGED use subsampling for finding the maximum 
//         	for ( u=-delayu; u<=delayu; u+=2 ) 
//             	for ( v=-delayv; v<=delayv; v+=2 ) { 
                	NCC_map[(u+delayu)*(2*delayv+1)+(v+delayv)] = compute_NCC(MIP_1,MIP_2,dimu,dimv, u, v, ps1, ps2);
                }
                
#endif // USE_OPENMP
        	
        	if(ps1) 
        		delete ps1;
        	if(ps2) 
        		delete ps2;
        		
        	TIMER_STOP; 
        	simpletimer+=TIMER_ELAPSED;
        } // endif defined(USECUDA)
        
#if defined(DUMPNCC)

        for ( int u=-delayu; u<=delayu; u++ ) {
                for ( int v=-delayv; v<=delayv; v++ ) {
                   fprintf(stderr,"NCC_map[%d]=%f\n", (u+delayu)*(2*delayv+1)+(v+delayv),
                                        NCC_map[(u+delayu)*(2*delayv+1)+(v+delayv)]);
                }
        }
        
#endif // defined(DUMPNCC)

# else //!defined(_PAR_VERSION) WARNING: the code enabled by symbol _PAR_VERSION is obsolete and not fully checked

        int npu, npv, nu1, nu2, nv1, nv2;
        int pu, u, pv, v, t, du, dv;

        // assume the processors are a power of 2
        int pow_2 = (int) floor(LOG2(par_degree));
        int par_degree_2 = (int) pow(2.0,pow_2);

        HANDLE *workerHandles = new HANDLE[par_degree_2];
        compute_NCC_map_params_t *compute_NCC_map_params = new compute_NCC_map_params_t[par_degree_2];

        /*
         *  work decomposition is done by partitioning the NCC maps
         *  maps are halved repeatedly until floor(log2(#procs)) blocks are generated
         *  each thread is given the indices of one block and computes the NCCs corresponding to that block
         */

        npv = pow_2 / 2;
        npu = npv + (pow_2 % 2);
        npv = (int) pow(2.0,npv); // partitions along v
        npu = (int) pow(2.0,npu); // partitions along u

        nu1 = (2*delayu+1) / npu; // minimum number of rows per partition along u
        nu2 = (2*delayu+1) % npu; // number of partitions along u that have one more row
        nv1 = (2*delayv+1) / npv; // minimum number of rows per partition along v
        nv2 = (2*delayu+1) % npv; // number of partitions along v that have one more column

        for ( pu=0, u=-delayu, t=0; pu<npu; pu++, u+=du ) {
                du = (pu<nu2) ? (nu1+1) : nu1;
                for ( pv=0, v=-delayv ; pv<npv; pv++, v+=dv, t++ ) {
                        dv = (pv<nv2) ? (nv1+1) : nv1;
                        compute_NCC_map_params[t].u_start = u;
                        compute_NCC_map_params[t].v_start = v;
                        compute_NCC_map_params[t].u_end = u + (du-1);
                        compute_NCC_map_params[t].v_end = v + (dv-1);
                }
        }

        for ( t=0; t<par_degree_2; t++ ) {
                compute_NCC_map_params[t].MIP_1 = MIP_1;
                compute_NCC_map_params[t].MIP_2 = MIP_2;
                compute_NCC_map_params[t].dimu = dimu;
                compute_NCC_map_params[t].dimv = dimv;
                compute_NCC_map_params[t].delayu = delayu;
                compute_NCC_map_params[t].delayv = delayv;
                compute_NCC_map_params[t].NCC_map = NCC_map;

                workerHandles[t] = CreateThread( NULL, 0, compute_NCC_map_worker, (compute_NCC_map_params+t), 0, NULL);
        }

        WaitForMultipleObjects(par_degree_2,workerHandles,TRUE,INFINITE);

        for ( t=0; t<par_degree_2; t++ )
                CloseHandle(workerHandles[t]);

        // deallocation
        delete compute_NCC_map_params;
        delete workerHandles;

# endif // !defined(_PAR_VERSION)

}


iom::real_t compute_NCC ( iom::real_t *MIP_1, iom::real_t *MIP_2, int dimu, int dimv, int u, int v, iom::real_t *ps1, iom::real_t *ps2) {
	// parallelization of compute_NCC_map makes parallelization of this operation pointless
    // 2014-10-31 Giulio. @CHANGED       iom::real_t f_mean, t_mean, f_prime, t_prime, numerator, factor1, factor2;
    double f_mean, t_mean, f_prime, t_prime, numerator, factor1, factor2;
    iom::real_t *pxl1, *pxl2;
    int i, j;

	int dimi = dimu - abs(u); // number of rows of the overlapping region
	int dimj = dimv - abs(v); // number of columns of the overlapping region
	int stride = abs(v);

	// these variables are needed by the code that is executed when optimization is not applied
	iom::real_t* im1 = MIP_1 + START_IND(u*dimv) + START_IND(v);
	iom::real_t* im2 = MIP_2 + START_IND(-u*dimv) + START_IND(-v);

    // computes means
    f_mean = t_mean = 0;
        
	// if auxiliary matrices have not been allocated do nothing
	if ( ps1 && ps2 && TILE_SIDE <= dimu && TILE_SIDE <= dimv ){ 
		// 2018.08-27. Giulio. FIXED ps1 and ps2 must be allocated AND both dimu and dimv must be larger than TILE_SIDE
		// ps1 and ps2 have been allocated and at least the sums of some tiles have already been computed and stored in ps1 and ps2

		// START_IND(i) = (((i)>0) ? (i) : 0)
		int pw = (dimv / TILE_SIDE); // width of matrix of precomputed sums
		int a_u = START_IND(u);  // u index over MIP_1 of overlapping region
		int a_v = START_IND(v);  // v index over MIP_1 of overlapping region
		int b_u = START_IND(-u); // u index over MIP_2 of overlapping region
		int b_v = START_IND(-v); // v index over MIP_2 of overlapping region


		// Indexes defining the tiled region inside the overlapping region of MIP_1
		int start_tiled_block_a_u = a_u - (a_u % TILE_SIDE); // u index of first tile intersecting the overlapping region
		int start_tiled_block_a_v = a_v - (a_v % TILE_SIDE); // v index of first tile intersecting the overlapping region
		// if indices are outside the overlapping region they are advanced by TILE_SIDE
		start_tiled_block_a_u = start_tiled_block_a_u == a_u ? start_tiled_block_a_u : start_tiled_block_a_u + TILE_SIDE;
		start_tiled_block_a_v = start_tiled_block_a_v == a_v ? start_tiled_block_a_v : start_tiled_block_a_v + TILE_SIDE;
		int end_tiled_block_a_u = a_u + dimi - ((a_u + dimi) % TILE_SIDE); // u index following the last one of a tile internal to the overlapping region
		int end_tiled_block_a_v = a_v + dimj - ((a_v + dimj) % TILE_SIDE); // v index following the last one of a tile internal to the overlapping region

		// Indexes defining the tiled region inside the overlapping region of MIP_2
		// the same comments of MIP_1 hold
		int start_tiled_block_b_u = b_u - (b_u % TILE_SIDE);
		int start_tiled_block_b_v = b_v - (b_v % TILE_SIDE);
		start_tiled_block_b_u = start_tiled_block_b_u == b_u ? start_tiled_block_b_u : start_tiled_block_b_u + TILE_SIDE;
		start_tiled_block_b_v = start_tiled_block_b_v == b_v ? start_tiled_block_b_v : start_tiled_block_b_v + TILE_SIDE;
		int end_tiled_block_b_u = b_u + dimi - ((b_u + dimi) % TILE_SIDE);
		int end_tiled_block_b_v = b_v + dimj - ((b_v + dimj) % TILE_SIDE);

		// First, sum all the values in the tiled region. We will move the threads only on special pixels of 
		// the overlapping region, that is, the ones with indexes that are multiples of TILE_SIDE and are inside the
		// tiled region of the overlapping region. (In a local coordinates of a tile, they are at position (0, 0))
		
		// first for image MIP_1
		// i, j move to the first vertex of tiles internal to the overlapping region
		for(i = start_tiled_block_a_u; i  < end_tiled_block_a_u; i+= TILE_SIDE){
			int row = (i/TILE_SIDE)*pw;
			for(j = start_tiled_block_a_v; j < end_tiled_block_a_v; j+= TILE_SIDE){
				f_mean += ps1[ row + (j/TILE_SIDE)];
			}
		}
		// then for image MIP_2
		// i, j move to the first vertex of tiles internal to the overlapping region
		for(i = start_tiled_block_b_u; i < end_tiled_block_b_u; i+= TILE_SIDE){
			int row = (i/TILE_SIDE)*pw;
			for(j = start_tiled_block_b_v; j < end_tiled_block_b_v; j+= TILE_SIDE){
				t_mean += ps2[ row + (j/TILE_SIDE)];
			}
		}

		// first for MIP_!
		// now sum the pixels at the border, that is, all the ones that are not inside the tiled region
		// each pixel of the overlapping region is processed at most once (it is not processed if it is in the tiled region
		for(i = a_u; i < dimi + a_u; i+= 1){
			j = a_v;
			while(j < dimj + a_v){
				if(j < start_tiled_block_a_v || j >= end_tiled_block_a_v || i < start_tiled_block_a_u || i >= end_tiled_block_a_u){
					f_mean += MIP_1[i*dimv + j];
					j+= 1;
				}else{ // index j is just entered in the tiled region: skip the tiled region
					int pixels_to_jump = end_tiled_block_a_v - j;
					j += pixels_to_jump;
				}
			}
		}
		
		// then for MIP_2
		// the same comments of MIP_1 hold
		for(i = b_u; i < dimi + b_u; i+= 1){
			j = b_v;
			while(j < dimj + b_v) {
				if(j < start_tiled_block_b_v || j >= end_tiled_block_b_v || i < start_tiled_block_b_u || i >= end_tiled_block_b_u){
					t_mean += MIP_2[i*dimv + j];
					j+= 1;
				}else{
					int pixels_to_jump = end_tiled_block_b_v - j;
					j += pixels_to_jump;
				}
			}
		}
	}else{
		// code without optimization
		for ( i=0, pxl1=im1, pxl2=im2; i<dimi; i++, pxl1+=stride, pxl2+=stride ) {
			for ( j=0; j<dimj; j++, pxl1++, pxl2++ ) {
				f_mean += *pxl1;
				t_mean += *pxl2;
			}
		}
	}

    f_mean /= (dimi*dimj);
    t_mean /= (dimi*dimj);

	// applies the optimization at the beginning of section 5 of Lewis article (t_prime has zero mean)
	numerator = 0;
	factor1 = 0;
	factor2 = 0;
	// 2019-04-08. Giulio. using two moving pointers optimizes the double loop
	for ( i=0, pxl1=im1, pxl2=im2; i<dimi; i++, pxl1+=stride, pxl2+=stride )
		for ( j=0; j<dimj; j++, pxl1++, pxl2++ ) {
			f_prime = *pxl1 - f_mean;
			t_prime = *pxl2 - t_mean;
			numerator += *pxl1 * t_prime;
			factor1 += f_prime * f_prime;
			factor2 += t_prime * t_prime;
		}

    return ((float) (numerator / sqrt(factor1*factor2))); // the result is converted to single precision

}

int compute_MAX_ind ( iom::real_t *vect, int len ) {
// actual len values are too small to deserve parallelization
        int i;
        iom::real_t val_max = vect[0];
        int ind_max = 0;
        for ( i=0; i<len; i++ )
                if ( vect[i] > val_max ) {
                        val_max = vect[i];
                        ind_max = i;
                }
        return ind_max;
}


/* Returns an NCC map of of size (2*newu+1) x (2*newv+1) centered around the NCC maximum or returns filed = true
 * if this map cannot be found (it is: newu = newv = NCC_params->wRangeThr)
 *
 * Parameters:
 *   NCC_params     : INPUT        : parameters of the MIP-NCC algorithm
 *   NCC            : INPUT        : initial NCC map of size (2*delayu+1) x (2*delayv+1), centered around the initial alignment
 *   delayu, delayv : INPUT        : vertical and horizontal half extensions of the initial NCC map 
 *   newu, newv     : INPUT        : vertical and horizontal half extensions of the output NCC map (NCC_new)
 *   ind_max        : INPUT        : linear index of maximum in NCC
 *   MIP_1, MIP_2   : INPUT        : MIPs of size dimu x dimv
 *   dimu, dimv     : INPUT        : vertical and horizontal extensions of MIP_1 and MIP_2
 *   NCCnew         : INPUT/OUTPUT : output (initially empty) NCC map of size (2*newu+1) x (2*newv+1); it extends NCC
 *   du, dv         : OUTPUT       : vertical and horizontal relative positions of the maximum of NCCnew with respect to the initial alignment
 *                                   (they are undefined if 'failed' is true)
 *   failed         : INPUT/OUTPUT : initialized to false, it is changed to true if the maximum is not centered in NCCnew
 */
void compute_Neighborhood ( NCC_parms_t *NCC_params, iom::real_t *NCC, int delayu, int delayv, int newu, int newv, int ind_max,
                                                   iom::real_t *MIP_1, iom::real_t *MIP_2, int dimu, int dimv, iom::real_t *NCCnew, int &du, int &dv, bool &failed) throw (iom::exception){

	// --- CRISTIAN MOD START --- (to optimize NCC computation)
	int ph = dimu / TILE_SIDE;
	int pw = dimv / TILE_SIDE;
	int psmatrix_size = ph*pw;

	iom::real_t *ps1 = (iom::real_t *) 0;
	iom::real_t *ps2 = (iom::real_t *) 0;

	if(psmatrix_size > 0){
		ps1 = (iom::real_t *) malloc(sizeof(iom::real_t) * ph * pw);
		ps2 = (iom::real_t *) malloc(sizeof(iom::real_t) * ph * pw);
		if(!ps1 || !ps2){
			// TODO: handle malloc failures according to policy
		}
	}

	seq_cpu_compute_partial_sums(MIP_1, MIP_2, dimu, dimv, ps1, ps2);
	// --- CRISTIAN MOD END ---	
	
	// suffixes u and v denote the vertical and the horizontal dimensions, respectively
	// suffix i denotes linear indices

	int u, v, i, d; // for variables

	// 2015-03-20. Giulio. @CHANGED newu and newv are moved as parameters
	//int newu = NCC_params->wRangeThr; // vertical half dimension of NCCnew
	//int newv = NCC_params->wRangeThr; // horizontal half dimension of NCCnew

	int ind_ref; // index of the center of NCCnew

	int initu; // vertical index of first pixel of subregion of NCC to be used to initially fill NCCnew
	int initv; // horizontal index of first pixel of subregion of NCC to be used to initially fill NCCnew
	int initi; // linear index of first pixel of subregion of NCC to be used to initially fill NCCnew

	int srcStartu; // vertical index of first pixel of the subregion of NCCnew to be reused when current maximum is moved to the center of NCCnew
	int srcStartv; // horizontal index of first pixel of the subregion of NCCnew to be reused when current maximum is moved to the center of NCCnew
	int srcStarti; // linear index of first pixel of the subregion of NCCnew to be reused when current maximum is moved to the center of NCCnew

	int dstStartu; // vertical index of first pixel of the subregion of NCCnew where the subregion to be used has to be copied
	int dstStartv; // horizontal index of first pixel of the subregion of NCCnew where the subregion to be used has to be copied
	int dstStarti; // linear index of first pixel of the subregion of NCCnew where the subregion to be used has to be copied

	int deltau; // vertical displacement of current maximum from the center of NCCnew
	int deltav; // horizontal displacement of current maximum from the center of NCCnew

	int firstv; // first horizontal index for copying elements of NCCnew to be reused
	int lastv;  // last horizontal index for copying elements of NCCnew to be reused

	int n_miss; // number of NCC to be computed to fill NCCnew
	int *missu = new int[(2*newu+1)*(2*newv+1)]; // list of vertical indices of NCC to be computed to fill NCCnew
	int *missv = new int[(2*newu+1)*(2*newv+1)]; // list of vertical indices of NCC to be computed to fill NCCnew

	// INITIALIZATION

	// fill NCCnew copying useful NCCs that have been already computed from NCC to NCCnew
	initu = MIN(MAX(0,ind_max/(2*delayv+1) - newu),2*(delayu - newu)); // initu is at least 2*(delayu - newu) to guarantee that NCCnew can be completely initialized
	initv = MIN(MAX(0,ind_max%(2*delayv+1) - newv),2*(delayv - newv)); // initv is at least 2*(delayv - newv) to guarantee that NCCnew can be completely initialized
	initi = initu * (2*delayv+1) + initv;
	if(initi < 0)
		throw iom::exception("CrossMIPs: negative index detected (initi)"); // Alessandro - 23/03/2013 - throw exception if initi is negative
	for ( u=0, i=0, d=0; u<(2*newu+1); u++, d+=2*(delayv-newv) ) // when row changes 2*(delayv-newv) values have to be skipped
		for ( v=0; v<(2*newv+1); v++ , i++)
			NCCnew[i] = NCC[i + initi + d];
	// compute displacement of the center of NCCnew with respect to the initial alignment (center of NCC)
	du = initu - delayu + newu; // displacement of first row + half dimension of NCCnew
	dv = initv - delayv + newv; // displacement of first column + half dimension of NCCnew
	// update ind_max with respect to NCCnew
	//      contribution due to rows (integer division is not distributive)               contribution due to columns
	ind_max = (2*newv+1) * (ind_max/(2*delayv+1) - initi/(2*delayv+1))   +   (ind_max%(2*delayv+1)) - (initi%(2*delayv+1));
	// index of the center of the new NCC
	ind_ref = (2*newv+1) * newu + newv;

	// UPDATE NEIGHBORHOOD AND SEARCH MAXIMUM

	int c=0; // NCC_params->maxIter iterations are allowed
	while ( c < NCC_params->maxIter && ind_max != ind_ref ) {
		// update NCCnew
		srcStartu = MAX(0,ind_max/(2*newv+1) - newu);
		srcStartv = MAX(0,ind_max%(2*newv+1) - newv);
		srcStarti = srcStartu * (2*newv+1) + srcStartv;
		deltau = ind_max/(2*newv+1) - ind_ref/(2*newv+1);
		deltav = ind_max%(2*newv+1) - ind_ref%(2*newv+1);
		dstStartu = srcStartu - deltau;
		dstStartv = srcStartv - deltav;
		dstStarti = dstStartu * (2*newv+1) + dstStartv;
		if ( srcStartu > 0 ) {     // forward copy of rows is safe
			if ( srcStartv > 0 ) { // forward copy of columns is safe
				i = 0; // first index of first row to be moved
				for ( u=0; u<((2*newu+1)-abs(deltau)); u++, i+=abs(deltav) ) { // when row changes |deltav| values have to be skipped forward
					for ( v=0; v<((2*newv+1)-abs(deltav)); v++ , i++) {
						NCCnew[i + dstStarti] = NCCnew[i + srcStarti];
					}
				}
			}
			else { // srcStartv == 0: backward copy of columns is safe
				i = (2*newv+1) - abs(deltav) - 1;  // last index of first row to be moved  
				for ( u=0; u<((2*newu+1)-abs(deltau)); u++, i+=(2*(2*newv+1) - abs(deltav)) ) { // when row changes two rows - |deltav| have to be skipped forward
					for ( v=0; v<((2*newv+1)-abs(deltav)); v++ , i--) {
						NCCnew[i + dstStarti] = NCCnew[i + srcStarti];
					}
				}
			}
		}
		else { // srcStartu == 0:     backward copy of rows is safe
			if ( srcStartv > 0 ) { // forward copy of columns is safe
				i = ((2*newu+1)-abs(deltau)-1)*(2*newv+1);     // first index last row to be moved
				for ( u=0; u<((2*newu+1)-abs(deltau)); u++, i-=(2*(2*newv+1) - abs(deltav)) ) { // when row changes two rows - |deltav| values have to be skipped backward
					for ( v=0; v<((2*newv+1)-abs(deltav)); v++ , i++) {
						NCCnew[i + dstStarti] = NCCnew[i + srcStarti];
					}
				}
			}
			else { // srcStartv == 0: backward copy of columns is safe
				i = ((2*newu+1)-abs(deltau))*(2*newv+1) - abs(deltav) - 1;     // last index of last rows to be moved
				for ( u=0; u<((2*newu+1)-abs(deltau)); u++, i-=abs(deltav) ) { // when row changes |deltav| values have to be skipped backward
					for ( v=0; v<((2*newv+1)-abs(deltav)); v++ , i--) {
						NCCnew[i + dstStarti] = NCCnew[i + srcStarti];
					}
				}
			}
		}

		// update position of the new center (current maximum)
		du += deltau;
		dv += deltav;

		// generate list of missing NCCs
		n_miss = 0;
		for ( u=0; u<(2*newu+1); u++ ) {
			if ( u+deltau < 0 || u+deltau >= (2*newu+1) ) { // all this row has to be computed
					firstv = 0;
					lastv  = 2*newv+1;
			}
			else { // only a fraction of the row has to be computed
					if ( deltav < 0 ) { // the initial part of the row has to be computed
							firstv = 0;
							lastv  = -deltav;
					}
					else if ( deltav > 0 ) { // the final part of the row has to be computed
							firstv = (2*newv+1) - deltav;
							lastv  = 2*newv+1;
					}
					else { // deltav == 0: no NCCs has to be computed
							firstv = 0;
							lastv  = 0;
					}
			}
			for ( v=firstv; v<lastv; v++ ) {
					missu[n_miss] = u;
					missv[n_miss] = v;
					n_miss++;
			}
		}
		// CHECK, MUST BE: n_miss == ((2*newu+1)*(2*newv+1) - ((2*newu+1)-abs(deltau))*((2*newv+1)-abs(deltav)))
		if ( n_miss != ((2*newu+1)*(2*newv+1) - ((2*newu+1)-abs(deltau))*((2*newv+1)-abs(deltav))) )
			throw iom::exception("CrossMIPs: incomplete NCC map in compute_Neighborhood");

		// compute missing NCCs
#if defined(USECUDA)
		if(usecuda>0) { // use the CUDA card
			unsigned int nthreads=NTHREADS;
			unsigned int nblocks=n_miss;
			// computes means
			if((2*newu+1)*(2*newv+1) > sizerv){
				// 2018-06-05. Giulio. @FIXED reallocation and initialization of d_missrv moved from here because dimu and/or dimv could have been changed
				if(sizerv > 0){
					MY_CUDA_CHECK(cudaFree(d_missrv));
				}
				sizerv = (2*newu+1)*(2*newv+1);
				MY_CUDA_CHECK (
					cudaMalloc (( void **) &d_missrv, sizerv*sizeof ( iom::real_t ) ));
				MY_CUDA_CHECK ( cudaMemset (d_missrv, 0, sizeof(iom::real_t)*sizerv));
			}
			if(sizemiss<n_miss) {
			   if(sizemiss>0) {
				   MY_CUDA_CHECK(cudaFree(d_missu));
				   MY_CUDA_CHECK(cudaFree(d_missv));
			   }
			   sizemiss=n_miss;
		   
			   MY_CUDA_CHECK (cudaMalloc (( void **) &d_missu, sizemiss*sizeof(int)));
			   MY_CUDA_CHECK (cudaMalloc (( void **) &d_missv, sizemiss*sizeof(int)));
			}
		
			MY_CUDA_CHECK( cudaMemcpy(d_missu, missu,
						  sizeof(int)*n_miss,cudaMemcpyHostToDevice) );
			MY_CUDA_CHECK( cudaMemcpy(d_missv, missv,
						  sizeof(int)*n_miss,cudaMemcpyHostToDevice) );
					  
			if(sizeimg<(dimu*dimv)) {
				fprintf(stderr,"Reallocating memory before gpu_NCC_miss: %d, now is %d\n",
				sizeimg,(dimu*dimv));
				assert(sizeimg);
				MY_CUDA_CHECK(cudaFree(d_im1));
				MY_CUDA_CHECK(cudaFree(d_im2));
				sizeimg=(dimu)*(dimv);
				MY_CUDA_CHECK (
					cudaMalloc (( void **) &d_im1, sizeimg*sizeof ( iom::real_t ) ));
				MY_CUDA_CHECK (
					cudaMalloc (( void **) &d_im2, sizeimg*sizeof ( iom::real_t ) ));
			}
			if(s_im1!=MIP_1) {
				MY_CUDA_CHECK( cudaMemcpy(d_im1, MIP_1,
						   sizeof(iom::real_t)*(dimu*dimv), cudaMemcpyHostToDevice) );
				s_im1=MIP_1;
			}
			if(s_im2!=MIP_2) {
				MY_CUDA_CHECK( cudaMemcpy(d_im2, MIP_2,
						   sizeof(iom::real_t)*(dimu*dimv), cudaMemcpyHostToDevice) );
				s_im2=MIP_2;
			}
			//fprintf(stderr,"----> %p %lu %d %d\n",d_missrv,sizeof(iom::real_t),newu,newv); // IANNELLO
			MY_CUDA_CHECK( cudaMemcpy(d_missrv, NCCnew,
				   sizeof(iom::real_t)*(2*newu+1)*(2*newv+1), cudaMemcpyHostToDevice) );
			gpu_NCC_miss<<<nblocks, nthreads>>>
							(d_im1, d_im2, d_missrv, dimu, dimv,
							 du, dv, newu, newv, d_missu, d_missv);
			MY_CUDA_CHECK( cudaMemcpy(NCCnew, d_missrv,
				   sizeof(iom::real_t)*(2*newu+1)*(2*newv+1), cudaMemcpyDeviceToHost) );
		} 
		else 
#endif // defined(USECUDA)
		{ // 2019-04-25. Giulio. This code should always be executed if CUDA code is not enabled
			//fprintf(stderr,"----> %lu %d %d\n",sizeof(iom::real_t),newu,newv); // IANNELLO
			for ( i=0; i<n_miss; i++ ) {
				// indices over MIPs have to be shifted to take into account their relative position with respecto to the center of NCCnew
				// and the relative position of the center with respect to the initial initial alignment (center of NCC)
				u = missu[i] - newu + du;
				v = missv[i] - newv + dv;
				NCCnew[missu[i]*(2*newv+1)+missv[i]] = compute_NCC(MIP_1, MIP_2,dimu,dimv, u, v, ps1, ps2);
			}
		}
#if defined(DUMPNCC)
		for ( i=0; i<n_miss; i++ ) {
			fprintf(stderr,"NCCnew[%d]=%f\n",missu[i]*(2*newv+1)+missv[i],
										NCCnew[missu[i]*(2*newv+1)+missv[i]]);
		}
#endif
		// find maximum
		ind_max = compute_MAX_ind(NCCnew,(2*newu+1)*(2*newv+1));

		c++;
	}

	if(ind_ref != ind_max)
	{
		deltau = ind_max/(2*newv+1) - ind_ref/(2*newv+1);
		deltav = ind_max%(2*newv+1) - ind_ref%(2*newv+1);

		// update position of the new center (current maximum)
		du += deltau;
		dv += deltav;

		failed=true;
	}

	// --- CRISTIAN MOD START ---
	if(ps1) 
		delete ps1;
	if(ps2) 
		delete ps2;
	// --- CRISTIAN MOD END ---

	delete[] missu;
	delete[] missv;
}


//void compute_Alignment( NCC_parms_t *NCC_params, REAL_T *NCC_xy, REAL_T *NCC_xz, REAL_T *NCC_yz,
//                                          int dimi, int dimj, int dimk, int ind_xy, int ind_xz, int ind_yz, NCC_descr_t *result) {
void compute_Alignment( NCC_parms_t *NCC_params, iom::real_t *NCC_xy, iom::real_t *NCC_xz, iom::real_t *NCC_yz,
                                            int dimi, int dimj, int dimk, int dx1, int dx2, int dy1, int dy2, int dz1, int dz2, bool failed_xy, bool failed_xz, bool failed_yz, NCC_descr_t *result) {

	int w1x, w2x, w1y, w2y, w1z, w2z;

	compute_NCC_width(NCC_params,NCC_xy,(2*dimj+1),(dimi*(2*dimj+1)+dimj),NCC_params->wRangeThr_i,NCC_params->wRangeThr_j,failed_xy, w1x,w1y);
	compute_NCC_width(NCC_params,NCC_xz,(2*dimk+1),(dimi*(2*dimk+1)+dimk),NCC_params->wRangeThr_i,NCC_params->wRangeThr_k,failed_xz, w2x,w1z);
	compute_NCC_width(NCC_params,NCC_yz,(2*dimk+1),(dimj*(2*dimk+1)+dimk),NCC_params->wRangeThr_j,NCC_params->wRangeThr_k,failed_yz, w2y,w2z);

	compute_NCC_alignment(NCC_params,result,0,dx1,NCC_xy[(dimi*(2*dimj+1)+dimj)],w1x,dx2,NCC_xz[(dimi*(2*dimk+1)+dimk)],w2x);
	compute_NCC_alignment(NCC_params,result,1,dy1,NCC_xy[(dimi*(2*dimj+1)+dimj)],w1y,dy2,NCC_yz[(dimj*(2*dimk+1)+dimk)],w2y);
	compute_NCC_alignment(NCC_params,result,2,dz1,NCC_xz[(dimi*(2*dimk+1)+dimk)],w1z,dz2,NCC_yz[(dimj*(2*dimk+1)+dimk)],w2z);
}


void enhance ( iom::real_t *im, int imLen, int graylevels, NCC_parms_t *NCC_params ) {
/*
 * the enhancement transformation is a multi-linear curve with n_transforms rescaled linear
 * tranformations
 * for i=0, ..., (n_transforms-1), percentiles[i] contains the fraction of pixels to which the
 * (i+1)-th rescaled linear transformation has to be applied; percentiles[n_transforms-1] must be 1.00
 * for i=1, ..., n_transforms, c[i-1] amd c[i] contain the lowest and highest value corresponding to
 * the the i-th rescaled linear transformation; c[0] must be 0.00 and c[n_transforms] must be 1.00
 */

	int n_transforms = NCC_params->n_transforms;
	iom::real_t *percentiles = NCC_params->percents;
	iom::real_t *c = NCC_params->c; // tranformed values of thresholds
	iom::real_t *thresholds = new iom::real_t[n_transforms+1];
	iom::real_t *a = new iom::real_t[n_transforms+1];
	iom::real_t *b = new iom::real_t[n_transforms+1];
	int i, j;
	bool found;

	thresholds[0] = (iom::real_t)0.00;

	stack_percentiles(im,imLen,graylevels,percentiles,thresholds,n_transforms);

	for ( i=1; i<=n_transforms; i++ ) {
			a[i] = (c[i] - c[i-1]) / (thresholds[i] - thresholds[i-1]);
			b[i] = c[i] - a[i]*thresholds[i];
	}

	for ( i=0; i<imLen; i++ ) {
			binary_search(thresholds,n_transforms,im[i],found,j);
			im[i] = a[j]*im[i] + b[j];
	}

	delete[] thresholds;
	delete[] a;
	delete[] b;
}


void stack_percentiles ( iom::real_t *im, int imLen, int graylevels,
                                                 iom::real_t *percentiles, iom::real_t *thresholds, int n_percentiles ) {

	iom::real_t d = (iom::real_t)1.0 / (iom::real_t)graylevels;
	iom::real_t *cumsum = new iom::real_t[graylevels];
	int i, j;

	for ( i=0; i<graylevels; i++ )
		cumsum[i] = 0;

	for ( i=0; i<imLen; i++ ) {
		j = MIN((int)floor(im[i]/d),graylevels-1); // guarantees that index is within limits
		cumsum[j]++;
	}

	cumsum[0] /= imLen;
	for ( i=1; i<graylevels; i++ )
		cumsum[i] = cumsum[i]/imLen + cumsum[i-1];

	// guarantees that last cumulative fraction is 1.0
	cumsum[graylevels-1] = 1.0;

	for ( i=0, j=1; j<n_percentiles; ) {
		if ( i == graylevels ) {
			DISPLAY_ERROR("i out of limits");
			exit(1);
		}
		if ( cumsum[i] >= percentiles[j-1] ) {
			thresholds[j] = d * i;
			j++;
		}
		else
			i++;
	}
	thresholds[n_percentiles] = (iom::real_t)1; // to avoid round off errors

	delete[] cumsum;
}
