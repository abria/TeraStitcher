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
* 2018-08-29. Giulio.     @CHECKED checked the code implementing the optimization
* 2018-08-26. Giulio.     @ADDED code and comments concerning the optimization in computing multiple NCCs on the same images with different overlap
* 2015-03-20. Giulio.     @CHANGED newu and newv have been moved as parameters in compute_Neighborhood
*/

/*
 * compute_funcs.h
 *
 *  Created on: September 2010
 *      Author: iannello
 *
 *  Last revision: May, 31 2011
 */


#ifndef COMPUTE_FUNCS_H
#define COMPUTE_FUNCS_H


# include "my_defs.h"

# include "CrossMIPs.h"


# define GRAY_LEVELS    65536


// Constants and functions to support efficient computation of the NCC.
// The images are tiled, i.e. partitioned in submatrices (tiles) of the same size
// (TILE_SIDE x TILE_SIDE), and the sum of each tile is computed once and stored
// in support matrices (named 'ps1' and 'ps2'). 
// These precomputed sums avoid the computation of the same sums many times.

// must be > 0
#define TILE_SIDE 32

void seq_cpu_compute_partial_sums(iom::real_t *image1, iom::real_t *image2, int height, int width, 
								  iom::real_t *ps1 = (iom::real_t *) 0, iom::real_t *ps2 = (iom::real_t *) 0 );
/* Computes the sums of the tiles which 2D matrices 'image1' and 'image2' are divided
 *
 * Parameters:
 *    image1, image2: 2D matrices of size height x width
 *    ps1, ps2:       pre-allocated matices of size (height/TILE_SIDE) x (width/TILE_SIDE)
 *                    containing te sum of each tile
 *
 * WARNING: if either height or width are not multiple of TILE_SIDE, the sum of some smaller 
 *          tiles at the borders are not computed (i.e. divisions are between integers)
 */

void compute_3_MIPs ( iom::real_t *A, iom::real_t *B,
					  iom::real_t *MIP_xy1, iom::real_t *MIP_xz1, iom::real_t *MIP_yz1, 
					  iom::real_t *MIP_xy2, iom::real_t *MIP_xz2, iom::real_t *MIP_yz2, 
					  int dimi_v, int dimj_v, int dimk_v, int stridei, int stridek );

void compute_NCC_map ( iom::real_t *NCC_map, iom::real_t *MIP_1, iom::real_t *MIP_2, 
					       int dimu, int dimv, int delayu, int delayv );
/* Parameters:
 *    NCC_map:      externaly allocated buffer of size (2*delayu+1) rows x (2*delayv+1)
 *    MIP_1, MIP_2: matrices of size dimu x dimv
 *
 * returns in NCC_map a map of Normalized Cross Correlations (NCC) 
 * each pixel of the map is the NCC computed on two 2D sub-matrices of identical size of 
 * matrices MIP_1 and MIP_2
 * more specifically, pixel (u,v) of the map is the NCC of the sub-matrices (using a 
 * Python-like notation:
 *
 *    MIP_1[ max(0,u-delayu) : (dimu - 1 + min(0,u-delayu)) , max(0,v-delayv) : (dimv - 1 + min(0,v-delayv)) ]
 *    MIP_2[ max(0,delayu-u) : (dimu - 1 - min(0,delayu-u)) , max(0,delayv-v) : (dimv - 1 + min(0,delayv-v)) ]
 *
 * these matrices have size:
 *
 *    (dimu - |u - delayu|) x (dimv - |v - delayv|)
 */

iom::real_t compute_NCC ( iom::real_t *MIP_1, iom::real_t *MIP_2, int dimu, int dimv, int u, int v, 
						  iom::real_t *ps1 = (iom::real_t *) 0, iom::real_t *ps2 = (iom::real_t *) 0);
/* Parameters:
 *    MIP_1, MIP_2: images of size dimu X dimv for which it is computed the NCC assuming an overlap defined by displacements u and v 
 *    (that indicate how much move MIP_1 with respect to MIP2).
 *    ps1, ps2: matrices of size (dimu / TILE_SIDE) * (dimv / TILE_SIDE) that contain precomputed sums of square submatrices of MIP_1, MIP_2
 *    of size TILE_SIDE x TILE_SIDE.
 *    Each square submatrix is referred to as "tile" and the original images are partitioned in "tiles". These precomputed sums avoid 
 *    to repeat in different calls of this in this function the same computations.
 *
 * returns the Normalized Cross Correlations (NCC) of the two matrices assuming an overlap defined by displacements u and v
 */

int compute_MAX_ind ( iom::real_t *vect, int len );

void compute_Neighborhood ( NCC_parms_t *NCC_params, iom::real_t *NCC, int delayu, int delayv, int newu, int newv, int ind, 
						    iom::real_t *MIP_1, iom::real_t *MIP_2, int dimu, int dimv, iom::real_t *NCCnew, int &du, int &dv, bool &failed)  ;

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
 * 
 * DETAILED EXPLANATION
 * Given an NCC map 'NCC' with dimensions 'delayu' x 'delayv' around the initial alignment, extensions 'newu' and 'newu'
 * of the new to be returned NCC map, the index 'ind' of its maximum and an empty NCC map 'NCCnew' with dimensions 
 * NCC_params->wRangeThr x NCC_params->wRangeThr, returns in 'NCCnew' the NCC map with center in the the NCC maximum of 
 * an NCC map between MIPs 'MIP_1' and 'MIP_2', that have dimensions 'dimu' x 'dimv';
 * also the position (du,dv) of this maximum, relative to the initial alignment is returned
 */

//void compute_Alignment( NCC_parms_t *NCC_params, REAL_T *NCC_xy, REAL_T *NCC_xz, REAL_T *NCC_yz, 
//					    int dimi, int dimj, int dimk, int ind_xy, int ind_xz, int ind_yz, NCC_descr_t *result);
void compute_Alignment( NCC_parms_t *NCC_params, iom::real_t *NCC_xy, iom::real_t *NCC_xz, iom::real_t *NCC_yz,
					    int dimi, int dimj, int dimk, int dx1, int dx2, int dy1, int dy2, int dz1, int dz2,  bool failed_xy, bool failed_xz, bool failed_yz, NCC_descr_t *result);
/* given the three NCCs, the corresponding indices of NCC maxima and their displacements with respect to the 
 * initial alignment, returns in the parameter result the three alignments and the corresponding reliability indices;
 * alignments represent offsets of the second stack with respect to the first one
 */

void enhance ( iom::real_t *im, int imLen, int graylevels, NCC_parms_t *NCC_params );
/* enhance contrast of an image
 * INPUT parameters:
 *    im is a pointer to actual pixels, represented as a linear sequence of length imLen
 *    graylevels is the number of discrete gray levels to be used to discriminate pixels for
 *    enhancement
 *    NCC_params contains data describing the transformation used for enhancement (see 
 *    CrossMIPs.h for details)
 */

void stack_percentiles ( iom::real_t *im, int imLen, int graylevels, 
						 iom::real_t *percentiles, iom::real_t *thresholds, int n_percentiles );
/* returns thresholds on pixel values to implement an arbitrary grayscale transformation
 * INPUT parameters:
 *    im is a pointer to actual pixels represented as a linear sequence of length imLen
 *    percentiles is a pointer to an array of n_percentiles real values representing
 *    the percentiles of pixels to which the i-th rescaled linear transformation has
 *    to be applied (i=1,...,n_percentiles)
 * INPUT/OUTPUT parmeter:
 *    thresholds is a pointer to an array of (n_percentiles+1) real values representing the
 *    pixel values interval to which the i-th linear transformation has to be applied
 *    (i=1,...,n_percentiles)
 * PRE: the value of percentiles[n_percentiles-1] must be (REAL_T)1.0
 *      the value of threshold[0] must be 0
 */

# ifdef _PAR_VERSION
void init_configuration ( );
# endif

# endif
