//------------------------------------------------------------------------------------------------
// Copyright (c) 2012  Alessandro Bria and Giulio Iannello (University Campus Bio-Medico of Rome).  
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
* 2019-04-25. Giulio.     @ADDED   using subsampling to compute the first NCC map
* 2017-04-20. Giulio.     @FIXED   a very old bug: in the third call to 'compute_Neighborhood' eas used 'NCC_params->wRangeThr_i' instead of 'NCC_params->wRangeThr_j' 
* 2016-01-27. Giulio.     @ADDED   checks to not consider NCC maps when they have a dimension that it is too small
* 2016-01-27. Giulio.     @CHANGED the way initial checks are managed; the search can now be as much lurge as allowed by NCC parameter minDim_NCCsrc
* 2016-01-27. Giulio.     @CHANGED the way initial checks are managed; the search are is now dynamically resized if overlap is too small
* 2015-03-20. Giulio.     @CHANGED different dimensions for the new NCC to be computed are passed to compute_Neighborhood
* 2015-03-20. Giulio.     @CHANGED newu and newv have been moved as parameters in compute_Neighborhood
*/

/*
 * libcrossmips.cpp
 *
 *  Created on: September 2010
 *      Author: iannello
 *
 *  Last revision: May, 31 2011
 *
 * --- revision of May, 31 2011 -------------------------------------------------------------------
 * Parameters can be controlled through a structure defined in CrossMIPs.h which can be 
 * passed as a parameter to norm_cross_corr_mips
 * Default values are set in implementation of norm_cross_corr_mips
 *
 * A structure is returned containing the thre offsets of the second 3d stack with respect to 
 * the first one. The structure is defined in CrossMIPs.h and it contains also information
 * about the reliability and precision of the alignement. Reliability of each offset is measured 
 * by the value of the NCC maxima used to evaluate it. Precision is measured by the distance
 * between NCC maxima and pixels equal to a given fraction of the maxima in the direction of
 * each offset. The algorithm used to determine all these quantities is in function 
 * compute_NCC_alignment, which is internal to the file compute_funcs.cpp
 * ----------------------------------------------------------------------------------------------
 */


# include <math.h>
# include <stdio.h>
# include <stdlib.h>

# include "CrossMIPs.h"

# include "compute_funcs.h"

// 2014-11-04 Giulio. @ADDED improved the conditionally compiled code that can be used to save images involved in NCC

/* Flag _WRITE_IMGS enables the saving of MIPs and NCCs
 * Flag _WRITE_STKS enables the saving of the input substacks too
 */

//#define _WRITE_IMGS
//#define _WRITE_NCC_ONLY

#if defined(_WRITE_IMGS) || defined(_WRITE_NCC_ONLY)
# include <string.h>
//# include "../imagemanager/IM_config.h"
//# include "../imagemanager/Tiff3DMngr.h"

# define _MAX_FNAME_LEN   1000
# define _PREFIX_LEN   8   // digits used to generate unique file names
# define _MAX_COUNTER  100000000 // max counter value must be 10^(_PREFIX_LEN) - 1

// used to generate unique identifiers for image files
static int _counter = 0; 
# endif


void fill_NCC_map (iom::real_t *MIP_1, iom::real_t  *MIP_2,  int dimu,  int dimv, 
					int delayu, int delayv, int wRangeThr_u, int wRangeThr_v, int ind_uv, iom::real_t *NCC_uv );

bool write_3D_stack ( char *fname, iom::real_t *stck, int dimi, int dimj, int dimk );


NCC_descr_t *norm_cross_corr_mips ( iom::real_t *A, iom::real_t *B,
						    int dimk, int dimi, int dimj, int nk, int ni, int nj,
							int delayk, int delayi, int delayj, int side, NCC_parms_t *NCC_params ) {
#if CM_VERBOSE > 1
	printf("\nin libcrossmips::norm_cross_corr_mips(A[%.6f-%.6f], B[%.6f-%.6f], dimk[%d], dimi[%d], dimj[%d], nk[%d], ni[%d], nj[%d], delayk[%d], delayi[%d], delayj[%d], side[%d]\n",
		A[0], A[(dimk-1)*dimj*dimi + (dimi-1)*dimj + dimj -1], B[0], B[(dimk-1)*dimj*dimi + (dimi-1)*dimj + dimj -1], dimk, dimi, dimj, nk, ni, nj, delayk, delayi, delayj, side );
#endif

#if defined(_WRITE_IMGS) || defined(_WRITE_NCC_ONLY)
	char _fname[_MAX_FNAME_LEN]; // used ti assemble file manes
	if ( _counter < _MAX_COUNTER )
		sprintf(_fname,"%08d",_counter); // WARNING: check the format string; the number of digits must be _PREFIX_LEN
	else
        throw iom::exception("CrossMIPs (_WRITE_IMGS enabled): too much substacks");
	_counter++;

#ifdef _WRITE_STKS
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_STCK_A.dat")-_PREFIX_LEN,A,dimi,dimj,dimk);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_STCK_B.dat")-_PREFIX_LEN,B,dimi,dimj,dimk);
#endif
#endif

	NCC_descr_t *result = new NCC_descr_t;
	int dimk_v, dimi_v, dimj_v;
	int stridei, stridek; // arrays are stored along j dimension, hence it does not need a stride
	iom::real_t *A1; // starting pixels of volume A to be scanned

	iom::real_t *MIP_xy1, *MIP_xz1, *MIP_yz1;
	iom::real_t *MIP_xy2, *MIP_xz2, *MIP_yz2;
	int i;

	iom::real_t *NCC_xy, *NCC_xz, *NCC_yz;

	iom::real_t *temp;

	int ind_xy, ind_xz, ind_yz;
	int dx1, dx2, dy1, dy2, dz1, dz2;

	bool failed_xy=false, failed_xz=false, failed_yz=false;

	//bool NCC_xy_valid = true, NCC_xz_valid = true, NCC_yz_valid = true;

	bool allocated = !NCC_params;
	if ( !NCC_params ) {

		// Alessandro - 31/05/2013 - parameter MUST be passed (and controlled) by the caller
        throw iom::exception("CrossMIPs: missing configuration parameters");

		//// ************************** SET DEFAULT PARAMETERS ****************
		//NCC_params = new NCC_parms_t;
		//// to change default behavior change the following parameters
		//NCC_params->enhance      = false;
		//NCC_params->maxIter		 = 2;
		//NCC_params->maxThr       = CM_DEF_MAX_THR;
		//NCC_params->widthThr     = CM_DEF_WIDTH_THR;
		//NCC_params->wRangeThr    = CM_DEF_W_RANGE_THR;
		//NCC_params->UNR_NCC      = CM_DEF_UNR_NCC;
		//NCC_params->INF_W        = CM_DEF_INF_W;
		//NCC_params->INV_COORD    = 0;

		///* Example of settings for fields 'percents' and 'c'
		// * percents = < 0.1, 0.9, 1.0 > and c = < 0.0, 0,0, 1.0, 1.0 >
		// * means that:
		// * - 10% of pixels with lower values are set to 0 (first transformation maps pixels from c[0]=0.0 to c[1]=0.0
		// * - 10% of pixels with higher values are set to 1 (third transformation maps pixels from c[2]=1.0 to c[3]=1.0
		// * - 80% of pixels with intermediate values are mapped to the interval [0,1] (second transformation maps pixels
		// *   from c[1]=0.0 to c[2]01.0
		// */

		//if ( enhance ) {
		//	//// two transformations defined by pairs (percentiles,value): {(0.0,0.0),(0.8,0.2),(1.0,1.0)}
		//	//NCC_params->n_transforms = 2;
		//	//NCC_params->percents     = new REAL_T[NCC_params->n_transforms];
		//	//NCC_params->c            = new REAL_T[NCC_params->n_transforms+1];
		//	//// percents must have n_transforms elements and percents[n_transforms-1] must be 1.00
		//	//NCC_params->percents[0]  = (REAL_T) 0.80;
		//	//NCC_params->percents[1]  = (REAL_T) 1.00;
		//	//// c must have (n_transforms+1) elements, c[0] must be 0.00 and c[n_transforms] must be 1.00
		//	//NCC_params->c[0]         = (REAL_T) 0.00;
		//	//NCC_params->c[1]         = (REAL_T) 0.20;
		//	//NCC_params->c[2]         = (REAL_T) 1.00;

		//	// three transformations defined by pairs (percentiles,value): {(0.0,0.0),(0.1,0.0),(0.95,1.0),(1.0,1.0)}
		//	NCC_params->n_transforms = 3;
		//	NCC_params->percents     = new iom::real_t[NCC_params->n_transforms];
		//	NCC_params->c            = new iom::real_t[NCC_params->n_transforms+1];
		//	// percents must have n_transforms elements and percents[n_transforms-1] must be 1.00
		//	NCC_params->percents[0]  = (iom::real_t) 0.10;
		//	NCC_params->percents[1]  = (iom::real_t) 0.99;
		//	NCC_params->percents[2]  = (iom::real_t) 1.00;
		//	// c must have (n_transforms+1) elements, c[0] must be 0.00 and c[n_transforms] must be 1.00
		//	NCC_params->c[0]         = (iom::real_t) 0.00;
		//	NCC_params->c[1]         = (iom::real_t) 0.00;
		//	NCC_params->c[2]         = (iom::real_t) 1.00;
		//	NCC_params->c[3]         = (iom::real_t) 1.00;
		//}
	}

	// Alessandro - 23/03/2013 - this is the old check, but is seems wrong and it does not throw any exception. The WARNING written into
	// CrossMIPs.h says: "moreover controlling parameter wRangeThr is supposed to be not greater than MAX(delayi,delayj,delayk)"
	//if ( NCC_params->wRangeThr > 2*MAX(delayi,MAX(delayj,delayk))+1 )
	//{
	//	char err_msg[500];
	//	sprintf(err_msg, "parameter wRangeThr[=%d] is too large with respect to 2*delayi/j/k +1 [=%d]", NCC_params->wRangeThr, 2*MAX(delayi,MAX(delayj,delayk))+1);
	//	DISPLAY_ERROR(err_msg);
	//}
	
	// These checks are required by the lines (in compute_funcs.cpp):
	// 	  initu = MIN(MAX(0,ind_max/(2*delayv+1) - newu),2*(delayu - newu));
	//    initv = MIN(MAX(0,ind_max%(2*delayv+1) - newv),2*(delayv - newv));
	// where newu = newv = NCC_params->wRangeThr and initu, initv must be positive integers
	if ( NCC_params->wRangeThr_i > delayi  || NCC_params->wRangeThr_j > delayj || NCC_params->wRangeThr_k > delayk )
	{
		// Alessandro - 31/05/2013 - throwing an exception instead of automatically correcting parameters
		char err_msg[1000];
		sprintf(err_msg, "CrossMIPs: one or more parameters: wRangeThr_i[=%d], wRangeThr_j[=%d], wRangeThr_k[=%d] are too large with respect to: delayi[=%d], delayj[=%d], delayik[=%d]", 
			NCC_params->wRangeThr_i, NCC_params->wRangeThr_j, NCC_params->wRangeThr_k, delayi, delayj ,delayk);
		throw iom::exception(err_msg);
	}

	/************************************************************************************************************************************************
	 * 2017-01-27. Giulio. OLD CODE - MAINTAIN UNTIL TESTS OF NEW CODE ARE COMPLETED
	 ***********************************************************************************************************************************************/
    //// skipping check for 2D images: see Alessandro's comment in PDAlgoMIPNCC.cpp on 21/08/2013
    //if(dimk != 1)
    //{
    //    // Alessandro - 23/03/2013 - added check to verify precondition written into CrossMIPs.h that says:
    //    // "in practice the dimensions of the MIPS (depending on dimi, dimj, dimk, ni, nj, nk) have to be large enough with respect to delayi, delayj, delayk"
    //    if(side == NORTH_SOUTH && ((dimi - ni < 2*delayi+1) || (dimj - nj < 2*delayj+1) || (dimk - nk < 2*delayk+1)))
    //        throw iom::exception("CrossMIPs: the search region is too big with respect to the overlapping region. "
    //                          "Overlapping extent should be > 2*delay+1 for each direction where delay is the "
    //                          "search region extent along that direction.");
    //    if(side == WEST_EAST   && ((dimj - nj < 2*delayi+1) || (dimi - ni < 2*delayj+1) || (dimk - nk < 2*delayk+1)))
    //        throw iom::exception("CrossMIPs: the search region is too big with respect to the overlapping region. "
    //                          "Overlapping extent should be > 2*delay+1 for each direction where delay is the "
    //                          "search region extent along that direction.");
    //}

 	/************************************************************************************************************************************************
	 * 2017-01-27. Giulio. NEW CODE - IMPLEMENT A MORE FLEXIBLE STRATEGY 
	 ***********************************************************************************************************************************************/
    // Alessandro - 23/03/2013 - added check to verify precondition written into CrossMIPs.h that says:
    // "in practice the dimensions of the MIPS (depending on dimi, dimj, dimk, ni, nj, nk) have to be large enough with respect to delayi, delayj, delayk"

	// 2017-02-26. Giulio. Old code. I do not remember because I distinguish these two cases. Now they look wrong
	//if ( side == NORTH_SOUTH ) { 
	//	delayi = MIN(delayi,(dimi - ni - 1) / 2);
	//	delayj = MIN(delayj,(dimj - nj - 1) / 2);
	//	delayk = MIN(delayk,(dimk - nk - 1) / 2);
	//}
	//else if ( side == WEST_EAST ) {
	//	delayi = MIN(delayj,(dimj - nj - 1) / 2);
	//	delayj = MIN(delayi,(dimi - ni - 1) / 2);
	//	delayk = MIN(delayk,(dimk - nk - 1) / 2);
	//}
	//else
	//	throw iom::exception("CrossMIPs: undefined side for tile alignment.");
	
	// 2017-02-26. Giulio. Check if search areas are too big with respect to overlap
	delayi = MIN(delayi,MAX(0,dimi - ni - NCC_params->minDim_NCCsrc));
	delayj = MIN(delayj,MAX(0,dimj - nj - NCC_params->minDim_NCCsrc));
	delayk = MIN(delayk,MAX(0,dimk - nk - NCC_params->minDim_NCCsrc));

	// 2017-02-26. Giulio. set flags controlling the computation of NCC maps
	//if ( delayi < NCC_params->minDim_NCCmap || delayj < NCC_params->minDim_NCCmap )
	//	NCC_xy_valid = false;
	//if ( delayi < NCC_params->minDim_NCCmap || delayk < NCC_params->minDim_NCCmap )
	//	NCC_xz_valid = false;
	//if ( delayj < NCC_params->minDim_NCCmap || delayk < NCC_params->minDim_NCCmap )
	//	NCC_yz_valid = false;

	// 2017-01-27. Giulio. Check if wRangeThr fields have to be reduced 
	// Note that if these values are very small it is likely that alignment will be unreliable 
	// because it is unlikely that values in the NCC map can quickly decrease if the search area is too narrow 
	NCC_params->wRangeThr_i = MIN(NCC_params->wRangeThr_i,delayi);
	NCC_params->wRangeThr_j = MIN(NCC_params->wRangeThr_j,delayj);
	NCC_params->wRangeThr_k = MIN(NCC_params->wRangeThr_k,delayk);

 	/************************************************************************************************************************************************
	 * 2017-01-27. Giulio. END NEW CODE 
	 ***********************************************************************************************************************************************/


	/*
	 * the following conde assumes that:
	 * - 2D matrices are stored along their second dimension
	 *   for instance, in a MIP projected along y, which includes dimensions x and z, the values of z are stored in sequence
	 * - 3D stacks are stored as a sequence of horizontal planes and each plane is stored as a 2D matrix (along its second dimension)
	 *   this assumption has effect only in the computation of the three MIP projections
	 */

# ifdef _PAR_VERSION

	init_configuration();

# endif
			
	// compute parameters to scan the volumes
	if ( side == NORTH_SOUTH ) {
		dimk_v = dimk;
		dimi_v = dimi - ni;
		dimj_v = dimj;
		stridei = 0; // rows are entirely scanned
		stridek = ni*dimj; // pixels to be skipped when changing slice are all contiguous
		A1 = A + stridek; // a block of contiguous pixels of volume A have to be skipped
	}
	else if ( side == WEST_EAST ) {
		dimk_v = dimk;
		dimi_v = dimi;
		dimj_v = dimj - nj;
		stridei = nj; // rows are partially scanned and nj pixels have to be skipped when changing row
		stridek = 0; // no more pixels have to be skipped when changing slice 
		A1 = A + stridei; // a partial row of pixels of volume A have to be skipped
	}
	else
		throw iom::exception("CrossMIPs: unexpected alignment configuration");

	// alloca le 6 immagini per i MIP
	MIP_xy1 = new iom::real_t[dimi_v*dimj_v];
	for ( i=0; i<(dimi_v*dimj_v); i++ )
		MIP_xy1[i] = 0;
	MIP_xz1 = new iom::real_t[dimi_v*dimk_v];
	for ( i=0; i<(dimi_v*dimk_v); i++ )
		MIP_xz1[i] = 0;
	MIP_yz1 = new iom::real_t[dimj_v*dimk_v];
	for ( i=0; i<(dimj_v*dimk_v); i++ )
		MIP_yz1[i] = 0;

	MIP_xy2 = new iom::real_t[dimi_v*dimj_v];
	for ( i=0; i<(dimi_v*dimj_v); i++ )
		MIP_xy2[i] = 0;
	MIP_xz2 = new iom::real_t[dimi_v*dimk_v];
	for ( i=0; i<(dimi_v*dimk_v); i++ )
		MIP_xz2[i] = 0;
	MIP_yz2 = new iom::real_t[dimj_v*dimk_v];
	for ( i=0; i<(dimj_v*dimk_v); i++ )
		MIP_yz2[i] = 0;

	compute_3_MIPs(A1,B,MIP_xy1,MIP_xz1,MIP_yz1,MIP_xy2,MIP_xz2,MIP_yz2,
									dimi_v,dimj_v,dimk_v,stridei,stridek);

#ifdef _WRITE_IMGS
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xy1.dat")-_PREFIX_LEN,MIP_xy1,dimi_v,dimj_v,1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xy2.dat")-_PREFIX_LEN,MIP_xy2,dimi_v,dimj_v,1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xz1.dat")-_PREFIX_LEN,MIP_xz1,dimi_v,dimk_v,1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xz2.dat")-_PREFIX_LEN,MIP_xz2,dimi_v,dimk_v,1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_yz1.dat")-_PREFIX_LEN,MIP_yz1,dimj_v,dimk_v,1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_yz2.dat")-_PREFIX_LEN,MIP_yz2,dimj_v,dimk_v,1);
#endif

	// calcola NCC su xy
	NCC_xy = new iom::real_t[(2*delayi+1)*(2*delayj+1)];
	for ( i=0; i<((2*delayi+1)*(2*delayj+1)); i++ )
		NCC_xy[i] = 0;

	if ( NCC_params->enhance ) {
		enhance(MIP_xy1,(dimi_v*dimj_v),GRAY_LEVELS,NCC_params);
		enhance(MIP_xy2,(dimi_v*dimj_v),GRAY_LEVELS,NCC_params);
	}

	//if ( NCC_xy_valid )
		compute_NCC_map(NCC_xy,MIP_xy1,MIP_xy2,dimi_v,dimj_v,delayi,delayj);

	// calcola NCC su xz
	NCC_xz = new iom::real_t[(2*delayi+1)*(2*delayk+1)];
	for ( i=0; i<((2*delayi+1)*(2*delayk+1)); i++ )
		NCC_xz[i] = 0;

	if ( NCC_params->enhance ) {
		enhance(MIP_xz1,(dimi_v*dimk_v),GRAY_LEVELS,NCC_params);
		enhance(MIP_xz2,(dimi_v*dimk_v),GRAY_LEVELS,NCC_params);
	}

	//if ( NCC_xz_valid )
		compute_NCC_map(NCC_xz,MIP_xz1,MIP_xz2,dimi_v,dimk_v,delayi,delayk);

	// calcola NCC su yz
	NCC_yz = new iom::real_t[(2*delayj+1)*(2*delayk+1)];
	for ( i=0; i<((2*delayj+1)*(2*delayk+1)); i++ )
		NCC_yz[i] = 0;

	if ( NCC_params->enhance ) {
		enhance(MIP_yz1,(dimj_v*dimk_v),GRAY_LEVELS,NCC_params);
		enhance(MIP_yz2,(dimj_v*dimk_v),GRAY_LEVELS,NCC_params);
	}

	//if ( NCC_yz_valid )
		compute_NCC_map(NCC_yz,MIP_yz1,MIP_yz2,dimj_v,dimk_v,delayj,delayk);

#ifdef _WRITE_IMGS
	if ( NCC_params->enhance ) {
		write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xy1_en.dat")-_PREFIX_LEN,MIP_xy1,dimi_v,dimj_v,1);
		write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xy2_en.dat")-_PREFIX_LEN,MIP_xy2,dimi_v,dimj_v,1);
		write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xz1_en.dat")-_PREFIX_LEN,MIP_xz1,dimi_v,dimk_v,1);
		write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_xz2_en.dat")-_PREFIX_LEN,MIP_xz2,dimi_v,dimk_v,1);
		write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_yz1_en.dat")-_PREFIX_LEN,MIP_yz1,dimj_v,dimk_v,1);
		write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_MIP_yz2_en.dat")-_PREFIX_LEN,MIP_yz2,dimj_v,dimk_v,1);
	}
#endif

#if defined(_WRITE_IMGS) || defined(_WRITE_NCC_ONLY)
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_NCC_xy.dat")-_PREFIX_LEN,NCC_xy,(2*delayi+1),(2*delayj+1),1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_NCC_xz.dat")-_PREFIX_LEN,NCC_xz,(2*delayi+1),(2*delayk+1),1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_NCC_yz.dat")-_PREFIX_LEN,NCC_yz,(2*delayj+1),(2*delayk+1),1);
#endif

	// max su xy, xz, yz
	ind_xy = compute_MAX_ind(NCC_xy,((2*delayi+1)*(2*delayj+1)));
	ind_xz = compute_MAX_ind(NCC_xz,((2*delayi+1)*(2*delayk+1)));
	ind_yz = compute_MAX_ind(NCC_yz,((2*delayj+1)*(2*delayk+1)));

// #if defined(USECUDA)
// 	if ( !(getenv("USECUDA_X_NCC")?1:0) )
// #endif
// 	{ // 2019-04-25. Giulio. This code must be executed only if NCC map has not be computed by the GPU
// 		// compute the missing points of the NCC map		
// 		fill_NCC_map(MIP_xy1,MIP_xy2,dimi_v,dimj_v,delayi,delayj,NCC_params->wRangeThr_i,NCC_params->wRangeThr_j,ind_xy,NCC_xy);		
// 		fill_NCC_map(MIP_xz1,MIP_xz2,dimi_v,dimk_v,delayi,delayk,NCC_params->wRangeThr_i,NCC_params->wRangeThr_k,ind_xz,NCC_xz);		
// 		fill_NCC_map(MIP_yz1,MIP_yz2,dimj_v,dimk_v,delayj,delayk,NCC_params->wRangeThr_j,NCC_params->wRangeThr_k,ind_yz,NCC_yz);
// 		// recompute the position of the maximum		
// 		ind_xy = compute_MAX_ind(NCC_xy,((2*delayi+1)*(2*delayj+1)));
// 		ind_xz = compute_MAX_ind(NCC_xz,((2*delayi+1)*(2*delayk+1)));
// 		ind_yz = compute_MAX_ind(NCC_yz,((2*delayj+1)*(2*delayk+1)));
// 	} 

	//if ( NCC_xy_valid ) {
		// NCC_xy: check neighborhood of maxima and search for better maxima if any
		temp = new iom::real_t[(2*NCC_params->wRangeThr_i+1)*(2*NCC_params->wRangeThr_j+1)];
		for ( i=0; i<((2*NCC_params->wRangeThr_i+1)*(2*NCC_params->wRangeThr_j+1)); i++ )
			temp[i] = 0;
		compute_Neighborhood(NCC_params,NCC_xy,delayi,delayj,NCC_params->wRangeThr_i,NCC_params->wRangeThr_j,ind_xy,MIP_xy1,MIP_xy2,dimi_v,dimj_v,temp,dx1,dy1, failed_xy);
		// substitute NCC map and delete the old one
		delete NCC_xy;
		NCC_xy = temp;
		temp = (iom::real_t *)0;
	//}
	//else {
	//	dx1 = dy1 = 0;
	//}

	//if ( NCC_xz_valid ) {
		// NCC_xz: check neighborhood of maxima and search for better maxima if any
		temp = new iom::real_t[(2*NCC_params->wRangeThr_i+1)*(2*NCC_params->wRangeThr_k+1)];
		for ( i=0; i<((2*NCC_params->wRangeThr_i+1)*(2*NCC_params->wRangeThr_k+1)); i++ )
			temp[i] = 0;
		compute_Neighborhood(NCC_params,NCC_xz,delayi,delayk,NCC_params->wRangeThr_i,NCC_params->wRangeThr_k,ind_xz,MIP_xz1,MIP_xz2,dimi_v,dimk_v,temp,dx2,dz1, failed_xz);
		// substitute NCC map and delete the old one
		delete NCC_xz;
		NCC_xz = temp;
		temp = (iom::real_t *)0;
	//} 
	//else {
	//	dx2 = dz1 = 0;
	//}


	//if ( NCC_yz_valid ) {
		// NCC_yz: check neighborhood of maxima and search for better maxima if any
		temp = new iom::real_t[(2*NCC_params->wRangeThr_j+1)*(2*NCC_params->wRangeThr_k+1)];
		for ( i=0; i<((2*NCC_params->wRangeThr_j+1)*(2*NCC_params->wRangeThr_k+1)); i++ )
			temp[i] = 0;
		compute_Neighborhood(NCC_params,NCC_yz,delayj,delayk,NCC_params->wRangeThr_j,NCC_params->wRangeThr_k,ind_yz,MIP_yz1,MIP_yz2,dimj_v,dimk_v,temp,dy2,dz2, failed_yz);
		// substitute NCC map and delete the old one
		delete NCC_yz;
		NCC_yz = temp;
		temp = (iom::real_t *)0;
	//} 
	//else {
	//	dy2 = dz2 = 0;
	//}

#if defined(_WRITE_IMGS) || defined(_WRITE_NCC_ONLY)
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_newNCC_xy.dat")-_PREFIX_LEN,NCC_xy,(2*NCC_params->wRangeThr_i+1),(2*NCC_params->wRangeThr_j+1),1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_newNCC_xz.dat")-_PREFIX_LEN,NCC_xz,(2*NCC_params->wRangeThr_i+1),(2*NCC_params->wRangeThr_k+1),1);
	write_3D_stack(strcpy(_fname+_PREFIX_LEN,"_newNCC_yz.dat")-_PREFIX_LEN,NCC_yz,(2*NCC_params->wRangeThr_j+1),(2*NCC_params->wRangeThr_k+1),1);
#endif

	//compute_Alignment(NCC_params,NCC_xy,NCC_xz,NCC_yz,delayi,delayj,delayk,ind_xy,ind_xz,ind_yz,result);
	
	compute_Alignment(NCC_params,NCC_xy,NCC_xz,NCC_yz,
		NCC_params->wRangeThr_i,NCC_params->wRangeThr_j,NCC_params->wRangeThr_k,dx1,dx2,dy1,dy2,dz1,dz2,failed_xy, failed_xz, failed_yz, result);

	if ( side == NORTH_SOUTH ) 
		result->coord[0] += ni;
	else if ( side == WEST_EAST ) 
		result->coord[1] += nj;
	else
		throw iom::exception("CrossMIPs: unexpected alignment configuration");

	if ( allocated ) {
		delete NCC_params->percents;
		delete NCC_params->c;
		delete NCC_params;
	}

	delete MIP_xy1;
	delete MIP_xz1;
	delete MIP_yz1;
	delete MIP_xy2;
	delete MIP_xz2;
	delete MIP_yz2;

	delete NCC_xy;
	delete NCC_xz;
	delete NCC_yz;

	// temp must not be deleted
#if CM_VERBOSE > 1
	printf("\tReturning\n\t\tV[%d, %.6f, %d]\n\t\tH[%d, %.6f, %d]\n\t\tD[%d, %.6f, %d]\n\n", 
		result->coord[0], result->NCC_maxs[0], result->NCC_widths[0],
		result->coord[1], result->NCC_maxs[1], result->NCC_widths[1],
		result->coord[2], result->NCC_maxs[2], result->NCC_widths[2]);
#endif
	return result;
}


void fill_NCC_map (iom::real_t *MIP_1, iom::real_t  *MIP_2,  int dimu,  int dimv, 
					int delayu, int delayv, int wRangeThr_u, int wRangeThr_v, int ind_uv, iom::real_t *NCC_uv ) {

	int r, c; // indices of max
	int rf1, rf2, rl1, rl2, cf1, cf2, cl1, cl2; // bound: '1' variables are for the first loop, '2' variables are for the second loop
	int ir, ic; // loop indices
//	iom::real_t temp;
	
	// xy
	r = ind_uv/(2*delayv+1); // row of max
	rf1 = rf2 = std::max<int>(0,r-wRangeThr_u); // first row of map
	if ( (r - rf1)%2 == 0 ) // adjust one of the two bounds
		rf1++;
	else 
		rf2++;
	rl1 = rl2 = std::min<int>(2*delayu,r+wRangeThr_u); // last row of map
	if ( (rl1 - r)%2 == 0 ) // adjust one of the two bounds
		rl1--;
	else
		rl2--;
	c = ind_uv%(2*delayv+1);  // column of max
	cf1 = cf2 = std::max<int>(0,c-wRangeThr_v); // first column of map
	if ( (c - cf2)%2 == 0 ) // adjust bound '2'
		cf2++;
	cl1 = cl2 = std::min<int>(2*delayv,c+wRangeThr_v); // last column of map
	if ( (c - cl2)%2 == 0 ) // adjust bound '2'
		cl2--;
	// compute all points of empty lines
// 	printf("%d %d %d %d %d %d \n", r, c, delayu, delayv, wRangeThr_u, wRangeThr_v);
// 	printf("%d %d %d %d %d %d %d %d \n", rf1, rf2, rl1, rl2, cf1, cf2, cl1, cl2);
	for ( ir=rf1; ir<=rl1; ir+=2 ) {
		for ( ic=cf1; ic<=cl1; ic++ ) { 
			NCC_uv[ir*(2*delayv+1)+ic] = compute_NCC(MIP_1,MIP_2,dimu,dimv,(ir-delayu),(ic-delayv),(iom::real_t *)0,(iom::real_t *)0);
// 			temp = compute_NCC(MIP_1,MIP_2,dimu,dimv,(ir-delayu),(ic-delayv),(iom::real_t *)0,(iom::real_t *)0);
// 			if ( NCC_uv[ir*(2*delayv+1)+ic] != temp )
// 				printf("loop 1: %d %d %f %f\n", ir, ic, NCC_uv[ir*(2*delayv+1)+ic], temp);
		}
	}
	// compute missing points in non empty lines
	for ( ir=rf2; ir<=rl2; ir+=2 ) {
		for ( ic=cf2; ic<=cl2; ic+=2 ) { 
			NCC_uv[r*(2*delayv+1)+c] = compute_NCC(MIP_1,MIP_2,dimu,dimv,(ir-delayu),(ic-delayv),(iom::real_t *)0,(iom::real_t *)0);
// 			temp = compute_NCC(MIP_1,MIP_2,dimu,dimv,(ir-delayu),(ic-delayv),(iom::real_t *)0,(iom::real_t *)0);
// 			if ( NCC_uv[ir*(2*delayv+1)+ic] != temp )
// 				printf("loop2: %d %d %f %f\n", ir, ic, NCC_uv[ir*(2*delayv+1)+ic], temp);
		}
	}
}



#if defined(_WRITE_IMGS) || defined(_WRITE_NCC_ONLY)

bool write_3D_stack ( char *fname, iom::real_t *stck, int dimi, int dimj, int dimk ) {
	FILE *fout;
	int i, j, k;

	if ( (fout = fopen(fname,"wb")) == NULL ) return false;
	
	fwrite(&dimi,sizeof(int),1,fout);
	fwrite(&dimj,sizeof(int),1,fout);
	fwrite(&dimk,sizeof(int),1,fout);

	for ( k=0; k<dimk; k++ )
		for ( j=0; j<dimj; j++ )
			for ( i=0; i<dimi; i++ )
				fwrite((stck+j+i*dimj+k*dimj*dimi),sizeof(iom::real_t),1,fout);

	fclose(fout);

	//char *err_tiff_fmt;
	//int img_depth = 16;

	//unsigned char *buffer = new unsigned char [dimi*dimj*dimk*(img_depth/8)];

	//for ( k=0; k<dimk; k++ )
	//	for ( j=0; j<dimj; j++ )
	//		for ( i=0; i<dimi; i++ ) {
	//			if ( img_depth == 8 )
	//				*(buffer + j+i*dimj+k*dimj*dimi) = static_cast<iim::uint8>(*(stck+j+i*dimj+k*dimj*dimi) * 255.0f + 0.5f);
	//			else // img_depth = 16
	//				*(((iim::uint16 *)buffer) + j+i*dimj+k*dimj*dimi) = static_cast<iim::uint16>(*(stck+j+i*dimj+k*dimj*dimi) * 65535.0f + 0.5f);
	//		}

	//// creates the file (2D image: depth is 1)
	//if ( (err_tiff_fmt = initTiff3DFile(fname,dimj,dimi,1,1,img_depth/8)) != 0 ) {
	//	throw iom::exception(iom::strprintf("unable to create tiff file (%s)",err_tiff_fmt), __iom__current__function__);
	//}

	//if ( (err_tiff_fmt = appendSlice2Tiff3DFile(fname,0,buffer,dimj,dimi)) != 0 ) {
	//	throw iom::exception(iom::strprintf("error in saving 2D image (%lld x %lld) in file %s (appendSlice2Tiff3DFile: %s)",dimj,dimi,fname,err_tiff_fmt), __iom__current__function__);
	//}

	//delete buffer;

	return true;
}

#endif

