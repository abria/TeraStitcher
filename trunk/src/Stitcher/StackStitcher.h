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
*       Bria, A., et al., (2012) "Stitching Terabyte-sized 3D Images Acquired in Confocal Ultramicroscopy", Proceedings of the 9th IEEE International Symposium on Biomedical Imaging.
*       Bria, A., Iannello, G., "TeraStitcher - A Tool for Fast 3D Automatic Stitching of Teravoxel-sized Microscopy Images", submitted for publication, 2012.
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

#ifndef STACK_STITCHER_H
#define STACK_STITCHER_H

#include <math.h>
#include "S_config.h"
#include "IOManager_defs.h"
#include "MyException.h"


class StackRestorer;
#ifndef STACKED_VOLUME_H
class StackedVolume;
#endif
#ifndef STACK_H_
class Stack;
#endif

class StackStitcher
{
	private:

		/******OBJECT MEMBERS******/
                StackedVolume *volume;					//pointer to the <StackedVolume> object to be stitched
                int V0, V1, H0, H1, D0, D1;				//voxel intervals that identify the final stitched volume
                int ROW_START, COL_START, ROW_END, COL_END;             //stack indexes that identify the stacks involved in stitching

		/******CLASS MEMBERS******/
		static double time_displ_comp;				//time employed for pairwise displacements computation
                static double time_merging;				//time employed to merge stacks
		static double time_stack_desc;				//time employed to compute stacks descriptions
		static double time_stack_restore;			//time employed to restore stacks
		static double time_multiresolution;			//time employed to obtain stitched volume at different resolutions


		/***OBJECT PRIVATE METHODS****/

		//default constructor will not be accessible
		StackStitcher(void){};


		/*************************************************************************************************************
		* Merges all slices of the given row at the given depth index, so obtaining the stripe that is returned.
		* Uses [...]_blending() functions to blend pixels in  overlapping zones.  The appropriate blending function is
		* selected by the [blending_algo] parameter. If a  <StackRestorer>  object has been passed,  each slice is re-
		* stored before it is combined into the final stripe.
		**************************************************************************************************************/
		real_t* getStripe(short row_index, short d_index, int restore_direction=-1, StackRestorer* stk_rst=NULL,
						  int blending_algo=S_SINUSOIDAL_BLENDING)    							   throw (MyException);

		/*************************************************************************************************************
		* Returns the (up = true -> TOP, up = false -> BOTTOM) V coordinate of the virtual stripe at <row_index> row. 
		**************************************************************************************************************/
		int getStripeABS_V(int row_index, bool up);



		/***CLASS PRIVATE METHODS****/

		/*************************************************************************************************************
		* Blending functions that returns the value at <angle> blending <pixel1> and <pixel2>.
		* <angle> is in [0, PI] where 0 occurs at the first pixel of the overlapping zone and PI at last pixel of  the 
		* overlapping zone. IMPORTANT: due to efficiency reasons, it is better to handle  different  types of blending 
		* functions with function pointers instead of using polymorphism of OOP.
		**************************************************************************************************************/
		static inline real_t sinusoidal_blending(double& angle, real_t& pixel1, real_t& pixel2){
			return (real_t)(  ((cos(angle)+1.0F)*0.5F)*pixel1 + ( 1.0F - ((cos(angle)+1.0F)*0.5F))*pixel2  );
		}

		static inline real_t no_blending(double& angle, real_t& pixel1, real_t& pixel2){
			return (angle <= S_PI/2 ? pixel1 : pixel2);		
		}

                /*************************************************************************************************************
                * This is a special blending function (together with the necessary static variables) which draws blank lines
                * along stacks borders without performing any blending. This enables easy checking of motorized stages coordi-
                * nates precision.
                **************************************************************************************************************/
                static double stack_marging_old_val;
                static bool blank_line_drawn;
                static inline real_t stack_margin(double& angle, real_t& pixel1, real_t& pixel2)
                {
                    if(angle == stack_marging_old_val && blank_line_drawn){
                        stack_marging_old_val = angle;
                        return 1;
                    }
                    else if(angle > S_PI/2.0 && stack_marging_old_val < S_PI/2.0){
                        stack_marging_old_val = angle;
                        blank_line_drawn = true;
                        return 1;
                    }
                    else{
                        stack_marging_old_val = angle;
                        blank_line_drawn = false;
                        return ( angle <= S_PI/2 ? pixel1 : pixel2);
                    }
		}
		
		/*************************************************************************************************************
		* Performs downsampling at a halved frequency on the given 3D image.  The given image is overwritten in order
		* to store its halvesampled version without allocating any additional resources.
		**************************************************************************************************************/
		static void halveSample(real_t* img, int height, int width, int depth);

	public:

                StackStitcher(StackedVolume* _volume);

		/*************************************************************************************************************
		* Method to be called for displacement computation. <> parameters are mandatory, while [] are optional.
		* <algorithm_type>		: ID of the pairwise displacement algorithm to be used.
		* [start/end_...]		: rows/columns intervals that possible identify the portion of volume to be processed.
		*						  If not given, all stacks will be processed.
		* [overlap_...]			: expected overlaps between the given stacks along V and H directions.These values can
		*   					  be used to determine the region of interest where the overlapping occurs. If not gi-
		*   					  ven,  default  values  are  assigned  by  computing the expected  overlaps using the 
		*   					  <MEC_...> members of the <StackedVolume> object.
		* [displ_max_...]		: maximum displacements along VHD between two  adjacent stacks  taking the given over-
		*						  lap as reference. These parameters, together with <overlap_...>,can be used to iden-
		*						  tify the region of interest where the correspondence between the given stacks has to
		*						  be found. When used, these parameters have to be tuned with respect to the precision 
		*						  of the motorized stages. If not given, value S_DISPL_SEARCH_RADIUS_DEF is assigned.
		* [subvol_DIM_D]		: desired subvolumes dimensions along D axis.  Each pair  of stacks is split into sub-
		*						  volumes along D axis in order to use memory efficiently.   Hence, multiple displace-
		*						  ments for each pair of adjacent stacks are computed. 
		*						  If not given, value S_SUBVOL_DIM_D_DEFAULT is assigned.
		* [restoreSPIM]			: enables SPIM artifacts removal (zebrated patterns) along the given direction.
		* [restore_direction]	: direction of SPIM zebrated patterns to be removed.
		* [show_progress_bar]	: enables/disables progress bar with estimated time remaining.
		**************************************************************************************************************/
		void computeDisplacements(int algorithm_type, int start_row = -1, int start_col = -1, int end_row = -1, 
								  int end_col = -1, int overlap_V = -1, int overlap_H =-1, 
								  int displ_max_V=S_DISPL_SEARCH_RADIUS_DEF, int displ_max_H=S_DISPL_SEARCH_RADIUS_DEF, 
								  int displ_max_D=S_DISPL_SEARCH_RADIUS_DEF, int subvol_DIM_D = S_SUBVOL_DIM_D_DEFAULT, 
								  bool restoreSPIM=false, int restore_direction=-1, bool show_progress_bar=true)
																								   throw (MyException);


		/*************************************************************************************************************
		* For each stack, the vector of redundant displacements along D is projected into the displacement which embe-
		* ds the most reliable parameters. After this operation, such vector will contain only the projected displace-
                * ment. Where for a pair of adjacent stacks no displacement is available,  a displacement  is generated using
                * nominal stage coordinates.
		**************************************************************************************************************/
		void projectDisplacements()																  throw (MyException);

		/*************************************************************************************************************
		* Assuming that for each pair of adjacent stacks  exists one  and only one displacement,  this displacement is 
		* thresholded according to the given <reliability_threshold>. When a displacement is not reliable enough,  its
		* parameters are set to default values (i.e. nominal motorized stage coordinates).
		* Moreover, stacks which do not have any reliable single-direction displacements with all 4 neighbors are mar-
		* ked as NON STITCHABLE.
		**************************************************************************************************************/
		void thresholdDisplacements(float reliability_threshold)								  throw (MyException);


		/*************************************************************************************************************
		* Executes the compute tiles placement algorithm associated to the given ID <algorithm_type>
		**************************************************************************************************************/
		void computeTilesPlacement(int algorithm_type)											  throw (MyException);


                /*************************************************************************************************************
                * Computes final stitched volume dimensions assuming that current <StackedVolume> object contains  the correct
                * stack coordinates. The given parameters identify the possible VOI (Volume Of Interest). If these are not us-
                * ed, the whole volume is  considered and the parameter <exclude_nonstitchable_stacks> is used to discard rows
                * or columns with no stitchable stacks
                **************************************************************************************************************/
                void computeVolumeDims(bool exclude_nonstitchable_stacks = true, int _ROW_START = -1,	   int _ROW_END = -1,
                                                           int _COL_START = -1, int _COL_END = -1, int _D0 = -1, int _D1 = -1) throw (MyException);

		/*************************************************************************************************************
		* Method to be called for tile merging. <> parameters are mandatory, while [] are optional.
		* <output_path>			: absolute directory path where merged tiles have to be stored.
		* [slice_height/width]	: desired dimensions of tiles  slices after merging.  It is actually an upper-bound of
		*						  the actual slice dimensions, which will be computed in such a way that all tiles di-
		*						  mensions can differ by 1 pixel only along both directions. If not given, the maximum
		*						  allowed dimensions will be set, which will result in a volume composed by  one large 
		*						  tile only.
		* [resolutions]			: pointer to an array of S_MAX_MULTIRES  size which boolean entries identify the acti-
		*						  vaction/deactivation of the i-th resolution.  If not given, all resolutions will  be
		*						  activated.
		* [exclude_nonstitc...] 
		* [_...START/END]		
		* [_D0/_D1]				: identify the possible VOI (Volume Of Interest). If these are not used, the whole vo-
		*						  lume is  considered and the parameter <exclude_nonstitchable_stacks> is used to dis-
		*						  card rows or columns with no stitchable stacks.
		* [restoreSPIM]			: enables SPIM artifacts removal (zebrated patterns) along the given direction.
		* [restore_direction]	: direction of SPIM zebrated patterns to be removed.
		* [blending_algo]		: ID of the blending algorithm to be used in the overlapping regions.		 
		* [test_mode]			: if enabled, the middle slice of the whole volume will be stitched and and  saved lo-
		*						  cally. Stage coordinates will be used, s o this can be used to test  their precision
		*						  as well as the selected reference system.
		* [show_progress_bar]	: enables/disables progress bar with estimated time remaining.
		* [saved_img_format]	: determines saved images format ("png","tif","jpeg", etc.).
		* [saved_img_depth]		: determines saved images bitdepth (16 or 8).
		**************************************************************************************************************/
		void mergeTiles(std::string output_path, int slice_height = -1, int slice_width = -1, bool* resolutions = NULL, 
						bool exclude_nonstitchable_stacks =true, int _ROW_START=-1, int _ROW_END=-1, int _COL_START=-1,
						int _COL_END=-1, int _D0=-1, int _D1=-1,	bool restoreSPIM=false,	  int restore_direction=-1,
						int blending_algo=S_SINUSOIDAL_BLENDING,	bool test_mode=false, bool show_progress_bar= true,
						const char* saved_img_format=IO_DEF_IMG_FORMAT, int saved_img_depth=IO_DEF_IMG_DEPTH) throw (MyException);
		

		/*************************************************************************************************************
		* Functions used to save single phase time performances
		**************************************************************************************************************/
		static void saveComputationTimes(const char *filename, StackedVolume &stk_org, double total_time=-1);
		static void resetComputationTimes();

                /*************************************************************************************************************
                * Get methods
                **************************************************************************************************************/
                int getV0(){return V0;}
                int getV1(){return V1;}
                int getH0(){return H0;}
                int getH1(){return H1;}
                int getD0(){return D0;}
                int getD1(){return D1;}
                int getROW0(){return ROW_START;}
                int getROW1(){return ROW_END;}
                int getCOL0(){return COL_START;}
                int getCOL1(){return COL_END;}

                /*************************************************************************************************************
                * Functions used to obtain absolute coordinates at different resolutions from relative coordinates
                **************************************************************************************************************/
                int getMultiresABS_V(int res, int REL_V);
                std::string getMultiresABS_V_string(int res, int REL_V);
                int getMultiresABS_H(int res, int REL_H);
                std::string getMultiresABS_H_string(int res, int REL_H);
                int getMultiresABS_D(int res, int REL_D);
                std::string getMultiresABS_D_string(int res, int REL_D);
};

#endif


