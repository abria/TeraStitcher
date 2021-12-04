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
* 2019-10-18 Giulio.      @CREATED
*/

#ifndef _COMPOSED_VOLUME_H
#define _COMPOSED_VOLUME_H

#include "VirtualVolume.h" // ADDED

//every object of this class has the default (1,2,3) reference system
class iim::ComposedVolume : public iim::VirtualVolume
{
	protected:	
		//******OBJECT ATTRIBUTES******

        iim::ref_sys reference_system;       //reference system of the stored volume
        float  VXL_1, VXL_2, VXL_3;         //voxel dimensions of the stored volume

		//***OBJECT PRIVATE METHODS****

		//Given the reference system, initializes all object's members using stack's directories hierarchy
        void init() ;

		//rotate stacks matrix around D axis (accepted values are theta=0,90,180,270)
		void rotate(int theta);

		//mirror stacks matrix along mrr_axis (accepted values are mrr_axis=1,2,3)
        void mirror(iim::axis mrr_axis);

		//extract spatial coordinates (in millimeters) of given Stack object reading directory and filenames as spatial coordinates
        void extractCoordinates(iim::Block* stk, int z, int* crd_1, int* crd_2, int* crd_3);
        
    	// generalized copy of strided data
    	void copy_strided_data ( iim::uint8 *dst, int dstOffs_V, int dstOffs_H, int dstOffs_D, int dstOffs_C,    // offsets on destination buffer
    								int dstSz_V, int dstSz_H, int dstSz_D, int dstSz_C,                          // size of destination buffer
    								iim::uint8 *src, int srcOffs_V, int srcOffs_H, int srcOffs_D, int srcOffs_C, // offsets on source buffer
    								int srcSz_V, int srcSz_H, int srcSz_D, int srcSz_C,                          // size of source buffer
									int height, int width, int depth, int chans );                               // size of data to be copied

	public:
		//CONSTRUCTORS-DESTRUCTOR
        // objects passed to constructors will not be deallocated by the destructor
		ComposedVolume(void);

        virtual ~ComposedVolume(void) ;

		//GET methods
        float  getVXL_1(){return VXL_1;}
        float  getVXL_2(){return VXL_2;}
        float  getVXL_3(){return VXL_3;}
        iim::axis   getAXS_1(){return reference_system.first;}
        iim::axis   getAXS_2(){return reference_system.second;}
        iim::axis   getAXS_3(){return reference_system.third;}
		iim::ref_sys getREF_SYS(){return reference_system;}

		// return 'volume_format' attribute of <TeraStitcher> node from the given xml. 
        static std::string getVolumeFormat(const std::string& xml_path) ;

        iim::real32 *loadSubvolume_to_real32(int V0=-1,int V1=-1, int H0=-1, int H1=-1, int D0=-1, int D1=-1)  ;

		//loads given subvolume in a 1-D array of iim::uint8 while releasing stacks slices memory when they are no longer needed
		iim::uint8 *loadSubvolume_to_UINT8(int V0=-1,int V1=-1, int H0=-1, int H1=-1, int D0=-1, int D1=-1,
                                               int *channels=0, int ret_type=iim::DEF_IMG_DEPTH) ;

        //loads given subvolume in a 1-D array of iim::uint8 and copy it into 'buffer' starting from  offsets 
        virtual iim::uint8 *loadSubvolume(
        		int V0=-1,int V1=-1, int H0=-1, int H1=-1, int D0=-1, int D1=-1, int *n_chans=0, int ret_type=iim::DEF_IMG_DEPTH,
        		iim::uint8 *buffer=0, int bufSize_V=0, int bufSize_H=0, int bufSize_D=0, int bufSize_C=0,
        		int bufOffs_V=0, int bufOffs_H=0, int bufOffs_D=0, int bufOffs_C=0 
                /*ret_type=iim::DEF_IMG_DEPTH*/
        )  = 0;
        
        // return true if the given dimension is tiled
        virtual bool isTiled(iim::dimension d) {return false;}
        // return vector of tiles along x-y-z (empty vector if the volume is not tiled)
        virtual std::vector< iim::voi3D<size_t> > tilesXYZ() {return std::vector< iim::voi3D<size_t> >();}

	// needed to enable the detection by the factory of volume format through use of the default constructor
        //friend class iim::VirtualVolume; 

};

#endif //_COMPOSED_VOLUME_H
