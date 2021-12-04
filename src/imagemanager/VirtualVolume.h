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
* 2018-06-30. Giulio.     @ADDED attribute 'depth_conv_algo' for specifying the ID of conversion algorithm to be used to convert from arbitrary depth to 8 bits
* 2017-10-21. Giulio.     @ADDED method to compact active channels when not all channels are active (n_active < DIM_C)
* 2016-04-27. Giulio.     @ADDED method to rearrange indices in order to meet pre-conditions of loadSubVolume methods
* 2016-04-27. Giulio.     @ADDED data members and methods to convert indices to coordinates and viceversa
* 2015-04-14. Alessandro. @ADDED 'instance_format' method with inputs = {path, format}.
* 2015-02-28. Giulio.     @FIXED added deallocation of data member 'active' in the destructor
* 2015-02-18. Giulio.     @CHANGED modified defalut values of parameters of loadSubvolume methods
* 2015-01-06. Giulio.     @ADDED changed interface of saveImage_from_UINT8_to_Tiff3D to introduce optimizations to reduce opend/close in append operations 
*/

#ifndef _IIM_VIRTUAL_VOLUME_H
#define _IIM_VIRTUAL_VOLUME_H

# define HALVE_BY_MEAN 1
# define HALVE_BY_MAX  2

#include <string>
#include <algorithm>

#include "IM_config.h"
#include "iomanager.config.h"

class iim::VirtualVolume {

protected:
	//******OBJECT ATTRIBUTES******
    char*  root_dir;				// C-string that contains the directory path of stacks matrix
    float  VXL_V, VXL_H, VXL_D;		// [microns]: voxel dimensions (in microns) along V(Vertical), H(horizontal) and D(Depth) axes
    float  ORG_V, ORG_H, ORG_D;		// [millimeters]: origin spatial coordinates (in millimeters) along VHD axes
    iim::uint32 DIM_V, DIM_H, DIM_D;// volume dimensions (in voxels) along VHD axes
    int    DIM_C;					// number of channels        (@ADDED by Iannello   on ..........)
    int    BYTESxCHAN;              // number of bytes per channel
    iim::uint32   *active;          // array  of active channels (@MOVED from "TiledMCVolume" by Alessandro on 2014-02-20)
    int           n_active;         // number of active channels (@MOVED from "TiledMCVolume" by Alessandro on 2014-02-20)

    int DIM_T;					    // number of time frames         (@ADDED by Alessandro on 2014-02-20)
    int t0, t1;                     // active frames are in [t0, t1] (@ADDED by Alessandro on 2014-02-20)
    
    int depth_conv_algo;             // ID of the algorithm used to convert from 16 to 8 bits the image color depth
    
    // this method should be used to rearrange coordinates in order to match the pre-conditions of loadSubVolume methods
    void rearrange_indices ( int &i0, int &i1 ) { 
    	if ( i0 < i1 )
    		;
    	else if ( i0 > i1 ) {
    		int i = i0; i0 = i1; i1 = i;
    	}
    	else // if indices are equal there is a potential error 
			throw IOException(strprintf("indices i0=%d and i1=%d are equal", i0, i1), __iim__current__function__);
    }

	// this method compacts the active channels starting from a buffer in which are stored all channels
	// assume n_active < DIM_C, i.e. that not all channels are active, which imply that in list 'active' some valid channel index is missing
	void compact_active_chans ( sint64 totalChanSize, uint8 *subvol ) {
		int c;          // active channels index 
		iim::uint32 cc; // all channels index
		uint8 *src;     // pointer to source buffer
		uint8 *dst;     // pointer to destination buffer

		for ( c=0, cc=0; c<n_active && active[c]==cc; c++, cc++ );
		src = dst = subvol + c * totalChanSize;
		for ( ; c<n_active; c++, cc++, src+=totalChanSize, dst+=totalChanSize ) {
			while ( active[c]!=cc ) {
				src += totalChanSize;
				cc++;
			}
			memcpy(dst,src,(size_t)totalChanSize);
		}
	}

public:
	//CONSTRUCTORS-DESTRUCTOR
    VirtualVolume();
    VirtualVolume(const char* _root_dir, float VXL_1=0, float VXL_2=0, float VXL_3=0)
    {

		this->root_dir = new char[strlen(_root_dir)+1];
		strcpy(this->root_dir, _root_dir);

		VXL_V = VXL_1;
		VXL_H = VXL_2;
		VXL_D = VXL_3;

		ORG_V = ORG_H = ORG_D = (float) 0.0;
		DIM_V = DIM_H = DIM_D = 0;

        DIM_C = 0;
		BYTESxCHAN = 0;
        active = (iim::uint32 *)0;
        n_active = 0;

        t0 = t1 = 0;
        DIM_T = 1;

    	depth_conv_algo = DEPTH_CONVERSION_DEFAULT;
    }

    virtual ~VirtualVolume()  {
		if(root_dir)
			delete[] root_dir;
		if(active)
			delete []active;
	}    

    virtual void initChannels ( )  = 0;


    //loads given subvolume in a 1-D array of iim::real32 while releasing stacks slices memory when they are no longer needed
    virtual iim::real32 *loadSubvolume_to_real32(int V0=-1,int V1=-1, int H0=-1, int H1=-1, int D0=-1, int D1=-1)   = 0;


    /* @ADDED by Giulio on 2015-12-06:
     * loads given subvolume in a 1-D array of iim::uint8 
     * data is returned in internal Vaa3D representation, i.e. voxels are stored as one 3D matrix per channel in C-like ordering
     * (dimension order is: H,V,D,channel)
     * 
     * V0-D1:    indices that define the subvolume; index intervals are open at right (e.g. [V0,V1))
     * channels: returns the number of channels loaded
     * ret_type: specifies the number of bits per voxel in each channel; the default is iim::DEF_IMG_DEPTH (see IM_config.h for definition)
     *           if the value of this parameter is iim::NATIVE_RTYPE (see IM_config.h for definition) the subvolume is stored in the 
     *           returned buffer using the number of bits per voxel in each channel of the current volume
     *           currently only iim::DEF_IMG_DEPTH and iim::NATIVE_RTYPE are permitted values for this parameter
     *
     * WARNING: the method returns a newly allocated buffer whose ownership is transferred to the caller which must provide to its
     *          deallocation; the size of the returned buffer is (V1-V0)*(H1-H0)*(D1-D0)*(ret_type/8)*n_active, where n_active is the 
     *          number of channels that have been previously set to be returned by method 'setActiveChannels'
     */
    virtual iim::uint8 *loadSubvolume_to_UINT8(int V0=-1,int V1=-1, int H0=-1, int H1=-1, int D0=-1, int D1=-1,
                                               int *channels=0, int ret_type=iim::DEF_IMG_DEPTH)  = 0;

    // ******GET METHODS******
    float   getORG_V() {return ORG_V;}
    float   getORG_H() {return ORG_H;}
    float   getORG_D() {return ORG_D;}
    int     getDIM_V() {return DIM_V;}
    int     getDIM_H() {return DIM_H;}
    int     getDIM_D() {return DIM_D;}
    int     getDIM(iim::axis dir)
    {
        if(dir == iim::vertical || dir == iim::inv_vertical)
            return DIM_V;
        else if(dir == iim::horizontal || dir == iim::inv_horizontal)
            return DIM_H;
        else if(dir == iim::depth || dir == iim::inv_depth)
            return DIM_D;
        else
            throw iim::IOException("VirtualVolume::getDIM(): axis invalid");
    }
    int     getABS_V(int ABS_PIXEL_V) {return iim::round( ORG_V * 1000 + ABS_PIXEL_V*this->getVXL_V());}
    int     getABS_H(int ABS_PIXEL_H) {return iim::round( ORG_H * 1000 + ABS_PIXEL_H*this->getVXL_H());}
    int     getABS_D(int ABS_PIXEL_D) {return iim::round( ORG_D * 1000 + ABS_PIXEL_D*this->getVXL_D());}
    float   getVXL_V() {return VXL_V;}
    float   getVXL_H() {return VXL_H;}
    float   getVXL_D() {return VXL_D;}
    int     getDIM_C() {return DIM_C;}
    int     getDIM_T() {return DIM_T;}  //@ADDED by Alessandro on 2014-02-18
    int     getT0()    {return t0;}
    int     getT1()    {return t1;}
    int     getBYTESxCHAN() {return BYTESxCHAN;}
    char*   getROOT_DIR() {return this->root_dir;}
    virtual float   getMVoxels(){return (DIM_V/1024.0f)*(DIM_H/1024.0f)*DIM_D*DIM_T;} // can be overriden
    int getNActiveFrames(){return t1 -t0 +1;}
    virtual int getNACtiveChannels() {return n_active;}
    virtual iim::uint32* getActiveChannels(){return active;}
    int     getDEPTH_CONV_ALGO(){return depth_conv_algo;}
    
    // @ADDED by Giulio. on 2016-04-27: methods to convert indices to coordinates and viceversa
    float ind2coord_V(int v) { return ORG_V + v*VXL_V; } 
    float ind2coord_H(int h) { return ORG_H + h*VXL_H; } 
    float ind2coord_D(int d) { return ORG_D + d*VXL_D; } 
    int coord2ind_V(float v) { return (int) ((v - ORG_V) / VXL_V); }    
    int coord2ind_H(float h) { return (int) ((h - ORG_H) / VXL_H); }    
    int coord2ind_D(float d) { return (int) ((d - ORG_D) / VXL_D); }   

    // @ADDED by Alessandro on 2014-02-18: returns a unique ID that identifies the volume format
    virtual std::string getPrintableFormat() = 0;

// @ADDED by Alessandro on 2014-02-18: additional info on the reference system (where available)
    virtual float getVXL_1() = 0;
    virtual float getVXL_2() = 0;
    virtual float getVXL_3() = 0;
    virtual iim::axis getAXS_1() = 0;
    virtual iim::axis getAXS_2() = 0;
    virtual iim::axis getAXS_3() = 0;

    // set active channels for 4D data (@MOVED from TileMCVolume.h by Alessandro on 2014-02-20)
	// WARNING: caller loses ownership of array '_active' 
    virtual void setActiveChannels ( iim::uint32 *_active, int _n_active );
    
    // set the remap or conversion algorithm used for visualization of subregions
    void setDEPTH_CONV_ALGO(int algoID);

    // set active frame for 5D data (@MOVED from TimeSeries.h by Alessandro on 2014-02-20)
    void setActiveFrames(int _t0, int _t1)
    {
        t0 = std::max(0, std::min(_t0,DIM_T-1));
        t1 = std::max(0, std::min(_t1,DIM_T-1));
        iim::debug(iim::LEV_MAX, iim::strprintf("asked to set [%d, %d], but set [%d, %d]", _t0, _t1, t0, t1).c_str(), __iim__current__function__);
    }

    // @ADDED by Alessandro on 2016-12-19
    // return true if the given dimension is tiled
    virtual bool isTiled(iim::dimension d) = 0;
    // return vector of tiles along x-y-z (empty vector if the volume is not tiled)
    virtual std::vector< iim::voi3D<size_t> > tilesXYZ() = 0;

	/*************************************************************************************************************
    * Save image method. <> parameters are mandatory, while [] are optional.
    * <img_path>                : absolute path of image to be saved. It DOES NOT include its extension, which is
    *                             provided by the [img_format] parameter.
    * <raw_img>                 : image to be saved. Raw data is in [0,1] and it is stored row-wise in a 1D array.
    * <raw_img_height/width>    : dimensions of raw_img.
    * [start/end_height/width]  : optional ROI (region of interest) to be set on the given image.
    * [img_format]              : image format extension to be used (e.g. "tif", "png", etc.)
    * [img_depth]               : image bitdepth to be used (8 or 16)
	**************************************************************************************************************/
    static void saveImage(std::string img_path,   iim::real32* raw_img,       int raw_img_height,   int   raw_img_width,
                              int start_height = 0,   int end_height = - 1,  int start_width = 0,  int end_width = - 1,
                              const char* img_format = iim::DEF_IMG_FORMAT.c_str(),    int img_depth = iim::DEF_IMG_DEPTH		 )
                                                                                                   ;

    /*************************************************************************************************************
    * Save image method from iim::uint8 raw data. <> parameters are mandatory, while [] are optional.
    * <img_path>                : absolute path of image to be saved. It DOES NOT include its extension, which is
    *                             provided by the [img_format] parameter.
    * <raw_ch1>                 : raw data of the first channel with values in [0,255].
    *                             For grayscale images this is the pointer to the raw image data.
    *                             For colour images this is the pointer to the raw image data of the RED channel.
    * <raw_ch2>                 : raw data of the second channel with values in [0,255].
    *                             For grayscale images this should be a null pointer.
    *                             For colour images this is the pointer to the raw image data of the GREEN channel.
    * <raw_ch3>                 : raw data of the second channel with values in [0,255].
    *                             For grayscale images this should be a null pointer.
    *                             For colour images this is the pointer to the raw image data of the BLUE channel.
    * <raw_img_height/width>    : dimensions of raw_img.        
    * [start/end_height/width]  : optional ROI (region of interest) to be set on the given image.
    * [img_format]              : image format extension to be used (e.g. "tif", "png", etc.)
    * [img_depth]               : image bitdepth to be used (8 or 16)
	*
	* WARNING: this method is intended to be used to save one slice into a file using a 2D format (i.e. a 2D plugin)
	* use saveImage_from_UINT8_to_Tiff3D to append one slice to a file storing images in 3D format
    **************************************************************************************************************/
    static void saveImage_from_UINT8 (std::string img_path, 
                                      iim::uint8* raw_ch1, iim::uint8* raw_ch2, iim::uint8* raw_ch3,
									  int raw_img_height, int raw_img_width,
                                      int start_height=0, int end_height =-1, int start_width=0, int end_width=-1,
                                      const char* img_format = iim::DEF_IMG_FORMAT.c_str(), int img_depth = iim::DEF_IMG_DEPTH ) ;


 	/*************************************************************************************************************
    * Save image method to Vaa3D raw format. <> parameters are mandatory, while [] are optional.
    * <img_path>                : absolute path of image to be saved. It DOES NOT include its extension, which is
    *                             provided by the [img_format] parameter.
    * <raw_img>                 : image to be saved. Raw data is in [0,1] and it is stored row-wise in a 1D array.
    * <raw_img_height/width>    : dimensions of raw_img.
    * [start/end_height/width]  : optional ROI (region of interest) to be set on the given image.
    * [img_format]              : image format extension to be used (e.g. "tif", "png", etc.)
    * [img_depth]               : image bitdepth to be used (8 or 16)
	**************************************************************************************************************/
    static void saveImage_to_Vaa3DRaw(int slice, std::string img_path, iim::real32* raw_img, int raw_img_height, int raw_img_width,
                              int start_height = 0, int end_height = - 1, int start_width = 0, int end_width = - 1,
                              const char* img_format = iim::DEF_IMG_FORMAT.c_str(), int img_depth = iim::DEF_IMG_DEPTH
							  )
                                                                                                   ;

   /*************************************************************************************************************
    * Save image method from iim::uint8 raw data to Vaa3D raw format. <> parameters are mandatory, while [] are optional.
    * <img_path>                : absolute path of image to be saved. It DOES NOT include its extension, which is
    *                             provided by the [img_format] parameter.
    * <raw_ch>                  : array of pointers to raw data of the channels with values in [0,255].
    *                             For grayscale images raw_ch[0] is the pointer to the raw image data.
    *                             For colour images raw_ch[0] is the pointer to the raw image data of the RED channel.
    * <n_chans>                 : number of channels (length of raw_ch).
	* <offset>                  : offset to be added to raw_ch[i] to get actual data
    * <raw_img_height/width>    : dimensions of raw_img.        
    * [start/end_height/width]  : optional ROI (region of interest) to be set on the given image.
    * [img_format]              : image format extension to be used (e.g. "tif", "png", etc.)
    * [img_depth]               : image bitdepth to be used (8 or 16)
    **************************************************************************************************************/
    static void saveImage_from_UINT8_to_Vaa3DRaw (int slice, std::string img_path, 
                                      iim::uint8** raw_ch, int n_chans, iim::sint64 offset,
									  int raw_img_height, int raw_img_width,
                                      int start_height=0, int end_height =-1, int start_width=0, int end_width=-1,
                                      const char* img_format = iim::DEF_IMG_FORMAT.c_str(), int img_depth = iim::DEF_IMG_DEPTH ) ;



	/*************************************************************************************************************
    * Save image method from iim::uint8 raw data to Tiff 3D (multipage) format. <> parameters are mandatory, while [] are optional.
    * <img_path>                : absolute path of image to be saved. It DOES NOT include its extension, which is
    *                             provided by the [img_format] parameter.
    * <raw_ch>                  : array of pointers to raw data of the channels with values in [0,255].
    *                             For grayscale images raw_ch[0] is the pointer to the raw image data.
    *                             For colour images raw_ch[0] is the pointer to the raw image data of the RED channel.
    * <n_chans>                 : number of channels (length of raw_ch).
	* <offset>                  : offset to be added to raw_ch[i] to get actual data
    * <raw_img_height/width>    : dimensions of raw_img.        
    * [start/end_height/width]  : optional ROI (region of interest) to be set on the given image.
    * [img_format]              : image format extension to be used (e.g. "tif", "png", etc.)
    * [img_depth]               : image bitdepth to be used (8 or 16)
	* [fhandle]                 : handle to the (open) file which image has to be appended (used only if do_open = false)
	* [n_pages]                 : total number of slice to be appended to the file (used only if do_open = false) 
	* [do_open]                 : if true img_path has to be used to open the file which is closed after the image has been appended
	*                             if false the image is appended without opening the file which is left open
    **************************************************************************************************************/
	// WARNING: current implementation could not work if the slice is not appended after the last one
	//          in other words the multipage file should have exactly (slice-1) pages
    static void saveImage_from_UINT8_to_Tiff3D (int slice, std::string img_path, 
                                      iim::uint8** raw_ch, int n_chans, iim::sint64 offset,
									  int raw_img_height, int raw_img_width,
                                      int start_height=0, int end_height =-1, int start_width=0, int end_width=-1,
                                      const char* img_format = iim::DEF_IMG_FORMAT.c_str(), int img_depth = iim::DEF_IMG_DEPTH, 
                                      void *fhandle = 0, int n_pages = -1, bool do_open = true ) ;

	/*************************************************************************************************************
	* Performs downsampling at a halved frequency on the given 3D image.  The given image is overwritten in order
	* to store its halvesampled version without allocating any additional resources.
	*
	* WARNING: Since the downsampling is carried out for more slices, a stride is introduced between downsampled
	* slices. The stride introduced is (height*width)
	**************************************************************************************************************/
    static void halveSample( iim::real32* img, int height, int width, int depth, int method = HALVE_BY_MEAN );

    static void halveSample_UINT8 ( iim::uint8** img, int height, int width, int depth, int channels, int method = HALVE_BY_MEAN, int bytes_chan = 1 );

    static void halveSample2D( iim::real32* img, int height, int width, int depth, int method = HALVE_BY_MEAN );

    static void halveSample2D_UINT8 ( iim::uint8** img, int height, int width, int depth, int channels, int method = HALVE_BY_MEAN, int bytes_chan = 1 );

	//utility function: returns true if "fullString" ends with "ending"
	inline static bool hasEnding (std::string const &fullString, std::string const &ending)
	{
	   if (fullString.length() >= ending.length())
		  return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
	   else
		  return false;
	}

    // tries to automatically detect the volume format and returns the imported volume if succeeds (otherwise returns 0)
    // WARNING: all metadata files (if needed by that format) are assumed to be present. Otherwise, that format will be skipped.
    static VirtualVolume* instance(const char* path) ;

	// should be used to instantiate BDV HDF5 volumes
    static VirtualVolume* instance(const char* fname, int res, void *descr, int timepoint = 0) ;

    // returns the imported volume if succeeds (otherwise returns 0)
    // WARNING: no assumption is made on metadata files, which are possibly (re-)generated using the additional informations provided.
    static VirtualVolume* instance(const char* path, std::string format,
                                   iim::axis AXS_1 = iim::axis_invalid, iim::axis AXS_2 = iim::axis_invalid, iim::axis AXS_3 = iim::axis_invalid,
                                   float VXL_1=0.0f, float VXL_2=0.0f, float VXL_3=0.0f) ;

    // 2015-04-14. Alessandro. @ADDED 'instance_format' method with inputs = {path, format}.
    static VirtualVolume* instance_format(const char* path, std::string format) ;
    
    // checks whether the volume stored in "path" can be imported directly (i.e., w/o additional metadata provided by the user)
    static bool isDirectlyImportable(const char* path)
    {
        /**/iim::debug(iim::LEV3, iim::strprintf("path = \"%s\"", path).c_str(), __iim__current__function__);

        VirtualVolume* vol = 0;
        try{vol = instance(path);}
        catch(iim::IOException &ex){/**/iim::debug(iim::LEV3, iim::strprintf("error = %s", ex.what()).c_str(), __iim__current__function__);}
        catch(...){}
        bool result = vol != 0;
        delete vol;
        return result;
    }

    // returns true if the given format is hierarchical, i.e. if it consists of nested folders (1 level at least)
    static bool isHierarchical(std::string format) ;
};

#endif
