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
* 2019-10-18. Giulio.      @CREATED
*/

//#include <iostream>
//#include <string>
#include "ComposedVolume.h"
//#include "iomanager.config.h"

#ifdef _WIN32
#include "dirent_win.h"
#else
#include <dirent.h>
#endif

#include "ProgressBar.h"
#include "tinyxml.h"


using namespace std;
using namespace iim;

// 2015-04-15. Alessandro. @ADDED definition for default constructor.
ComposedVolume::ComposedVolume(void) : VirtualVolume()
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);   

	reference_system.first  = iim::vertical;
	reference_system.second = iim::horizontal;
	reference_system.third  = iim::depth;
	VXL_1 = VXL_2 = VXL_3 = 0.0f;
}


ComposedVolume::~ComposedVolume(void) 
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);
}


void ComposedVolume::init()
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);
}


iim::real32 *ComposedVolume::loadSubvolume_to_real32(int V0,int V1, int H0, int H1, int D0, int D1)   {

	throw iom::exception(iom::strprintf("not implemented yet"), __iom__current__function__);
}


iim::uint8 *ComposedVolume::loadSubvolume_to_UINT8(int V0,int V1, int H0, int H1, int D0, int D1,
										   int *channels, int ret_type)  {

    if ( (ret_type != iim::NATIVE_RTYPE) && (ret_type != iim::DEF_IMG_DEPTH) ) {
		// return type should be converted, but not to 8 bits per channel
        char err_msg[STATIC_STRINGS_SIZE];
		sprintf(err_msg,"MultiSliceVolume::loadSubvolume_to_UINT8: non supported return type (%d bits) - native type is %d bits",ret_type, 8*this->BYTESxCHAN); 
        throw IOException(err_msg);
	}
	
	// load subvolume in a newly allocated buffer with all offsets set to 0
	iim::uint8 *subvol = loadSubvolume(V0,V1,H0,H1,D0,D1,channels,ret_type,(iim::uint8 *)0,0,0,0,0,0,0,0,0);
	
	return subvol;
}


// return 'volume_format' attribute of <TeraStitcher> node from the given xml. 
std::string ComposedVolume::getVolumeFormat(const std::string& xml_path) 
{
	// open xml
	TiXmlDocument xml;
	if(!xml.LoadFile(xml_path.c_str()))
		throw iom::exception(iom::strprintf("in ComposedVolume::getVolumeFormat(): cannot open xml file at \"%s\"", xml_path.c_str()));

	// get root node
	TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));

	// get 'volume_format' attribute
	const char *volformat = hRoot.ToElement()->Attribute("volume_format");
	if(!volformat)
		throw iom::exception(iom::strprintf("in ComposedVolume::getVolumeFormat(): cannot find 'volume_format' <TeraStitcher> attribute in xml file at \"%s\". Too old xml, please regenerate it.", xml_path.c_str()).c_str());

	return volformat;
}


void ComposedVolume::copy_strided_data ( iim::uint8 *dst, int dstOffs_V, int dstOffs_H, int dstOffs_D, int dstOffs_C, 
    						int dstSz_V, int dstSz_H, int dstSz_D, int dstSz_C,
    						iim::uint8 *src, int srcOffs_V, int srcOffs_H, int srcOffs_D, int srcOffs_C, 
    						int srcSz_V, int srcSz_H, int srcSz_D, int srcSz_C,
							int height, int width, int depth, int chans ) {

	int c, d, v;

	iim::uint64 src_stride_H   = srcSz_H;
	iim::uint64 src_stride_VH  = srcSz_V * src_stride_H;
	iim::uint64 src_stride_DVH = srcSz_D * src_stride_VH;

	iim::uint64 dst_stride_H   = dstSz_H;
	iim::uint64 dst_stride_VH  = dstSz_V * dst_stride_H;
	iim::uint64 dst_stride_DVH = dstSz_D * dst_stride_VH;

	iim::uint8 *srcPtr_V;
	iim::uint8 *srcPtr_D;
	iim::uint8 *srcPtr_C;

	iim::uint8 *dstPtr_V;
	iim::uint8 *dstPtr_D;
	iim::uint8 *dstPtr_C;

	for ( c=0, dstPtr_C=dst+dstOffs_C*dst_stride_DVH, srcPtr_C=src+srcOffs_C*src_stride_DVH; 
		  c<srcSz_C; 
		  c++, dstPtr_C+=dst_stride_DVH, srcPtr_C+=src_stride_DVH ) {

		for ( d=0, dstPtr_D=dstPtr_C+dstOffs_D*dst_stride_VH, srcPtr_D=srcPtr_C+srcOffs_D*src_stride_VH; 
			  d<srcSz_D; 
			  d++, dstPtr_D+=dst_stride_VH, srcPtr_D+=src_stride_VH ) {

			for ( v=0, dstPtr_V=dstPtr_D+dstOffs_V*dst_stride_H+dstOffs_H, srcPtr_V=srcPtr_D+srcOffs_V*src_stride_H+srcOffs_H; 
				  v<srcSz_V; 
				  v++, dstPtr_V+=dst_stride_H, srcPtr_V+=src_stride_H ) {

				memcpy(dstPtr_V,srcPtr_V,width);

			}
		}
	}
}          
