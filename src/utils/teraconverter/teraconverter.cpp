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

# include <stdio.h> 
# include <stdlib.h>

#include "VolumeConverter.h"
#include <CmdLine.h>
#include "TemplateCLI.h"
#include "iomanager.config.h"

/*
 * This is the driver program for the library class VolumeConverter
 *
 * User has to pass on the command line the following parameters:
 * - directory or file name of the source image (the file name if the image 
 *   is stored in a single file, e.g. in V3D raw format)
 * - directory where to store the output image
 * - the height of the slices of the substacks in the output image
 * - the width of the slices of the substacks in the output image
 * - format of the source image ("Stacked" for Terastitcher stacked format, 
 *   "Simple" for sequence of numbered .tif images in the same directory,
 *   "Raw" for V3D raw 4D format)
 * - format of the output image ("intensity" for real valued pixels in [0,1],
 *   "graylevel" for integer valued pixels in [0,255], "RGB" for pixel represented
 *   according to RGB format)
 * - which resolutions have to be generated; specified resolutions are generated
 *
 * If the source image is multi channel the format of the output image is 
 * automatically set to "RGB"
 *
 * Allowed suffixes for V3D raw 4D format are: .raw .RAW, .v3draw .V3DRAW
 *
 *
 * ****************** NEW FEATURE **********************
 *
 * It is possoble to generate the multiresolution image in different formats
 * Allowed formats:
 * - stacks of 2D tiff (default)
 * - blocks stored in Vaa3D raw format (-f Vaa3DRaw, --outFmt Vaa3DRaw)
 * - one image per channel and single channel blocks stored in Vaa3D raw format 
 *   (-f Vaa3DRawMC, --outFmt Vaa3DRawMC)
 *
 */


int main ( int argc, char *argv[] ) {

	//char err_msg[IM_STATIC_STRINGS_SIZE];
	try {

		//importing command-line arguments from <TeraStitcherCLI> object
		TemplateCLI cli;
		cli.readParams(argc, argv);
		cli.checkParams();
		
		// do what you have to do
		VolumeConverter vc;
		
		vc.setSrcVolume(cli.src_root_dir.c_str(),cli.src_format.c_str(),cli.dst_format.c_str());

		if ( cli.outFmt == "Tiff2DStck" )
			vc.generateTiles(cli.dst_root_dir.c_str(),cli.resolutions,
				cli.slice_height,cli.slice_width,cli.halving_method,
				cli.show_progress_bar,"tif",8*vc.getVolume()->getBYTESxCHAN());
		else if ( cli.outFmt == "Vaa3DRaw" )
			vc.generateTilesVaa3DRaw(cli.dst_root_dir.c_str(),cli.resolutions,
				cli.slice_height,cli.slice_width,cli.slice_depth,cli.halving_method,
				cli.show_progress_bar,"vaa3DRaw",8*vc.getVolume()->getBYTESxCHAN());
		else if ( cli.outFmt == "Tiff3D" )
			vc.generateTilesVaa3DRaw(cli.dst_root_dir.c_str(),cli.resolutions,
				cli.slice_height,cli.slice_width,cli.slice_depth,cli.halving_method,
				cli.show_progress_bar,"Tiff3D",8*vc.getVolume()->getBYTESxCHAN());
		else if ( cli.outFmt == "Vaa3DRawMC" )
			vc.generateTilesVaa3DRawMC(cli.dst_root_dir.c_str(),cli.resolutions,
				cli.slice_height,cli.slice_width,cli.slice_depth,cli.halving_method,
				cli.show_progress_bar,"Vaa3DRaw",8*vc.getVolume()->getBYTESxCHAN());
		else if ( cli.outFmt == "Tiff3DMC" )
			vc.generateTilesVaa3DRawMC(cli.dst_root_dir.c_str(),cli.resolutions,
				cli.slice_height,cli.slice_width,cli.slice_depth,cli.halving_method,
				cli.show_progress_bar,"Tiff3D",8*vc.getVolume()->getBYTESxCHAN());
		else if ( cli.outFmt == "Fiji_HDF5" )
			vc.generateTilesBDV_HDF5(cli.dst_root_dir.c_str(),cli.resolutions,
				cli.slice_height,cli.slice_width,cli.slice_depth,cli.halving_method,
				cli.show_progress_bar,"Tiff3D",8*vc.getVolume()->getBYTESxCHAN());
	}
	catch( iom::exception& exception){
		cout<<"ERROR: "<<exception.what()<<endl<<endl;
	}
	catch( iim::IOException& exception){
		cout<<"ERROR: "<<exception.what()<<endl<<endl;
	}
	catch(bad_exception& be){
		cout<<"GENERIC ERROR: "<<be.what()<<endl<<endl;
	}
	catch(char* error){
		cout<<"GENERIC ERROR: "<<error<<endl<<endl;
	}

	return EXIT_SUCCESS;
}