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

#include <iostream>
#include <stdio.h>
#include <CmdLine.h>
#include "StackedVolume.h"
#include "StackStitcher.h"
#include "TeraStitcherCLI.h"

using namespace std;

int main(int argc, char** argv)
{
	try
	{
		//importing command-line arguments from <TeraStitcherCLI> object
		TeraStitcherCLI cli;
		cli.readParams(argc, argv);
		cli.checkParams();
		//system("PAUSE");

		//importing volume from directory or XML file
		StackedVolume* volume = NULL;
		if(cli.import || (cli.computedisplacements && cli.projfile_load_path.compare("null")==0) || cli.stitch || cli.test)
			volume = new StackedVolume(cli.volume_load_path.c_str(), cli.reference_system, cli.VXL_1, cli.VXL_2, cli.VXL_3);
		if( (cli.computedisplacements && cli.volume_load_path.compare("null")==0) || cli.projdisplacements || cli.thresholddisplacements || cli.placetiles || cli.mergetiles)
			volume = new StackedVolume(cli.projfile_load_path.c_str());

		//processing volume		
		double total_time = -TIME(0);
		StackStitcher* stitcher = new StackStitcher(volume);
		if(cli.computedisplacements || cli.stitch)
			stitcher->computeDisplacements(cli.pd_algo, cli.start_stack_row, cli.start_stack_col, cli.end_stack_row, cli.end_stack_col, 
										   cli.overlap_V, cli.overlap_H, cli.search_radius_V, cli.search_radius_H, cli.search_radius_D, 
										   cli.subvol_dim_D, cli.enable_restore, cli.restoring_direction, cli.show_progress_bar);
		if(cli.projdisplacements || cli.stitch)
			stitcher->projectDisplacements();
		if(cli.thresholddisplacements || cli.stitch)
			stitcher->thresholdDisplacements(cli.reliability_threshold);
		if(cli.placetiles || cli.stitch)
			stitcher->computeTilesPlacement(cli.tp_algo);
		if(cli.mergetiles || cli.stitch)
			stitcher->mergeTiles(cli.volume_save_path, cli.slice_height, cli.slice_width, cli.resolutions, cli.exclude_nonstitchables, 
			                     cli.start_stack_row, cli.end_stack_row, cli.start_stack_col, cli.end_stack_col, cli.D0, cli.D1, 
								 cli.enable_restore, cli.restoring_direction, cli.tm_blending, false, cli.show_progress_bar, cli.img_format.c_str(), cli.img_depth);
		if(cli.test)
			stitcher->mergeTiles(cli.volume_save_path, -1, -1, NULL, false, -1, -1, -1, -1, volume->getN_SLICES()/2, volume->getN_SLICES()/2 +1, 
			                     false, false, S_NO_BLENDING, true, false, cli.img_format.c_str(), cli.img_depth);
		total_time += TIME(0);

		//saving project file and computation times
		if(cli.projfile_save_path.compare("null") != 0)
			volume->saveXML(cli.projfile_save_path.c_str());
		if(cli.save_execution_times)
			stitcher->saveComputationTimes(cli.execution_times_filename.c_str(), *volume, total_time);

		//releasing objects
		if(stitcher)
			delete stitcher;
		if(volume)
			delete volume;
	}
	catch( MyException& exception){
		cout<<"ERROR: "<<exception.what()<<endl<<endl;
	}
	catch(bad_exception& be){
		cout<<"GENERIC ERROR: "<<be.what()<<endl<<endl;
	}
	catch(char* error){
		cout<<"GENERIC ERROR: "<<error<<endl<<endl;
	}
	//system("PAUSE");
	return EXIT_SUCCESS;
}
