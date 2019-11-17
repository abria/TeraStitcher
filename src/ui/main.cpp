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
* 2019-11-17. Giulio.     @ADDED code to launch parastitcher using the script in the share/python directory of the distribution
* 2019-11-07. Giulio.     @ADDED command line option 'fixed_tiling' controlling the strategy used to partition the volume into tiles
* 2019-09-26. Giulio.     @CHANGED the construction of the parastitcher command to make it compatible with compilers that do not support C++11 extensions
* 2019-08-25. Giulio.     @ADDED possibility to launch Parastitcher
* 2018-08-22. Giulio.     @ADDED display of total execution time
* 2018-03-02. Giulio.     @ADDED the possibility to set a path and a name for the mdata.bin file generated when volumes are created from data
* 2018-01-23. Giulio.     @ADDED check to set always the input format of a MultiVolume dataset as a non-interleaved input format
* 2017-12-01. Giulio      @ADDED conditional code to activate the new merge step 
* 2017-06-30. Giulio.     @FIXED the coding of iom::CHANS was not coherent with the one used internally by vmVirtualVolume
* 2017-06-30. Giulio.     @ADDED control over displacement computation of last row and last column of tiles
* 2017-04-02. Giulio.     @ADDED support for creation of BigTiff files
* 2017-04-02. Giulio.     @ADDED setting of the input channel in case the plugin is not interleaved (allow more than 3 channels)
* 2016-09-04. Giulio.     @ADDED the setting of the configuration of the LibTIFF library
* 2015-08-16. Giulio.     @ADDED the 'method' and 'isotropi' parameters to calls to methods merging 3D tiled images
* 2015-06-12. Giulio.     @ADDED calla to check method on the imported/tested volume
* 2015-02-03. Alessandro. @FIXED bug in the call of mergeTiles (wrong D1 if D1 was not set).
* 2014-12-10. Alessandro. @FIXED not working --import step with XML file.
* 2014-12-06. Giulio      @ADDED a case in main to be executed when the makedirs flag is set.
*/

#include <iostream>
#include <stdio.h>
#include <CmdLine.h>
#include "StackStitcher.h"
#include "CLI.h"
#include "vmVirtualVolume.h"
#include "vmStackedVolume.h"
#include "vmBlockVolume.h"
#include "vmMCVolume.h"
#include "IOPluginAPI.h"
#include "ProgressBar.h"

#include "Tiff3DMngr.h"

#ifdef NEW_MERGE_ACTIVE
#include "VolumeConverter.h"
#include "UnstitchedVolume.h"
#endif

using namespace std;

//extracts the filename from the given path and stores it into <filename>
inline std::string getFileName(std::string const & path, bool save_ext = true){

	std::string filename = path;

	// Remove directory if present.
	// Do this before extension removal in case directory has a period character.
	const size_t last_slash_idx = filename.find_last_of("\\/");
	if (std::string::npos != last_slash_idx)
		filename.erase(0, last_slash_idx + 1);

	// Remove extension if present.
	if(!save_ext)
	{
		const size_t period_idx = filename.rfind('.');
		if (std::string::npos != period_idx)
			filename.erase(period_idx);
	}

	return filename;
}

//returns true if the given string <fullString> ends with <ending>
inline bool hasEnding (std::string const &fullString, std::string const &ending){
	if (fullString.length() >= ending.length()) {
		return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
	} else {
		return false;
	}
}

//returns file extension, if any (otherwise returns "")
inline std::string getFileExtension(const std::string& FileName){
	if(FileName.find_last_of(".") != std::string::npos)
		return FileName.substr(FileName.find_last_of(".")+1);
	return "";
}

//returns the absolute path (if present) of the given filepath, i.e. its path without the file name
//e.g. "C:/somedir/somefile.ext" -> "C:/somedir"
inline std::string getAbsolutePath(const std::string& filepath){

	string file_path_string = filepath;
	string file_name;
	string ris;

	//loading file name from 'file_path_string' into 'file_name' by strtoking with "/" OR "\" character
	char * tmp;
	tmp = strtok (&file_path_string[0],"/\\");
	while (tmp)
	{
		file_name = tmp;
		tmp = strtok (0, "/\\");
	}

	//restoring content of file_path_string due to possible alterations done by strtok
	file_path_string = filepath;

	//loading file name substring index in 'file_path_string'
	int index_of_filename= (int) file_path_string.find(file_name);

	//extract substring that contains absolute path but the filename
	ris = file_path_string.substr(0,index_of_filename);
	return ris;
}

inline string fillPath(const std::string & incompletePath, string defPath, string defName, string defExt){
	
	string ext = getFileExtension(incompletePath);
	string name = getFileName(incompletePath, false);
	string path = getAbsolutePath(incompletePath);

	string fullpath;
	if(path.empty())
		fullpath += defPath;
	else
		fullpath += path;

	if(name.empty())
		fullpath += "/" + defName;
	else
		fullpath += "/" + name;

	if(ext.empty())
		fullpath += "." + defExt;
	else
		fullpath += "." + ext;

	return fullpath;
}

int main(int argc, char** argv)
{
	try
	{
		// import command-line arguments from <TeraStitcherCLI> object
		TeraStitcherCLI cli;
		cli.readParams(argc, argv);
		cli.checkParams();
		string defaultOutputFileName = "null";
		
		// 2019-04-25. Giulio. check first id Parastitcher has to be launced
		if ( cli.launch_pyscript ) {

			std::string python_path;
			char *temp;
			if ( (temp = getenv("__PYSCRIPTS_PATH__")) ) {
				python_path = temp;
			}
			else {
				python_path = getAbsolutePath(cli.arg0);
				if ( python_path.find('\\') != std::string::npos )
					python_path += "..\\share\\python";
				else 
					python_path += "../share/python";
				if ( !vm::isDirectory(python_path) )
					throw iom::exception(vm::strprintf("Path of Parastitcher not available: check environment variable __PYSCRIPTS_PATH__").c_str());
				//else
				//	throw iom::exception(vm::strprintf("Path of Parastitcher is: %s",python_path.c_str()).c_str());
			}
			
			if ( python_path[python_path.size()-1] != '/' && python_path[python_path.size()-1] != '\\' ) {
				if ( python_path.find('\\') != std::string::npos )
					python_path += "\\";
				else 
					python_path += "/";
			}

			std::string python_cmd;
			// 2019-09-26. Giulio. @CHANGED to make it compatible with compilers that do not support C++11 extensions
			char numprocs[5];
			sprintf(numprocs,"%d",cli.num_procs);
			std::string np = (const char *)numprocs;
			//python_cmd = "mpirun -np " + std::to_string(cli.num_procs) + " python " + python_path + "Parastitcher.py";
			python_cmd = "mpirun -np " + np + " python " + python_path + "Parastitcher.py";
			for ( int opt=1; opt<argc; opt++ ) {
				if ( strstr(argv[opt],"--parastitcher") == 0 )
					python_cmd += std::string(" ") + argv[opt];
			}
			//printf("---> %s\n",python_cmd.c_str());
			if ( system(python_cmd.c_str()) != 0 ) 
				throw iom::exception(iom::strprintf("in main(): cannot execute python command \"%s\"", python_cmd.c_str()));
			
			return EXIT_SUCCESS;
		}

		// start timer
		double proctime = -TIME(0);

		setLibTIFFcfg(!cli.libtiff_uncompressed,cli.libtiff_bigtiff,cli.libtiff_rowsPerStrip);

		// redirect progress bar output to command line
		terastitcher::ProgressBar::instance()->setToGUI(false);

		// perform simple tasks
		if(cli.pluginsinfo)
		{
			std::cout << iom::IOPluginFactory::pluginsInfo();
			return EXIT_SUCCESS;
		}
		if(cli.dumpMData)
		{
			vm::StackedVolume::dumpMData(cli.volume_load_path.c_str());
			return EXIT_SUCCESS;
		}

		// import volume from directory or XML file
		volumemanager::VirtualVolume* volume = 0;
		// 2014-12-10. Alessandro. @FIXED not working --import step with XML file.
		if((cli.test || cli.import) && cli.projfile_load_path.compare("null")!=0) {
            volume = volumemanager::VirtualVolumeFactory::createFromXML(cli.projfile_load_path.c_str(), cli.rescanFiles || !vm::IMG_FILTER_REGEX.empty());
			if ( !volume->check() )
				throw iom::exception(vm::strprintf("Volume \"%s\" is incomplete or not coherent", volume->getSTACKS_DIR()).c_str());
		}
		else if(cli.import || (cli.computedisplacements && cli.projfile_load_path.compare("null")==0) || cli.stitch || cli.test) {
		volume = volumemanager::VirtualVolumeFactory::createFromData(vm::VOLUME_INPUT_FORMAT_PLUGIN, cli.volume_load_path.c_str(), cli.reference_system, cli.VXL_1, cli.VXL_2, cli.VXL_3, true, cli.mdata_fname);
			if ( !volume->check(fillPath(cli.errlogfile_path, volume->getSTACKS_DIR(), "err_log_file", "txt").c_str()) ) {
				// save xml for further checking
				volume->saveXML(0, fillPath(cli.projfile_save_path, volume->getSTACKS_DIR(), "xml_import", "xml").c_str());
				throw iom::exception(vm::strprintf("Volume \"%s\" is incomplete or not coherent: look at error log file \"%s\"",
					volume->getSTACKS_DIR(),fillPath(cli.errlogfile_path, volume->getSTACKS_DIR(), "err_log_file", "txt").c_str()).c_str());
			}
			if ( vm::VOLUME_INPUT_FORMAT_PLUGIN.compare(vm::MCVolume::id) == 0 ) { // it is a Multivolume: set the enabled channel
				// channels are not interleaved 'active' channel has to be set after checking that the channel exists
				if ( iom::CHANS == iom::NONE ) { // a channel has been specified: the enabled subvolume should be reset
					if ( iom::CHANS_no < volume->getDIM_C() ) {
						((vm::MCVolume *) volume)->setENABLEDSV(iom::CHANS_no);
						if ( !((vm::MCVolume *) volume)->getALIGNED() )
							((vm::MCVolume *) volume)->resetTilePositions();
					}
					else
						throw iom::exception(vm::strprintf("Channel %d does not exist", iom::CHANS_no).c_str());
				}
			}
		}
		else if( (cli.computedisplacements && cli.volume_load_path.compare("null")==0) || cli.projdisplacements || cli.thresholddisplacements || cli.placetiles || cli.mergetiles || cli.makeDirs || cli.metaData) {
            volume = volumemanager::VirtualVolumeFactory::createFromXML(cli.projfile_load_path.c_str(), cli.rescanFiles || !vm::IMG_FILTER_REGEX.empty());
			if ( !volume->check() )
				throw iom::exception(vm::strprintf("Volume \"%s\" is incomplete or not coherent", volume->getSTACKS_DIR()).c_str());
		}

		// set input channel if image data has to be read
		if ( cli.test || cli.import || cli.computedisplacements || cli.mergetiles || cli.stitch ) {
			// check if the input plugin support more than three channels
			bool flag = false;
			std::string plugin_type;
			try{
				plugin_type = "3D image-based I/O plugin";
				flag = iom::IOPluginFactory::getPlugin3D(iom::IMIN_PLUGIN)->desc().find(plugin_type) != std::string::npos;
			}
			catch (...) {
				plugin_type = "2D image-based I/O plugin";
				flag = iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->desc().find(plugin_type) != std::string::npos;
			}
			if ( !flag )
 				throw iom::exception(iomanager::strprintf("cannot determine the type of the input plugin"), __iom__current__function__);
			if ( (plugin_type.compare("3D image-based I/O plugin") == 0) ?
									iomanager::IOPluginFactory::getPlugin3D(iom::IMIN_PLUGIN)->isChansInterleaved() :
									iomanager::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->isChansInterleaved() && 
															vm::VOLUME_INPUT_FORMAT_PLUGIN.compare(vm::MCVolume::id) != 0 ) // 2018-01-23. Giulio. if is an MCVolume then channels are not interleaved

  				; // channels are interleaved: do nothing
			else {
				// channels are not interleaved 'active' channel has to be set after checking that the channel exists
				if ( iom::CHANS == iom::NONE ) { // no R,G, or B have been specified
					if ( iom::CHANS_no < volume->getDIM_C() ) 
						volume->setACTIVE_CHAN(iom::CHANS_no);
					else
						throw iom::exception(vm::strprintf("Channel %d does not exist", iom::CHANS_no).c_str());
				}
				else {
					if ( iom::CHANS == iom::ALL ) {
						if ( volume->getDIM_C() > 1 && vm::VOLUME_INPUT_FORMAT_PLUGIN.compare(vm::MCVolume::id) != 0 ) 
							// the exception is not raised when it is a MultiVolume; in this case both the enabled subvolume and the active channel have been set (possibily to 0 which is the default value)
							throw iom::exception(vm::strprintf("Conversion from multi-channel to intensity images not supported").c_str());
					}
					else if ( iom::CHANS == iom::G ) {
						if ( volume->getDIM_C() < 2 )
							throw iom::exception(vm::strprintf("Green channel does not exist").c_str());
					}
					else if ( iom::CHANS == iom::B ) {
						if ( volume->getDIM_C() < 3 )
							throw iom::exception(vm::strprintf("Blue channel does not exist").c_str());
					}
					// iom::CHANS is an enumerated type and it represents RGB channels with 1(R), 2(G), 3(B)
					// conversely, in a vmVirtualVolume channels are coded starting form 0
					// hence the value of iom::CHANS must be decremented by 1
					volume->setACTIVE_CHAN(iom::CHANS-1); 
				}
			}
		}

		// process volume		
		double total_time = -TIME(0);
		StackStitcher* stitcher = new StackStitcher(volume);
		if(cli.import)
			defaultOutputFileName = "xml_import";
		if(cli.computedisplacements || cli.stitch)
		{
			stitcher->computeDisplacements(cli.pd_algo, cli.start_stack_row, cli.start_stack_col, cli.end_stack_row, cli.end_stack_col, 
										   cli.overlap_V, cli.overlap_H, cli.search_radius_V, cli.search_radius_H, cli.search_radius_D, 
										   cli.subvol_dim_D, cli.enable_restore, cli.restoring_direction, cli.show_progress_bar, cli.D0, cli.D1, 
										   cli.disable_last_row, cli.disable_last_col);
			defaultOutputFileName = "xml_displcomp";
		}
		if(cli.projdisplacements || cli.stitch)
		{
			stitcher->projectDisplacements();
			defaultOutputFileName = "xml_displproj";
		}
		if(cli.thresholddisplacements || cli.stitch)
		{
			stitcher->thresholdDisplacements(cli.reliability_threshold);
			defaultOutputFileName = "xml_displthres";
		}
		if(cli.placetiles || cli.stitch)
		{
			stitcher->computeTilesPlacement(cli.tp_algo);
			defaultOutputFileName = "xml_merging";
		}
		if(cli.makeDirs) {
			stitcher->createDirectoryHierarchy(cli.volume_save_path, cli.slice_height, cli.slice_width, cli.slice_depth, cli.resolutions, cli.exclude_nonstitchables, 
									 cli.start_stack_row, cli.end_stack_row, cli.start_stack_col, cli.end_stack_col, cli.D0, cli.D1, 
									 cli.enable_restore, cli.restoring_direction, cli.tm_blending, cli.halving_method, cli.isotropic, false, cli.show_progress_bar, cli.img_format.c_str(), cli.img_depth);
		}
		if(cli.metaData) {
			stitcher->mdataGenerator(cli.volume_save_path, cli.slice_height, cli.slice_width, cli.slice_depth, cli.resolutions, cli.exclude_nonstitchables, 
									 cli.start_stack_row, cli.end_stack_row, cli.start_stack_col, cli.end_stack_col, cli.D0, cli.D1, 
									 cli.enable_restore, cli.restoring_direction, cli.tm_blending, cli.halving_method, cli.isotropic, false, cli.show_progress_bar, cli.img_format.c_str(), cli.img_depth);
		}
		if(cli.mergetiles || cli.stitch)
		{
#ifdef NEW_MERGE_ACTIVE

			// 2017-12-01. Giulio. @ADDED step 6 is substituted by TeraConverter

			string dst_format;

			if ( vm::VOLUME_OUTPUT_FORMAT_PLUGIN.compare(vm::BlockVolume::id)==0 ) {
				if ( cli.img_format == "raw" )
					dst_format = iim::TILED_FORMAT;
				else if ( cli.img_format == "tif" )
					dst_format = iim::TILED_TIF3D_FORMAT;
				else
					throw iom::exception(iom::strprintf("unsupported file format %s",cli.img_format.c_str()), __iom__current__function__);
			}
			else if ( vm::VOLUME_OUTPUT_FORMAT_PLUGIN.compare(vm::StackedVolume::id)==0 ) {
				if ( cli.img_format == "raw" )
					dst_format = iim::STACKED_RAW_FORMAT;
				else if ( cli.img_format == "tif" )
					dst_format = iim::STACKED_FORMAT;
				else
					throw iom::exception(iom::strprintf("unsupported file format %s",cli.img_format.c_str()), __iom__current__function__);
			}
			else
				throw iom::exception(vm::strprintf("Unsupported output volume format plugin \"%s\"", vm::VOLUME_OUTPUT_FORMAT_PLUGIN.c_str()).c_str());

			UnstitchedVolume *uvolume = new UnstitchedVolume(stitcher->getVOLUME());  
			
			vcDriver(
				uvolume, // the volume has not been imported, the source path has to be passed
				"", // not used
				cli.volume_save_path, 
				"", // not used
				dst_format, 
				0, // invalid value: imout_depth is set inside
				cli.resolutions, 
				"", // select all channels
				"", // not used
				"", // not used
				cli.slice_depth, 
				cli.slice_height, 
				cli.slice_width, 
				1, // no down sampling
				cli.halving_method, 
				cli.libtiff_rowsPerStrip, 
				cli.libtiff_uncompressed, 
				cli.libtiff_bigtiff, 
				cli.show_progress_bar, 
				cli.isotropic, 
				-1, -1, -1, -1, cli.D0, (cli.D1 == -1 ? -1 : cli.D1+1), // all tiles are always selected, if cli.D1 is not default, it is the last slice to be merged (included)
				false, // no time series allowed
				cli.makeDirs, 
				cli.metaData, 
				cli.parallel,
				cli.fixed_tiling,
				"", // not used
				0); // no bit compression allowed
				
				// WARNING: uvolume is not deallocated to avoid deallocation of the volume internal to stitcher which is still needed
#else

			// 2015-02-03. Alessandro. @FIXED bug in the call of mergeTiles (wrong D1 if D1 was not set).
			if ( vm::VOLUME_OUTPUT_FORMAT_PLUGIN.compare(vm::BlockVolume::id)==0 )
				stitcher->mergeTilesVaa3DRaw(cli.volume_save_path, cli.slice_height, cli.slice_width, cli.slice_depth, cli.resolutions, cli.exclude_nonstitchables, 
									 cli.start_stack_row, cli.end_stack_row, cli.start_stack_col, cli.end_stack_col, cli.D0, (cli.D1 == -1 ? -1 : cli.D1+1), 
									 cli.enable_restore, cli.restoring_direction, cli.tm_blending, cli.halving_method, cli.isotropic, false, cli.show_progress_bar, cli.img_format.c_str(), cli.img_depth,cli.parallel);
			else if ( vm::VOLUME_OUTPUT_FORMAT_PLUGIN.compare(vm::StackedVolume::id)==0 )
				stitcher->mergeTiles(cli.volume_save_path, cli.slice_height, cli.slice_width, cli.resolutions, cli.exclude_nonstitchables, 
									 cli.start_stack_row, cli.end_stack_row, cli.start_stack_col, cli.end_stack_col, cli.D0, (cli.D1 == -1 ? -1 : cli.D1+1), 
									 cli.enable_restore, cli.restoring_direction, cli.tm_blending, cli.halving_method, false, cli.show_progress_bar, cli.img_format.c_str(), cli.img_depth);
			else
				throw iom::exception(vm::strprintf("Unsupported output volume format plugin \"%s\"", vm::VOLUME_OUTPUT_FORMAT_PLUGIN.c_str()).c_str());

#endif //  _NEW_MERGE_ACTIVE
		}
		if(cli.test || cli.import)
		{
			if ( cli.test )
				stitcher->mergeTiles(cli.volume_save_path, -1, -1, NULL, false, -1, -1, -1, -1, volume->getN_SLICES()/2, volume->getN_SLICES()/2 +1, 
									 false, false, S_SHOW_STACK_MARGIN, cli.halving_method, true, false, cli.img_format.c_str(), cli.img_depth);
			defaultOutputFileName = "xml_import";
		}

		// save project file and computation times
		if(volume)
			volume->saveXML(0, fillPath(cli.projfile_save_path, volume->getSTACKS_DIR(), defaultOutputFileName, "xml").c_str());
		if(stitcher && cli.save_execution_times)
			stitcher->saveComputationTimes(cli.execution_times_filename.c_str(), *volume, total_time);

		// release objects
		if(stitcher)
			delete stitcher;
		if(volume)
			delete volume;

		// display elapsed time
		printf("\nTime elapsed: %.1f seconds\n\n", proctime + TIME(0));
	}
	catch( iom::exception& exception){
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
