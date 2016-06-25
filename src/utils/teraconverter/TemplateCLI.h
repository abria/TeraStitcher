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
* 2016-06-18  Giulio.     @ADDED option for downsampling the reading of data
* 2016-04-13  Giulio.     @ADDED options for parallelizing teraconverter
*/

#ifndef _TEMPLATE_COMMAND_LINE_INTERFACE_H
#define _TEMPLATE_COMMAND_LINE_INTERFACE_H

#include <string>
#include "iomanager.config.h"
#include "GUI_config.h"

using namespace std;

class TemplateCLI
{
	public:

		// switch parameters
		//bool highest_resolution;						// generate highest resolution (default: false)
        bool makeDirs;                          //creates the directory hiererchy
        bool metaData;                          //creates the mdata.bin file of the output volume
        bool parallel;                          //parallel mode: does not perform side-effect operations during merge
        bool isotropic;                         //generate lowest resolutiona with voxels as much isotropic as possible

		bool pluginsinfo;						//display plugins information

		// other parameters
		// int/float/double/string XXXX;	// description
		string src_root_dir;
		string dst_root_dir;
		int slice_depth;
		int slice_height;
		int slice_width;
		string src_format;
		string dst_format;
		bool resolutions[S_MAX_MULTIRES];
		int halving_method;
		bool show_progress_bar;					//enables/disables progress bar with estimated time remaining

		string outFmt;
		string infofile_path;					//file path of the info log file to be saved
		int downsamplingFactor;                 //downsampling factor to be used to read source volume (only if it is a serie of 2D slices)

		// vertices defining the subvolume to be converted
		int V0;
		int V1;
		int H0;
		int H1;
		int D0;
		int D1;

		int tm_blending;						//tiles merging blending type

		//constructor - deconstructor
		TemplateCLI(void);					//set default params
		~TemplateCLI(void){};

		//reads options and parameters from command line
		void readParams(int argc, char** argv) throw (iom::exception);

		//checks parameters correctness
		void checkParams() throw (iom::exception);

		//returns help text
		string getHelpText();

		//print all arguments
		void print();
};

#endif /* _TERASTITCHER_COMMAND_LINE_INTERFACE_H */


