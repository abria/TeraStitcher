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
* 2021-02-12. Giulio.     @FIXED empty constructor 'vm::StackedVolume' must initialize 'STACKS'
* 2019-11-02. Giulio.     @ADDED 'mdata_fname' parameter to constructor from xml
* 2018-03-02. Giulio.     @ADDED the possibility to set a path and a name for the mdata.bin file
* 2017-04-12. Giulio.     @ADDED method to release all buffers allocated in VirtualStack
* 2015-06-12. Giulio      @ADDED 'check' method to check completeness and coherence of a volume
* 2015-02-26. Giulio.     @ADDED initChannels private method to initialize fields DIM_C and BYTESxCHAN
* 2014-09-20. Alessandro. @ADDED overwrite_mdata flag to the XML-based constructor.
* 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
* 2014-09-05. Alessandro. @ADDED 'normalize_stacks_attributes()' method to normalize stacks attributes (width, height, etc.)
*/

#ifndef _VM_STACKED_VOLUME_H
#define _VM_STACKED_VOLUME_H

#include <string>
#include "volumemanager.config.h"
#include <sstream>
#include "iomanager.config.h"
#include "vmVirtualVolume.h" 
#include <cstdarg>
#include <vector>
#include <sstream>
#include <limits>
#include <cstring>


class vm::StackedVolume : public vm::VirtualVolume
{
	private:

		// 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
		static const std::string creator_id1, creator_id2;							
        static vm::VirtualVolume* createFromXML(const char* xml_path, bool ow_mdata) { return new StackedVolume(xml_path, ow_mdata); }
		static vm::VirtualVolume* createFromData(const char* data_path, vm::ref_sys ref, float vxl1, float vxl2, float vxl3, bool ow_mdata, std::string mdata_fname) { 
			return new StackedVolume(data_path, ref, vxl1, vxl2, vxl3, ow_mdata, mdata_fname); 
		}


		//******OBJECT ATTRIBUTES******
		vm::Stack ***STACKS;					//2-D array of <Stack*>	

		//initialization methods
		void init() ;
        void initChannels() ;
		void applyReferenceSystem(vm::ref_sys reference_system, float VXL_1, float VXL_2, float VXL_3) ;

		//binary metadata load/save methods
		void saveBinaryMetadata(char *metadata_filepath) ;
		void loadBinaryMetadata(char *metadata_filepath) ;

		//rotates stacks matrix around D vm::axis (accepted values are theta=0,90,180,270)
		void rotate(int theta);

		//mirrors stacks matrix along mrr_axis (accepted values are mrr_axis=1,2,3)
		void mirror(vm::axis mrr_axis) ;

		// 2014-09-05. Alessandro. @ADDED 'normalize_stacks_attributes()' method to normalize stacks attributes (width, height, etc.)
		void normalize_stacks_attributes() ;

	public:

		// 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
		static const std::string id;	

		//CONSTRUCTORS-DECONSTRUCTOR
		StackedVolume() : vm::VirtualVolume(){ STACKS = (vm::Stack ***) 0; }
        StackedVolume(const char* _stacks_dir, vm::ref_sys reference_system, float VXL_1=0, float VXL_2=0, float VXL_3=0, bool overwrite_mdata=false, std::string mdata_fname="") ;
        StackedVolume(const char *xml_filepath, bool overwrite_mdata=false, std::string mdata_fname="") ;
		~StackedVolume();

		// ******GET METHODS******
		int		 getStacksHeight();
		int		 getStacksWidth();
		vm::VirtualStack*** getSTACKS();

		//loads/saves metadata from/in the given xml filename
		void loadXML(const char *xml_filename) ;
		void initFromXML(const char *xml_filename) ;
        void saveXML(const char *xml_filename=0, const char *xml_filepath=0) ;

		void releaseBuffers() ;

        /**********************************************************************************
        * UTILITY methods
        ***********************************************************************************/

        //check if volume is complete and coherent; return true if the volume is ok, false otherwise
		//if a file name is passed and thevolume is not ok an error log file is generated
		bool check(const char *errlogFileName = 0) ;

        //counts the total number of displacements and the number of displacements per pair of adjacent stacks
        void countDisplacements(int& total, float& per_stack_pair);

        //counts the number of single-direction displacements having a reliability measure above the given threshold
        void countReliableSingleDirectionDisplacements(float threshold, int& total, int& reliable);

        //counts the number of stitchable stacks given the reliability threshold
        int countStitchableStacks(float threshold);

		// print mdata.bin content to stdout
		static void dumpMData(const char* volumePath) ;
};

namespace{																
	const vm::StackedVolume* objectStackedVolume = new vm::StackedVolume();
} 


#endif /* STACKED_VOLUME_H */

