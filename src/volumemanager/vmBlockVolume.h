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
* 2021-02-12. Giulio.     @FIXED empty constructor 'vm::BlockedVolume' must initialize 'BLOCKS'
* 2019-11-02. Giulio.     @ADDED 'mdata_fname' parameter to constructor from xml
* 2018-03-02. Giulio.     @ADDED the possibility to set a path and a name for the mdata.bin file
* 2017-06-26. Giulio.     @ADDED methods to set active resolution and active timepoint
* 2017-04-12. Giulio.     @ADDED method to release all buffers allocated in VirtualStack
* 2016-10-27. Giulio.     @ADDED string fields for control over the subimage to be exposed through the xml import file  
* 2015-06-12. Giulio      @ADDED 'check' method to check completeness and coherence of a volume
* 2015-02-26. Giulio.     @ADDED initChannels private method to initialize fields DIM_C and BYTESxCHAN
* 2015-01-17. Alessandro. @FIXED missing throw(iom::exception) declaration in loadXML and initFromXML methods.
* 2014-09-20. Alessandro. @ADDED overwrite_mdata flag to the XML-based constructor.
* 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
* 2014-09-05. Alessandro. @ADDED 'normalize_stacks_attributes()' method to normalize stacks attributes (width, height, etc.)
*/

#ifndef _VM_BLOCK_VOLUME_H
#define _VM_BLOCK_VOLUME_H

#include <string>
#include <sstream>
#include "volumemanager.config.h"
#include "iomanager.config.h"
#include "vmVirtualVolume.h" 
#include "vmBlock.h"


class vm::BlockVolume : public vm::VirtualVolume
{
	private:

		// 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
		static const std::string creator_id1, creator_id2;							
        static vm::VirtualVolume* createFromXML(const char* xml_path, bool ow_mdata) { return new BlockVolume(xml_path, ow_mdata); }
		static vm::VirtualVolume* createFromData(const char* data_path, vm::ref_sys ref, float vxl1, float vxl2, float vxl3, bool ow_mdata, std::string mdata_fname) { 
			return new BlockVolume(data_path, ref, vxl1, vxl2, vxl3, ow_mdata,mdata_fname); 
		}

		vm::Block ***BLOCKS;			    //2-D array of <Block*>

		// these private fields are initialized internally from the xml import file if specified
		// there are setters for active_res and acrive_tp, but they are intended to be used by teraconverter only when it has to stitch timeseries from unstitched datasets
		std::string active_res;   // active resolution (default. 0)
		std::string active_tp;    // active timepoint (default: 0)
		bool series_no;
		bool additionalIOPluginParams;  // to avoid passing unnecessary additional parameters to ioplugins

		//Given the reference system, initializes all object's members using stack's directories hierarchy
        void init() ;
        void initChannels() ;

		void applyReferenceSystem(vm::ref_sys reference_system, float VXL_1, float VXL_2, float VXL_3) ;

		//binary metadata load/save methods
		void saveBinaryMetadata(char *metadata_filepath) ;
		void loadBinaryMetadata(char *metadata_filepath) ;

		//rotate stacks matrix around D vm::axis (accepted values are theta=0,90,180,270)
		void rotate(int theta);

		//mirror stacks matrix along mrr_axis (accepted values are mrr_axis=1,2,3)
		void mirror(vm::axis mrr_axis);

		// iannello returns the number of channels of images composing the volume
		//void initChannels ( ) throw (iom::iom::exception);

		// 2014-09-05. Alessandro. @ADDED 'normalize_stacks_attributes()' method to normalize stacks attributes (width, height, etc.)
		void normalize_stacks_attributes() ;

	public:

		// 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
		static const std::string id;		

		//CONSTRUCTORS-DECONSTRUCTOR
		BlockVolume() : vm::VirtualVolume(){ BLOCKS = (vm::Block ***) 0; }
        BlockVolume(const char* _stacks_dir, vm::ref_sys reference_system, float VXL_1=0, float VXL_2=0, float VXL_3=0, bool overwrite_mdata=false, std::string mdata_fname="") ;
        BlockVolume(const char *xml_filepath, bool overwrite_mdata=false, std::string mdata_fname="") ;
		~BlockVolume();

		// ******GET METHODS******
		int		 getStacksHeight()   {return BLOCKS[0][0]->getHEIGHT();}
		int		 getStacksWidth()    {return BLOCKS[0][0]->getWIDTH();}
		vm::VirtualStack*** getSTACKS()  {return (vm::VirtualStack***)this->BLOCKS;}

		std::string getACTIVE_RES ( ) { return active_res; }
		void        setACTIVE_REC (int r) { active_res = r; }

		std::string getACTIVE_TP  ( ) { return active_tp; }
		void        setACTIVE_TP  (int t) { active_tp = t; }

		bool getSERIES_NO ( ) { return series_no; }
		bool getADDITIONAL_IOPLUGIN_PARAMETERS ( ) { return additionalIOPluginParams; }

		void setADDITIONAL_IOPLUGIN_PARAMETERS ( bool _flag) { additionalIOPluginParams = _flag; }

		//loads/saves metadata from/in the given xml filename
		void loadXML(const char *xml_filename) ;
		void initFromXML(const char *xml_filename) ;
        void saveXML(const char *xml_filename=0, const char *xml_filepath=0) ;

		void releaseBuffers() ;

        /**********************************************************************************
        * UTILITY methods
        ***********************************************************************************/

        //check if volume is complete and coherent; return true if the volume is ok, false otherwise
		//if a file name is passed and the volume is not ok an error log file is generated
		bool check(const char *errlogFileName = 0) ;

        //counts the total number of displacements and the number of displacements per pair of adjacent stacks
        void countDisplacements(int& total, float& per_stack_pair);

        //counts the number of single-direction displacements having a reliability measure above the given threshold
        void countReliableSingleDirectionDisplacements(float threshold, int& total, int& reliable);

        //counts the number of stitchable stacks given the reliability threshold
        int countStitchableStacks(float threshold);
};

namespace{																
	const vm::BlockVolume* objectBlockVolume = new vm::BlockVolume();
} 


#endif /* BLOCK_VOLUME_H */
