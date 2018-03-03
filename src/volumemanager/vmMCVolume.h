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
* 2018-03-02. Giulio. @ADDED a parameter to set a path and a name for the mdata.bin file for compatibility, but it is not used
* 2018-01-20. Giulio. @CREATED
*/

#ifndef _VM_MC_VOLUME_H
#define _VM_MC_VOLUME_H

#include <string>
#include <sstream>
#include "volumemanager.config.h"
#include "iomanager.config.h"
#include "vmVirtualVolume.h"
#include "vmBlock.h"


class vm::MCVolume : public vm::VirtualVolume
{
	private:

		// 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
		static const std::string creator_id1, creator_id2;							
        static vm::VirtualVolume* createFromXML(const char* xml_path, bool ow_mdata = false) { return new MCVolume(xml_path, ow_mdata); }
		// 2018-01-21. Giulio @ADDED default parameters to support creation from data (folder containing xml import files)
		static vm::VirtualVolume* createFromData(const char* data_path, vm::ref_sys ref = vm::ref_sys(vm::axis(0),vm::axis(0),vm::axis(0)), float vxl1 = 0, float vxl2 = 0, float vxl3 = 0, bool ow_mdata = false, std::string mdata_fname = "") { 
			return new MCVolume(data_path, ref, vxl1, vxl2, vxl3, ow_mdata); // the mdata_fname parameter is not used for MCVolumes
		}

		vm::VirtualVolume **subvolumes;
		std::string *sv_format;
		std::string *xml_file_names;
		int N_SUBVOLS;
		int enabledSV;

		bool checked; // to avoid multiple checking
		bool aligned; // to signal that tiles do not need to be reset (re-aligned)

		// these private fields are initialized internally from the xml import file if specified
		// there are setters for active_res and acrive_tp, but they are intended to be used by teraconverter only when it has to stitch timeseries from unstitched datasets
		std::string active_res;   // active resolution (default. 0)
		std::string active_tp;    // active timepoint (default: 0)
		bool series_no;
		bool additionalIOPluginParams;  // to avoid passing unnecessary additional parameters to ioplugins

		//Given the reference system, initializes all object's members using stack's directories hierarchy
        void init() throw (iom::exception);
        void initChannels() throw (iom::exception);

		void applyReferenceSystem(vm::ref_sys reference_system, float VXL_1, float VXL_2, float VXL_3) throw (iom::exception);

		//binary metadata load/save methods
		void saveBinaryMetadata(char *metadata_filepath) throw (iom::exception);
		void loadBinaryMetadata(char *metadata_filepath) throw (iom::exception);

		//rotate stacks matrix around D vm::axis (accepted values are theta=0,90,180,270)
		void rotate(int theta);

		//mirror stacks matrix along mrr_axis (accepted values are mrr_axis=1,2,3)
		void mirror(vm::axis mrr_axis);

		// 2014-09-05. Alessandro. @ADDED 'normalize_stacks_attributes()' method to normalize stacks attributes (width, height, etc.)
		//void normalize_stacks_attributes() throw (iom::exception);

	public:

		// 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'StackedVolume' a volume format plugin.
		static const std::string id;		

		//CONSTRUCTORS-DECONSTRUCTOR
		MCVolume() : vm::VirtualVolume(){}
        MCVolume(const char* _stacks_dir, vm::ref_sys reference_system = vm::ref_sys(vm::axis(0),vm::axis(0),vm::axis(0)), float VXL_1=0, float VXL_2=0, float VXL_3=0, bool overwrite_mdata=false) throw (iom::exception);
        MCVolume(const char *xml_filepath, bool overwrite_mdata=false) throw (iom::exception);
		~MCVolume();

		// ******GET METHODS******
		int		 getStacksHeight()   {return subvolumes[0]->getSTACKS()[0][0]->getHEIGHT();}
		int		 getStacksWidth()    {return subvolumes[0]->getSTACKS()[0][0]->getWIDTH();}
		vm::VirtualStack*** getSTACKS()  {return subvolumes[active_channel]->getSTACKS();}

		std::string getACTIVE_RES ( ) { return active_res; }
		void        setACTIVE_REC (int r) { active_res = r; }

		std::string getACTIVE_TP  ( ) { return active_tp; }
		void        setACTIVE_TP  (int t) { active_tp = t; }

		int         getN_ROWS() { return subvolumes[active_channel]->getN_ROWS(); }
		int	        getN_COLS() { return subvolumes[active_channel]->getN_COLS(); }
		int         getN_SLICES() { return subvolumes[active_channel]->getN_SLICES(); }
		int	        getOVERLAP_V() { return subvolumes[active_channel]->getOVERLAP_V(); }
		int	        getOVERLAP_H() { return subvolumes[active_channel]->getOVERLAP_H(); }
		int         getDEFAULT_DISPLACEMENT_V() { return subvolumes[active_channel]->getDEFAULT_DISPLACEMENT_V(); }
		int         getDEFAULT_DISPLACEMENT_H() { return subvolumes[active_channel]->getDEFAULT_DISPLACEMENT_H(); }
		virtual int getDEFAULT_DISPLACEMENT_D() { return subvolumes[active_channel]->getDEFAULT_DISPLACEMENT_D(); }
		vm::ref_sys getREF_SYS() {return subvolumes[active_channel]->getREF_SYS(); }


		bool getSERIES_NO ( ) { return series_no; }
		bool getADDITIONAL_IOPLUGIN_PARAMETERS ( ) { return additionalIOPluginParams; }

		bool getALIGNED () { return aligned; }

		void setADDITIONAL_IOPLUGIN_PARAMETERS ( bool _flag) { additionalIOPluginParams = _flag; }
		void setENABLEDSV ( int sv ) { enabledSV = sv; }

		//loads/saves metadata from/in the given xml filename
		void loadXML(const char *xml_filename) throw (iom::exception);
		void initFromXML(const char *xml_filename) throw (iom::exception);
        void saveXML(const char *xml_filename=0, const char *xml_filepath=0) throw (iom::exception);

		void releaseBuffers();

		void resetTilePositions ( );

        /**********************************************************************************
        * UTILITY methods
        ***********************************************************************************/

        //check if volume is complete and coherent; return true if the volume is ok, false otherwise
		//if a file name is passed and the volume is not ok an error log file is generated
		bool check(const char *errlogFileName = 0) throw (iom::exception);

        //counts the total number of displacements and the number of displacements per pair of adjacent stacks
        void countDisplacements(int& total, float& per_stack_pair);

        //counts the number of single-direction displacements having a reliability measure above the given threshold
        void countReliableSingleDirectionDisplacements(float threshold, int& total, int& reliable);

        //counts the number of stitchable stacks given the reliability threshold
        int countStitchableStacks(float threshold);
};

namespace{																
	const vm::MCVolume* objectMCVolume = new vm::MCVolume();
} 


#endif /* _VM_MC_VOLUME_H */
