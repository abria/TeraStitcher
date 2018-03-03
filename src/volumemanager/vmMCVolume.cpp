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
* 2018-03-01. Giulio. @FIXED global variable 'vm::VOLUME_INPUT_FORMAT_PLUGIN' is set to vm::MCVolume::id in the constructor to signal that it is a non interleaved channel representation
* 2018-01-20. Giulio. @CREATED
*/

/****************************
* Management of sparse data *
*****************************
*
* If --sparse_data flag is set, in the import step, after reading files names missing blocks are detected
*
* Each tile is a list of blocks some of which may be empty
* empty blocks have a null pointer as file name, but have first index and size correctly set
* this structure is binarized into the mdata.bin file
* the z_ranges variable store for each tile the intervals corresponding to exisiting blocks
* thi information is stored into the xml file
* both internal tile structure and z_ranges field are set every time the volume is created 
*/

#ifdef _WIN32
#include "dirent_win.h"
#else
#include <dirent.h>
#endif
#include <list>
#include <fstream>
#include <sstream>
#include <set>
#include "vmMCVolume.h"
#include "S_config.h"
//#include <string>
#include "../imagemanager/IM_config.h"
#include "vmStackedVolume.h"
#include "vmBlockVolume.h"


using namespace std;
using namespace iom;
using namespace vm;

// 2014-09-10. Alessandro. @ADDED plugin creation/registration functions to make 'MCVolume' a volume format plugin.
const std::string MCVolume::id = "MultiVolume";
const std::string MCVolume::creator_id1 = volumemanager::VirtualVolumeFactory::registerPluginCreatorXML(&createFromXML, MCVolume::id);
const std::string MCVolume::creator_id2 = volumemanager::VirtualVolumeFactory::registerPluginCreatorData(&createFromData, MCVolume::id);


MCVolume::MCVolume(const char* _stacks_dir, vm::ref_sys _reference_system, float VXL_1, float VXL_2, float VXL_3, bool overwrite_mdata) throw (iom::exception)
	: VirtualVolume(_stacks_dir, _reference_system, VXL_1, VXL_2, VXL_3)
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin MCVolume::MCVolume(_stacks_dir=%s, reference_system = {%d,%d,%d}, VXL_1 = %.2f, VXL_2 = %.2f, VXL_3 = %.2f)\n", 
		  _stacks_dir,reference_system.first, reference_system.second, reference_system.third, VXL_1, VXL_2, VXL_3);
	#endif

	// 2018-03-01. Giulio. @FIXED global variable 'vm::VOLUME_INPUT_FORMAT_PLUGIN' has to be set to vm::MCVolume::id to signal that it is a non interleaved channel representation
	vm::VOLUME_INPUT_FORMAT_PLUGIN = vm::MCVolume::id;
	
	checked = false;
	aligned = false; // no guarantee that subvolumes are aligned

	// in this case active resolution and timepoint are always 0
	active_res = active_tp = "0";
	series_no = additionalIOPluginParams = false;

	// subvolume 0 is the active one by default
	enabledSV = 0;

	init();
	initChannels();
	check();

	// no guarantee that subvolumes are aligned: check tile matices and align homologous tiles
	resetTilePositions();

	// all subvolumes refer to the same space: use the origin of the enabled subvolume
	ORG_V = subvolumes[enabledSV]->getORG_V();
	ORG_H = subvolumes[enabledSV]->getORG_H();
	ORG_D = subvolumes[enabledSV]->getORG_D();

	// set to invalid value: each subvolume may have its own mechanical displacement
	MEC_V = -1;
	MEC_H = -1;

	// set to invalid value: each subvolume may have a different tile and a different number of slices
	N_ROWS   = -1;
	N_COLS   = -1;
	N_SLICES = -1;

	// he active subvolume is by default 0
	active_channel = 0;
}

MCVolume::MCVolume(const char *xml_filepath, bool overwrite_mdata) throw (iom::exception)
	: VirtualVolume(xml_filepath)
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin MCVolume::MCVolume(xml_filepath=%s)\n", xml_filepath);
	#endif

	// 2018-03-01. Giulio. @FIXED global variable 'vm::VOLUME_INPUT_FORMAT_PLUGIN' has to be set to vm::MCVolume::id to signal that it is a non interleaved channel representation
	vm::VOLUME_INPUT_FORMAT_PLUGIN = vm::MCVolume::id;
	
 	checked = false;
	aligned = false; // no guarantee that subvolumes are aligned

	// in this case active resolution and timepoint set to the dafault, but can be changed when xml import file is loaded
	active_res = active_tp = "0";
	series_no = additionalIOPluginParams = false;

   //extracting <stacks_dir> field from XML
    TiXmlDocument xml;
    if(!xml.LoadFile(xml_filepath))
    {
        char errMsg[2000];
        sprintf(errMsg,"in MCVolume::MCVolume(xml_filepath = \"%s\") : unable to load xml", xml_filepath);
        throw iom::exception(errMsg);
    }
    TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));
    TiXmlElement * pelem = hRoot.FirstChildElement("subvolumes_dir").Element();
    this->stacks_dir = new char[strlen(pelem->Attribute("value"))+1];
    strcpy(this->stacks_dir, pelem->Attribute("value"));

	// load xml content
	initFromXML(xml_filepath);

	initChannels();
	check();

	if ( !aligned ) 
		resetTilePositions();
}

MCVolume::~MCVolume()
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin MCVolume::~MCVolume(void)\n");
	#endif

	if(stacks_dir)
		delete[] stacks_dir;

	if(subvolumes)
	{
		for(int sv=0; sv<N_SUBVOLS; sv++)
		{
			delete subvolumes[sv];
		}
		delete[] subvolumes;
	}
	if(sv_format)
	{
		delete[] sv_format;
	}
	if(xml_file_names)
	{
		delete[] xml_file_names;
	}
}



void MCVolume::init() throw (iom::exception)
{
	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	//LOCAL VARIABLES
	string tmp_path;				//string that contains temp paths during computation
    string tmp;						//string that contains temp data during computation
	DIR *cur_dir_lev1;				//pointer to DIR, the data structure that represents a DIRECTORY (level 1 of hierarchical structure)
	dirent *entry_lev1;				//pointer to DIRENT, the data structure that represents a DIRECTORY ENTRY inside a directory (level 1)
	int i=0, j=0;					//for counting of N_ROWS, N_COLS
 //   list<Block*> blocks_list;                       //each stack found in the hierarchy is pushed into this list
    list<string> entries_lev1;                      //list of entries of first level of hierarchy
    list<string>::iterator entry_i;                 //iterator for list 'entries_lev1'
    list<string> entries_lev2;                      //list of entries of second level of hierarchy
    list<string>::iterator entry_j;                 //iterator for list 'entries_lev2'

	//obtaining DIR pointer to root_dir (=NULL if directory doesn't exist)
	if (!(cur_dir_lev1=opendir(stacks_dir)))
	{
        char msg[iim::STATIC_STRINGS_SIZE];
        sprintf(msg,"in MCVolume::init(...): Unable to open directory \"%s\"", stacks_dir);
        throw iim::IOException(msg);
	}

	//scanning first level of hierarchy which entries need to be ordered alphabetically. This is done using STL.
	while ((entry_lev1=readdir(cur_dir_lev1)))
	{
        tmp=entry_lev1->d_name;
        if((tmp.find(".") == string::npos && tmp.find(" ") == string::npos) || (tmp.find(".xml") != string::npos))
                entries_lev1.push_front(entry_lev1->d_name);
	}
	closedir(cur_dir_lev1);
	entries_lev1.sort();

	N_SUBVOLS = (int) entries_lev1.size();
	subvolumes = new VirtualVolume *[N_SUBVOLS];
	sv_format = new string[N_SUBVOLS];
	xml_file_names = new string[N_SUBVOLS];

	//for each entry creates a VirtualVolume
	for(entry_i = entries_lev1.begin(), i=0; entry_i!= entries_lev1.end(); entry_i++, i++) {
		//building absolute path of first level entry
		tmp_path=stacks_dir;
		tmp_path.append("/");
		tmp_path.append(*entry_i);
		subvolumes[i] = volumemanager::VirtualVolumeFactory::createFromXML(tmp_path.c_str(),false);
		sv_format[i] = subvolumes[i]->getVolumeFormat(tmp_path.c_str());
		xml_file_names[i] = tmp_path.c_str();
	}
	reference_system = subvolumes[enabledSV]->getREF_SYS();

	entries_lev1.clear();
}

void MCVolume::initChannels()  throw (iom::exception) 
{
	DIM_C = N_SUBVOLS;
	BYTESxCHAN = subvolumes[enabledSV]->getSTACKS()[0][0]->getN_BYTESxCHAN();
}

void MCVolume::resetTilePositions ( ) {
	int n_rows = subvolumes[enabledSV]->getN_ROWS();
	int n_cols = subvolumes[enabledSV]->getN_COLS();

	for ( int sv=0; sv<N_SUBVOLS; sv++ ) {
		if ( sv != enabledSV ) { // enabledSV is skipped
			if ( n_rows != subvolumes[sv]->getN_ROWS() || n_cols != subvolumes[sv]->getN_COLS() )
			{
				char msg[iim::STATIC_STRINGS_SIZE];
				sprintf(msg,"in MCVolume::resetTilePositions(...): subvolume %d from xml file \"%s\" has a tile matrix (%d,%d) instead of (%d,%d)", 
													sv, xml_file_names[sv].c_str(), subvolumes[sv]->getN_ROWS(), subvolumes[sv]->getN_COLS(), n_rows, n_cols);
				throw iim::IOException(msg);
			}
			else { // the tile positions of the subvolume can be reset
				for ( int i=0; i<n_rows; i++ ) {
					for ( int j=0; j<n_cols; j++ ) {
						subvolumes[sv]->getSTACKS()[i][j]->setABS_V(subvolumes[enabledSV]->getSTACKS()[i][j]->getABS_V());
						subvolumes[sv]->getSTACKS()[i][j]->setABS_H(subvolumes[enabledSV]->getSTACKS()[i][j]->getABS_H());
						subvolumes[sv]->getSTACKS()[i][j]->setABS_D(subvolumes[enabledSV]->getSTACKS()[i][j]->getABS_D());
					}
				}
			}
		}
	}
}


void MCVolume::applyReferenceSystem(vm::ref_sys reference_system, float VXL_1, float VXL_2, float VXL_3) throw (iom::exception)
{
	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);
}

void MCVolume::saveBinaryMetadata(char *metadata_filepath) throw (iom::exception)
{
	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);
}

void MCVolume::loadBinaryMetadata(char *metadata_filepath) throw (iom::exception)
{
	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);
}

//rotate stacks matrix around D vm::axis (accepted values are theta=0,90,180,270)
void MCVolume::rotate(int theta)
{
	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);
}

//mirror stacks matrix along mrr_axis (accepted values are mrr_axis=1,2,3)
void MCVolume::mirror(vm::axis mrr_axis)
{
	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);
}


//check if volume is complete and coherent
bool MCVolume::check(const char *errlogFileName) throw (iom::exception)
{
	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	bool ok = true;

	// check that subvolumens are intensity images
	for ( int sv=0; sv<N_SUBVOLS; sv++ ) {
		if ( subvolumes[sv]->getDIM_C() > 1 ) 
		{
 			iom::warning(iom::strprintf("in MCVolume::check(...): subvolumes with more than one channel not supported (subvolume %d from xml file \"%s\" has %d channels)",
																				sv, xml_file_names[sv].c_str(), subvolumes[sv]->getDIM_C()).c_str(), __iom__current__function__);
			ok = false;
		}
		// check that all subvolumes have the same input plugin (there is a global variable)
		if ( subvolumes[sv]->getInputPlugin(xml_file_names[sv]).compare(iom::IMIN_PLUGIN) != 0 )
		{
 			iom::warning(iom::strprintf("in MCVolume::check(...): subvolume %d from xml file \"%s\" has input plugin %s instead of %s",
				sv, xml_file_names[sv].c_str(), subvolumes[sv]->getInputPlugin(xml_file_names[sv]).c_str(), iom::IMIN_PLUGIN.c_str()).c_str(), __iom__current__function__);
			ok = false;
		}
		// check that all subvolumes have the same bytes per channel 
		if ( subvolumes[sv]->getBYTESxCHAN() != BYTESxCHAN )
		{
 			iom::warning(iom::strprintf("in MCVolume::check(...): subvolume %d from xml file \"%s\" has %d bytes per channel instead of %d",
															sv, xml_file_names[sv].c_str(), subvolumes[sv]->getBYTESxCHAN(), BYTESxCHAN).c_str(), __iom__current__function__);
			ok = false;
		}
	}

	if ( ok ) 
		checked = true; // the volume has been successfully checked

	return ok;
}


void MCVolume::loadXML(const char *xml_filepath) throw (iom::exception)
{
//	#if VM_VERBOSE > 3
//	printf("\t\t\t\tin MCVolume::loadXML(char *xml_filepath = %s)\n", xml_filepath);
//	#endif
//
//	TiXmlDocument xml;
//	if(!xml.LoadFile(xml_filepath))
//	{
//		char errMsg[2000];
//		sprintf(errMsg,"in MCVolume::loadXML(xml_filepath = \"%s\") : unable to load xml", xml_filepath);
//		throw iom::exception(errMsg);
//	}
//
//	//setting ROOT element (that is the first child, i.e. <TeraStitcher> node)
//	TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));
//
//	// 2014-09-10. Alessandro. @ADDED 'volume_format' attribute to <TeraStitcher> XML node
//	const char *volformat = hRoot.ToElement()->Attribute("volume_format");
//	if(volformat && strcmp(volformat, id.c_str()) != 0)
//		throw iom::exception(vm::strprintf("in MCVolume::initFromXML(): unsupported volume_format = \"%s\" (current format is \"%s\")", volformat, id.c_str()).c_str());
//
//	// 2017-04-27. Giulio. ADDED 'input_plugin' attribute to <TeraStitcher> XML node
//	const char *inplugin = hRoot.ToElement()->Attribute("input_plugin"); 
//	if(inplugin)
//		iom::IMIN_PLUGIN = inplugin;
//	
//	//reading fields and checking coherence with metadata previously read from VM_BIN_METADATA_FILE_NAME
//	TiXmlElement * pelem = hRoot.FirstChildElement("stacks_dir").Element();
//	if(strcmp(pelem->Attribute("value"), stacks_dir) != 0)
//	{
//		char errMsg[2000];
//		sprintf(errMsg, "in MCVolume::loadXML(...): Mismatch in <stacks_dir> field between xml file (=\"%s\") and %s (=\"%s\").", pelem->Attribute("value"), vm::BINARY_METADATA_FILENAME.c_str(), stacks_dir);
//		throw iom::exception(errMsg);
//	}
//	// 2014-11-06. Giulio. @ADDED saved reference system into XML file
//	vm::ref_sys reference_system_read;
//	if ( (pelem = hRoot.FirstChildElement("ref_sys").Element()) != 0 ) { // skip if not present (for compatibility with previous versions)
//		pelem->QueryIntAttribute("ref1", (int *) &reference_system_read.first);
//		pelem->QueryIntAttribute("ref2", (int *) &reference_system_read.second);
//		pelem->QueryIntAttribute("ref3", (int *) &reference_system_read.third);
//		if (reference_system_read.first != reference_system.first || reference_system_read.second != reference_system.second || reference_system_read.third != reference_system.third ) 
//		{
//			char errMsg[2000];
//			sprintf(errMsg, "in MCVolume::loadXML(...): Mismatch in <ref_sys> field between xml file (= (%d,%d,%d) ) and %s (= (%d,%d,%d) ).", 
//				reference_system_read.first, reference_system_read.second, reference_system_read.third, vm::BINARY_METADATA_FILENAME.c_str(), reference_system.first, reference_system.second, reference_system.third);
//			throw iom::exception(errMsg);
//		}
//	}
//	pelem = hRoot.FirstChildElement("voxel_dims").Element();
//	float VXL_V_read=0.0f, VXL_H_read=0.0f, VXL_D_read=0.0f;
//	pelem->QueryFloatAttribute("V", &VXL_V_read);
//	pelem->QueryFloatAttribute("H", &VXL_H_read);
//	pelem->QueryFloatAttribute("D", &VXL_D_read);
//	if(VXL_V_read != VXL_V || VXL_H_read != VXL_H || VXL_D_read != VXL_D)
//	{
//		char errMsg[2000];
//		sprintf(errMsg, "in MCVolume::loadXML(...): Mismatch in <voxel_dims> field between xml file (= %.2f x %.2f x %.2f ) and %s (= %.2f x %.2f x %.2f ).", VXL_V_read, VXL_H_read, VXL_D_read, vm::BINARY_METADATA_FILENAME.c_str(), VXL_V, VXL_H, VXL_D);
//		throw iom::exception(errMsg);
//	}
//	pelem = hRoot.FirstChildElement("origin").Element();
//	float ORG_V_read=0.0f, ORG_H_read=0.0f, ORG_D_read=0.0f;
//	pelem->QueryFloatAttribute("V", &ORG_V_read);
//	pelem->QueryFloatAttribute("H", &ORG_H_read);
//	pelem->QueryFloatAttribute("D", &ORG_D_read);
//	/*if(ORG_V_read != ORG_V || ORG_H_read != ORG_H || ORG_D_read != ORG_D)
//	{
//		char errMsg[2000];
//		sprintf(errMsg, "in MCVolume::loadXML(...): Mismatch in <origin> field between xml file (= {%.7f, %.7f, %.7f} ) and %s (= {%.7f, %.7f, %.7f} ).", ORG_V_read, ORG_H_read, ORG_D_read, VM_BIN_METADATA_FILE_NAME, ORG_V, ORG_H, ORG_D);
//		throw iom::iom::exception(errMsg);
//	} @TODO: bug with float precision causes often mismatch */ 
//
//	// 2016-10-27. Giulio. New field in the xml import file to select a subimage (resolution, timepoint, series_no)
//	if ( (pelem = hRoot.FirstChildElement("subimage").Element()) != 0 ) { // skip if not present (for compatibility with previous versions)
//		int value;
//		std::stringstream str;
//		if ( pelem->QueryIntAttribute("resolution", &value) == TIXML_SUCCESS ) {
//			// 2017-06-27. Giulio. If the attribute is present additional parameters are always enabled 
//			//if ( value ) {// additional parameters are not needed if resolution is zero
//				additionalIOPluginParams = true;
//				str << value;
//				active_res = str.str();
//			//}
//		}
//		if ( pelem->QueryIntAttribute("timepoint", &value) == TIXML_SUCCESS ) {
//			// 2017-06-27. Giulio. If the attribute is present additional parameters are always enabled 
//			//if ( value ) { // additional parameters are not needed if resolution is zero
//				additionalIOPluginParams = true;
//				str.str("");
//				str << value;
//				active_tp = str.str();
//			//}
//		}
//		const char *series_no_flag=pelem->Attribute("series_no");
//		if ( series_no_flag ) {
//			if ( strcmp(series_no_flag,"true") == 0 )
//				series_no = additionalIOPluginParams = true;
//		}
//	}
//
//	pelem = hRoot.FirstChildElement("mechanical_displacements").Element();
//	float MEC_V_read=0.0f, MEC_H_read=0.0f;
//	pelem->QueryFloatAttribute("V", &MEC_V_read);
//	pelem->QueryFloatAttribute("H", &MEC_H_read);
//	if(fabs(MEC_V_read - MEC_V) > MECH_MISMATCH || fabs(MEC_H_read - MEC_H) > MECH_MISMATCH)
//	{
//		char errMsg[2000];
//		sprintf(errMsg, "in MCVolume::loadXML(...): Mismatch in <mechanical_displacements> field between xml file (= %.1f x %.1f ) and %s (= %.1f x %.1f ).", MEC_V_read, MEC_H_read, vm::BINARY_METADATA_FILENAME.c_str(), MEC_V, MEC_H);
//		throw iom::exception(errMsg);
//	}
//	pelem = hRoot.FirstChildElement("dimensions").Element();
//	int N_ROWS_read, N_COLS_read, N_SLICES_read;
//	pelem->QueryIntAttribute("stack_rows", &N_ROWS_read);
//	pelem->QueryIntAttribute("stack_columns", &N_COLS_read);
//	pelem->QueryIntAttribute("stack_slices", &N_SLICES_read);
//	if ( N_SLICES_read <= 0 ) // the xml file is an incomplete one: data loaded from mdata.bin wins
//		N_SLICES_read = N_SLICES;
//	if(N_ROWS_read != N_ROWS || N_COLS_read != N_COLS || N_SLICES_read != N_SLICES)
//	{
//		char errMsg[2000];
//		sprintf(errMsg, "in MCVolume::loadXML(...): Mismatch between in <dimensions> field xml file (= %d x %d x %d ) and %s (= %d x %d x %d).", N_ROWS_read, N_COLS_read, N_SLICES_read, vm::BINARY_METADATA_FILENAME.c_str(), N_ROWS, N_COLS, N_SLICES);
//		throw iom::exception(errMsg);
//	}
//
//	//pelem = hRoot.FirstChildElement("STACKS").Element()->FirstChildElement();
//	//int i,j;
//	//for(i=0; i<N_ROWS; i++) {
//	//	for(j=0; j<N_COLS; j++, pelem = pelem->NextSiblingElement()) {
//	//		subvolumes[i]->loadXML(pelem);
//	//		// 2016-11-27. Giulio. @ADDED 'SERIES_NO' attribute in the xml node to identify stack into multi-stack files
//	//		if ( series_no ) {
//	//			const char* series_no_str = pelem->Attribute("SERIES_NO");
//	//			if ( series_no_str ) {
//	//				if ( atoi(series_no_str) >= 0 ) // series_no id is a valid value
//	//					BLOCKS[i][j]->series_no = series_no_str;
//	//			}
//	//		}
//	//	}
//	//}
}

void MCVolume::initFromXML(const char *xml_filepath) throw (iom::exception)
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin MCVolume::initFromXML(char *xml_filename = %s)\n", xml_filename);
		#endif

	TiXmlDocument xml;
	if(!xml.LoadFile(xml_filepath))
	{
		char errMsg[2000];
		sprintf(errMsg,"in MCVolume::initFromXML(xml_filepath = \"%s\") : unable to load xml", xml_filepath);
		throw iom::exception(errMsg);
	}

	//setting ROOT element (that is the first child, i.e. <TeraStitcher> node)
	TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));

	// 2014-09-10. Alessandro. @ADDED 'volume_format' attribute to <TeraStitcher> XML node
	const char *volformat = hRoot.ToElement()->Attribute("volume_format");
	if(volformat && strcmp(volformat, id.c_str()) != 0)
		throw iom::exception(vm::strprintf("in MCVolume::initFromXML(): unsupported volume_format = \"%s\" (current format is \"%s\")", volformat, id.c_str()).c_str());

	// 2017-04-27. Giulio. ADDED 'input_plugin' attribute to <TeraStitcher> XML node
	const char *inplugin = hRoot.ToElement()->Attribute("input_plugin"); 
	if(inplugin)
		iom::IMIN_PLUGIN = inplugin;
	
	//reading fields
	TiXmlElement * pelem = hRoot.FirstChildElement("stacks_dir").Element();
	// 2014-11-06. Giulio. @ADDED saved reference system into XML file
	if ( (pelem = hRoot.FirstChildElement("ref_sys").Element()) != 0 ) { // skip if not present (for compatibility with previous versions)
		pelem->QueryIntAttribute("ref1", (int *) &reference_system.first);
		pelem->QueryIntAttribute("ref2", (int *) &reference_system.second);
		pelem->QueryIntAttribute("ref3", (int *) &reference_system.third);
	}
	else {
		// 2014-11-06. Giulio. @MOVED in case XML is old
		// 2014-09-09. Alessandro. @FIXED. Added default reference system if volume is imported from xml.
		reference_system = vm::ref_sys(vm::vertical,vm::horizontal,vm::depth);
	}
	pelem = hRoot.FirstChildElement("voxel_dims").Element();
	pelem->QueryFloatAttribute("V", &VXL_V);
	pelem->QueryFloatAttribute("H", &VXL_H);
	pelem->QueryFloatAttribute("D", &VXL_D);
	pelem = hRoot.FirstChildElement("origin").Element();
	pelem->QueryFloatAttribute("V", &ORG_V);
	pelem->QueryFloatAttribute("H", &ORG_H);
	pelem->QueryFloatAttribute("D", &ORG_D);

	//// 2016-10-27. Giulio. New field in the xml import file to select a subimage (resolution, timepoint, series_no)
	//if ( (pelem = hRoot.FirstChildElement("subimage").Element()) != 0 ) { // skip if not present (for compatibility with previous versions)
	//	int value;
	//	std::stringstream str;
	//	if ( pelem->QueryIntAttribute("resolution", &value) == TIXML_SUCCESS ) {
	//		if ( value ) {// additional parameters are not needed if resolution is zero
	//			additionalIOPluginParams = true;
	//			str << value;
	//			active_res = str.str();
	//		}
	//	}
	//	if ( pelem->QueryIntAttribute("timepoint", &value) == TIXML_SUCCESS ) {
	//		if ( value ) { // additional parameters are not needed if resolution is zero
	//			additionalIOPluginParams = true;
	//			str.str("");
	//			str << value;
	//			active_tp = str.str();
	//		}
	//	}
	//	const char *series_no_flag=pelem->Attribute("series_no");
	//	if ( series_no_flag ) {
	//		if ( strcmp(series_no_flag,"true") == 0 ) 
	//			series_no = additionalIOPluginParams = true;
	//	}
	//}

	pelem = hRoot.FirstChildElement("mechanical_displacements").Element();
	pelem->QueryFloatAttribute("V", &MEC_V);
	pelem->QueryFloatAttribute("H", &MEC_H);
	pelem = hRoot.FirstChildElement("dimensions").Element();
	int nrows, ncols, nslices;
	pelem->QueryIntAttribute("stack_rows", &nrows);
	pelem->QueryIntAttribute("stack_columns", &ncols);
	N_ROWS = nrows;
	N_COLS = ncols;
	pelem->QueryIntAttribute("stack_slices", &nslices);
	N_SLICES = nslices;

	pelem = hRoot.FirstChildElement("SUBVOLUMES").Element();
	pelem->QueryIntAttribute("N_SUBVOLUMES", &N_SUBVOLS);
	pelem->QueryIntAttribute("ENABLED_SUBVOLUME", &enabledSV);

	const char *strvalue = pelem->Attribute("ALIGNED");
	if ( strcmp(strvalue,"true") == 0 )
		aligned = true;
	else
		aligned = false;

	subvolumes = new vm::VirtualVolume *[N_SUBVOLS];
	sv_format = new string[N_SUBVOLS];
	xml_file_names = new string[N_SUBVOLS];

	//for each entry creates a VirtualVolume
	TiXmlElement * pelem2;
	pelem2 = pelem->FirstChildElement("Subvolume");
	for( int i=0; i<N_SUBVOLS; i++ ) {
		xml_file_names[i] = pelem2->Attribute("xml_fname");
		subvolumes[i] = volumemanager::VirtualVolumeFactory::createFromXML(xml_file_names[i].c_str(),false);
		sv_format[i] = subvolumes[i]->getVolumeFormat(xml_file_names[i].c_str());
		pelem2 = pelem2->NextSiblingElement();
	}
}

void MCVolume::saveXML(const char *xml_filename, const char *xml_filepath) throw (iom::exception)
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin MCVolume::saveXML(char *xml_filename = %s)\n", xml_filename);
	#endif

	//LOCAL VARIABLES
	char xml_abs_path[S_STATIC_STRINGS_SIZE];
	TiXmlDocument xml;
	TiXmlElement * root;
	TiXmlElement * pelem;

    //obtaining XML absolute path
    if(xml_filename)
        sprintf(xml_abs_path, "%s/%s.xml", stacks_dir, xml_filename);
    else if(xml_filepath)
        strcpy(xml_abs_path, xml_filepath);
    else
        throw iom::exception("in MCVolume::saveXML(...): no xml path provided");

	//initializing XML file with DTD declaration
	fstream XML_FILE(xml_abs_path, ios::out);
	XML_FILE<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>"<<endl;
	XML_FILE<<"<!DOCTYPE TeraStitcher SYSTEM \"TeraStitcher.DTD\">"<<endl;
	XML_FILE.close();

	//loading previously initialized XML file 
    if(!xml.LoadFile(xml_abs_path))
	{
		char errMsg[5000];
                sprintf(errMsg, "in MCVolume::saveToXML(...) : unable to load xml file at \"%s\"", xml_abs_path);
		throw iom::exception(errMsg);
	}

	//inserting root node <TeraStitcher> and children nodes
	root = new TiXmlElement("TeraStitcher");  
	xml.LinkEndChild( root );  

	// 2014-09-10. Alessandro. @ADDED 'volume_format' attribute to <TeraStitcher> XML node
	root->SetAttribute("volume_format", id.c_str());

	// 2017-04-27. Giulio. ADDED 'input_plugin' attribute to <TeraStitcher> XML node
	root->SetAttribute("input_plugin", "MultiVolume");

	pelem = new TiXmlElement("subvolumes_dir");
	pelem->SetAttribute("value", stacks_dir);
	root->LinkEndChild(pelem);
	// 2014-11-06. Giulio. @ADDED saved reference system into XML file
	pelem = new TiXmlElement("ref_sys");
	pelem->SetAttribute("ref1", reference_system.first);
	pelem->SetAttribute("ref2", reference_system.second);
	pelem->SetAttribute("ref3", reference_system.third);
	root->LinkEndChild(pelem);
	pelem = new TiXmlElement("voxel_dims");
	pelem->SetDoubleAttribute("V", VXL_V);
	pelem->SetDoubleAttribute("H", VXL_H);
	pelem->SetDoubleAttribute("D", VXL_D);
	root->LinkEndChild(pelem);
	pelem = new TiXmlElement("origin");
	pelem->SetDoubleAttribute("V", ORG_V);
	pelem->SetDoubleAttribute("H", ORG_H);
	pelem->SetDoubleAttribute("D", ORG_D);
	root->LinkEndChild(pelem);

	// 2016-11-27. Giulio. @ADDED additional parameters attributes must be saved in the xml
	if ( additionalIOPluginParams ) {
		pelem = new TiXmlElement("subimage");
		pelem->SetAttribute("resolution", active_res.c_str());
		pelem->SetAttribute("timepoint", active_tp.c_str());
		pelem->SetAttribute("series_no", (series_no ? "true" : "false"));
		root->LinkEndChild(pelem);
	}

	pelem = new TiXmlElement("mechanical_displacements");
	pelem->SetDoubleAttribute("V", MEC_V);
	pelem->SetDoubleAttribute("H", MEC_H);
	root->LinkEndChild(pelem);
	pelem = new TiXmlElement("dimensions");
	pelem->SetAttribute("stack_rows", N_ROWS);
	pelem->SetAttribute("stack_columns", N_COLS);
	pelem->SetAttribute("stack_slices", N_SLICES);
	root->LinkEndChild(pelem);

	//inserting subvolumes
	pelem = new TiXmlElement("SUBVOLUMES");
	pelem->SetAttribute("N_SUBVOLUMES", N_SUBVOLS);
	pelem->SetAttribute("ENABLED_SUBVOLUME", enabledSV);
	pelem->SetAttribute("ALIGNED", aligned ? "true" : "false");
	for( int i=0; i<N_SUBVOLS; i++) {
		TiXmlElement *pelem2 = new TiXmlElement("Subvolume");
		pelem2->SetAttribute("xml_fname", xml_file_names[i].c_str());
		pelem->LinkEndChild(pelem2);
	}
	root->LinkEndChild(pelem);
	
	//saving the file
	xml.SaveFile();
}

//counts the total number of displacements and the number of displacements per stack
void MCVolume::countDisplacements(int& total, float& per_stack_pair)
{
 //   /* PRECONDITIONS: none */
 //   total = 0;
	//per_stack_pair = 0.0f;
 //   for(int i=0; i<N_ROWS; i++)
 //       for(int j=0; j<N_COLS; j++)
 //       {
 //           total+= static_cast<int>(BLOCKS[i][j]->getEAST().size());
 //           total+= static_cast<int>(BLOCKS[i][j]->getSOUTH().size());
 //           per_stack_pair += static_cast<int>(BLOCKS[i][j]->getEAST().size());
 //           per_stack_pair += static_cast<int>(BLOCKS[i][j]->getSOUTH().size());
 //       }
 //   per_stack_pair /= 2*(N_ROWS*N_COLS) - N_ROWS - N_COLS;
}

//counts the number of single-direction displacements having a reliability measure above the given threshold
void MCVolume::countReliableSingleDirectionDisplacements(float threshold, int& total, int& reliable)
{
    /* PRECONDITIONS:
     *   - for each pair of adjacent stacks one and only one displacement exists (CHECKED) */

    //total = reliable = 0;
    //for(int i=0; i<N_ROWS; i++)
    //    for(int j=0; j<N_COLS; j++)
    //    {
    //        if(j != (N_COLS-1) && BLOCKS[i][j]->getEAST().size()==1)
    //        {
    //            total+=3;
    //            reliable += BLOCKS[i][j]->getEAST()[0]->getReliability(dir_vertical) >= threshold;
    //            reliable += BLOCKS[i][j]->getEAST()[0]->getReliability(dir_horizontal) >= threshold;
    //            reliable += BLOCKS[i][j]->getEAST()[0]->getReliability(dir_depth) >= threshold;
    //        }
    //        if(i != (N_ROWS-1) && BLOCKS[i][j]->getSOUTH().size()==1)
    //        {
    //            total+=3;
    //            reliable += BLOCKS[i][j]->getSOUTH()[0]->getReliability(dir_vertical) >= threshold;
    //            reliable += BLOCKS[i][j]->getSOUTH()[0]->getReliability(dir_horizontal) >= threshold;
    //            reliable += BLOCKS[i][j]->getSOUTH()[0]->getReliability(dir_depth) >= threshold;
    //        }
    //    }
}

//counts the number of stitchable stacks given the reliability threshold
int MCVolume::countStitchableStacks(float threshold)
{
    /* PRECONDITIONS:
     *   - for each pair of adjacent stacks one and only one displacement exists (CHECKED) */

    //stitchable stacks are stacks that have at least one reliable single-direction displacement
    int stitchables = 0;
    //bool stitchable;
    //for(int i=0; i<N_ROWS; i++)
    //    for(int j=0; j<N_COLS; j++)
    //    {
    //        stitchable = false;
    //        Block* stk = BLOCKS[i][j];
    //        if(i!= 0 && BLOCKS[i][j]->getNORTH().size()==1)
    //            for(int k=0; k<3; k++)
    //                stitchable = stitchable || (stk->getNORTH()[0]->getReliability(direction(k)) >= threshold);
    //        if(j!= (N_COLS -1) && BLOCKS[i][j]->getEAST().size()==1)
    //            for(int k=0; k<3; k++)
    //                stitchable = stitchable || (stk->getEAST()[0]->getReliability(direction(k)) >= threshold);
    //        if(i!= (N_ROWS -1) && BLOCKS[i][j]->getSOUTH().size()==1)
    //            for(int k=0; k<3; k++)
    //                stitchable = stitchable || (stk->getSOUTH()[0]->getReliability(direction(k)) >= threshold);
    //        if(j!= 0 && BLOCKS[i][j]->getWEST().size()==1)
    //            for(int k=0; k<3; k++)
    //                stitchable = stitchable || (stk->getWEST()[0]->getReliability(direction(k)) >= threshold);
    //        stitchables += stitchable;
    //    }
    return stitchables;
}

void MCVolume::releaseBuffers() {
	for ( int sv=0; sv<N_SUBVOLS; sv++ ) 
		subvolumes[sv]->releaseBuffers();
}
