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
* 2019-10-22.  Giulio.     @CREATED 
*/


#include "MultiCycleVolume.h"
//#include "VirtualVolume.h"
#include "TiledVolume.h"
#include "UnstitchedVolume.h"
#include "MultiSliceVolume.h"

//#include "../stitcher/DisplacementMIPNCC.h"

#include "volumemanager.config.h"

#include <fstream>
#include <sstream>
#include <list>
#include <string.h>

#ifdef _WIN32
#include "dirent_win.h"
#else
#include <dirent.h>
#endif

using namespace std;
using namespace iim;

MultiCycleVolume::MultiCycleVolume ( string _cycles_dir, float _norm_factor_D )  : 
ComposedVolume() {

	root_dir = new char[_cycles_dir.size() + 1];
	strcpy(root_dir,_cycles_dir.c_str());

	// in case init does not initialize them
	cycle_formats         = (std::string *) 0;
	cycles_new_xml_fnames = (std::string *) 0;

	init();

    VXL_V = VXL_1 = CYCLES[0]->getVXL_V();
	VXL_H = VXL_2 = CYCLES[0]->getVXL_H();
	VXL_D = VXL_3 = CYCLES[0]->getVXL_D();

    ORG_V = CYCLES[0]->getORG_V();
	ORG_H = CYCLES[0]->getORG_H();
	ORG_D = CYCLES[0]->getORG_D();	

	DIM_V = CYCLES[0]->getDIM_V();
	for ( int i=1; i<N_CYCLES; i++ ) {
		if ( CYCLES[i]->getDIM_V() != DIM_V) 
			throw iom::exception(vm::strprintf("in MultiCycleVolume::MultiCycleVolume(): cycles have different V size [%d (cycle 0) <> %d (cycle %d)]", DIM_V, CYCLES[i]->getDIM_V(),i).c_str());
	}
	DIM_H = CYCLES[0]->getDIM_H();
	for ( int i=1; i<N_CYCLES; i++ ) {
		if ( CYCLES[i]->getDIM_H() != DIM_H) 
			throw iom::exception(vm::strprintf("in MultiCycleVolume::MultiCycleVolume(): cycles have different H size [%d (cycle 0) <> %d (cycle %d)]", DIM_H, CYCLES[i]->getDIM_H(),i).c_str());
	}
	DIM_D = CYCLES[0]->getDIM_D();
	for ( int i=1; i<N_CYCLES; i++ ) {
		if ( CYCLES[i]->getDIM_D() != DIM_D) 
			throw iom::exception(vm::strprintf("in MultiCycleVolume::MultiCycleVolume(): cycles have different D size [%d (cycle 0) <> %d (cycle %d)]", DIM_D, CYCLES[i]->getDIM_D(),i).c_str());
	}

	normal_factor_D = _norm_factor_D;
	//cut_depth = _cut_depth;

	// compute the nominal dimension and nominal displacements along D
	// compute also nominal coords of layers
	cycles_coords = new VHD_coords[N_CYCLES];
	memset(cycles_coords,0,N_CYCLES*sizeof(VHD_coords));

	disps = new vector<Displacement *> *[N_CYCLES];
	for ( int i=0; i<N_CYCLES; i++ ) {
		disps[i] = (vector<Displacement *> *) 0;
	}
	
	DIM_C = 0;
	for ( int i=0; i<N_CYCLES; i++ ) {
		DIM_C += CYCLES[i]->getDIM_C();
	}

    BYTESxCHAN = CYCLES[0]->getBYTESxCHAN();  
	for ( int i=1; i<N_CYCLES; i++ ) {
		if ( CYCLES[i]->getBYTESxCHAN() != BYTESxCHAN) 
			throw iom::exception(vm::strprintf("in MultiCycleVolume::MultiCycleVolume(): cycles have different number of bytes per channel [%d (cycle 0) <> %d (cycle %d)]", BYTESxCHAN, CYCLES[i]->getBYTESxCHAN(),i).c_str());
	}

	initChannels();
}

MultiCycleVolume::MultiCycleVolume ( const char *xml_filepath ) 
: ComposedVolume() {
    //extracting <stacks_dir> field from XML
    TiXmlDocument xml;
    if(!xml.LoadFile(xml_filepath))
    {
        char errMsg[2000];
        sprintf(errMsg,"in MultiCycleVolume::MultiCycleVolume(_cycles_dir = \"%s\") : unable to load xml", xml_filepath);
        throw IOException(errMsg);
    }
    TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));
    TiXmlElement * pelem = hRoot.FirstChildElement("root_dir").Element();
    this->root_dir = new char[strlen(pelem->Attribute("value"))+1];
    strcpy(this->root_dir, pelem->Attribute("value"));
	xml.Clear();

	cycles_new_xml_fnames = (std::string *) 0; // in case InitFromXML does not initialize it

	// load xml content and generate mdata.bin
	initFromXML(xml_filepath);

	initChannels();
}

MultiCycleVolume::~MultiCycleVolume ( )  {
	if ( CYCLES ) {
		for ( int i=0; i<N_CYCLES; i++ ) 
			if ( CYCLES[i] )
				delete CYCLES[i];
		delete []CYCLES;
	}

	if ( cycle_formats ) {
		delete []cycle_formats;
	}

	if ( cycles_coords )
		delete []cycles_coords;

	//if ( nominal_D_overlap )
	//	delete []nominal_D_overlap;

	if ( disps ) {
		for ( int i=0; i<N_CYCLES; i++ ) 
			if ( disps[i] )
				delete disps[i];
		delete []disps;
	}
	
	if ( cycles_new_xml_fnames )
		delete []cycles_new_xml_fnames;
}

void MultiCycleVolume::init ( ) {

	/**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	//LOCAL VARIABLES
	string tmp_path;				//string that contains temp paths during computation
    string tmp;						//string that contains temp data during computation
	DIR *cur_dir_lev1;				//pointer to DIR, the data structure that represents a DIRECTORY (level 1 of hierarchical structure)
	DIR *cur_dir_lev2;				//pointer to DIR, the data structure that represents a DIRECTORY (level 2 of hierarchical structure)
	dirent *entry_lev1;				//pointer to DIRENT, the data structure that represents a DIRECTORY ENTRY inside a directory (level 1)
	dirent *entry_lev2;				//pointer to DIRENT, the data structure that represents a DIRECTORY ENTRY inside a directory (level 2)
	int i=0, j=0;					//for counting of N_ROWS, N_COLS
 //   list<Block*> blocks_list;                       //each stack found in the hierarchy is pushed into this list
    list<string> entries_lev1;                      //list of entries of first level of hierarchy
    list<string>::iterator entry_i;                 //iterator for list 'entries_lev1'
    list<string> entries_lev2;                      //list of entries of second level of hierarchy
    list<string>::iterator entry_j;                 //iterator for list 'entries_lev2'
    char block_i_j_path[STATIC_STRINGS_SIZE];

	//obtaining DIR pointer to root_dir (=NULL if directory doesn't exist)
	if (!(cur_dir_lev1=opendir(root_dir)))
	{
        char msg[STATIC_STRINGS_SIZE];
        sprintf(msg,"in MultiCycleVolume::init(...): Unable to open directory \"%s\"", root_dir);
        throw IOException(msg);
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

	N_CYCLES = (int) entries_lev1.size();
	CYCLES = new VirtualVolume *[N_CYCLES];
	cycle_formats = new std::string[N_CYCLES];
	cycles_new_xml_fnames = new std::string[N_CYCLES];


	// check if volumes are stitiched on unstitched
	if ( strstr(entries_lev1.front().c_str(),".xml") == 0 ) { // volume is stitched and not composed
		// WARNING: this code has not be checked - it is likely wrong!
		//for each entry creates a VirtualVolume
		for(entry_i = entries_lev1.begin(), i=0; entry_i!= entries_lev1.end(); entry_i++, i++)
		{
			//building absolute path of first level entry
			tmp_path=root_dir;
			tmp_path.append("/");
			tmp_path.append(*entry_i);
			cur_dir_lev2 = opendir(tmp_path.c_str());
			if (!cur_dir_lev2)
				throw IOException("in MultiCycleVolume::init(...): A problem occurred during scanning of subdirectories");

			//scanning second level of hierarchy, actually just one entry ("RES(...)" directory)
			while ((entry_lev2=readdir(cur_dir_lev2)))
			{
				tmp=entry_lev2->d_name;
				if(tmp.find(".") == string::npos && tmp.find(" ") == string::npos)
					entries_lev2.push_back(entry_lev2->d_name);
			}
			closedir(cur_dir_lev2);
			entries_lev2.sort();

			//for each entry of the second level, allocating a new Block
			for(entry_j = entries_lev2.begin(), j=0; entry_j!= entries_lev2.end(); entry_j++, j++)
			{
				//allocating new layer
				sprintf(block_i_j_path,"%s/%s/%s",root_dir,(*entry_i).c_str(), (*entry_j).c_str());
				CYCLES[i] = new TiledVolume(block_i_j_path); //,iim::ref_sys(iim::axis(1),iim::axis(2),iim::axis(3)),(float)1.0,(float)1.0,(float)1.0);
				cycle_formats[i] = TILED_TIF3D_FORMAT;
				reference_system = ((TiledVolume *) CYCLES[0])->getREF_SYS();
			}
			cycles_new_xml_fnames[i] = tmp_path;
			entries_lev2.clear();
		
		}
		reference_system = ((TiledVolume *) CYCLES[0])->getREF_SYS();
	}
	else { // volume is unstitched or composed

		//for each entry creates a VirtualVolume
		for(entry_i = entries_lev1.begin(), i=0; entry_i!= entries_lev1.end(); entry_i++, i++) {
			//building absolute path of first level entry
			tmp_path=root_dir;
			tmp_path.append("/");
			tmp_path.append(*entry_i);
			if ( strstr(getVolumeFormat(tmp_path).c_str(),"TiledXY") == 0 ) { // volume is composed
				if ( getVolumeFormat(tmp_path).compare(MULTISLICE_FORMAT) == 0 ) {
					CYCLES[i] = new MultiSliceVolume(tmp_path.c_str()); 
					cycle_formats[i] = MULTISLICE_FORMAT;
					reference_system = ((MultiSliceVolume *) CYCLES[0])->getREF_SYS(); // to be improved: assign multiple times the same reference system
				}
				else 
					throw iom::exception(vm::strprintf("in MultiCycleVolume::init(): unknown volume_format = \"%s\")", getVolumeFormat(tmp_path).c_str()).c_str());
			}
			else { // volume is unstitched
				CYCLES[i] = new UnstitchedVolume(tmp_path.c_str()); 
				cycle_formats[i] = UNST_TIF3D_FORMAT;
				reference_system = ((UnstitchedVolume *) CYCLES[0])->getREF_SYS(); // to be improved: assign multiple times the same reference system
			}
			cycles_new_xml_fnames[i] = tmp_path;
		}
	}

	entries_lev1.clear();
}


void MultiCycleVolume::initChannels ( ) 
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

    n_active = DIM_C;
    active = new iim::uint32[n_active];
    for ( int c=0; c<DIM_C; c++ )
        active[c] = c; // all cycles are assumed active
}


//void MultiCycleVolume::updateLayerCoords ( ) {
//
//	int offs;
//
//	DIM_V = CYCLES[0]->getDIM_V();
//	for ( int i=1; i<N_CYCLES; i++ ) {
//		DIM_V = (iim::uint32) ((CYCLES[i]->getDIM_V() > (int)DIM_V) ? CYCLES[i]->getDIM_V() : DIM_V);
//	}
//	DIM_H = CYCLES[0]->getDIM_H();
//	for ( int i=1; i<N_CYCLES; i++ ) {
//		DIM_H = (iim::uint32) ((CYCLES[i]->getDIM_H() > (int)DIM_H) ? CYCLES[i]->getDIM_H() : DIM_H);
//	}
//
//	DIM_D = 0;
//	for ( int i=0; i<(N_CYCLES - 1); i++ ) {
//		offs = (int) ROUND((CYCLES[i+1]->getORG_D() - CYCLES[i]->getORG_D()) * 1000.0F / VXL_D);
//		DIM_D += offs;
//		cycles_coords[i+1][2] = DIM_D; // other nominal coords are 0
//		//nominal_D_overlap[i] = CYCLES[i]->getDIM_D() - offs; // nominal overlap may become negative if offset is too large
//	}
//	DIM_D += CYCLES[N_CYCLES-1]->getDIM_D();
//}

int	MultiCycleVolume::getCYCLE_DIM(int i, int j) {
	if ( j==0 )
		return CYCLES[i]->getDIM_V();
	else if ( j==1 )
		return CYCLES[i]->getDIM_H();
	else if ( j==2 )
		return CYCLES[i]->getDIM_D();
	else
        throw IOException("in MultiCycleVolume::getCYCLE_DIM(...): invalid direction");
}


void MultiCycleVolume::initFromXML(const char *xml_filename)  {
	#if VM_VERBOSE > 3
    printf("\t\t\t\tin MultiCycleVolume::initFromXML(char *xml_filename = %s)\n", xml_filepath);
	#endif

	TiXmlDocument xml;
	if(!xml.LoadFile(xml_filename))
	{
		char errMsg[2000];
		sprintf(errMsg,"in MultiCycleVolume::initFromXML(xml_filepath = \"%s\") : unable to load xml", xml_filename);
		throw IOException(errMsg);
	}

	//setting ROOT element (that is the first child, i.e. <TeraStitcher> node)
	TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));

	// 2014-09-10. Alessandro. @ADDED 'volume_format' attribute to <TeraStitcher> XML node
	const char *volformat = hRoot.ToElement()->Attribute("volume_format");
	if(volformat && strcmp(volformat, getPrintableFormat().c_str()) != 0)
		throw iom::exception(vm::strprintf("in MultiCycleVolume::initFromXML(): unsupported volume_format = \"%s\" (current format is \"%s\")", volformat, getPrintableFormat().c_str()).c_str());

	//reading fields
	TiXmlElement * pelem = hRoot.FirstChildElement("root_dir").Element();
	pelem = hRoot.FirstChildElement("voxel_dims").Element();
	pelem->QueryFloatAttribute("V", &VXL_V);
	pelem->QueryFloatAttribute("H", &VXL_H);
	pelem->QueryFloatAttribute("D", &VXL_D);
	pelem->QueryFloatAttribute("norm_factor_D", &normal_factor_D);
	//pelem->QueryFloatAttribute("cut_depth", &cut_depth);
	pelem = hRoot.FirstChildElement("origin").Element();
	pelem->QueryFloatAttribute("V", &ORG_V);
	pelem->QueryFloatAttribute("H", &ORG_H);
	pelem->QueryFloatAttribute("D", &ORG_D);
	pelem = hRoot.FirstChildElement("dimensions").Element();
	int dim_v, dim_h, dim_d;
	pelem->QueryIntAttribute("DIM_V", &dim_v);
	pelem->QueryIntAttribute("DIM_H", &dim_h);
	pelem->QueryIntAttribute("DIM_D", &dim_d);
	DIM_V = dim_v;
	DIM_H = dim_h;
	DIM_D = dim_d;
	pelem->QueryIntAttribute("DIM_C", &DIM_C);
	pelem->QueryIntAttribute("BYTESxCHAN", &BYTESxCHAN);
	pelem = hRoot.FirstChildElement("ref_sys").Element();
	pelem->QueryIntAttribute("ref1", (int *) &reference_system.first);
	pelem->QueryIntAttribute("ref2", (int *) &reference_system.second);
	pelem->QueryIntAttribute("ref3", (int *) &reference_system.third);

	pelem = hRoot.FirstChildElement("CYCLES").Element();
	pelem->QueryIntAttribute("value", &N_CYCLES);
	CYCLES = new VirtualVolume *[N_CYCLES];
	cycle_formats = new std::string[N_CYCLES];
	cycles_new_xml_fnames = new std::string[N_CYCLES];
	cycles_coords = new VHD_coords[N_CYCLES];
	
	pelem = pelem->FirstChildElement();
	int index;
	//char *volume_dir;
	for(int i = 0; i < N_CYCLES; i++)
	{
		pelem->QueryIntAttribute("INDEX", &index);
		//volume_dir = new char[strlen(pelem->Attribute("value"))+1];
		//strcpy(volume_dir, pelem->Attribute("value"));
		cycles_new_xml_fnames[index] = pelem->Attribute("value");

		// get cycle format
		cycle_formats[index] = pelem->Attribute("format");

		// check if volumes are stitiched on unstitched
		// 2019-10-23. Giulio. @CHANGED //if ( strstr(volume_dir,".xml") == 0 ) { // volumes are stitched
		if ( cycle_formats[index].compare(TILED_TIF3D_FORMAT) == 0 ) {
			CYCLES[index] = new TiledVolume(cycles_new_xml_fnames[index].c_str());
			reference_system = ((TiledVolume *) CYCLES[0])->getREF_SYS();
		}
		else if ( cycle_formats[index].compare(UNST_TIF3D_FORMAT) == 0 ) { // volumes are unstitched or composed
			CYCLES[index] = new UnstitchedVolume(cycles_new_xml_fnames[index].c_str());
			reference_system = ((UnstitchedVolume *) CYCLES[0])->getREF_SYS();
		}
		else if ( cycle_formats[index].compare(MULTISLICE_FORMAT) == 0 ) {
			CYCLES[index] = new MultiSliceVolume(cycles_new_xml_fnames[index].c_str());
			reference_system = ((MultiSliceVolume *) CYCLES[0])->getREF_SYS();
		}
		pelem->QueryIntAttribute("coord_V", &cycles_coords[i][0]);
		pelem->QueryIntAttribute("coord_H", &cycles_coords[i][1]);
		pelem->QueryIntAttribute("coord_D", &cycles_coords[i][2]);
		pelem = pelem->NextSiblingElement();
		//delete[] volume_dir;
	}

	pelem = hRoot.FirstChildElement("INTER_CYCLES").Element();
	int dummy_val;
	pelem->QueryIntAttribute("value", &dummy_val); // just to consume the value
	//nominal_D_overlap = new int[N_CYCLES - 1];
	disps = new vector<Displacement *> *[N_CYCLES];
	for ( int i=0; i<N_CYCLES; i++ ) {
		disps[i] = (vector<Displacement *> *) 0;
	}

	pelem = pelem->FirstChildElement();
	for(int i = 0; i < N_CYCLES; i++)
	{
		pelem->QueryIntAttribute("INDEX", &index);
		if ( strcmp(pelem->Attribute("disps"),"yes") == 0 ) {
			disps[i] = new vector<Displacement *>(N_CYCLES,(Displacement *) 0);
			TiXmlElement *pelem2 = pelem->FirstChildElement("chan"); // skip indices
			for (int ii=0; ii<N_CYCLES; ii++) {
				TiXmlElement *pelem3 = pelem2->FirstChildElement("Displacement");
				if ( pelem3 ) // 2019. Giulio. @ADDED to deal with cycles that have no displacements computed
					disps[i]->at(ii) = Displacement::getDisplacementFromXML(pelem3);
				pelem2 = pelem2->NextSiblingElement("chan"); // next indices 
			}
		}

		pelem = pelem->NextSiblingElement();
	}
}


void MultiCycleVolume::saveXML(const char *xml_filename, const char *xml_filepath)  {
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin MultiCycleVolume::saveXML(char *xml_filename = %s)\n", xml_filename);
	#endif

	//LOCAL VARIABLES
    char xml_abs_path[STATIC_STRINGS_SIZE];
	TiXmlDocument xml;
	TiXmlElement * root;
	TiXmlElement * pelem;
	int i;
	int ii;

	//obtaining XML absolute path
	if(xml_filename)
		sprintf(xml_abs_path, "%s/%s.xml", root_dir, xml_filename);
	else if(xml_filepath)
		strcpy(xml_abs_path, xml_filepath);
	else
		throw IOException("in MultiCycleVolume::saveXML(...): no xml path provided");

	//initializing XML file with DTD declaration
	fstream XML_FILE(xml_abs_path, ios::out);
	XML_FILE<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>"<<endl;
	XML_FILE<<"<!DOCTYPE TeraStitcher SYSTEM \"TeraStitcher.DTD\">"<<endl;
	XML_FILE.close();

	//loading previously initialized XML file 
	if(!xml.LoadFile(xml_abs_path))
	{
		char errMsg[5000];
		sprintf(errMsg, "in MultiCycleVolume::saveToXML(...) : unable to load xml file at \"%s\"", xml_abs_path);
		throw IOException(errMsg);
	}

	//inserting root node <TeraStitcher> and children nodes
	root = new TiXmlElement("TeraStitcher");  
	// 2014-09-10. Alessandro. @ADDED 'volume_format' attribute to <TeraStitcher> XML node
	root->SetAttribute("volume_format", getPrintableFormat().c_str());
	xml.LinkEndChild( root );  

	pelem = new TiXmlElement("root_dir");
	pelem->SetAttribute("value", root_dir);
	root->LinkEndChild(pelem);
	pelem = new TiXmlElement("voxel_dims");
	pelem->SetDoubleAttribute("V", VXL_V);
	pelem->SetDoubleAttribute("H", VXL_H);
	pelem->SetDoubleAttribute("D", VXL_D);
	pelem->SetDoubleAttribute("norm_factor_D", normal_factor_D);
	//pelem->SetDoubleAttribute("cut_depth", cut_depth);
	root->LinkEndChild(pelem);
	pelem = new TiXmlElement("origin");
	pelem->SetDoubleAttribute("V", ORG_V);
	pelem->SetDoubleAttribute("H", ORG_H);
	pelem->SetDoubleAttribute("D", ORG_D);
	root->LinkEndChild(pelem);
	pelem = new TiXmlElement("dimensions");
	pelem->SetAttribute("DIM_V", DIM_V);
	pelem->SetAttribute("DIM_H", DIM_H);
	pelem->SetAttribute("DIM_D", DIM_D);
	pelem->SetAttribute("DIM_C", DIM_C);
	pelem->SetAttribute("BYTESxCHAN", BYTESxCHAN);
	root->LinkEndChild(pelem);
	pelem = new TiXmlElement("ref_sys");
	pelem->SetAttribute("ref1", reference_system.first);
	pelem->SetAttribute("ref2", reference_system.second);
	pelem->SetAttribute("ref3", reference_system.third);
	root->LinkEndChild(pelem);

	//inserting layers nodes
	pelem = new TiXmlElement("CYCLES");
	pelem->SetAttribute("value", N_CYCLES);
	for(i=0; i<N_CYCLES; i++) {
		TiXmlElement * pelem2 = new TiXmlElement("CYCLE");
		pelem2->SetAttribute("INDEX", i);
		pelem2->SetAttribute("value", cycles_new_xml_fnames[i].c_str());
		//if ( cycles_new_xml_fnames ) 
		//	// use new xml files associated to layers have been generated
		//	pelem2->SetAttribute("value", cycles_new_xml_fnames[i].c_str());
		//else 
		//	pelem2->SetAttribute("value", CYCLES[i]->getROOT_DIR());
		pelem2->SetAttribute("format", cycle_formats[i].c_str());
		pelem2->SetAttribute("coord_V", cycles_coords[i][0]);
		pelem2->SetAttribute("coord_H", cycles_coords[i][1]);
		pelem2->SetAttribute("coord_D", cycles_coords[i][2]);
		pelem->LinkEndChild(pelem2);
	}
	root->LinkEndChild(pelem);

	//inserting inter-layer information (overlap, displacements)
	pelem = new TiXmlElement("INTER_CYCLES");
	pelem->SetAttribute("value", N_CYCLES);
	for(i=0; i<N_CYCLES; i++) {
		TiXmlElement * pelem2 = new TiXmlElement("INTER_CYCLE_INFO");
		pelem2->SetAttribute("INDEX", i);

		if ( disps[i] ) {
			for (ii=0; ii<disps[i]->size(); ii++ ) {
				TiXmlElement * pelem3 = new TiXmlElement("chan");
				if ( disps[i]->at(ii) ) { // 2019-09-25. Giulio. @ADDED do not add displacements if they have not be computed
					pelem3->LinkEndChild(disps[i]->at(ii)->getXML());
					//pelem3->SetAttribute("disp_V", disps[i]->at(ii).at(jj)->getDisplacement(dir_vertical));
					//pelem3->SetAttribute("disp_H", disps[i]->at(ii).at(jj)->getDisplacement(dir_horizontal));
					//pelem3->SetAttribute("disp_D", disps[i]->at(ii).at(jj)->getDisplacement(dir_depth));
				}
				pelem2->LinkEndChild(pelem3);
			}
			pelem2->SetAttribute("disps", "yes");
		}
		else {
			pelem2->SetAttribute("disps", "no");
		}

		pelem->LinkEndChild(pelem2);
	}
	root->LinkEndChild(pelem);

	//saving the file
	xml.SaveFile();
}


//void MultiCycleVolume::saveCyclesXML(const char *xml_filename, const char *xml_filepath)  {
//	#if VM_VERBOSE > 3
//	printf("\t\t\t\tin MultiCycleVolume::saveLayersXML(char *xml_filename = %s)\n", xml_filename);
//	#endif
//
//	//LOCAL VARIABLES
//    std::string xml_base_abs_path = "";
//
//	//obtaining XML absolute path
//	if(xml_filename)
//		xml_base_abs_path = xml_base_abs_path + root_dir + xml_filename;
//	else if(xml_filepath)
//		xml_base_abs_path = xml_filepath;
//	else
//		throw IOException("in MultiCycleVolume::saveLayersXML(...): no xml path provided");
//
//	// eliminate suffix if any
//	if ( xml_base_abs_path.find(".xml") == (xml_base_abs_path.size() - 4) )
//		//for ( int i=0; i<4; i++ )
//		//	xml_base_abs_path.pop_back();
//		xml_base_abs_path = xml_base_abs_path.substr(0,xml_base_abs_path.size() - 4);
//
//	// computing cycle directory names
//	int n_digits = 1;
//	int _N_CYCLES = N_CYCLES / 10;	
//	while ( _N_CYCLES ) {
//		n_digits++;
//		_N_CYCLES /= 10;
//	}
//
//	// create list of new xml file names associated to layers
//	cycles_new_xml_fnames = new std::string [N_CYCLES];
//
//	for ( int i=0; i<N_CYCLES; i++ ) {
//		std::stringstream xmlfile_num;
//		xmlfile_num.width(n_digits);
//		xmlfile_num.fill('0');
//		xmlfile_num << i;
//		cycles_new_xml_fnames[i] = xml_base_abs_path + "L" + xmlfile_num.str() + ".xml";
//		((UnstitchedVolume *) CYCLES[i])->volume->saveXML(0,cycles_new_xml_fnames[i].c_str());
//	}
//}


void MultiCycleVolume::initDISPS(int i, int _DIM_V, int _DIM_H) {
	if ( disps[i] != ((vector<Displacement *> *) 0) )
		delete disps[i];

	disps[i] = new vector<Displacement *>(N_CYCLES,(Displacement *) 0);
}


void MultiCycleVolume::insertDisplacement(int i, int k, Displacement *displacement)  {

	displacement->evalReliability(dir_vertical);
	displacement->evalReliability(dir_horizontal);
	displacement->evalReliability(dir_depth);

	// we assume that by default corresponding cycles in adjacent layers are aligned with respect to motorized stages coordinates 
	displacement->setDefaultV(0);	
	displacement->setDefaultH(0);
	displacement->setDefaultD(0);

	if ( disps[k]->at(i) ) 
		delete disps[k]->at(i);

	disps[k]->at(i) = displacement;
}


iim::uint8 *MultiCycleVolume::loadSubvolume(
		int V0,int V1, int H0, int H1, int D0, int D1, int *n_chans, int ret_type,
		iim::uint8 *buffer, int bufSize_V, int bufSize_H, int bufSize_D, int bufSize_C,
		int bufOffs_V, int bufOffs_H, int bufOffs_D, int bufOffs_C 
)  {

	//throw IOException("in MultiCycleVolume::loadSubvolume(...): not implemented yet");
	
    if ( (ret_type != iim::NATIVE_RTYPE) && (ret_type != iim::DEF_IMG_DEPTH) ) {
		// return type should be converted, but not to 8 bits per cycle
        char err_msg[STATIC_STRINGS_SIZE];
		sprintf(err_msg,"MultiCycleVolume::loadSubvolume: non supported return type (%d bits) - native type is %d bits",ret_type, 8*this->BYTESxCHAN); 
        throw IOException(err_msg);
	}
	
	//initializations
	V0 = (V0 == -1 ? 0	     : V0);
	V1 = (V1 == -1 ? DIM_V   : V1);
	H0 = (H0 == -1 ? 0	     : H0);
	H1 = (H1 == -1 ? DIM_H   : H1);
	D0 = (D0 == -1 ? 0		 : D0);
	D1 = (D1 == -1 ? DIM_D	 : D1);

	sint64 sbv_height = V1 - V0;
	sint64 sbv_width  = H1 - H0;
	sint64 sbv_depth  = D1 - D0;

	if ( buffer==0 ) { // a buffer must be allocated
		//allocation		
		buffer = new iim::uint8[sbv_height * sbv_width * sbv_depth * BYTESxCHAN * n_active];
		memset(buffer,0,(sbv_height * sbv_width * sbv_depth * BYTESxCHAN * n_active));
		bufSize_V = (int) sbv_height;
		bufSize_H = (int) sbv_width;
		bufSize_D = (int) sbv_depth;
		bufSize_C = (int) n_active;
	}

	*n_chans = 0;
	
	int n_chans_cycle; // number of channels of cycle 'cy'
	int na;            // number of active channel of cycle 'cy'
	int nr;            // number of channels returned by cycle 'cy'
	int cc = 0;        // current channel to be considered
	for ( int cy=0; cy<N_CYCLES; cy++, cc+=n_chans_cycle ) {
		// cc is the number of channels of the cycles that have been already processed
		// set active channels of cycle 'cy'
		n_chans_cycle = CYCLES[cy]->getDIM_C(); // get the number of channels of cycle 'cy'
		iim::uint32 *active_cycle = new iim::uint32[n_chans_cycle]; // allocate active channels vector for cycle 'cy'
		memset(active_cycle,0,(n_chans_cycle*sizeof(iim::uint32))); // set all channels as inactive
		na = 0;
		for ( int c=0; c<n_chans_cycle; c++ ) { // check if channel cc+c is active
			bool found = false;
			int ic = 0;
			while ( ic<n_active && !found )
				if ( active[ic] == cc+c ) // the channel is active
					found = true;
				else
					ic++;
			if ( found ) {
				active_cycle[na++] = c;
			}
		}
		// set cycles of components to be loaded
		CYCLES[cy]->setActiveChannels(active_cycle,na); // ownership of array 'active_cycle' is passed to CYCLES[cy] 
		
		// compute component subvolume to be loaded
		int csV0 = std::max<int>(0,V0 - cycles_coords[cy][0]);
		int csV1 = std::min<int>(V1 - cycles_coords[cy][0],CYCLES[cy]->getDIM_V());
		int csH0 = std::max<int>(0,H0 - cycles_coords[cy][1]);
		int csH1 = std::min<int>(H1 - cycles_coords[cy][1],CYCLES[cy]->getDIM_H());
		int csD0 = std::max<int>(0,D0 - cycles_coords[cy][2]); 
		int csD1 = std::min<int>(D1 - cycles_coords[cy][2],CYCLES[cy]->getDIM_D());;
		 
		// get subvolume
		iim::uint8 *componentSubvol = CYCLES[cy]->loadSubvolume_to_UINT8(csV0,csV1,csH0,csH1,csD0,csD1,&nr,ret_type);
		if ( na != nr )
			throw iom::exception(vm::strprintf("in MultiCycleVolume::loadSubvolume(): the number of channel returned by cyle %d (%d) differ from the number of active channels (%d)", cy, nr).c_str());
		*n_chans += nr;

		// compute offsets
		int dstBufOffs_V = bufOffs_V + std::max<int>(0,cycles_coords[cy][0] - V0);
		int dstBufOffs_H = bufOffs_H + std::max<int>(0,cycles_coords[cy][1] - H0);
		int dstBufOffs_D = bufOffs_C + std::max<int>(0,cycles_coords[cy][2] - D0);
		int dstBufOffs_C = bufOffs_C + cc; // cycles to be loaded have been set: no additional offset is needed
		 
		// copy slice 'cy' to buffer
		// copy componentSubvol to buffer
		copy_strided_data(buffer,dstBufOffs_V,(dstBufOffs_H*BYTESxCHAN),dstBufOffs_D,dstBufOffs_C,
							bufSize_V, (bufSize_H*BYTESxCHAN), bufSize_D, bufSize_C,
							componentSubvol, 0, (0*BYTESxCHAN), 0, 0,
							(csV1-csV0),((csH1-csH0)*BYTESxCHAN),(csD1-csD0),nr,
							(csV1-csV0),((csH1-csH0)*BYTESxCHAN),(csD1-csD0),nr); // data dimensions coincide with source buffer dimensions	

		delete []componentSubvol;
	} 

	return (iim::uint8 *) buffer;
}

