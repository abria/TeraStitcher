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
* 2019-10-19.  Giulio.     @CREATED 
*/


#include "MultiSliceVolume.h"
#include "VirtualVolume.h"
#include "TiledVolume.h"
#include "UnstitchedVolume.h"

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

MultiSliceVolume::MultiSliceVolume ( string _layers_dir, float _cut_depth, float _norm_factor_D ) : ComposedVolume() {

	int offs;

	root_dir = new char[_layers_dir.size() + 1];
	strcpy(root_dir,_layers_dir.c_str());

	init();

    VXL_V = VXL_1 = LAYERS[0]->getVXL_V();
	VXL_H = VXL_2 = LAYERS[0]->getVXL_H();
	VXL_D = VXL_3 = LAYERS[0]->getVXL_D();

    ORG_V = LAYERS[0]->getORG_V();
	ORG_H = LAYERS[0]->getORG_H();
	ORG_D = LAYERS[0]->getORG_D();	

	DIM_V = LAYERS[0]->getDIM_V();
	for ( int i=1; i<N_LAYERS; i++ ) {
		DIM_V = ((LAYERS[i]->getDIM_V() > (int)DIM_V) ? LAYERS[i]->getDIM_V() : DIM_V);
	}
	DIM_H = LAYERS[0]->getDIM_H();
	for ( int i=1; i<N_LAYERS; i++ ) {
		DIM_H = (iim::uint32) ((LAYERS[i]->getDIM_H() > (int)DIM_H) ? LAYERS[i]->getDIM_H() : DIM_H);
	}

	normal_factor_D = _norm_factor_D;
	cut_depth = _cut_depth;

	// compute the nominal dimension and nominal displacements along D
	// compute also nominal coords of layers
	layers_coords = new VHD_coords[N_LAYERS];
	nominal_D_overlap = new int[N_LAYERS-1];

	memset(layers_coords,0,N_LAYERS*sizeof(VHD_coords));
	DIM_D = 0;
	for ( int i=0; i<(N_LAYERS - 1); i++ ) {
		offs = (int) ROUND((LAYERS[i+1]->getORG_D() - LAYERS[i]->getORG_D()) * 1000.0F / VXL_D);
		DIM_D += offs;
		layers_coords[i+1][2] = DIM_D; // other nominal coords are 0
		nominal_D_overlap[i] = LAYERS[i]->getDIM_D() - offs; // nominal overlap may become negative if offset is too large
	}
	DIM_D += LAYERS[N_LAYERS-1]->getDIM_D();

	disps = new vector< vector<Displacement *> > *[N_LAYERS-1];
	for ( int i=0; i<(N_LAYERS - 1); i++ ) {
		disps[i] = (vector< vector<Displacement *> > *) 0;
	}
	
	iBest = new int[N_LAYERS-1];
	jBest = new int[N_LAYERS-1];
	memset(iBest,0,sizeof(int)*(N_LAYERS-1));
	memset(jBest,0,sizeof(int)*(N_LAYERS-1));

    DIM_C = LAYERS[0]->getDIM_C();	
    BYTESxCHAN = LAYERS[0]->getBYTESxCHAN();  

	layers_new_xml_fnames = (std::string *) 0;

	initChannels();
}

MultiSliceVolume::MultiSliceVolume ( const char *xml_filepath ) : ComposedVolume() {
    //extracting <stacks_dir> field from XML
    TiXmlDocument xml;
    if(!xml.LoadFile(xml_filepath))
    {
        char errMsg[2000];
        sprintf(errMsg,"in MultiSliceVolume::MultiSliceVolume(_layers_dir = \"%s\") : unable to load xml", xml_filepath);
        throw IOException(errMsg);
    }
    TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));
    TiXmlElement * pelem = hRoot.FirstChildElement("root_dir").Element();
    this->root_dir = new char[strlen(pelem->Attribute("value"))+1];
    strcpy(this->root_dir, pelem->Attribute("value"));
	xml.Clear();

	// load xml content and generate mdata.bin
	initFromXML(xml_filepath);

	initChannels();

	layers_new_xml_fnames = (std::string *) 0;
}

MultiSliceVolume::~MultiSliceVolume ( )  {
	if ( LAYERS ) {
		for ( int i=0; i<N_LAYERS; i++ ) 
			if ( LAYERS[i] )
				delete LAYERS[i];
		delete []LAYERS;
	}

	if ( layers_coords )
		delete []layers_coords;

	if ( nominal_D_overlap )
		delete []nominal_D_overlap;

	if ( disps ) {
		for ( int i=0; i<(N_LAYERS - 1); i++ ) 
			if ( disps[i] )
				delete disps[i];
		delete disps;
	}
	
	if ( iBest )
		delete[] iBest;

	if ( jBest )
		delete[] jBest;

	if ( layers_new_xml_fnames ) 
		delete []layers_new_xml_fnames;
		
}

void MultiSliceVolume::init ( ) {

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
        sprintf(msg,"in MultiSliceVolume::init(...): Unable to open directory \"%s\"", root_dir);
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

	N_LAYERS = (int) entries_lev1.size();
	LAYERS = new VirtualVolume *[N_LAYERS];

	// check if volumes are stitiched on unstitched
	if ( strstr(entries_lev1.front().c_str(),".xml") == 0 ) { // volumes are stitched

		//for each entry creates a VirtualVolume
		for(entry_i = entries_lev1.begin(), i=0; entry_i!= entries_lev1.end(); entry_i++, i++)
		{
			//building absolute path of first level entry
			tmp_path=root_dir;
			tmp_path.append("/");
			tmp_path.append(*entry_i);
			cur_dir_lev2 = opendir(tmp_path.c_str());
			if (!cur_dir_lev2)
				throw IOException("in MultiSliceVolume::init(...): A problem occurred during scanning of subdirectories");

			//scanning second level of hierarchy, actuallt just one entry ("RES(...)" directory)
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
				LAYERS[i] = new TiledVolume(block_i_j_path); //,iim::ref_sys(iim::axis(1),iim::axis(2),iim::axis(3)),(float)1.0,(float)1.0,(float)1.0);
			}
			entries_lev2.clear();
		
		}
		reference_system = ((TiledVolume *) LAYERS[0])->getREF_SYS();
	}
	else { // volumes are unstitched
		//for each entry creates a VirtualVolume
		for(entry_i = entries_lev1.begin(), i=0; entry_i!= entries_lev1.end(); entry_i++, i++) {
			//building absolute path of first level entry
			tmp_path=root_dir;
			tmp_path.append("/");
			tmp_path.append(*entry_i);
			LAYERS[i] = new UnstitchedVolume(tmp_path.c_str()); 
		}
		reference_system = ((UnstitchedVolume *) LAYERS[0])->getREF_SYS();
	}

	entries_lev1.clear();
}


void MultiSliceVolume::initChannels ( ) 
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

    n_active = DIM_C;
    active = new iim::uint32[n_active];
    for ( int c=0; c<DIM_C; c++ )
        active[c] = c; // all channels are assumed active
}


void MultiSliceVolume::updateLayerCoords ( ) {

	int offs;

	DIM_V = LAYERS[0]->getDIM_V();
	for ( int i=1; i<N_LAYERS; i++ ) {
		DIM_V = (iim::uint32) ((LAYERS[i]->getDIM_V() > (int)DIM_V) ? LAYERS[i]->getDIM_V() : DIM_V);
	}
	DIM_H = LAYERS[0]->getDIM_H();
	for ( int i=1; i<N_LAYERS; i++ ) {
		DIM_H = (iim::uint32) ((LAYERS[i]->getDIM_H() > (int)DIM_H) ? LAYERS[i]->getDIM_H() : DIM_H);
	}

	DIM_D = 0;
	for ( int i=0; i<(N_LAYERS - 1); i++ ) {
		offs = (int) ROUND((LAYERS[i+1]->getORG_D() - LAYERS[i]->getORG_D()) * 1000.0F / VXL_D);
		DIM_D += offs;
		layers_coords[i+1][2] = DIM_D; // other nominal coords are 0
		nominal_D_overlap[i] = LAYERS[i]->getDIM_D() - offs; // nominal overlap may become negative if offset is too large
	}
	DIM_D += LAYERS[N_LAYERS-1]->getDIM_D();
}

int	MultiSliceVolume::getLAYER_DIM(int i, int j) {
	if ( j==0 )
		return LAYERS[i]->getDIM_V();
	else if ( j==1 )
		return LAYERS[i]->getDIM_H();
	else if ( j==2 )
		return LAYERS[i]->getDIM_D();
	else
        throw IOException("in MultiSliceVolume::getLAYER_DIM(...): invalid direction");
}


void MultiSliceVolume::initFromXML(const char *xml_filename)  {
	#if VM_VERBOSE > 3
    printf("\t\t\t\tin MultiSliceVolume::initFromXML(char *xml_filename = %s)\n", xml_filepath);
	#endif

	TiXmlDocument xml;
	if(!xml.LoadFile(xml_filename))
	{
		char errMsg[2000];
		sprintf(errMsg,"in MultiSliceVolume::initFromXML(xml_filepath = \"%s\") : unable to load xml", xml_filename);
		throw IOException(errMsg);
	}

	//setting ROOT element (that is the first child, i.e. <TeraStitcher> node)
	TiXmlHandle hRoot(xml.FirstChildElement("TeraStitcher"));

	// 2014-09-10. Alessandro. @ADDED 'volume_format' attribute to <TeraStitcher> XML node
	const char *volformat = hRoot.ToElement()->Attribute("volume_format");
	if(volformat && strcmp(volformat, getPrintableFormat().c_str()) != 0)
		throw iom::exception(vm::strprintf("in MultiChannelVolume::initFromXML(): unsupported volume_format = \"%s\" (current format is \"%s\")", volformat, getPrintableFormat().c_str()).c_str());

	//reading fields
	TiXmlElement * pelem = hRoot.FirstChildElement("root_dir").Element();
	pelem = hRoot.FirstChildElement("voxel_dims").Element();
	pelem->QueryFloatAttribute("V", &VXL_V);
	pelem->QueryFloatAttribute("H", &VXL_H);
	pelem->QueryFloatAttribute("D", &VXL_D);
	pelem->QueryFloatAttribute("norm_factor_D", &normal_factor_D);
	pelem->QueryFloatAttribute("cut_depth", &cut_depth);
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

	pelem = hRoot.FirstChildElement("LAYERS").Element();
	pelem->QueryIntAttribute("value", &N_LAYERS);
	LAYERS = new VirtualVolume *[N_LAYERS];
	layers_coords = new VHD_coords[N_LAYERS];
	
	pelem = pelem->FirstChildElement();
	int index;
	char *volume_dir;
	for(int i = 0; i < N_LAYERS; i++)
	{
		pelem->QueryIntAttribute("INDEX", &index);
		volume_dir = new char[strlen(pelem->Attribute("value"))+1];
		strcpy(volume_dir, pelem->Attribute("value"));

		// check if volumes are stitiched on unstitched
		if ( strstr(volume_dir,".xml") == 0 ) { // volumes are stitched
			LAYERS[index] = new TiledVolume(volume_dir);
			reference_system = ((TiledVolume *) LAYERS[0])->getREF_SYS();
		}
		else { // volumes are unstitched
			LAYERS[index] = new UnstitchedVolume(volume_dir);
			reference_system = ((UnstitchedVolume *) LAYERS[0])->getREF_SYS();
		}
		pelem->QueryIntAttribute("coord_V", &layers_coords[i][0]);
		pelem->QueryIntAttribute("coord_H", &layers_coords[i][1]);
		pelem->QueryIntAttribute("coord_D", &layers_coords[i][2]);
		pelem = pelem->NextSiblingElement();
		delete[] volume_dir;
	}

	pelem = hRoot.FirstChildElement("INTER_LAYERS").Element();
	int dummy_val;
	pelem->QueryIntAttribute("value", &dummy_val); // just to consume the value
	nominal_D_overlap = new int[N_LAYERS - 1];
	disps = new vector< vector<Displacement *> > *[N_LAYERS-1];
	for ( int i=0; i<(N_LAYERS - 1); i++ ) {
		disps[i] = (vector< vector<Displacement *> > *) 0;
	}
	iBest = new int[N_LAYERS-1];
	jBest = new int[N_LAYERS-1];

	pelem = pelem->FirstChildElement();
	for(int i = 0; i < (N_LAYERS - 1); i++)
	{
		pelem->QueryIntAttribute("INDEX", &index);
		pelem->QueryIntAttribute("nominal_overlap", &nominal_D_overlap[index]);
		if ( strcmp(pelem->Attribute("disps"),"yes") == 0 ) {
			int dim_i, dim_j;
			pelem->QueryIntAttribute("dim_i", &dim_i);
			pelem->QueryIntAttribute("dim_j", &dim_j);
			disps[i] = new vector< vector<Displacement *> >(dim_i,vector<Displacement *>(dim_j,(Displacement *) 0));
			TiXmlElement *pelem2 = pelem->FirstChildElement("tile"); // skip indices
			for (int ii=0; ii<dim_i; ii++) {
				for (int jj=0; jj<dim_j; jj++) {
					pelem2->QueryIntAttribute("i", &dummy_val);
					pelem2->QueryIntAttribute("j", &dummy_val);
					TiXmlElement *pelem3 = pelem2->FirstChildElement("Displacement");
					if ( pelem3 ) // 2019. Giulio. @ADDED to deal with tiles that have no displacements computed
						disps[i]->at(ii).at(jj) = Displacement::getDisplacementFromXML(pelem3);
					pelem2 = pelem2->NextSiblingElement("tile"); // next indices 
				}
			}
			pelem->QueryIntAttribute("iBest", &iBest[i]);
			pelem->QueryIntAttribute("jBest", &jBest[i]);
		}

		pelem = pelem->NextSiblingElement();
	}
}


void MultiSliceVolume::saveXML(const char *xml_filename, const char *xml_filepath)  {
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin MultiSliceVolume::saveXML(char *xml_filename = %s)\n", xml_filename);
	#endif

	//LOCAL VARIABLES
    char xml_abs_path[STATIC_STRINGS_SIZE];
	TiXmlDocument xml;
	TiXmlElement * root;
	TiXmlElement * pelem;
	int i;
	int ii,jj;

	//obtaining XML absolute path
	if(xml_filename)
		sprintf(xml_abs_path, "%s/%s.xml", root_dir, xml_filename);
	else if(xml_filepath)
		strcpy(xml_abs_path, xml_filepath);
	else
		throw IOException("in MultiSliceVolume::saveXML(...): no xml path provided");

	//initializing XML file with DTD declaration
	fstream XML_FILE(xml_abs_path, ios::out);
	XML_FILE<<"<?xml version=\"1.0\" encoding=\"UTF-8\"?>"<<endl;
	XML_FILE<<"<!DOCTYPE TeraStitcher SYSTEM \"TeraStitcher.DTD\">"<<endl;
	XML_FILE.close();

	//loading previously initialized XML file 
	if(!xml.LoadFile(xml_abs_path))
	{
		char errMsg[5000];
		sprintf(errMsg, "in MultiSliceVolume::saveToXML(...) : unable to load xml file at \"%s\"", xml_abs_path);
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
	pelem->SetDoubleAttribute("cut_depth", cut_depth);
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
	pelem = new TiXmlElement("LAYERS");
	pelem->SetAttribute("value", N_LAYERS);
	for(i=0; i<N_LAYERS; i++) {
		TiXmlElement * pelem2 = new TiXmlElement("LAYER");
		pelem2->SetAttribute("INDEX", i);
		if ( layers_new_xml_fnames ) 
			// use new xml files associated to layers have been generated
			pelem2->SetAttribute("value", layers_new_xml_fnames[i].c_str());
		else 
			pelem2->SetAttribute("value", LAYERS[i]->getROOT_DIR());
		pelem2->SetAttribute("coord_V", layers_coords[i][0]);
		pelem2->SetAttribute("coord_H", layers_coords[i][1]);
		pelem2->SetAttribute("coord_D", layers_coords[i][2]);
		pelem->LinkEndChild(pelem2);
	}
	root->LinkEndChild(pelem);

	//inserting inter-layer information (overlap, displacements)
	pelem = new TiXmlElement("INTER_LAYERS");
	pelem->SetAttribute("value", N_LAYERS-1);
	for(i=0; i<(N_LAYERS - 1); i++) {
		TiXmlElement * pelem2 = new TiXmlElement("INTER_LAYER_INFO");
		pelem2->SetAttribute("INDEX", i);
		pelem2->SetAttribute("nominal_overlap", nominal_D_overlap[i]);

		if ( disps[i] ) {
			for (ii=0; ii<disps[i]->size(); ii++ ) {
				for (jj=0; jj<disps[i]->at(ii).size(); jj++ ) {
					TiXmlElement * pelem3 = new TiXmlElement("tile");
					pelem3->SetAttribute("i", ii);
					pelem3->SetAttribute("j", jj);
					if ( disps[i]->at(ii).at(jj) ) { // 2019-09-25. Giulio. @ADDED do not add displacements if they have not be computed
						pelem3->LinkEndChild(disps[i]->at(ii).at(jj)->getXML());
						//pelem3->SetAttribute("disp_V", disps[i]->at(ii).at(jj)->getDisplacement(dir_vertical));
						//pelem3->SetAttribute("disp_H", disps[i]->at(ii).at(jj)->getDisplacement(dir_horizontal));
						//pelem3->SetAttribute("disp_D", disps[i]->at(ii).at(jj)->getDisplacement(dir_depth));
					}
					pelem2->LinkEndChild(pelem3);
				}
			}
			pelem2->SetAttribute("disps", "yes");
			pelem2->SetAttribute("dim_i", (int)disps[i]->size());
			pelem2->SetAttribute("dim_j", (int)disps[i]->at(0).size()); // all rows of tiles have the same number of tiles
			pelem2->SetAttribute("iBest", iBest[i]);
			pelem2->SetAttribute("jBest", jBest[i]); // all rows of tiles have the same number of tiles
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


//void MultiSliceVolume::saveLayersXML(const char *xml_filename, const char *xml_filepath)  {
//	#if VM_VERBOSE > 3
//	printf("\t\t\t\tin MultiSliceVolume::saveLayersXML(char *xml_filename = %s)\n", xml_filename);
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
//		throw IOException("in MultiSliceVolume::saveLayersXML(...): no xml path provided");
//
//	// eliminate suffix if any
//	if ( xml_base_abs_path.find(".xml") == (xml_base_abs_path.size() - 4) )
//		//for ( int i=0; i<4; i++ )
//		//	xml_base_abs_path.pop_back();
//		xml_base_abs_path = xml_base_abs_path.substr(0,xml_base_abs_path.size() - 4);
//
//	// computing channel directory names
//	int n_digits = 1;
//	int _N_LAYERS = N_LAYERS / 10;	
//	while ( _N_LAYERS ) {
//		n_digits++;
//		_N_LAYERS /= 10;
//	}
//
//	// create list of new xml file names associated to layers
//	layers_new_xml_fnames = new std::string [N_LAYERS];
//
//	for ( int i=0; i<N_LAYERS; i++ ) {
//		std::stringstream xmlfile_num;
//		xmlfile_num.width(n_digits);
//		xmlfile_num.fill('0');
//		xmlfile_num << i;
//		layers_new_xml_fnames[i] = xml_base_abs_path + "L" + xmlfile_num.str() + ".xml";
//		((UnstitchedVolume *) LAYERS[i])->volume->saveXML(0,layers_new_xml_fnames[i].c_str());
//	}
//}


void MultiSliceVolume::initDISPS(int i, int _DIM_V, int _DIM_H) {
	if ( disps[i] != ((vector< vector<Displacement *> > *) 0) )
		delete disps[i];

	disps[i] = new vector< vector<Displacement *> >(_DIM_V,vector<Displacement *>(_DIM_H,(Displacement *) 0));
}


void MultiSliceVolume::insertDisplacement(int i, int j, int k, Displacement *displacement)  {

	displacement->evalReliability(dir_vertical);
	displacement->evalReliability(dir_horizontal);
	displacement->evalReliability(dir_depth);

	// we assume that by default corresponding tiles in adjacent layers are aligned with respect to motorized stages coordinates 
	displacement->setDefaultV(0);	
	displacement->setDefaultH(0);
	displacement->setDefaultD(0);

	if ( disps[k]->at(i).at(j) ) 
		delete disps[k]->at(i).at(j);

	disps[k]->at(i).at(j) = displacement;
}


iim::uint8 *MultiSliceVolume::loadSubvolume(
		int V0,int V1, int H0, int H1, int D0, int D1, int *n_chans, int ret_type,
		iim::uint8 *buffer, int bufSize_V, int bufSize_H, int bufSize_D, int bufSize_C,
		int bufOffs_V, int bufOffs_H, int bufOffs_D, int bufOffs_C 
)  {

	//throw IOException("in MultiSliceVolume::loadSubvolume(...): not implemented yet");
	
    if ( (ret_type != iim::NATIVE_RTYPE) && (ret_type != iim::DEF_IMG_DEPTH) ) {
		// return type should be converted, but not to 8 bits per channel
        char err_msg[STATIC_STRINGS_SIZE];
		sprintf(err_msg,"MultiSliceVolume::loadSubvolume: non supported return type (%d bits) - native type is %d bits",ret_type, 8*this->BYTESxCHAN); 
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
	
	for ( int sl=D0; sl<D1; sl++ ) {
		// make a copy of array 'active'
		iim::uint32 *set_active = new iim::uint32[n_active];
		for ( int c=0; c<n_active; c++ )
			set_active[c] = active[c];
		// set channels of components to be loaded
		LAYERS[sl]->setActiveChannels(set_active,n_active); // ownership of array 'set_active' is passed to LAYERS[sl] 
		
		// compute component subvolume to be loaded
		int csV0 = std::max<int>(0,V0 - layers_coords[sl][0]);
		int csV1 = std::min<int>(V1 - layers_coords[sl][0],LAYERS[sl]->getDIM_V());
		int csH0 = std::max<int>(0,H0 - layers_coords[sl][1]);
		int csH1 = std::min<int>(H1 - layers_coords[sl][1],LAYERS[sl]->getDIM_H());
		int csD0 = std::max<int>(0,D0 - layers_coords[sl][2]); 
		int csD1 = std::min<int>(D1 - layers_coords[sl][2],LAYERS[sl]->getDIM_D());;
		 
		// get subvolume
		iim::uint8 *componentSubvol = LAYERS[sl]->loadSubvolume_to_UINT8(csV0,csV1,csH0,csH1,csD0,csD1,n_chans,ret_type);
		
		// compute offsets
		int dstBufOffs_V = bufOffs_V + std::max<int>(0,layers_coords[sl][0] - V0);
		int dstBufOffs_H = bufOffs_H + std::max<int>(0,layers_coords[sl][1] - H0);
		int dstBufOffs_D = bufOffs_C + std::max<int>(0,layers_coords[sl][2] - D0);
		int dstBufOffs_C = bufOffs_C; // channels to be loaded have been set: no additional offset is needed 
		 
		// copy slice 'sl' to buffer
		// copy componentSubvol to buffer
		copy_strided_data(buffer,dstBufOffs_V,(dstBufOffs_H*BYTESxCHAN),dstBufOffs_D,dstBufOffs_C,
							bufSize_V, (bufSize_H*BYTESxCHAN), bufSize_D, bufSize_C,
							componentSubvol, 0, (0*BYTESxCHAN), 0, 0,
							(csV1-csV0),((csH1-csH0)*BYTESxCHAN),(csD1-csD0),(*n_chans),
							(csV1-csV0),((csH1-csH0)*BYTESxCHAN),(csD1-csD0),(*n_chans)); // data dimensions coincide with source buffer dimensions	

		delete []componentSubvol;
	} 

	return (iim::uint8 *) buffer;
}

