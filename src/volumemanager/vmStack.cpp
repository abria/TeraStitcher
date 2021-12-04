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
* 2018-01-21. Giulio.     @ADDED in 'loadImageStack2' conditions to check if the requested region does not overlap with stack
* 2017-07-06. Giulio.     @ADDED method 'loadImageStack2' to enable selective reads of data
* 2017-04-26. Alessandro. @FIXED issue on Windows with Thumbs.db files that need to be ignored
* 2017-04-12. Giulio.     @ADDED check on a precondition of 'loadImageStack'
* 2016-11-14. Giulio.     @ADDED management of the case when z_end is invalid (i.e. when import is from an xml import file generated externally
* 2016-09-01. Giulio.     @ADDED cache management in loadImageStack (only if just one slice has to be loaded)
* 2015-08-24. Giulio.     @FIXED memory leak in loadImageStack
* 2015-02-28. Giulio.     @ADDED saving of fields N_CHANS and N_BYTESxCHAN in the xml files
* 2015-02-26. Giulio.     @ADDED initialization of fields N_CHANS and N_BYTESxCHAN in constructor; this information is NOT saved in mdata.bin and xml files
* 2015-01-17. Alessandro. @ADDED constructor for initialization from XML.
* 2015-01-17. Alessandro. @ADDED support for all-in-one-folder data (import from xml only).
* 2014-09-12. Alessandro. @FIXED 'init()' method to deal with non-empty tiles.
* 2014-09-05. Alessandro. @ADDED 'z_end' parameter in 'loadXML()' method to support sparse data feature.
* 2014-09-03. Alessandro. @ADDED 'Z_RANGES' attribute in the xml node.
* 2014-09-01. Alessandro. @FIXED 'binarizeInto()' and 'unbinarizeFrom()' to support sparse tiles.
* 2014-08-30. Alessandro. @ADDED 'hasValidImageData()' method.
* 2014-08-25. Alessandro. @REMOVED unused 'entry_k' variable declared in 'init()'.
*/

#ifdef _WIN32
    #include "dirent_win.h"
#else
    #include <dirent.h>
#endif
#include <iostream>
#include <string>
#include <list>
#include <sstream>
#include <boost/xpressive/xpressive.hpp>
#include "vmStackedVolume.h"
#include "vmStack.h"
#include "StackStitcher.h"
#include "Displacement.h"
#include "IOPluginAPI.h"

#include "vmCacheManager.h"

using namespace std;
using namespace vm;
using namespace iom;

//CONSTRUCTOR WITH ARGUMENTS
Stack::Stack(StackedVolume* _CONTAINER, int _ROW_INDEX, int _COL_INDEX, const char* _DIR_NAME) 
	: VirtualStack()
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack::Stack(StackedVolume* _CONTAINER, int _ROW_INDEX=%d, int _COL_INDEX=%d, char* _DIR_NAME=%s)\n",
		_ROW_INDEX, _COL_INDEX, _DIR_NAME);
	#endif

    CONTAINER = _CONTAINER;
    DIR_NAME = new char[strlen(_DIR_NAME)+1];
    strcpy(DIR_NAME, _DIR_NAME);
    ROW_INDEX = _ROW_INDEX;
    COL_INDEX = _COL_INDEX;

	init();
}

// 2015-01-17. Alessandro. @ADDED constructor for initialization from XML.
Stack::Stack(StackedVolume* _CONTAINER, int _ROW_INDEX, int _COL_INDEX, TiXmlElement* stack_node, int &z_end) 
	: VirtualStack()
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack::Stack(StackedVolume* _CONTAINER, int _ROW_INDEX=%d, int _COL_INDEX=%d, TiXmlElement*, int z_end=%d)\n",
		_ROW_INDEX, _COL_INDEX, z_end);
	#endif

	// check for valid stack node
	if(!stack_node)
		throw iom::exception("not an xml node", __iom__current__function__);
	if( strcmp(stack_node->ToElement()->Value(), "Stack") != 0)
		throw iom::exception(iom::strprintf("invalid xml node name: expected \"Stack\", found \"%s\"", stack_node->ToElement()->Value()), __iom__current__function__);

	CONTAINER = _CONTAINER;
	DIR_NAME = new char[strlen(stack_node->Attribute("DIR_NAME"))+1];
	strcpy(DIR_NAME, stack_node->Attribute("DIR_NAME"));
	ROW_INDEX = _ROW_INDEX;
	COL_INDEX = _COL_INDEX;  

	// first read image regex field (if any) from xml node
	readImgRegex(stack_node);

	// then scan folder for images
	init();

	// 2016-11-14. Giulio. @ADDED to manage the case when z_end has an invalid value and it must be initialized
	if ( z_end <= 0 ) { // invalid value
		z_end = DEPTH; 	// WARNING: this works if the volue is not sparse (check is not neeeded: it is done before calling the Stack constructor)
	}

	// finally load other xml attributes
	loadXML(stack_node, z_end);
}

Stack::Stack(StackedVolume* _CONTAINER, int _ROW_INDEX, int _COL_INDEX, FILE* bin_file) 
	: VirtualStack()
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack::Stack(StackedVolume* _CONTAINER, int _ROW_INDEX=%d, int _COL_INDEX=%d, FILE* bin_file)\n",
		_ROW_INDEX, _COL_INDEX);
	#endif

	CONTAINER = _CONTAINER;
	ROW_INDEX = _ROW_INDEX;
	COL_INDEX = _COL_INDEX;

	unBinarizeFrom(bin_file);
}

Stack::~Stack()
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack[%d,%d]::~Stack()\n",ROW_INDEX, COL_INDEX);
	#endif

	NORTH.clear();
	EAST.clear();
	SOUTH.clear();
	WEST.clear();

    if(FILENAMES)
    {
		for(int z=0; z<DEPTH; z++)
			if(FILENAMES[z])
				delete[] FILENAMES[z];
		delete[] FILENAMES;
	}
	if(STACKED_IMAGE)
		delete[] STACKED_IMAGE;
	if(DIR_NAME)
		delete[] DIR_NAME;
}

void Stack::init() 
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack[%d,%d]::init()\n",ROW_INDEX, COL_INDEX);
	#endif

	//LOCAL variables
	string tmp;
	DIR *cur_dir_lev3;
	dirent *entry_lev3;
    list<string> entries_lev3;
	string entry;

	//opening stack directory
	char abs_path[S_STATIC_STRINGS_SIZE];
	sprintf(abs_path,"%s/%s", CONTAINER->getSTACKS_DIR(), DIR_NAME);
	cur_dir_lev3 = opendir(abs_path);
	if (!cur_dir_lev3)
	{
		char errMsg[S_STATIC_STRINGS_SIZE];
		sprintf(errMsg, "in Stack::init(): can't open directory \"%s\"", abs_path);
		throw iom::exception(errMsg);
	}

	//scanning third level of hierarchy which entries need to be ordered alphabetically. This is done using STL.
	while ((entry_lev3=readdir(cur_dir_lev3)))
	{
		tmp = entry_lev3->d_name;

		// 2015-01-17. Alessandro. @ADDED support for all-in-one-folder data (import from xml only).
        if(img_regex.empty())
        {
			if(tmp.compare(".") != 0 && tmp.compare("..") != 0 && tmp.find(".") != string::npos 
				&& tmp != "Thumbs.db")	// @2017-04-26. Alessandro. @FIXED issue on Windows with Thumbs.db files that need to be ignored
                entries_lev3.push_back(tmp);
        }
        else
        {
            boost::xpressive::sregex rex = boost::xpressive::sregex::compile(img_regex.c_str());
            boost::xpressive::smatch what;
            if(boost::xpressive::regex_match(tmp, what, rex))
                entries_lev3.push_back(tmp);
        }
	}
	entries_lev3.sort();
    DEPTH = (int)entries_lev3.size();

    // close dir
	closedir(cur_dir_lev3);

	// 2014-08-30. Alessandro. @FIXED to deal with sparse data.
    // if stack is empty...
    if(DEPTH == 0)
	{
        // ...and SPARSE_DATA option is active, then exit
        if(vm::SPARSE_DATA)
            return;
        // ...otherwise throw an exception
        else
        {
            if(img_regex.empty())
                throw iom::exception(vm::strprintf("in Stack[%s]::init(): stack is empty. If that was your intent, please use the '--sparse_data' option", DIR_NAME).c_str());
            else
                throw iom::exception(vm::strprintf("in Stack[%s]::init(): no files found that match regular expression \"%s\"", DIR_NAME, img_regex.c_str()).c_str());
        }
	}

    // convert filenames_list (STL list of <string> objects) into FILENAMES (1-D array of C-strings)
    FILENAMES = new char*[DEPTH];
    for(int z=0; z<DEPTH; z++)
	{
		entry = entries_lev3.front();
		FILENAMES[z] = new char[entry.size()+1];
		strcpy(FILENAMES[z], entry.c_str());
		entries_lev3.pop_front();
	}
    entries_lev3.clear();

    // extract HEIGHT and WIDTH attributes from first slice
	iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->readMetadata(vm::strprintf("%s/%s/%s", CONTAINER->getSTACKS_DIR(), DIR_NAME, FILENAMES[0]), WIDTH, HEIGHT, N_BYTESxCHAN, N_CHANS);

	// 2014-09-12. Alessandro. @FIXED 'init()' method to deal with non-empty tiles.
	// add to 'z_ranges' the full range
	z_ranges.clear();
	z_ranges.push_back(vm::interval<int>(0, DEPTH));
}

//binarizing-unbinarizing methods
void Stack::binarizeInto(FILE* file) 
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack[%d,%d]::binarizeInto(...)\n",ROW_INDEX, COL_INDEX);
	#endif

	//LOCAL VARIABLES
	uint16 str_size=0;
	int i=0;

	fwrite(&HEIGHT, sizeof(int), 1, file);
	fwrite(&WIDTH, sizeof(int), 1, file);
	fwrite(&DEPTH, sizeof(int), 1, file);
	fwrite(&ABS_V, sizeof(int), 1, file);
	fwrite(&ABS_H, sizeof(int), 1, file);
	fwrite(&ABS_D, sizeof(int), 1, file);
	str_size = (uint16)strlen(DIR_NAME) + 1;
	fwrite(&str_size, sizeof(uint16), 1, file);
	fwrite(DIR_NAME, str_size, 1, file);
	for(i = 0; i < DEPTH; i++)
	{
		// 2014-09-01. Alessandro. @FIXED: added support for sparse tiles
		str_size = static_cast<uint16>( FILENAMES[i] ? strlen(FILENAMES[i]) + 1 : 0);
		fwrite(&str_size, sizeof(uint16), 1, file);
		if(FILENAMES[i])
			fwrite(FILENAMES[i], str_size, 1, file);
	}
}

void Stack::unBinarizeFrom(FILE* file) 
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack[%d,%d]::unBinarizeFrom(...)\n",ROW_INDEX, COL_INDEX);
	#endif

	//LOCAL VARIABLES
	uint16 str_size;
	int i;
	size_t fread_return_val;

	fread_return_val = fread(&HEIGHT, sizeof(int), 1, file);
	if(fread_return_val != 1) {
		fclose(file);
		throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
	}
	fread_return_val = fread(&WIDTH, sizeof(int), 1, file);
	if(fread_return_val != 1) {
		fclose(file);
		throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
	}
	fread_return_val = fread(&DEPTH, sizeof(int), 1, file);
	if(fread_return_val != 1) {
		fclose(file);
		throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
	}
	fread_return_val = fread(&ABS_V, sizeof(int), 1, file);
	if(fread_return_val != 1) {
		fclose(file);
		throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
	}
	fread_return_val = fread(&ABS_H, sizeof(int), 1, file);
	if(fread_return_val != 1) {
		fclose(file);
		throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
	}
	fread_return_val = fread(&ABS_D, sizeof(int), 1, file);
	if(fread_return_val != 1) {
		fclose(file);
		throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
	}
	fread_return_val = fread(&str_size, sizeof(uint16), 1, file);
	if(fread_return_val != 1) {
		fclose(file);
		throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
	}
	DIR_NAME = new char[str_size];
	fread_return_val = fread(DIR_NAME, str_size, 1, file);
	FILENAMES = new char*[DEPTH];
	for(i = 0; i < DEPTH; i++)
	{
		fread_return_val = fread(&str_size, sizeof(uint16), 1, file);
		if(fread_return_val != 1) {
			fclose(file);
			throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
		}		

		// 2014-09-01. Alessandro. @FIXED: added support for sparse tiles.
		FILENAMES[i] = str_size ? new char[str_size] : 0;		
		if(str_size)
		{
			fread_return_val = fread(FILENAMES[i], str_size, 1, file);
			if(fread_return_val != 1) {
				fclose(file);
				throw iom::exception("in Stack::unBinarizeFrom(...): error while reading binary metadata file");
			}
		}
	}

	// 2014-09-01. Alessandro. @ADDED support for sparse data.
	compute_z_ranges();
}

//XML methods
TiXmlElement* Stack::getXML()
{
	#if VM_VERBOSE > 3
	printf("......in Stack[%d,%d]::getXML()\n",ROW_INDEX, COL_INDEX);
	#endif

	TiXmlElement *xml_representation = new TiXmlElement("Stack");
	xml_representation->SetAttribute("N_CHANS",N_CHANS);
	xml_representation->SetAttribute("N_BYTESxCHAN",N_BYTESxCHAN);
	xml_representation->SetAttribute("ROW",ROW_INDEX);
	xml_representation->SetAttribute("COL",COL_INDEX);
	xml_representation->SetAttribute("ABS_V",ABS_V);
	xml_representation->SetAttribute("ABS_H",ABS_H);
	xml_representation->SetAttribute("ABS_D",ABS_D);
	xml_representation->SetAttribute("STITCHABLE",stitchable ? "yes" : "no");
	xml_representation->SetAttribute("DIR_NAME",DIR_NAME);

	// 2014-09-03. Alessandro. @ADDED 'Z_RANGES' attribute in the xml node
	std::string z_ranges_string;
	for(int k=0; k<z_ranges.size(); k++)
		z_ranges_string += vm::strprintf("[%d,%d)%s", z_ranges[k].start, z_ranges[k].end, k == z_ranges.size()-1 ? "" : ";");
	xml_representation->SetAttribute("Z_RANGES",z_ranges_string.c_str());
	writeImgRegex(xml_representation);

	vector<Displacement*>::iterator i;
	TiXmlElement *NORTH_displacements = new TiXmlElement("NORTH_displacements");
	for(i = NORTH.begin(); i != NORTH.end(); i++)
		NORTH_displacements->LinkEndChild((*i)->getXML());
	TiXmlElement *EAST_displacements = new TiXmlElement("EAST_displacements");
	for(i = EAST.begin(); i != EAST.end(); i++)
		EAST_displacements->LinkEndChild((*i)->getXML());
	TiXmlElement *SOUTH_displacements = new TiXmlElement("SOUTH_displacements");
	for(i = SOUTH.begin(); i != SOUTH.end(); i++)
		SOUTH_displacements->LinkEndChild((*i)->getXML());
	TiXmlElement *WEST_displacements = new TiXmlElement("WEST_displacements");
	for(i = WEST.begin(); i != WEST.end(); i++)
		WEST_displacements->LinkEndChild((*i)->getXML());
	xml_representation->LinkEndChild(NORTH_displacements);
	xml_representation->LinkEndChild(EAST_displacements);
	xml_representation->LinkEndChild(SOUTH_displacements);
	xml_representation->LinkEndChild(WEST_displacements);

	return xml_representation;
}

void Stack::loadXML(
	TiXmlElement *stack_node, 
	int z_end)					// 2014-09-05. Alessandro. @ADDED 'z_end' parameter to support sparse data feature
	//										   Here 'z_end' identifies the range [0, z_end) that slices can span

{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack[%d,%d]::loadXML(TiXmlElement *stack_node)\n",ROW_INDEX, COL_INDEX);
	#endif

	stack_node->QueryIntAttribute("N_CHANS",&N_CHANS);
	stack_node->QueryIntAttribute("N_BYTESxCHAN",&N_BYTESxCHAN);
	if ( N_CHANS == -1 ||N_BYTESxCHAN == -1 ) {
		iom::warning(iom::strprintf("old xml file: missing N_CHANS and N_BYTESxCHAN tags; class data members set to 1.").c_str(), __iom__current__function__);
		// data members set to 1 assuming that there is only one 8bits channel
		N_CHANS = 1;
		N_BYTESxCHAN = 1;
	}

	stack_node->QueryIntAttribute("ABS_V", &ABS_V);
	stack_node->QueryIntAttribute("ABS_H", &ABS_H);
	stack_node->QueryIntAttribute("ABS_D", &ABS_D);
	stitchable = strcmp(stack_node->Attribute("STITCHABLE"),"yes")==0 ? true : false;
	if(strcmp(stack_node->Attribute("DIR_NAME"), DIR_NAME) != 0)
	{
		char errMsg[2000];
		sprintf(errMsg, "in Stack[%s]::loadXML(...): Mismatch between xml file and %s in <DIR_NAME> field.", DIR_NAME, vm::BINARY_METADATA_FILENAME.c_str());
		throw iom::exception(errMsg);
	}

	// 2015-01-17. Alessandro. @ADDED support for all-in-one-folder data (import from xml only).
	readImgRegex(stack_node);

	// 2014-09-03. Alessandro. @ADDED 'Z_RANGES' attribute in the xml node
	const char* z_ranges_c = stack_node->Attribute("Z_RANGES");
	if( z_ranges_c )	// field is present: we can run 'sparse data' code
	{
		// nonempty tile
		if ( strlen(z_ranges_c) != 0 )
		{
			// parse 'Z_RANGES' field
			z_ranges.clear();
			std::string z_ranges_c_str = z_ranges_c;
			std::string z_ranges_string = vm::cls(z_ranges_c_str);
			std::vector<std::string> tokens;
			vm::split(z_ranges_string, ";", tokens);
			for(int i=0; i<tokens.size(); i++)
			{
				// ignore empty tokens
				if(tokens[i].empty())
					continue;

				// parse [start, end) interval
				std::vector<std::string> extremes;
				vm::split(tokens[i].substr(1, tokens[i].size()-2), ",", extremes);

				// check correct parsing
				if(extremes.size() != 2)
					throw iom::exception(vm::strprintf("in Stack(%s)::loadXML(): cannot parse 'Z_RANGES' subentry \"%s\"", DIR_NAME, tokens[i].c_str()).c_str());

				// get integral range
				int start = vm::str2num<int>(extremes[0]);
				int end   = vm::str2num<int>(extremes[1]);

				// check valid range
				if(start < 0 ||  start >= end ||  end > z_end)
					throw iom::exception(vm::strprintf("in Stack(%s)::loadXML(): 'Z_RANGES' subentry \"%s\" is out of range [%d,%d) ", DIR_NAME, tokens[i].c_str(), 0, z_end).c_str());

				// push range
				z_ranges.push_back(vm::interval<int>(start, end));
			}

			// check precondition: z_ranges should contain consecutive but not contiguous intervals [a_1,b_1), [a_2,b_2), ... such that a_n > b_(n-1)
			for(int i=1; i<z_ranges.size(); i++)
				if(z_ranges[i].start <= z_ranges[i-1].end)
					throw iom::exception(vm::strprintf("in Stack(%s)::loadXML(): wrong sequence in 'Z_RANGES' attribute. 'Z-RANGES' should contain "
					"consecutive but not contiguous intervals [a_1,b_1), [a_2,b_2), ... such that a_n > b_(n-1). "
					"Found a_%d(%d) <= b_%d(%d)", DIR_NAME, i, z_ranges[i].start, i-1, z_ranges[i-1].end).c_str());

			// if 'FILENAMES' is not a sparse list...
			if(DEPTH != z_end)
			{
				// make 'FILENAMES' a sparse list
				char **FILENAMES_sparse = new char*[z_end];
				for(int z=0; z<z_end; z++)
					FILENAMES_sparse[z] = 0;
				int i=0;
				for(int k=0; k<z_ranges.size(); k++)
					for(int z=z_ranges[k].start; z<z_ranges[k].end; z++)
					{
						// if we have used all the available slices and the current 'z_range' asks for others, throw an exception
						if(i >= DEPTH)
							throw iom::exception(vm::strprintf("in Stack(%s)::loadXML(): no more slices available to cover z-range [%d,%d)", 
							DIR_NAME, z_ranges[k].start, z_ranges[k].end).c_str());

						FILENAMES_sparse[z] = FILENAMES[i++];
					}
				
				// substitute list of filenames with its sparse version
				if(FILENAMES)
					delete[] FILENAMES;
				FILENAMES = FILENAMES_sparse;
				DEPTH = z_end;
			}
			// ...otherwise check if the sparse list matches with z_ranges
			else
			{
				for(int i=0; i<DEPTH; i++)
					if(FILENAMES[i] && isEmpty(i,i))
						throw iom::exception(vm::strprintf("in Stack(%s)::loadXML(): mismatch between sparse list loaded from mdata.bin and 'Z_RANGES' loaded from xml. Please delete mdata.bin and re-load the xml.", DIR_NAME));
					else if(FILENAMES[i] == 0 && isComplete(i,i))
						throw iom::exception(vm::strprintf("in Stack(%s)::loadXML(): mismatch between sparse list loaded from mdata.bin and 'Z_RANGES' loaded from xml. Please delete mdata.bin and re-load the xml.", DIR_NAME));
			}
		}
		// 'Z_RANGES' field is present, but is empty: assume empty tile
		else
		{
			// make 'FILENAMES' a sparse list
			if(FILENAMES)
				delete[] FILENAMES;
			FILENAMES = new char*[z_end];
			for(int z=0; z<z_end; z++)
				FILENAMES[z] = 0;
			DEPTH = z_end;
		}
	}
	// no 'Z_RANGES' field (old or incomplete XML)
	else
	{
		// if 'sparse data' option is active, 'Z_RANGES' must be provided (then throw an exception)
		if(SPARSE_DATA)
			throw iom::exception(vm::strprintf("in Stack::loadXML(): cannot find 'Z_RANGES' attribute in stack \"%s\" xml node ('sparse data' option is active). "
			"Please edit your XML file accordingly or disable 'sparse data' option.", DIR_NAME).c_str());
		
		// 'sparse data' option is not active and 'Z_RANGES' field is not present, but DEPTH does not match with z_end: throw an exception
		if(DEPTH != z_end)
			throw iom::exception(vm::strprintf("in Stack::loadXML(): stack %s has %d slices, but volume has %d. If that was your intent, please use "
			"'sparse data' option.", DIR_NAME, DEPTH, z_end).c_str());

		// all checks were successful: assume tile is not sparse
        z_ranges.clear();
        z_ranges.push_back(vm::interval<int>(0, DEPTH));
	}
	

	TiXmlElement *NORTH_displacements = stack_node->FirstChildElement("NORTH_displacements");
	for(TiXmlElement *displ_node = NORTH_displacements->FirstChildElement("Displacement"); displ_node; displ_node = displ_node->NextSiblingElement())
		NORTH.push_back(Displacement::getDisplacementFromXML(displ_node));
	TiXmlElement *EAST_displacements = stack_node->FirstChildElement("EAST_displacements");
	for(TiXmlElement *displ_node = EAST_displacements->FirstChildElement("Displacement"); displ_node; displ_node = displ_node->NextSiblingElement())
		EAST.push_back(Displacement::getDisplacementFromXML(displ_node));
	TiXmlElement *SOUTH_displacements = stack_node->FirstChildElement("SOUTH_displacements");
	for(TiXmlElement *displ_node = SOUTH_displacements->FirstChildElement("Displacement"); displ_node; displ_node = displ_node->NextSiblingElement())
		SOUTH.push_back(Displacement::getDisplacementFromXML(displ_node));
	TiXmlElement *WEST_displacements = stack_node->FirstChildElement("WEST_displacements");
	for(TiXmlElement *displ_node = WEST_displacements->FirstChildElement("Displacement"); displ_node; displ_node = displ_node->NextSiblingElement())
		WEST.push_back(Displacement::getDisplacementFromXML(displ_node));

}


//loads image stack from <first_file> to <last_file> extremes included, if not specified loads entire Stack
iom::real_t* Stack::loadImageStack(int first_file, int last_file) 
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack[%d,%d](%s, empty = %s)::loadImageStack(first_file = %d, last_file = %d)\n",ROW_INDEX, COL_INDEX, DIR_NAME, isEmpty() ? "true": "false", first_file, last_file);
	#endif

	uint64 image_stack_width=0, image_stack_height=0, z=0;
 	iom::real_t *raw_data;

	unsigned char *data = (unsigned char *) 0;
	int img_width;
	int img_height;
	int img_bytes_x_chan;
	int img_chans;

	if ( STACKED_IMAGE ) {
		delete[] STACKED_IMAGE;
    	STACKED_IMAGE=0;
		throw iom::exception(iom::strprintf("precondition violated in Stack [%d,%d] loading slices [%d,%d]: buffer is already allocated", ROW_INDEX, COL_INDEX,first_file, last_file), __iom__current__function__);
	}

	// check and adjust file selection
	first_file = (first_file == -1 ? 0 : first_file);
	last_file  = (last_file  == -1 ? DEPTH - 1 : last_file );
	first_file = min(first_file, DEPTH-1);
	last_file = min(last_file, DEPTH-1);

	if(isEmpty(first_file, last_file))
	{
		// allocate and initialize a black stack
		uint64 image_size = static_cast<uint64>(WIDTH) * static_cast<uint64>(HEIGHT) * static_cast<uint64>(last_file-first_file+1);
		STACKED_IMAGE = new iom::real_t[image_size];
		for(uint64 k=0; k<image_size; k++)
			STACKED_IMAGE[k] = 0;

		return STACKED_IMAGE;
	}

	if ( first_file == last_file ) { // only one slice has been requested
		if ( CONTAINER->getCACHEBUFFER()->getSlice(ROW_INDEX,COL_INDEX,first_file,data) ) {
			img_width        = WIDTH;
			img_height       = HEIGHT;
			img_bytes_x_chan = N_BYTESxCHAN;
			img_chans        = N_CHANS;

			/**** DUPLICATED CODE (BEGIN) ***/
			// initialize output data and image dimensions
			if(!STACKED_IMAGE)
			{
				image_stack_height = (uint64) img_height;
				image_stack_width  = (uint64) img_width;
				uint64 image_stack_size = image_stack_width * image_stack_height * (last_file - first_file + 1);
				STACKED_IMAGE = new iom::real_t[image_stack_size];
				for(uint64 j=0; j < image_stack_size; j++)
					STACKED_IMAGE[j] = 0;		// default is 0 (black)
			}
			else if(image_stack_width != (uint64)img_width || image_stack_height != (uint64)img_height)
				throw iom::exception(iom::strprintf("images in stack have not the same dimensions (%d x %d != %d x %d)", img_width, img_height, image_stack_width, image_stack_height), __iom__current__function__);

			int offset;

			// convert image to [0.0,1.0]-valued array
			if ( img_chans == 1 ) {
				raw_data = &STACKED_IMAGE[z*img_height*img_width];
				if ( img_bytes_x_chan == 1 ) {
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) data[i]/255.0f;
				}
				else { // img_bytes_x_chan == 2
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) ((uint16 *)data)[i]/65535.0f; // data must be interpreted as a uint16 array
				}
			}
			else { // conversion to an intensity image
				// test how channel are stored in the returned buffer 'data'
				if ( iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->isChansInterleaved() ) {
					if ( CHANS == iom::ALL ) {
						throw iom::exception("conversion from multi-channel to intensity images not supported.", __iom__current__function__);
					}
					else if ( CHANS == iom::R ) {
						offset = 0;
					}
					else if ( CHANS == iom::G ) {
						offset = 1;
					}
					else if ( CHANS == iom::B ) {
						offset = 2;
					}
					else {
						throw iom::exception("wrong value for parameter iom::CHANNEL_SELECTION.", __iom__current__function__);
					}
					raw_data = &STACKED_IMAGE[z*img_height*img_width];
					if ( img_bytes_x_chan == 1 )
						for(int i = 0; i < (img_height * img_width); i++)
							raw_data[i] = (iom::real_t) data[3*i + offset]/255.0f;
					else // img_bytes_x_chan == 2
						for(int i = 0; i < (img_height * img_width); i++)
							raw_data[i] = (iom::real_t) ((uint16 *)data)[3*i + offset]/65535.0f; // data must be interpreted as a uint16 array
				}
				else {
					throw iom::exception("2D image formats with non-interleaved channels are not supported yet.");
				}
			}
			/**** DUPLICATED CODE (END) ***/

			return STACKED_IMAGE;
		}
	}

	// build absolute image path of first file in the stack
	int i = first_file;
	while ( isSparse() && i<=last_file && FILENAMES[i]==0 ) // look for existing file
		i++;
	if ( i > last_file )
		throw iom::exception("unexpected exception: substack is entirely empty", __iom__current__function__);

	std::string image_path = std::string(CONTAINER->getSTACKS_DIR()) + "/" + DIR_NAME + "/" + FILENAMES[i];

	iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->readMetadata(image_path,img_width,img_height,img_bytes_x_chan,img_chans);

	// 2016-08-29. Giulio. @ADDED check of coherence
	if ( (img_width != WIDTH) ||  (img_height != HEIGHT) || (img_bytes_x_chan != N_BYTESxCHAN) || (img_chans != N_CHANS) )
		throw iom::exception(iom::strprintf("inconsistent metadata in file %s",image_path.c_str()).c_str(), __iom__current__function__);

	// 2016-04-29. Giulio. patch to test opencv2D when using jpeg
	//img_chans = 1;

	// allocate buffer
	data = new unsigned char[img_width * img_height * img_bytes_x_chan * img_chans];

	// loop over files
	for(int file_i = first_file; file_i <= last_file; file_i++, z++) {

		// skip missing slices if stack is sparse
		if(isSparse() && !FILENAMES[file_i])
			continue;

		// if stack is not sparse, a missing slice must throw an iom::exception
		if(!FILENAMES[file_i])
			throw iom::exception("invalid slice filename in non-sparse tile", __iom__current__function__);

		// build absolute image path
		image_path = std::string(CONTAINER->getSTACKS_DIR()) + "/" + DIR_NAME + "/" + FILENAMES[file_i];

		iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->readData(image_path,img_width,img_height,img_bytes_x_chan,img_chans,data);

		// initialize output data and image dimensions
		if(!STACKED_IMAGE)
		{
			image_stack_height = (uint64) img_height;
			image_stack_width  = (uint64) img_width;
			uint64 image_stack_size = image_stack_width * image_stack_height * (last_file - first_file + 1);
			STACKED_IMAGE = new iom::real_t[image_stack_size];
			for(uint64 j=0; j < image_stack_size; j++)
				STACKED_IMAGE[j] = 0;		// default is 0 (black)
		}
		else if(image_stack_width != (uint64)img_width || image_stack_height != (uint64)img_height)
			throw iom::exception(iom::strprintf("images in stack have not the same dimensions (%d x %d != %d x %d)", img_width, img_height, image_stack_width, image_stack_height), __iom__current__function__);

		int offset;

		// convert image to [0.0,1.0]-valued array
		if ( img_chans == 1 ) {
			raw_data = &STACKED_IMAGE[z*img_height*img_width];
			if ( img_bytes_x_chan == 1 ) {
				for(int i = 0; i < (img_height * img_width); i++)
					raw_data[i] = (iom::real_t) data[i]/255.0f;
			}
			else { // img_bytes_x_chan == 2
				for(int i = 0; i < (img_height * img_width); i++)
					raw_data[i] = (iom::real_t) ((uint16 *)data)[i]/65535.0f; // data must be interpreted as a uint16 array
			}
		}
		else { // conversion to an intensity image
			// test how channel are stored in the returned buffer 'data'
			if ( iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->isChansInterleaved() ) {
				if ( CHANS == iom::ALL ) {
					throw iom::exception("conversion from multi-channel to intensity images not supported.", __iom__current__function__);
				}
				else if ( CHANS == iom::R ) {
					offset = 0;
				}
				else if ( CHANS == iom::G ) {
					offset = 1;
				}
				else if ( CHANS == iom::B ) {
					offset = 2;
				}
				else {
					throw iom::exception("wrong value for parameter iom::CHANNEL_SELECTION.", __iom__current__function__);
				}
				raw_data = &STACKED_IMAGE[z*img_height*img_width];
				if ( img_bytes_x_chan == 1 )
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) data[3*i + offset]/255.0f;
				else // img_bytes_x_chan == 2
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) ((uint16 *)data)[3*i + offset]/65535.0f; // data must be interpreted as a uint16 array
			}
			else {
				throw iom::exception("2D image formats with non-interleaved channels are not supported yet.");
			}
		}
	}

	if ( first_file == last_file ) {
		if ( !CONTAINER->getCACHEBUFFER()->cacheSlice(ROW_INDEX,COL_INDEX,first_file,data) ) { // slice cannot be cached
			// 2015-08-24. Giulio. data must be released
			delete [] data;
		}
	}
	else {
		// 2015-08-24. Giulio. data must be released
		delete [] data;
	}

	return STACKED_IMAGE;
}


//deallocates memory used by STACKED_IMAGE
void Stack::releaseImageStack()
{
	if(STACKED_IMAGE)
		delete[] STACKED_IMAGE;
    STACKED_IMAGE=0;
}

// compute 'z_ranges'
void 
	Stack::compute_z_ranges(
	std::set<std::string> const * z_coords /*= 0*/)	// set of z-coordinates where at least one slice (of a certain stack) is available
									// if null, 'z_ranges' will be compute based on 'FILENAMES' vector
{
	// if 'z_coords' has been provided, we use it to associate each file in 'FILENAMES' to the correspondent z-coordinate
	if(z_coords)
	{
		// check for nonempty 'z_coords'
		if(z_coords->empty())
			throw iom::exception("in Stack::compute_z_ranges(): z_coords is empty");

		// associate each file to the correspondent z-coordinate
		char **FILENAMES_sparse = new char*[z_coords->size()];
		int DEPTH_sparse = 0;
		int k=0;
        for(std::set<std::string>::const_iterator it = z_coords->begin(); it != z_coords->end(); it++)
		{
			if(k<DEPTH && it->compare(VirtualVolume::name2coordZ(FILENAMES[k])) == 0)
				FILENAMES_sparse[DEPTH_sparse++] = FILENAMES[k++];
			else
				FILENAMES_sparse[DEPTH_sparse++] = 0;
		}


		// substitute list of filenames with its sparse version
		if(FILENAMES)
			delete[] FILENAMES;
		FILENAMES = FILENAMES_sparse;
		DEPTH = DEPTH_sparse;
	}

	// compute 'z_range' from 'FILENAMES'
	bool interval_start = true;
	bool interval_stop = false;
    z_ranges.clear();
	for(int k=0; k<DEPTH; k++)
	{
		// valid slice found
		if(FILENAMES[k] && interval_start)
		{
			interval_start = false;
			interval_stop  = true;
			z_ranges.push_back(vm::interval<int>(k, -1));
		}
		if(FILENAMES[k] == 0 && interval_stop)
		{
			interval_start = true;
			interval_stop  = false;
			z_ranges.back().end = k;
		}
	}
	if(!z_ranges.empty() && z_ranges.back().end == -1)
		z_ranges.back().end = DEPTH;
}


//loads image stack from <first_file> to <last_file> extremes included, if not specified loads entire Stack
iom::real_t* Stack::loadImageStack2(int first_file, int last_file, int V0, int V1, int H0, int H1) 
{
	#if VM_VERBOSE > 3
	printf("\t\t\t\tin Stack[%d,%d](%s, empty = %s)::loadImageStack(first_file = %d, last_file = %d)\n",ROW_INDEX, COL_INDEX, DIR_NAME, isEmpty() ? "true": "false", first_file, last_file);
	#endif

	uint64 image_stack_width=0, image_stack_height=0, z=0;
 	iom::real_t *raw_data;

	unsigned char *data = (unsigned char *) 0;
	int img_width;
	int img_height;
	int img_bytes_x_chan;
	int img_chans;

	if ( STACKED_IMAGE ) {
		delete[] STACKED_IMAGE;
    	STACKED_IMAGE=0;
		throw iom::exception(iom::strprintf("precondition violated in Stack [%d,%d] loading slices [%d,%d]: buffer is already allocated", ROW_INDEX, COL_INDEX,first_file, last_file), __iom__current__function__);
	}

	// check and adjust file selection
	first_file = (first_file == -1 ? 0 : first_file);
	last_file  = (last_file  == -1 ? DEPTH - 1 : last_file );
	first_file = min(first_file, DEPTH-1);
	last_file = min(last_file, DEPTH-1);

	if(isEmpty(first_file, last_file) ||
		V0 >= ABS_V+HEIGHT || V1 <= ABS_V || H0 >= ABS_H+WIDTH || H1 <= ABS_H ) //2018-01-21. Giulio. @ADDED conditions to check if the requested region does not overlap with stack
	{
		// allocate and initialize a black stack
		uint64 image_size = static_cast<uint64>(WIDTH) * static_cast<uint64>(HEIGHT) * static_cast<uint64>(last_file-first_file+1);
		STACKED_IMAGE = new iom::real_t[image_size];
		for(uint64 k=0; k<image_size; k++)
			STACKED_IMAGE[k] = 0;

		return STACKED_IMAGE;
	}

	if ( first_file == last_file ) { // only one slice has been requested
		if ( CONTAINER->getCACHEBUFFER()->getSlice(ROW_INDEX,COL_INDEX,first_file,data) ) {
			img_width        = WIDTH;
			img_height       = HEIGHT;
			img_bytes_x_chan = N_BYTESxCHAN;
			img_chans        = N_CHANS;

			/**** DUPLICATED CODE (BEGIN) ***/
			// initialize output data and image dimensions
			if(!STACKED_IMAGE)
			{
				image_stack_height = (uint64) img_height;
				image_stack_width  = (uint64) img_width;
				uint64 image_stack_size = image_stack_width * image_stack_height * (last_file - first_file + 1);
				STACKED_IMAGE = new iom::real_t[image_stack_size];
				for(uint64 j=0; j < image_stack_size; j++)
					STACKED_IMAGE[j] = 0;		// default is 0 (black)
			}
			else if(image_stack_width != (uint64)img_width || image_stack_height != (uint64)img_height)
				throw iom::exception(iom::strprintf("images in stack have not the same dimensions (%d x %d != %d x %d)", img_width, img_height, image_stack_width, image_stack_height), __iom__current__function__);

			int offset;

			// convert image to [0.0,1.0]-valued array
			if ( img_chans == 1 ) {
				raw_data = &STACKED_IMAGE[z*img_height*img_width];
				if ( img_bytes_x_chan == 1 ) {
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) data[i]/255.0f;
				}
				else { // img_bytes_x_chan == 2
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) ((uint16 *)data)[i]/65535.0f; // data must be interpreted as a uint16 array
				}
			}
			else { // conversion to an intensity image
				// test how channel are stored in the returned buffer 'data'
				if ( iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->isChansInterleaved() ) {
					if ( CHANS == iom::ALL ) {
						throw iom::exception("conversion from multi-channel to intensity images not supported.", __iom__current__function__);
					}
					else if ( CHANS == iom::R ) {
						offset = 0;
					}
					else if ( CHANS == iom::G ) {
						offset = 1;
					}
					else if ( CHANS == iom::B ) {
						offset = 2;
					}
					else {
						throw iom::exception("wrong value for parameter iom::CHANNEL_SELECTION.", __iom__current__function__);
					}
					raw_data = &STACKED_IMAGE[z*img_height*img_width];
					if ( img_bytes_x_chan == 1 )
						for(int i = 0; i < (img_height * img_width); i++)
							raw_data[i] = (iom::real_t) data[3*i + offset]/255.0f;
					else // img_bytes_x_chan == 2
						for(int i = 0; i < (img_height * img_width); i++)
							raw_data[i] = (iom::real_t) ((uint16 *)data)[3*i + offset]/65535.0f; // data must be interpreted as a uint16 array
				}
				else {
					throw iom::exception("2D image formats with non-interleaved channels are not supported yet.");
				}
			}
			/**** DUPLICATED CODE (END) ***/

			return STACKED_IMAGE;
		}
	}

	// build absolute image path of first file in the stack
	int i = first_file;
	while ( isSparse() && i<=last_file && FILENAMES[i]==0 ) // look for existing file
		i++;
	if ( i > last_file )
		throw iom::exception("unexpected exception: substack is entirely empty", __iom__current__function__);

	std::string image_path = std::string(CONTAINER->getSTACKS_DIR()) + "/" + DIR_NAME + "/" + FILENAMES[i];

	iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->readMetadata(image_path,img_width,img_height,img_bytes_x_chan,img_chans);

	// 2016-08-29. Giulio. @ADDED check of coherence
	if ( (img_width != WIDTH) ||  (img_height != HEIGHT) || (img_bytes_x_chan != N_BYTESxCHAN) || (img_chans != N_CHANS) )
		throw iom::exception(iom::strprintf("inconsistent metadata in file %s",image_path.c_str()).c_str(), __iom__current__function__);

	// 2016-04-29. Giulio. patch to test opencv2D when using jpeg
	//img_chans = 1;

	// allocate buffer
	data = new unsigned char[img_width * img_height * img_bytes_x_chan * img_chans];

	// loop over files
	for(int file_i = first_file; file_i <= last_file; file_i++, z++) {

		// skip missing slices if stack is sparse
		if(isSparse() && !FILENAMES[file_i])
			continue;

		// if stack is not sparse, a missing slice must throw an iom::exception
		if(!FILENAMES[file_i])
			throw iom::exception("invalid slice filename in non-sparse tile", __iom__current__function__);

		// build absolute image path
		image_path = std::string(CONTAINER->getSTACKS_DIR()) + "/" + DIR_NAME + "/" + FILENAMES[file_i];

		iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->readData(image_path,img_width,img_height,img_bytes_x_chan,img_chans,data);

		// initialize output data and image dimensions
		if(!STACKED_IMAGE)
		{
			image_stack_height = (uint64) img_height;
			image_stack_width  = (uint64) img_width;
			uint64 image_stack_size = image_stack_width * image_stack_height * (last_file - first_file + 1);
			STACKED_IMAGE = new iom::real_t[image_stack_size];
			for(uint64 j=0; j < image_stack_size; j++)
				STACKED_IMAGE[j] = 0;		// default is 0 (black)
		}
		else if(image_stack_width != (uint64)img_width || image_stack_height != (uint64)img_height)
			throw iom::exception(iom::strprintf("images in stack have not the same dimensions (%d x %d != %d x %d)", img_width, img_height, image_stack_width, image_stack_height), __iom__current__function__);

		int offset;

		// convert image to [0.0,1.0]-valued array
		if ( img_chans == 1 ) {
			raw_data = &STACKED_IMAGE[z*img_height*img_width];
			if ( img_bytes_x_chan == 1 ) {
				for(int i = 0; i < (img_height * img_width); i++)
					raw_data[i] = (iom::real_t) data[i]/255.0f;
			}
			else { // img_bytes_x_chan == 2
				for(int i = 0; i < (img_height * img_width); i++)
					raw_data[i] = (iom::real_t) ((uint16 *)data)[i]/65535.0f; // data must be interpreted as a uint16 array
			}
		}
		else { // conversion to an intensity image
			// test how channel are stored in the returned buffer 'data'
			if ( iom::IOPluginFactory::getPlugin2D(iom::IMIN_PLUGIN)->isChansInterleaved() ) {
				if ( CHANS == iom::ALL ) {
					throw iom::exception("conversion from multi-channel to intensity images not supported.", __iom__current__function__);
				}
				else if ( CHANS == iom::R ) {
					offset = 0;
				}
				else if ( CHANS == iom::G ) {
					offset = 1;
				}
				else if ( CHANS == iom::B ) {
					offset = 2;
				}
				else {
					throw iom::exception("wrong value for parameter iom::CHANNEL_SELECTION.", __iom__current__function__);
				}
				raw_data = &STACKED_IMAGE[z*img_height*img_width];
				if ( img_bytes_x_chan == 1 )
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) data[3*i + offset]/255.0f;
				else // img_bytes_x_chan == 2
					for(int i = 0; i < (img_height * img_width); i++)
						raw_data[i] = (iom::real_t) ((uint16 *)data)[3*i + offset]/65535.0f; // data must be interpreted as a uint16 array
			}
			else {
				throw iom::exception("2D image formats with non-interleaved channels are not supported yet.");
			}
		}
	}

	if ( first_file == last_file ) {
		if ( !CONTAINER->getCACHEBUFFER()->cacheSlice(ROW_INDEX,COL_INDEX,first_file,data) ) { // slice cannot be cached
			// 2015-08-24. Giulio. data must be released
			delete [] data;
		}
	}
	else {
		// 2015-08-24. Giulio. data must be released
		delete [] data;
	}

	return STACKED_IMAGE;
}
