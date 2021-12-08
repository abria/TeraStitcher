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
* 2018-12-28. Giulio.     @FIXED in 'loadSubvolume_to_UINT8' in call to 'compact_active_chans' the channel size depends also on sbv_depth 
* 2018-10-23. Yang.       @ADDED support partial data, YY 10/23/2018
* 2018-07-05. Giulio.     @ADDED remapping of 8 bits images for better visualization
* 2018-06-30. Giulio.     @ADDED parameter for specifying the conversion algorithm to be used to convert from arbitrary depth to 8 bits
* 2017-10-21. Giulio.     @ADDED compact active channels if not all channels are active in 'loadSubvolume_to_UINT8'
* 2016-09-01. Giulio.     @DELETED old commented code
* 2016-04-10. Giulio.     @FIXED the set of the right input plugin when input volume contains TIFF 3D files (in init or load, not in the constructor)
* 2015-04-15. Alessandro. @FIXED bad/missing exception handling in loadSubvolume_to_UINT8.
* 2015-04-15. Alessandro. @ADDED definition for default constructor.
* 2015-04-06. Giulio.     @CHANGED Modified prunt method: printing stacks information is now off by default
* 2015-03-03. Giulio.     @ADDED check that ioplugin interleaves channels in loadSubvolume_to_UINT8.
* 2015-03-03. Giulio.     @ADDED selection of IO plugin in the constructors if not provided.
* 2015-02-27. Alessandro. @ADDED automated selection of IO plugin if not provided.
* 2014-11-22  Giulio.     @CHANGED code using OpenCV has been commente. It can be found searching comments containing 'Giulio_CV'
*/

#include<iostream>
#include<string>
#include<thread>
#include<chrono>
#include "TiledVolume.h"
#include "imBlock.h"
#include "VirtualFmtMngr.h"
#include "Tiff3DMngr.h"
#include "RawFmtMngr.h"
#include "iomanager.config.h"

#ifdef _WIN32
#include "dirent_win.h"
#else
#include <dirent.h>
#endif

// Giulio_CV #include <cxcore.h>
// Giulio_CV #include <cv.h>
// Giulio_CV #include <highgui.h>

#include <list>
#include <fstream>
#include "ProgressBar.h"

#include "IOPluginAPI.h" // 2015-03-03. Giulio.

using namespace std;
using namespace iim;

// 2015-04-15. Alessandro. @ADDED definition for default constructor.
TiledVolume::TiledVolume(void) : VirtualVolume()
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);
    
    N_ROWS = N_COLS = 0;
    BLOCKS = 0;
    reference_system.first = reference_system.second = reference_system.third = iim::axis_invalid;
    VXL_1 = VXL_2 = VXL_3 = 0.0f;
    fmtMngr = 0;
}

TiledVolume::TiledVolume(const char* _root_dir)  throw (IOException)
: VirtualVolume(_root_dir) // iannello ADDED
{
    /**/iim::debug(iim::LEV3, strprintf("_root_dir=%s", _root_dir).c_str(), __iim__current__function__);

	// 2016-04-10. Giulio. @REMOVED because a check is done in both init and load methods and the right input plugin is set
	// 2015-03-03. Giulio. @ADDED selection of IO plugin if not provided.
	//if(iom::IMIN_PLUGIN.compare("empty") == 0)
	//{
	//	iom::IMIN_PLUGIN = "tiff3D";
	//}

	//DIM_V = DIM_H = DIM_D = 0;
	VXL_1 = VXL_2 = VXL_3 = 0;
	N_ROWS = N_COLS = 0;
	BLOCKS = NULL;
	reference_system.first = reference_system.second = reference_system.third = axis_invalid;

	ffmt = "";
	fmtMngr = 0;

	//without any configuration parameter, volume import must be done from the metadata file stored in the root directory, if it exists
    char mdata_filepath[STATIC_STRINGS_SIZE];
    sprintf(mdata_filepath, "%s/%s", root_dir, iim::MDATA_BIN_FILE_NAME.c_str());
    if(iim::isFile(mdata_filepath))
	{
		load(mdata_filepath);
		initChannels();
	}
	else
	{
        char errMsg[STATIC_STRINGS_SIZE];
		sprintf(errMsg, "in TiledVolume::TiledVolume(...): unable to find metadata file at %s", mdata_filepath);
        throw IOException(errMsg);
	}
}

TiledVolume::TiledVolume(const char* _root_dir, ref_sys _reference_system, float _VXL_1, float _VXL_2, float _VXL_3, bool overwrite_mdata, bool save_mdata)  throw (IOException)
: VirtualVolume(_root_dir) // iannello ADDED
{
    /**/iim::debug(iim::LEV3, strprintf("_root_dir=%s, ref_sys reference_system={%d,%d,%d}, VXL_1=%.4f, VXL_2=%.4f, VXL_3=%.4f",
                                        _root_dir, _reference_system.first, _reference_system.second, _reference_system.third, _VXL_1, _VXL_2, _VXL_3).c_str(), __iim__current__function__);

	// 2016-04-10. Giulio. @REMOVED because a check is done in both init and load methods and the right input plugin is set
	// 2015-03-03. Giulio. @ADDED selection of IO plugin if not provided.
	//if(iom::IMIN_PLUGIN.compare("empty") == 0)
	//{
	//	iom::IMIN_PLUGIN = "tiff3D";
	//}

	//DIM_V = DIM_H = DIM_D = 0;
	VXL_1 = VXL_2 = VXL_3 = 0;
	N_ROWS = N_COLS = 0;
	BLOCKS = NULL;
	reference_system.first = reference_system.second = reference_system.third = axis_invalid;

	ffmt = "";
	fmtMngr = 0;

	//trying to unserialize an already existing metadata file, if it doesn't exist the full initialization procedure is performed and metadata is saved
    char mdata_filepath[STATIC_STRINGS_SIZE];
    sprintf(mdata_filepath, "%s/%s", root_dir, iim::MDATA_BIN_FILE_NAME.c_str());
    if(iim::isFile(mdata_filepath) && !overwrite_mdata)
		load(mdata_filepath);
	else
	{
        if(_reference_system.first == axis_invalid ||  _reference_system.second == axis_invalid ||
          _reference_system.third == axis_invalid || _VXL_1 == 0 || _VXL_2 == 0 || _VXL_3 == 0)
            throw IOException("in TiledVolume::TiledVolume(...): invalid importing parameters");

        reference_system.first  = _reference_system.first;
        reference_system.second = _reference_system.second;
        reference_system.third  = _reference_system.third;
        VXL_1 = _VXL_1;
        VXL_2 = _VXL_2;
        VXL_3 = _VXL_3;
        init();
        if(save_mdata)
            save(mdata_filepath);
	}
	initChannels();
}

TiledVolume::~TiledVolume(void) throw (iim::IOException)
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	// iannello if(root_dir)
	// iannello 	delete[] root_dir;

	if(BLOCKS)
	{
		for(int row=0; row<N_ROWS; row++)
		{
			for(int col=0; col<N_COLS; col++)
				delete BLOCKS[row][col];
			delete[] BLOCKS[row];
		}
		delete[] BLOCKS;
	}

	if ( fmtMngr )
		delete fmtMngr;
}


int TiledVolume::getStacksHeight(){return BLOCKS[0][0]->getHEIGHT();}
int TiledVolume::getStacksWidth(){return BLOCKS[0][0]->getWIDTH();}

void TiledVolume::save(char* metadata_filepath) throw (IOException)
{
    /**/iim::debug(iim::LEV3, strprintf("metadata_filepath=%s", metadata_filepath).c_str(), __iim__current__function__);

	//LOCAL VARIABLES
	FILE *file;
	int i,j;

    file = fopen(metadata_filepath, "wb");

    // --- Alessandro 2013-04-23: added exception when file can't be opened in write mode
    if(!file)
    {
        char errMsg[STATIC_STRINGS_SIZE];
        sprintf(errMsg, "in TiledVolume::save(): cannot write metadata binary file at \"%s\".\n\nPlease check write permissions on this storage.", metadata_filepath);
        throw IOException(errMsg);
    }

    float mdata_version = static_cast<float>(iim::MDATA_BIN_FILE_VERSION);
    fwrite(&mdata_version, sizeof(float), 1, file); // --- Alessandro 2012-12-31: added field for metadata file version
    fwrite(&reference_system.first, sizeof(axis), 1, file);
    fwrite(&reference_system.second, sizeof(axis), 1, file); // iannello CORRECTED
    fwrite(&reference_system.third, sizeof(axis), 1, file);  // iannello CORRECTED
    fwrite(&VXL_1, sizeof(float), 1, file);
    fwrite(&VXL_2, sizeof(float), 1, file);
    fwrite(&VXL_3, sizeof(float), 1, file);
	fwrite(&VXL_V, sizeof(float), 1, file);
	fwrite(&VXL_H, sizeof(float), 1, file);
	fwrite(&VXL_D, sizeof(float), 1, file);
	fwrite(&ORG_V, sizeof(float), 1, file);
	fwrite(&ORG_H, sizeof(float), 1, file);
	fwrite(&ORG_D, sizeof(float), 1, file);
	fwrite(&DIM_V, sizeof(iim::uint32), 1, file);
	fwrite(&DIM_H, sizeof(iim::uint32), 1, file);
	fwrite(&DIM_D, sizeof(iim::uint32), 1, file);
	fwrite(&N_ROWS, sizeof(iim::uint16), 1, file);
	fwrite(&N_COLS, sizeof(iim::uint16), 1, file);

	for(i = 0; i < N_ROWS; i++)
		for(j = 0; j < N_COLS; j++)
			BLOCKS[i][j]->binarizeInto(file);

	fclose(file);
}

void TiledVolume::load(char* metadata_filepath) throw (IOException)
{
    /**/iim::debug(iim::LEV3, strprintf("metadata_filepath=%s", metadata_filepath).c_str(), __iim__current__function__);

	//LOCAL VARIABLES
	FILE *file;
	int i,j;
	size_t fread_return_val;

	file = fopen(metadata_filepath, "rb");

    // --- Alessandro 2013-04-23: added exception when file can't be opened in read mode
    if(!file)
    {
        char errMsg[STATIC_STRINGS_SIZE];
        sprintf(errMsg, "in TiledVolume::load(): cannot read metadata binary file at \"%s\".\n\nPlease check read permissions on this storage.", metadata_filepath);
        throw IOException(errMsg);
    }

    // --- Alessandro 2012-12-31: added field for metadata file version
    float mdata_version_read = 0;
    float mdata_version = static_cast<float>(iim::MDATA_BIN_FILE_VERSION);
    fread_return_val = fread(&mdata_version_read, sizeof(float), 1, file);
    if(fread_return_val != 1 || mdata_version_read != mdata_version)
    {
        // --- Alessandro 2013-01-06: instead of throwing an exception, it is better to mantain compatibility
//            char errMsg[STATIC_STRINGS_SIZE];
//            sprintf(errMsg, "in Block::unBinarizeFrom(...): metadata file version (%.2f) is different from the supported one (%.2f). "
//                    "Please re-import the current volume.", mdata_version_read, mdata_version);
//            throw iim::IOException(errMsg);

        fclose(file);
        file = fopen(metadata_filepath, "rb");
        iim::uint16 str_size;
        fread_return_val = fread(&str_size, sizeof(iim::uint16), 1, file);
        if(fread_return_val != 1)
        {
            fclose(file);
            throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
        }
        char stored_root_dir[STATIC_STRINGS_SIZE];
        fread_return_val = fread(stored_root_dir, str_size, 1, file);
        if(fread_return_val != 1)
        {
            fclose(file);
            throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
        }
    }

    fread_return_val = fread(&reference_system.first, sizeof(axis), 1, file);
    if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

    fread_return_val = fread(&reference_system.second, sizeof(axis), 1, file);
    if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

    fread_return_val = fread(&reference_system.third, sizeof(axis), 1, file);
    if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

    fread_return_val = fread(&VXL_1, sizeof(float), 1, file);
    if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

    fread_return_val = fread(&VXL_2, sizeof(float), 1, file);
    if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

    fread_return_val = fread(&VXL_3, sizeof(float), 1, file);
    if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&VXL_V, sizeof(float), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&VXL_H, sizeof(float), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&VXL_D, sizeof(float), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&ORG_V, sizeof(float), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&ORG_H, sizeof(float), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&ORG_D, sizeof(float), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&DIM_V, sizeof(iim::uint32), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&DIM_H, sizeof(iim::uint32), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&DIM_D, sizeof(iim::uint32), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&N_ROWS, sizeof(iim::uint16), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	fread_return_val = fread(&N_COLS, sizeof(iim::uint16), 1, file);
	if(fread_return_val != 1)
    {
        fclose(file);
        throw IOException("in Block::unBinarizeFrom(...): error while reading binary metadata file");
    }

	BLOCKS = new Block **[N_ROWS];
	for(i = 0; i < N_ROWS; i++)
	{
		BLOCKS[i] = new Block *[N_COLS];
		for(j = 0; j < N_COLS; j++)
			BLOCKS[i][j] = new Block(this, i, j, file);
	}

	// get the file format of the blocks
	ffmt = BLOCKS[0][0]->getFMT();
	if ( ffmt == "Tiff3D" )
	{
		// 2015-02-27. Alessandro. @ADDED automated selection of IO plugin if not provided.
		if(iom::IMIN_PLUGIN.compare("tiff3D") != 0)
			iom::IMIN_PLUGIN = "tiff3D";
		fmtMngr = new Tiff3DFmtMngr();
	}
	else if ( ffmt == "Vaa3DRaw" )
	{
		fmtMngr = new Vaa3DRawFmtMngr();
	}
	else {
        char msg[STATIC_STRINGS_SIZE];
        sprintf(msg,"in TiledVolume::unBinarizeFrom(...): Unknown file format \"%s\"", ffmt.c_str());
        throw IOException(msg);
	}

	fclose(file);
}

void TiledVolume::init() throw (IOException)
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	/************************* 1) LOADING STRUCTURE *************************    
	*************************************************************************/

	//LOCAL VARIABLES
	string tmp_path;				//string that contains temp paths during computation
    string tmp;						//string that contains temp data during computation
	DIR *cur_dir_lev1;				//pointer to DIR, the data structure that represents a DIRECTORY (level 1 of hierarchical structure)
	DIR *cur_dir_lev2;				//pointer to DIR, the data structure that represents a DIRECTORY (level 2 of hierarchical structure)
	dirent *entry_lev1;				//pointer to DIRENT, the data structure that represents a DIRECTORY ENTRY inside a directory (level 1)
	dirent *entry_lev2;				//pointer to DIRENT, the data structure that represents a DIRECTORY ENTRY inside a directory (level 2)
	int i=0,j=0;					//for counting of N_ROWS, N_COLS
    list<Block*> blocks_list;                       //each stack found in the hierarchy is pushed into this list
    list<string> entries_lev1;                      //list of entries of first level of hierarchy
    list<string>::iterator entry_i;                 //iterator for list 'entries_lev1'
    list<string> entries_lev2;                      //list of entries of second level of hierarchy
    list<string>::iterator entry_j;                 //iterator for list 'entries_lev2'
    char block_i_j_path[STATIC_STRINGS_SIZE];

	//obtaining DIR pointer to root_dir (=NULL if directory doesn't exist)
    const unsigned int num_retries = 40;
    bool error_open_dir = false;
#ifdef _MSC_VER
#pragma loop( no_vector )
#elif __GNUC__
#pragma GCC unroll 1
#elif __MINGW32__
#pragma GCC unroll 1
#elif __clang__
#pragma clang loop unroll(disable)
#elif __BORLANDC__
#pragma nounroll
#elif __INTEL_COMPILER
#pragma nounroll
#endif
    for (int attempt = 0; attempt < num_retries; attempt++) {
        try {
            cur_dir_lev1 = opendir(root_dir);
            if (!cur_dir_lev1) throw std::exception();
            error_open_dir = false;
        }
        catch (const std::exception& exc) {
            (void)exc;
            error_open_dir = true;
            std::this_thread::sleep_for(std::chrono::milliseconds(100)); //ms
            continue;
        }
        break;
    }
    if (error_open_dir)
	{
        char msg[STATIC_STRINGS_SIZE];
        sprintf(msg,"in TiledVolume::init(...): Unable to open directory \"%s\"", root_dir);
        throw IOException(msg);
	}

	//scanning first level of hierarchy which entries need to be ordered alphabetically. This is done using STL.
	while ((entry_lev1=readdir(cur_dir_lev1)))
	{
        tmp=entry_lev1->d_name;
        if(tmp.find(".") == string::npos && tmp.find(" ") == string::npos)
                entries_lev1.push_front(entry_lev1->d_name);
	}
	closedir(cur_dir_lev1);
	entries_lev1.sort();
	N_ROWS = (iim::uint16) entries_lev1.size();
	N_COLS = 0;

	//for each entry of first level, scanning second level
	for(entry_i = entries_lev1.begin(), i=0; entry_i!= entries_lev1.end(); entry_i++, i++)
	{
		//building absolute path of first level entry to be used for "opendir(...)"
		tmp_path=root_dir;
		tmp_path.append("/");
		tmp_path.append(*entry_i);
		cur_dir_lev2 = opendir(tmp_path.c_str());
		if (!cur_dir_lev2)
            throw IOException("in TiledVolume::init(...): A problem occurred during scanning of subdirectories");

		//scanning second level of hierarchy which entries need to be ordered alphabetically. This is done using STL.
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
            //allocating new stack
            sprintf(block_i_j_path,"%s/%s",(*entry_i).c_str(), (*entry_j).c_str());
            Block *new_stk = new Block(this,i,j,block_i_j_path);
            blocks_list.push_back(new_stk);
		}
		entries_lev2.clear();
		if (N_COLS == 0)
			N_COLS = j;
		else if (j != N_COLS)
            throw IOException("in TiledVolume::init(...): Number of second-level directories is not the same for all first-level directories!");
	}
	entries_lev1.clear();

	//intermediate check
	if(N_ROWS == 0 || N_COLS == 0)
            throw IOException("in TiledVolume::init(...): Unable to find stacks in the given directory");

	//converting stacks_list (STL list of Block*) into BLOCKS (2-D array of Block*)
	BLOCKS = new Block**[N_ROWS];
	for(int row=0; row < N_ROWS; row++)
            BLOCKS[row] = new Block*[N_COLS];
	for(list<Block*>::iterator i = blocks_list.begin(); i != blocks_list.end(); i++)
            BLOCKS[(*i)->getROW_INDEX()][(*i)->getCOL_INDEX()] = (*i);

	// get the file format of the blocks
	ffmt = BLOCKS[0][0]->getFMT();
	if ( ffmt == "Tiff3D" )
	{
		// 2015-02-27. Alessandro. @ADDED automated selection of IO plugin if not provided.
		if(iom::IMIN_PLUGIN.compare("tiff3D") != 0)
			iom::IMIN_PLUGIN = "tiff3D";
		fmtMngr = new Tiff3DFmtMngr();
	}
	else if ( ffmt == "Vaa3DRaw" )
		fmtMngr = new Vaa3DRawFmtMngr();
	else {
        char msg[STATIC_STRINGS_SIZE];
        sprintf(msg,"in TiledVolume::init(...): Unknown file format \"%s\"", ffmt.c_str());
        throw IOException(msg);
	}


	/******************* 2) SETTING THE REFERENCE SYSTEM ********************
	The entire application uses a vertical-horizontal reference system, so
	it is necessary to fit the original reference system into the new one.     
	*************************************************************************/

	//adjusting possible sign mismatch betwwen reference system and VXL
	//in these cases VXL is adjusted to match with reference system
    if(sgn(reference_system.first) != sgn(VXL_1))
            VXL_1*=-1.0f;
    if(sgn(reference_system.second) != sgn(VXL_2))
            VXL_2*=-1.0f;
    if(sgn(reference_system.third) != sgn(VXL_3))
            VXL_3*=-1.0f;

	//HVD --> VHD
        if  (abs(reference_system.first)==2 && abs(reference_system.second)==1  && reference_system.third==3)
	{
        this->rotate(90);
        this->mirror(axis(2));

        if(reference_system.first == -2)
                this->mirror(axis(2));
        if(reference_system.second == -1)
                this->mirror(axis(1));

        int computed_ORG_1, computed_ORG_2, computed_ORG_3;
        extractCoordinates(BLOCKS[0][0], 0, &computed_ORG_1, &computed_ORG_2, &computed_ORG_3);
        ORG_V = computed_ORG_2/10000.0F;
        ORG_H = computed_ORG_1/10000.0F;
        ORG_D = computed_ORG_3/10000.0F;
        VXL_V = VXL_2 ;
        VXL_H = VXL_1 ;
        VXL_D = VXL_3 ;
	}
	//VHD --> VHD
	else if (abs(reference_system.first)==1 && abs(reference_system.second)==2 && reference_system.third==3)
	{		
        if(reference_system.first == -1)
                this->mirror(axis(1));
        if(reference_system.second == -2)
                this->mirror(axis(2));

        int computed_ORG_1, computed_ORG_2, computed_ORG_3;
        extractCoordinates(BLOCKS[0][0], 0, &computed_ORG_1, &computed_ORG_2, &computed_ORG_3);
        ORG_V = computed_ORG_1/10000.0F;
        ORG_H = computed_ORG_2/10000.0F;
        ORG_D = computed_ORG_3/10000.0F;
        VXL_V = VXL_1;
        VXL_H = VXL_2;
        VXL_D = VXL_3;
	}
	//unsupported reference system
	else
	{
        char msg[STATIC_STRINGS_SIZE];
        sprintf(msg, "in TiledVolume::init(...): the reference system {%d,%d,%d} is not supported.",
                reference_system.first, reference_system.second, reference_system.third);
        throw IOException(msg);
	}

	// GI_140628: this adjustment should not be performed by the converter since it operates on volumes already adjusted
	//some little adjustments of the origin
	//if(VXL_V < 0)
	//	ORG_V -= (BLOCKS[0][0]->getHEIGHT()-1)* VXL_V/1000.0f;

	//if(VXL_H < 0)
	//	ORG_H -= (BLOCKS[0][0]->getWIDTH() -1)* VXL_H/1000.0f;

	/******************** 3) COMPUTING VOLUME DIMENSIONS ********************  
	*************************************************************************/
	for(int row=0; row < N_ROWS; row++) {
        for(int col=0; col < N_COLS; col++)
        {
            if(row==0)
				DIM_H+=BLOCKS[row][col]->getWIDTH();
            if(col==0)
				DIM_V+=BLOCKS[row][col]->getHEIGHT();
            DIM_D = BLOCKS[row][col]->getDEPTH() > DIM_D ? BLOCKS[row][col]->getDEPTH() : DIM_D;
        }
	}

	/**************** 4) COMPUTING STACKS ABSOLUTE POSITIONS ****************  
	*************************************************************************/
	for(int row=0; row < N_ROWS; row++) {
        for(int col=0; col < N_COLS; col++)
        {
            if(row)
				BLOCKS[row][col]->setABS_V(BLOCKS[row-1][col]->getABS_V()+BLOCKS[row-1][col]->getHEIGHT());
            else
				BLOCKS[row][col]->setABS_V(0);

            if(col)
				BLOCKS[row][col]->setABS_H(BLOCKS[row][col-1]->getABS_H()+BLOCKS[row][col-1]->getWIDTH());
            else
				BLOCKS[row][col]->setABS_H(0);
        }
	}
}

void TiledVolume::initChannels ( ) throw (IOException)
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

    DIM_C = BLOCKS[0][0]->getN_CHANS();
	BYTESxCHAN = (int)BLOCKS[0][0]->getN_BYTESxCHAN();

    n_active = DIM_C;
    active = new iim::uint32[n_active];
    for ( int c=0; c<DIM_C; c++ )
        active[c] = c; // all channels are assumed active
}

//PRINT method
void TiledVolume::print( bool print_stacks )
{
	printf("*** Begin printing TiledVolume object...\n\n");
	printf("\tDirectory:\t%s\n", root_dir);
	printf("\tDimensions:\t\t%d(V) x %d(H) x %d(D)\n", DIM_V, DIM_H, DIM_D);
	printf("\tVoxels:\t\t\t%.4f(V) x %.4f(H) x %.4f(D)\n", VXL_V, VXL_H, VXL_D);
	printf("\tOrigin:\t\t\t%.4f(V) x %.4f(H) x %.4f(D)\n", ORG_V, ORG_H, ORG_D);
	printf("\tChannels:\t\t%d\n", DIM_C);
	printf("\tBytes per channel:\t%d\n", BYTESxCHAN);
	printf("\tReference system:\tref1=%d, ref2=%d, ref3=%d\n",reference_system.first,reference_system.second,reference_system.third);
	printf("\tStacks matrix:\t\t%d(V) x %d(H)\n", N_ROWS, N_COLS);
	if ( print_stacks ) {
		printf("\t |\n");
		for(int row=0; row<N_ROWS; row++)
			for(int col=0; col<N_COLS; col++)
				BLOCKS[row][col]->print();
	}
	printf("\n*** END printing TiledVolume object...\n\n");
}

//rotate stacks matrix around D axis (accepted values are theta=0,90,180,270)
void TiledVolume::rotate(int theta)
{
    /**/iim::debug(iim::LEV3, strprintf("theta=%d", theta).c_str(), __iim__current__function__);

    //PRECONDITIONS:
    //	1) current TiledVolume object has been initialized (init() method has been called)
    //	2) accepted values for 'theta' are 0,90,180,270

    //POSTCONDITIONS:
    //  1) a new 2D-array of Block* objects is created considering a rotation of 'theta' angle of current TiledVolume object

    Block*** new_Block_2D_ARRAY = NULL;
    int new_N_ROWS = 0, new_N_COLS = 0;

	switch(theta)
	{
		case(0): break;

		case(90):
        {
            new_N_COLS = N_ROWS;
            new_N_ROWS = N_COLS;

            //allocating new_Block_2D_ARRAY
            new_Block_2D_ARRAY = new Block**[new_N_ROWS];
            for(int i=0; i<new_N_ROWS; i++)
                new_Block_2D_ARRAY[i] = new Block*[new_N_COLS];

            //populating new_Block_2D_ARRAY
            for(int i=0; i<new_N_ROWS; i++)
                for(int j=0; j<new_N_COLS; j++){
                    new_Block_2D_ARRAY[i][j] = BLOCKS[N_ROWS-1-j][i];
                    new_Block_2D_ARRAY[i][j]->setROW_INDEX(i);
                    new_Block_2D_ARRAY[i][j]->setCOL_INDEX(j);
                }
			break;
        }
		case(180):
        {
            new_N_COLS=N_COLS;
            new_N_ROWS=N_ROWS;

            //allocating new_Block_2D_ARRAY
            new_Block_2D_ARRAY = new Block**[new_N_ROWS];
            for(int i=0; i<new_N_ROWS; i++)
                new_Block_2D_ARRAY[i] = new Block*[new_N_COLS];

            //populating new_Block_2D_ARRAY
            for(int i=0; i<new_N_ROWS; i++)
                for(int j=0; j<new_N_COLS; j++){
                    new_Block_2D_ARRAY[i][j] = BLOCKS[N_ROWS-1-i][N_COLS-1-j];
                    new_Block_2D_ARRAY[i][j]->setROW_INDEX(i);
                    new_Block_2D_ARRAY[i][j]->setCOL_INDEX(j);
                }
			break;
        }
		case(270):
        {
            new_N_COLS=N_ROWS;
            new_N_ROWS=N_COLS;

            //allocating new_Block_2D_ARRAY
            new_Block_2D_ARRAY = new Block**[new_N_ROWS];
            for(int i=0; i<new_N_ROWS; i++)
                new_Block_2D_ARRAY[i] = new Block*[new_N_COLS];

            //populating new_Block_2D_ARRAY
            for(int i=0; i<new_N_ROWS; i++)
                for(int j=0; j<new_N_COLS; j++){
                    new_Block_2D_ARRAY[i][j] = BLOCKS[j][N_COLS-1-i];
                    new_Block_2D_ARRAY[i][j]->setROW_INDEX(i);
                    new_Block_2D_ARRAY[i][j]->setCOL_INDEX(j);
                }
			break;
        }
	}


    //deallocating current Block_2DARRAY object
    for(int row=0; row<N_ROWS; row++)
        delete[] BLOCKS[row];
    delete[] BLOCKS;

    BLOCKS = new_Block_2D_ARRAY;
    N_COLS = new_N_COLS;
    N_ROWS = new_N_ROWS;
}

//mirror stacks matrix along mrr_axis (accepted values are mrr_axis=1,2,3)
void TiledVolume::mirror(axis mrr_axis)
{
    /**/iim::debug(iim::LEV3, strprintf("mrr_axis=%d", mrr_axis).c_str(), __iim__current__function__);

    //PRECONDITIONS:
    //	1) current TiledVolume object has been initialized (init() method has been called)
    //	2) accepted values for 'mrr_axis' are 1(V axis), 2(H axis) or 3(D axis)

    //POSTCONDITIONS:
    //  1) a new 2D-array of Stack* objects is created considering a mirrorization along 'axis' of current TiledVolume object

    if(mrr_axis!= 1 && mrr_axis != 2)
    {
        char msg[1000];
        sprintf(msg,"in TiledVolume::mirror(axis mrr_axis=%d): unsupported axis mirroring", mrr_axis);
        throw IOException(msg);
    }

    Block*** new_BLOCK_2D_ARRAY;

    switch(mrr_axis)
    {
        case(1):
        {
            //allocating new_STACK_2D_ARRAY
            new_BLOCK_2D_ARRAY = new Block**[N_ROWS];
            for(int i=0; i<N_ROWS; i++)
                new_BLOCK_2D_ARRAY[i] = new Block*[N_COLS];

            //populating new_STACK_2D_ARRAY
            for(int i=0; i<N_ROWS; i++)
                for(int j=0; j<N_COLS; j++){
                    new_BLOCK_2D_ARRAY[i][j]=BLOCKS[N_ROWS-1-i][j];
                    new_BLOCK_2D_ARRAY[i][j]->setROW_INDEX(i);
                    new_BLOCK_2D_ARRAY[i][j]->setCOL_INDEX(j);
                }
            break;
        }
        case(2):
        {
            //allocating new_STACK_2D_ARRAY
            new_BLOCK_2D_ARRAY = new Block**[N_ROWS];
            for(int i=0; i<N_ROWS; i++)
                new_BLOCK_2D_ARRAY[i] = new Block*[N_COLS];

            //populating new_STACK_2D_ARRAY
            for(int i=0; i<N_ROWS; i++)
                for(int j=0; j<N_COLS; j++){
                    new_BLOCK_2D_ARRAY[i][j] = BLOCKS[i][N_COLS-1-j];
                    new_BLOCK_2D_ARRAY[i][j]->setROW_INDEX(i);
                    new_BLOCK_2D_ARRAY[i][j]->setCOL_INDEX(j);
                }
            break;
        }
        default: break;
    }

    //deallocating current STACK_2DARRAY object
    for(int row=0; row<N_ROWS; row++)
        delete[] BLOCKS[row];
    delete[] BLOCKS;

    BLOCKS = new_BLOCK_2D_ARRAY;
}

//extract spatial coordinates (in millimeters) of given Stack object
void TiledVolume::extractCoordinates(Block* blk, int z, int* crd_1, int* crd_2, int* crd_3)
{
    bool found_ABS_X=false;
    bool found_ABS_Y=false;

    //loading estimations for absolute X and Y stack positions
    char * pch;
    char buffer[100];
    strcpy(buffer,&(blk->getDIR_NAME()[0]));
    pch = strtok (buffer,"/_");
    pch = strtok (NULL, "/_");

    while (pch != NULL)
    {
        if(!found_ABS_X)
        {
            if(sscanf(pch, "%d", crd_1) == 1)
                found_ABS_X=true;
        }
        else if(!found_ABS_Y)
        {
            if(sscanf(pch, "%d", crd_2) == 1)
                found_ABS_Y=true;
        }
        else
            break;

        pch = strtok (NULL, "/_");
    }

    if(!found_ABS_X || !found_ABS_Y)
    {
        char msg[200];
        sprintf(msg,"in TiledVolume::extractCoordinates(directory_name=\"%s\"): format 000000_000000 or X_000000_X_000000 not found", blk->getDIR_NAME());
        throw msg;
    }

    //loading estimation for absolute Z stack position
    if(crd_3!= NULL)
    {
        char* first_file_name = blk->getFILENAMES()[z];

        char * pch;
        char lastTokenized[100];
        char buffer[500];
        strcpy(buffer,&(first_file_name[0]));

        pch = strtok (buffer,"_");
        while (pch != NULL)
        {
            strcpy(lastTokenized,pch);
            pch = strtok (NULL, "_");
        }

        pch = strtok (lastTokenized,".");
        strcpy(lastTokenized,pch);

        if(sscanf(lastTokenized, "%d", crd_3) != 1)
        {
            char msg[200];
            sprintf(msg,"in TiledVolume::extractCoordinates(...): unable to extract Z position from filename %s", first_file_name);
            throw msg;
        }
    }
}

//loads given subvolume in a 1-D array of float
real32* TiledVolume::loadSubvolume(int V0,int V1, int H0, int H1, int D0, int D1, list<Block*> *involved_blocks, bool release_blocks) throw (IOException)
{
    /**/iim::debug(iim::LEV3, strprintf("V0=%d, V1=%d, H0=%d, H1=%d, D0=%d, D1=%d, %s", V0, V1, H0, H1, D0, D1, (involved_blocks? ", involved_stacks" : "")).c_str(), __iim__current__function__);

 //   char msg[STATIC_STRINGS_SIZE];
	//sprintf(msg,"in TiledVolume::loadSubvolume: not completed yet");
 //   throw IOException(msg);

	int channels;
	float scale_factor;

    //initializations
    V0 = V0 < 0 ? 0 : V0;
    H0 = H0 < 0 ? 0 : H0;
    D0 = D0 < 0 ? 0 : D0;
    V1 = (V1 < 0 || V1 > (int)DIM_V) ? DIM_V : V1; // iannello MODIFIED
    H1 = (H1 < 0 || H1 > (int)DIM_H) ? DIM_H : H1; // iannello MODIFIED
    D1 = (D1 < 0 || D1 > (int)DIM_D) ? DIM_D : D1; // iannello MODIFIED

	iim::uint8 *data = loadSubvolume_to_UINT8(V0,V1,H0,H1,D0,D1,&channels,BYTESxCHAN*8);

	//conversion from unsigned char to iom::real_t

	if (DIM_C == 2 || DIM_C > 3) // only monocromatic or RGB images are supported
	{
		char errMsg[2000];
		sprintf(errMsg, "in TiledVolume::loadSubvolume(...): %d channels are not supported.", DIM_C);
		throw IOException(errMsg);
	}

	if ( BYTESxCHAN == 1)
		scale_factor  = 255.0F;
	else if (BYTESxCHAN == 2)
		scale_factor = 65535.0F;
	else
	{
		char errMsg[2000];
		sprintf(errMsg, "in TiledVolume::loadSubvolume(...): Too many bytes per channel (%d).", BYTESxCHAN);
		throw IOException(errMsg);
	}

	real32 *subvol = new real32[(V1-V0) * (H1-H0) * (D1-D0)]; // this image is an intensity image (only one channel)

	int offset;

	if ( DIM_C == 1 ) {
		for(int i = 0; i < ((V1-V0) * (H1-H0) * (D1-D0)); i++)
			subvol[i] = (real32) data[i]/scale_factor;
	}
	else { // conversion to an intensity image
		if ( iim::CHANNEL_SELECTION == iim::ALL ) {
			char errMsg[2000];
			sprintf(errMsg, "in TiledVolume::loadSubvolume(...): conversion from multi-channel to intensity images not supported.");
			throw IOException(errMsg);
		}
		else if ( iim::CHANNEL_SELECTION == iim::R ) {
			offset = 0;
		}
		else if ( iim::CHANNEL_SELECTION == iim::G ) {
			offset = 1;
		}
		else if ( iim::CHANNEL_SELECTION == iim::B ) {
			offset = 2;
		}
		else {
			char errMsg[2000];
			sprintf(errMsg, "in TiledVolume::loadSubvolume(...): wrong value for parameter iom::CHANNEL_SELECTION.");
			throw IOException(errMsg);
		}
		for(int i = 0; i < ((V1-V0) * (H1-H0) * (D1-D0)); i++)
			subvol[i] = (real32) data[3*i + offset]/scale_factor;
	}

	delete [] data;


	return subvol;
}

//loads given subvolume in a 1-D array of iim::uint8 while releasing stacks slices memory when they are no longer needed
//---03 nov 2011: added color support
iim::uint8* TiledVolume::loadSubvolume_to_UINT8(int V0,int V1, int H0, int H1, int D0, int D1, int *channels, int ret_type ) throw (IOException, iom::exception)
{
    /**/iim::debug(iim::LEV3, strprintf("V0=%d, V1=%d, H0=%d, H1=%d, D0=%d, D1=%d, *channels=%d, ret_type=%d", V0, V1, H0, H1, D0, D1, channels ? *channels : -1, ret_type).c_str(), __iim__current__function__);

    //checking for non implemented features
	//if( this->BYTESxCHAN != 1 ) {
    //	char err_msg[STATIC_STRINGS_SIZE];
	//	sprintf(err_msg,"TiledVolume::loadSubvolume_to_UINT8: invalid number of bytes per channel (%d)",this->BYTESxCHAN); 
	//	throw iim::IOException(err_msg);
	//}

    //if ( (ret_type == iim::DEF_IMG_DEPTH) && ((8 * this->BYTESxCHAN) != iim::DEF_IMG_DEPTH) ) {
		// does not support depth conversion: 
		// return type is 8 bits, but native depth is not 8 bits
    if ( (ret_type != iim::NATIVE_RTYPE) && (ret_type != iim::DEF_IMG_DEPTH) ) {
		// return type should be converted, but not to 8 bits per channel
        char err_msg[STATIC_STRINGS_SIZE];
		sprintf(err_msg,"TiledVolume::loadSubvolume_to_UINT8: non supported return type (%d bits) - native type is %d bits",ret_type, 8*this->BYTESxCHAN); 
        throw IOException(err_msg);
	}

	// reduction factor to be applied to the loaded buffer
    int red_factor = (ret_type == iim::NATIVE_RTYPE) ? 1 : ((8 * this->BYTESxCHAN) / ret_type);

	char *err_rawfmt;

    //initializations
    V0 = V0 < 0 ? 0 : V0;
    H0 = H0 < 0 ? 0 : H0;
    D0 = D0 < 0 ? 0 : D0;
    V1 = (V1 < 0 || V1 > (int)DIM_V) ? DIM_V : V1; // iannello MODIFIED
    H1 = (H1 < 0 || H1 > (int)DIM_H) ? DIM_H : H1; // iannello MODIFIED
    D1 = (D1 < 0 || D1 > (int)DIM_D) ? DIM_D : D1; // iannello MODIFIED
    iim::uint8 *subvol = 0;

    //checking that the interval is valid
    if(V1-V0 <=0 || H1-H0 <= 0 || D1-D0 <= 0)
        throw IOException("in TiledVolume::loadSubvolume_to_UINT8: invalid subvolume intervals");

    //computing dimensions
    iim::sint64 sbv_height = V1 - V0;
    iim::sint64 sbv_width  = H1 - H0;
    iim::sint64 sbv_depth  = D1 - D0;

    //initializing the number of channels with an undefined value (it will be detected from the first slice read)
    iim::sint64 sbv_channels = -1;
	iim::sint64 sbv_bytes_chan = -1;

    //scanning of stacks matrix for data loading and storing into subvol
    Rect_t subvol_area;
    subvol_area.H0 = H0;
    subvol_area.V0 = V0;
    subvol_area.H1 = H1;
    subvol_area.V1 = V1;

    char slice_fullpath[STATIC_STRINGS_SIZE];

    bool first_time = true;

	Segm_t *intersect_segm = BLOCKS[0][0]->Intersects(D0,D1);

	if (intersect_segm) // there is intersection
	{
		for(int row=0; row<N_ROWS; row++) 
		{
			for(int col=0; col<N_COLS; col++)
			{
				Rect_t *intersect_area = BLOCKS[row][col]->Intersects(subvol_area);
				if(intersect_area) // there is intersection
				{
                    //printf("\t\t\t\tin TiledVolume::loadSubvolume_to_UINT8(): using STACK[%d,%d] for area %d-%d(V) x %d-%d(H)\n", row, col, intersect_area->V0-V0, intersect_area->V1-V0, intersect_area->H0-H0, intersect_area->H1-H0);

					for(int k=intersect_segm->ind0; k<=intersect_segm->ind1; k++)
					{
					    //if this is the first time a block is loaded: safely allocating memory for data
					    if(first_time)
					    {
					        first_time = false;
							sbv_channels = this->DIM_C;
							sbv_bytes_chan = this->BYTESxCHAN;

							// 2016-04-20. Alessandro. @FIXED getPlugin3D call when ffmt == "Vaa3DRaw": undefined behaviour
							// 2015-04-15. Alessandro. @FIXED bad/missing exception handling in loadSubvolume_to_UINT8.
                            try
                            {
                                if ( ffmt != "Vaa3DRaw" && sbv_channels >1 && !(iom::IOPluginFactory::getPlugin3D(iom::IMIN_PLUGIN)->isChansInterleaved()) ) {
                                    throw iim::IOException("the plugin do not store channels in interleaved mode: more than one channel not supported yet.");
                                }
                            }
                            catch(iom::exception &ex)
                            {
                                throw iim::IOException(ex.what());
                            }

							try
					        {
                                iim::sint64 sizesubvol = sbv_height * sbv_width * sbv_depth * sbv_channels * sbv_bytes_chan;
                                subvol = new iim::uint8[sizesubvol];
                                memset(subvol, 0, sizesubvol); // 2018-10-23. Yang. @ADDED init by zeros, YY 10/23/2018
					            //if ( !subvol )
								//	throw iim::IOException("in TiledVolume::loadSubvolume_to_UINT8: unable to allocate memory");
					        }
                            catch(...){throw IOException("in TiledVolume::loadSubvolume_to_UINT8: unable to allocate memory");}
					    }
						//loading region
						sprintf(slice_fullpath, "%s/%s/%s", root_dir, BLOCKS[row][col]->getDIR_NAME(), BLOCKS[row][col]->getFILENAMES()[k]);
						
						/* rationale (for V dimension, the same for H)
						 *
						 * V0, ABS_V are the first indices, V1 is the last index + 1
						 *
						 * case V0<intersect_V0 && V1=intersect_V1
						 *
						 *   indices                0   (V1 - ABS_V)
						 *                          |      | 
						 *                          ABS_V  |   ABS_V+SIZE               
						 *   source (file)          |----------|
						 *                          |------|
						 *   buffer             |---|------|
						 *                      V0  |      V1
						 *                          |      |
						 *                intesect_V0      intersect_V1
						 *                          |      |
						 *   indices  (intersect_V0 - V0)  (intersect_V1 - V0)=(V1 - V0)=sbv_height 
						 *
						 *
						 * V0=intersect_V0 && case V1>intersect_V1 
						 *
						 *   indices   (V0 - ABS_V)  SIZE
						 *                     |      |
						 *               ABS_V |   ABS_V+SIZE               
						 *   source (file) |----------|           
						 *                     |------|
						 *   buffer            |------|---|
						 *                     V0     |   V1
						 *                     |      |
						 *           intesect_V0      intersect_V1
						 *                     |      |
						 *   indices           0    (intersect_V1 - V0) 
						 */

						// vertices of file block
						int sV0 = (V0<intersect_area->V0) ? 0 : (V0 - BLOCKS[row][col]->getABS_V());
						int sV1 = (V1>intersect_area->V1) ? BLOCKS[row][col]->getHEIGHT() : (V1 - BLOCKS[row][col]->getABS_V());
						int sH0 = (H0<intersect_area->H0) ? 0 : (H0 - BLOCKS[row][col]->getABS_H());
						int sH1 = (H1>intersect_area->H1) ? BLOCKS[row][col]->getWIDTH() : (H1 - BLOCKS[row][col]->getABS_H());
						int sD0 = (D0<BLOCKS[row][col]->getBLOCK_ABS_D()[k]) ? 0 : (D0 - BLOCKS[row][col]->getBLOCK_ABS_D()[k]);
						int sD1 = (D1>(int)(BLOCKS[row][col]->getBLOCK_ABS_D()[k]+BLOCKS[row][col]->getBLOCK_SIZE()[k])) ? (int)BLOCKS[row][col]->getBLOCK_SIZE()[k] : (int)(D1 - BLOCKS[row][col]->getBLOCK_ABS_D()[k]);

						// vertices of buffer block
						int bV0 = (V0>intersect_area->V0) ? 0 : (int)(intersect_area->V0 - V0);
						//int bV1 = (V1<intersect_area->V1) ? (int)sbv_height : (intersect_area->V1 - V0); // unused
						int bH0 = (H0>intersect_area->H0) ? 0 : (intersect_area->H0 - H0);
						//int bH1 = (H1<intersect_area->H1) ? (int)sbv_width : (int)(intersect_area->H1 - H0); // unused
						int bD0 = (D0>BLOCKS[row][col]->getBLOCK_ABS_D()[k]) ? 0 : (BLOCKS[row][col]->getBLOCK_ABS_D()[k] - D0);
						//int bD1 = (D1<(int)(BLOCKS[row][col]->getBLOCK_ABS_D()[k]+BLOCKS[row][col]->getBLOCK_SIZE()[k])) ? (int)sbv_depth : (BLOCKS[row][col]->getBLOCK_ABS_D()[k]+BLOCKS[row][col]->getBLOCK_SIZE()[k] - D0); // unused

                        // 2018-10-23. Yang. @ADDED support partial data, YY 10/23/2018
                        if(string(slice_fullpath).find("NULL.tif") != string::npos)
                        {
                            continue;
                        }
                        else
                        {
                            if ( (err_rawfmt = fmtMngr->copyFileBlock2Buffer(
                                    slice_fullpath,
                                    sV0,sV1,sH0,sH1,sD0,sD1,
                                    (unsigned char *)subvol,
                                    (int)sbv_bytes_chan, // this is native rtype, it has substituted sizeof(iim::uint8)
                                    bH0+bV0*sbv_width+bD0*sbv_width*sbv_height,
                                    sbv_width,
                                    sbv_width*sbv_height,
                                    sbv_width*sbv_height*sbv_depth) ) != 0 )
                            {
                                char err_msg[STATIC_STRINGS_SIZE];
                                sprintf(err_msg,
                                    "TiledVolume::loadSubvolume_to_UINT8: error in extracting a block from file %s (%s)",
                                    BLOCKS[row][col]->getFILENAMES()[k], err_rawfmt);
                                throw IOException(err_msg);
                            }
                        }
					}
					delete intersect_area;
				}
			}
		}
		delete intersect_segm;
	}
	else
        throw IOException("in TiledVolume::loadSubvolume_to_UINT8: depth interval out of range");
	
	if ( n_active < DIM_C ) { // not all channels are active
		compact_active_chans((sbv_height * sbv_width * sbv_depth * sbv_bytes_chan),subvol);
		sbv_channels = n_active;
	}

    //returning outputs
    if(channels)
        *channels = (int)sbv_channels;

	if ( red_factor > 1 ) { // the buffer has to be reduced

		if ( (err_rawfmt = convert2depth8bits(red_factor,(sbv_height*sbv_width*sbv_depth),sbv_channels,subvol,(iim::uint8 *)0,depth_conv_algo)) != 0  ) {
            char err_msg[STATIC_STRINGS_SIZE];
			sprintf(err_msg,"TiledVolume::loadSubvolume_to_UINT8: %s", err_rawfmt);
            throw IOException(err_msg);
		}
	}
	else if ( (BYTESxCHAN == 1) && (depth_conv_algo & MASK_REMAP_ALGORITHM) ) { // a remap algorithm must be used
		
		if ( (err_rawfmt = remap2depth8bits((sbv_height*sbv_width*sbv_depth),sbv_channels,subvol,depth_conv_algo)) != 0  ) {
            char err_msg[STATIC_STRINGS_SIZE];
			sprintf(err_msg,"TiledVolume::loadSubvolume_to_UINT8: %s", err_rawfmt);
            throw IOException(err_msg);
		}
	}

	return subvol;
}
		

//releases allocated memory of stacks
void TiledVolume::releaseStacks(int first_file, int last_file)
{
    /**/iim::debug(iim::LEV3, strprintf("first_file=%d, last_file=%d", first_file, last_file).c_str(), __iim__current__function__);

    char msg[STATIC_STRINGS_SIZE];
	sprintf(msg,"in TiledVolume::releaseStacks: not implemented yet");
    throw IOException(msg);

	first_file = (first_file == -1 ? 0		: first_file);
	last_file  = (last_file  == -1 ? DIM_D	: last_file);
	//for(int row_index=0; row_index<N_ROWS; row_index++)
	//	for(int col_index=0; col_index<N_COLS; col_index++)
	//		STACKS[row_index][col_index]->releaseStack(first_file,last_file);
}



// OPERATIONS FOR STREAMED SUBVOLUME LOAD 

void *TiledVolume::streamedLoadSubvolume_open ( int steps, iim::uint8 *buf, int V0,int V1, int H0, int H1, int D0, int D1 )
{
    /**/iim::debug(iim::LEV3, strprintf("steps=%d, V0=%d, V1=%d, H0=%d, H1=%d, D0=%d, D1=%d", steps, V0, V1, H0, H1, D0, D1).c_str(), __iim__current__function__);

    //checking for non implemented features
	if( this->BYTESxCHAN != 1 ) {
        char err_msg[STATIC_STRINGS_SIZE];
		sprintf(err_msg,"TiledVolume::streamedLoadSubvolume_open: invalid number of bytes per channel (%d)",this->BYTESxCHAN); 
        throw IOException(err_msg);
	}

    //checks
	if ( V0>=V1 || H0>=H1 || D0>=D1 || V0<0 || H0<0 || D0<0 || V1>(int)DIM_V || H1>(int)DIM_H || D1>(int)DIM_D ) {
		char msg[1000];
		sprintf(msg,"in TiledVolume::streamedLoadSubvolume_open: invalid sub-block vertices");
        throw IOException(msg);
	}

	char *err_rawfmt;

	iim::sint64 stridex = H1 - H0;
	iim::sint64 stridexy = stridex * (V1 - V0);
	iim::sint64 stridexyz = stridexy * (D1 - D0);
	Streamer_Descr_t *stream_descr = new Streamer_Descr_t(buf,sizeof(unsigned char),stridex,stridexy,stridexyz,this->DIM_C,steps);
	if ( !stream_descr ) {
		char msg[1000];
		sprintf(msg,"in TiledVolume::streamedLoadSubvolume_open: Unable to allocate stream descriptor");
        throw IOException(msg);
	}

    // iim::sint64 sbv_channels = this->CHANS; // unused

    //scanning of stacks matrix for data loading and storing into subvol
    Rect_t subvol_area;
    subvol_area.H0 = H0;
    subvol_area.V0 = V0;
    subvol_area.H1 = H1;
    subvol_area.V1 = V1;

    char slice_fullpath[STATIC_STRINGS_SIZE];

	Segm_t *intersect_segm = BLOCKS[0][0]->Intersects(D0,D1);

	if (intersect_segm) // there is intersection
	{
		for(int row=0; row<N_ROWS; row++) 
		{
			for(int col=0; col<N_COLS; col++)
			{
				Rect_t *intersect_area = BLOCKS[row][col]->Intersects(subvol_area);
				if(intersect_area) // there is intersection
				{
                    //printf("\t\t\t\tin TiledVolume::streamedLoadSubvolume_open: using STACK[%d,%d] for area %d-%d(V) x %d-%d(H)\n", row, col, intersect_area->V0-V0, intersect_area->V1-V0, intersect_area->H0-H0, intersect_area->H1-H0);

					for(int k=intersect_segm->ind0; k<=intersect_segm->ind1; k++)
					{
						//loading region
						sprintf(slice_fullpath, "%s/%s/%s", root_dir, BLOCKS[row][col]->getDIR_NAME(), BLOCKS[row][col]->getFILENAMES()[k]);
						
						// vertices of file block
						int sV0 = (V0<intersect_area->V0) ? 0 : (V0 - BLOCKS[row][col]->getABS_V());
						int sV1 = (V1>intersect_area->V1) ? BLOCKS[row][col]->getHEIGHT() : (V1 - BLOCKS[row][col]->getABS_V());
						int sH0 = (H0<intersect_area->H0) ? 0 : (H0 - BLOCKS[row][col]->getABS_H());
						int sH1 = (H1>intersect_area->H1) ? BLOCKS[row][col]->getWIDTH() : (H1 - BLOCKS[row][col]->getABS_H());
						int sD0 = (D0<BLOCKS[row][col]->getBLOCK_ABS_D()[k]) ? 0 : (D0 - BLOCKS[row][col]->getBLOCK_ABS_D()[k]);
						int sD1 = (D1>(int)(BLOCKS[row][col]->getBLOCK_ABS_D()[k]+BLOCKS[row][col]->getBLOCK_SIZE()[k])) ? (int)BLOCKS[row][col]->getBLOCK_SIZE()[k] : (int)(D1 - BLOCKS[row][col]->getBLOCK_ABS_D()[k]);

						// vertices of buffer block
						int bV0 = (V0>intersect_area->V0) ? 0 : (int)(intersect_area->V0 - V0);
						int bH0 = (H0>intersect_area->H0) ? 0 : (intersect_area->H0 - H0);
						int bD0 = (D0>BLOCKS[row][col]->getBLOCK_ABS_D()[k]) ? 0 : (BLOCKS[row][col]->getBLOCK_ABS_D()[k] - D0);

						if ( (err_rawfmt = stream_descr->addSubBlock(
								slice_fullpath,bH0+bV0*stridex+bD0*stridexy,sV0,sV1,sH0,sH1,sD0,sD1)) != 0 ) 
						{
                            char err_msg[STATIC_STRINGS_SIZE];
							sprintf(err_msg,
								"TiledVolume::streamedLoadSubvolume_open: error in adding the block of file %s (%s)", 
								BLOCKS[row][col]->getFILENAMES()[k], err_rawfmt);
                            throw IOException(err_msg);
						}
					}
					delete intersect_area;
				}
			}
		}
		delete intersect_segm;
	}

	if ( (err_rawfmt = streamer_open(stream_descr)) != 0 ) {
		return err_rawfmt;
	}

	return stream_descr;
}

iim::uint8 *TiledVolume::streamedLoadSubvolume_dostep ( void *stream_descr, unsigned char *buffer2 )
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	char *err_rawfmt;
	
	if ( (err_rawfmt = streamer_dostep((Streamer_Descr_t *)stream_descr,buffer2)) != 0 ) {
		char msg[1000];
		sprintf(msg,"in TiledVolume::streamedLoadSubvolume_dostep: Unable to perform next step (%s)",err_rawfmt);
        throw IOException(msg);
	}

	return ((Streamer_Descr_t *)stream_descr)->buf;
}

void TiledVolume::streamedLoadSubvolume_cpydata ( void *stream_descr, unsigned char *buffer, unsigned char *buffer2 )
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	char *err_rawfmt;
	
	if ( (err_rawfmt = streamer_cpydata((Streamer_Descr_t *)stream_descr,buffer,buffer2)) != 0 ) {
		char msg[1000];
		sprintf(msg,"in TiledVolume::streamedLoadSubvolume_dostep: Unable to perform next step (%s)",err_rawfmt);
        throw IOException(msg);
	}
}

iim::uint8 *TiledVolume::streamedLoadSubvolume_close ( void *stream_descr, bool return_buffer )
{
    /**/iim::debug(iim::LEV3, 0, __iim__current__function__);

	iim::uint8 *temp = ((Streamer_Descr_t *)stream_descr)->buf;

	streamer_close((Streamer_Descr_t *)stream_descr);

	delete ((Streamer_Descr_t *)stream_descr);

	if ( return_buffer )
		return temp;
	else {
		delete[] ((Streamer_Descr_t *)stream_descr)->buf;
		return 0;
	}
}

// return vector of tiles along x-y-z (empty vector if the volume is not tiled)
std::vector< iim::voi3D<size_t> > TiledVolume::tilesXYZ()
{
    std::vector< iim::voi3D<size_t> > tiles;

    for(int row_index=0; row_index<N_ROWS; row_index++)
        for(int col_index=0; col_index<N_COLS; col_index++)
            for(int block_index = 0; block_index<BLOCKS[row_index][col_index]->getN_BLOCKS(); block_index++)
                tiles.push_back(
                    iim::voi3D<size_t> (
                        iim::xyz<size_t>(
                            BLOCKS[row_index][col_index]->getABS_H(),
                            BLOCKS[row_index][col_index]->getABS_V(),
                            BLOCKS[row_index][col_index]->getBLOCK_ABS_D()[block_index]
                        ),
                        iim::xyz<size_t>(
                            BLOCKS[row_index][col_index]->getABS_H() + BLOCKS[row_index][col_index]->getWIDTH(),
                            BLOCKS[row_index][col_index]->getABS_V() + BLOCKS[row_index][col_index]->getHEIGHT(),
                            BLOCKS[row_index][col_index]->getBLOCK_ABS_D()[block_index] + BLOCKS[row_index][col_index]->getBLOCK_SIZE()[block_index]
                        )
                    )
                );

    return tiles;
}
