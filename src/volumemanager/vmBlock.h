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
* 2017-07-06. Giulio.     @ADDED method 'loadImageStack2' to enable selective reads of data
* 2016-11-14. Giulio.     @CHANGED interface of constructor from xml to manage of the case when z_end is invalid (i.e. when import is from an xml import file generated externally
* 2015-07-22. Giluio.     @ADDED supporto for spase data.
* 2015-01-17. Alessandro. @ADDED constructor for initialization from XML.
* 2015-01-17. Alessandro. @FIXED missing throw(iom::exception) declaration in many methods.
* 2014-09-05. Alessandro. @ADDED 'z_end' parameter in 'loadXML()' method to support sparse data feature.
* 2014-08-25. Alessandro. @ADDED missing 'throw (iom::iom::exception)' statement in the 'loadImageStack()' method's signature
*/

#ifndef _VM_BLOCK_H
#define _VM_BLOCK_H

#include <set>
#include <string>

#include "volumemanager.config.h"
#include "iomanager.config.h"
#include "tinyxml.h"
#include "vmVirtualStack.h"



class Displacement;

//TYPE DEFINITIONS
//structure representing a substack
//D0: first slice, D1: last slice, ind0: index of 1st block (containing D0), ind1: index of last block (containing D1) 
typedef struct {int D0, D1, ind0, ind1;} Segm_t;

class vm::Block : public vm::VirtualStack
{
	private:

		vm::BlockVolume* CONTAINER;					//pointer to <VirtualVolume> object that contains the current object
		int          N_BLOCKS;                   //number of blocks along z
		int         *BLOCK_SIZE;                 //dimensions of blocks along z
		int         *BLOCK_ABS_D;                //absolute D voxel coordinates of blocks

		std::string series_no;  // index into multi-stack files
		/* default (invalid) value = -1; to be used when multiple stack are in the same  file (3D formats)
		 * is initialized and actually used only when import is performed from an xml import file
		 */

		//******** OBJECT PRIVATE METHODS *********
        Block(void){}

		//Initializes all object's members given DIR_NAME
		void init() ;

	    //binarizing-unbinarizing methods
		void binarizeInto(FILE* file) ;
		void unBinarizeFrom(FILE* file) ;

		// compute 'z_ranges'
		void 
			compute_z_ranges(
			std::pair<int,int> const * z_coords = 0)		// set of z-coordinates where at least one slice (of a certain stack) is available
		;								// if null, 'z_ranges' will be compute based on 'FILENAMES' vector

		//returns a pointer to the intersection segment (along D) if the given segment (D0,D1-1) intersects current stack, otherwise returns NULL
		//D0 first index of the segment
		//D1 last index of the segment + 1
		Segm_t* Intersects(int D0, int D1);

		//******** FRIEND CLASS DECLARATION *********
		//BlockVolume can access Block private members and methods
		friend class vm::BlockVolume;
	
	public:
		//CONSTRUCTORS
		Block(vm::BlockVolume* _CONTAINER, int _ROW_INDEX, int _COL_INDEX, const char* _DIR_NAME) ;				// build from scratch
		Block(vm::BlockVolume* _CONTAINER, int _ROW_INDEX, int _COL_INDEX, FILE* bin_file) ;						// build from mdata.bin
		Block(vm::BlockVolume* _CONTAINER, int _ROW_INDEX, int _COL_INDEX, TiXmlElement* stack_node, int &z_end) ;	// build from XML
		/* parameter z_end is passed by reference because it can contain an invalid value (when the xml import file is externally generated)
		 * in this case the constructor must initialize the parameter
		 */

		~Block(void);

		//GET methods
		int  getN_BLOCKS()		{return N_BLOCKS;}

		int  *getBLOCK_SIZE()   {return BLOCK_SIZE;}
		int  *getBLOCK_ABS_D()  {return BLOCK_ABS_D;}
		
        void *getCONTAINER()    {return CONTAINER;}

		//LOAD and RELEASE methods
        iom::real_t* loadImageStack(int first_file=-1, int last_file=-1) ;
        iom::real_t* loadImageStack2(int first_file=-1, int last_file=-1, int V0=-1, int V1=-1, int H0=-1, int H1=-1) ;
		void releaseImageStack();

		//XML methods
		TiXmlElement* getXML();
		void loadXML(
			TiXmlElement *stack_node,
			int z_end)					// 2014-09-05. Alessandro. @ADDED 'z_end' parameter to support sparse data feature
										//			   Here 'z_end' identifies the range [0, z_end) that slices can span
		;
};

#endif //_BLOCK_H
