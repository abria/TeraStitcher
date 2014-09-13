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

#ifndef DISPLACEMENT_MIP_NCC_H
#define DISPLACEMENT_MIP_NCC_H

#include "Displacement.h"
//#include "volumemanager.config.h"
#include "CrossMIPs.h"
#include <math.h>

class DisplacementMIPNCC : Displacement
{
	protected:

		//int TYPE;					//INHERITED from Displacement
		//int VHD_def_coords[3];    //INHERITED from Displacement
		int    VHD_coords[3];       //VHD displacements as offsets of the second stack with respect to the first one
		float  NCC_maxs[3];			//first reliability feature: measures the NCC peak, range in [0=unreliable,1=reliable]
		int    NCC_widths[3];		//second reliability feature: measures the NCC shape width sampled at the 75% of
									//the peak, range in [0=reliable, S_NCC_WIDTH_MAX=unreliable]
		float rel_factors[3];		//reliability factors of VHD displacements, range in [0=unreliable,1=reliable] 

		// Alessandro - 31/05/2013 - storing also MIPNCC parameters
		int delays[3];				//half range of NCC search around initial offset along VHD dimensions
        int wRangeThrs[3];			//ranges used to evaluate maximum NCC widths (when maximum width is greater or equal to
									//this value, width is set to invWidth
        int invWidths[3];			//see <wRangeThrs>

	public:

		//COSTRUCTORS - DECONSTRUCTOR
		DisplacementMIPNCC(void);
		DisplacementMIPNCC(const DisplacementMIPNCC &ex_instance);
		DisplacementMIPNCC(NCC_descr_t &mip_ncc_descr);
                DisplacementMIPNCC(int Vnominal, int Hnominal, int Dnominal);
		DisplacementMIPNCC(TiXmlElement *displ_node);
        ~DisplacementMIPNCC(void){}

		//evaluates displacement reliability possibly along the given direction. The result(s) should be stored
		//in one or more object members, so that they have to be computed once and then accessed by GET methods
		float			 evalReliability(direction _direction=invalid)	throw (iom::exception);

		//returns the reliability possibly along the given direction. An exception is thrown if the reliability
		//index(es) are not computed yet.
		float			 getReliability(direction _direction=invalid)	throw (iom::exception);

		//returns the displacement along the given direction
		int				 getDisplacement(direction _direction)			throw (iom::exception);

		//sets to default values the displacements with a reliability factor above the given threshold
		void			 threshold(float rel_threshold)					throw (iom::exception);

		//returns the displacement mirrored along the given direction.
		Displacement*	 getMirrored(direction _direction)				throw (iom::exception);

		//XML methods: convert/load displacement object into/from XML schema
		TiXmlElement*	 getXML() 								throw (iom::exception);
		void			 loadXML(TiXmlElement *displ_node)				throw (iom::exception);

		//combines the parameters of the current and the given displacement so that after this operation
		//the two displacements are more reliable (and are EQUAL).
		void			 combine(Displacement& displ)					throw (iom::exception);

		friend class PDAlgoMIPNCC;
};

#endif /* DISPLACEMENT_MIP_NCC_H */

