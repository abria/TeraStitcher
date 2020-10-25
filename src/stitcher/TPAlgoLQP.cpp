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
* 2019-04-25. Giulio      @ADDED the possibility to use the more general environment variable __PYSCRIPTS_PATH__ to find the Python script implementing LQP optimization
* 2019-04-25. Giulio      @FIXED find the maximum "delay" to be written in the auxiliary file to deal with initial empty stacks
* 2018-04-14. Giulio.     @CREATED global optimization algorithm based on Linear Quadtratic Programming + heuristics (S_FATPM_LQP_HE)
*/


#include <limits>
#include "TPAlgoLQP.h"
#include "S_config.h"
#include "volumemanager.config.h"
#include "vmStackedVolume.h"
#include "vmVirtualStack.h"
#include "Displacement.h"
#include "DisplacementMIPNCC.h"
#include "IM_config.h"

using namespace volumemanager;
using namespace iomanager;

//triplet data type definition
typedef struct 							
{
	float V;
	float H;
	float D;
}triplet;

static std::string python_path;

TPAlgoLQP::TPAlgoLQP(VirtualVolume * _volume) : TPAlgo(S_FATPM_LQP_HE, _volume)
{
	#if S_VERBOSE>4
	printf("........in TPAlgoLQP::TPAlgoLQP(VirtualVolume * _volume)");
	#endif
	
	// set optional parameters
	char *temp;
	if ( (temp = getenv("__LQP_PATH__")) ) {
		python_path = temp;
	}
	else if ( (temp = getenv("__PYSCRIPTS_PATH__")) ) {
		python_path = temp;
	}
	else
		python_path = "";

	if ( python_path != "" ) {
		if ( python_path[python_path.size()-1] != '/' && python_path[python_path.size()-1] != '\\' ) {
			if ( python_path.find('\\') != std::string::npos )
				python_path += "\\";
			else 
				python_path += "/";
		}
	}
}

/*************************************************************************************************************
* Finds the optimal tile placement on the <volume> object member via Linear Quadratic Programming + heuristics.
* PROs: 
* CONs: 
**************************************************************************************************************/
void TPAlgoLQP::execute() throw (iom::exception)
{
	#if S_VERBOSE > 2
	printf("....in TPAlgoLQP::execute()");
	#endif
	

	// generate the input file for python code

	FILE *fout_temp;
	std::string base_path;
	std::string fout_name;
	std::string fin_name;
	std::string slash = "";

	base_path = volume->getSTACKS_DIR();
	if ( base_path[base_path.size()-1] != '/' && base_path[base_path.size()-1] != '\\' ) {
		if ( base_path.find('\\') != std::string::npos )
			slash = "\\";
		else 
			slash = "/";
	}
	fout_name = base_path + slash +  S_FATPM_LQP_HE_FILE_TEMP_NAME;
	fin_name = base_path + slash + S_FATPM_LQP_HE_FILE_SOL_NAME;

	if ( (fout_temp = fopen(fout_name.c_str(),"w")) == 0 ) {
		throw iom::exception(iom::strprintf("in TPAlgoLQP::execute: cannot open temporary file \"%s\"", fout_name.c_str()));
	}

	int rows = volume->getN_ROWS();
	int cols = volume->getN_COLS();
	int n_constraints = (rows-1) * (cols-1);
	int n_vars = rows * (cols-1) + (rows-1) * cols;

	fprintf(fout_temp,"%d ",rows);
	fprintf(fout_temp,"%d ",cols);
	fprintf(fout_temp,"%d ",n_constraints);
	fprintf(fout_temp,"%d ",n_vars);
	fprintf(fout_temp,"\n");

	int *constr = new int[n_vars];
	for ( int i=0, v_ind=(rows * (cols-1) - 1); i<n_constraints; i++, v_ind++ ) {
		memset(constr,0,n_vars*sizeof(int));
		if ( i%(cols-1) == 0 )
			v_ind++;
		constr[i           ] =  1;
		constr[v_ind       ] = -1;
		constr[v_ind + 1   ] =  1;
		constr[i + (cols-1)] = -1;
		for ( int j=0; j<n_vars; j++ )
			fprintf(fout_temp,"%d ",constr[j]);
		fprintf(fout_temp,"\n");
	}
	delete []constr;

	for ( int d=0; d < 3; d++ ) {

		int r;
		int c;

		// 2019-04-25. Giulio. look for maximum delay for dealing with empty tiles that set the delay to -1
		int delay = 0;
		for ( int j=0; j<n_vars; j++ ) {
			r = (j < (rows * (cols-1))) ? (j / (cols-1)) : ((j-(rows * (cols-1))) / cols);
			c = (j < (rows * (cols-1))) ? (j % (cols-1)) : ((j-(rows * (cols-1))) % cols);
			// WARNING: eliminated multiplication by 2 that was in the original statement below because it should be wrong
			delay = std::max<int>(delay,((DisplacementMIPNCC *) volume->getSTACKS()[0][0]->getEAST()[0])->getDelays(direction(d)));
		}
		fprintf(fout_temp,"%d \n",delay);
		//fprintf(fout_temp,"%d \n",(2 * ((DisplacementMIPNCC *) volume->getSTACKS()[0][0]->getEAST()[0])->getDelays(direction(d))));

		for ( int j=0; j<n_vars; j++ ) {
			if ( j < (rows * (cols-1)) ) 
				fprintf(fout_temp,"%d ",volume->getSTACKS()[0][0]->getEAST()[0]->getDefaultDisplacement(direction(d)));
			else
				fprintf(fout_temp,"%d ",volume->getSTACKS()[0][0]->getSOUTH()[0]->getDefaultDisplacement(direction(d)));
		}
		fprintf(fout_temp,"\n");

		for ( int j=0; j<n_vars; j++ ) {
			r = (j < (rows * (cols-1))) ? (j / (cols-1)) : ((j-(rows * (cols-1))) / cols);
			c = (j < (rows * (cols-1))) ? (j % (cols-1)) : ((j-(rows * (cols-1))) % cols);
			if ( j < (rows * (cols-1)) ) 
				fprintf(fout_temp,"%d ",volume->getSTACKS()[r][c]->getEAST()[0]->getDisplacement(direction(d)));
			else
				fprintf(fout_temp,"%d ",volume->getSTACKS()[r][c]->getSOUTH()[0]->getDisplacement(direction(d)));
		}
		fprintf(fout_temp,"\n");

		for ( int j=0; j<n_vars; j++ ) {
			r = (j < (rows * (cols-1))) ? (j / (cols-1)) : ((j-(rows * (cols-1))) / cols);
			c = (j < (rows * (cols-1))) ? (j % (cols-1)) : ((j-(rows * (cols-1))) % cols);
			if ( j < (rows * (cols-1)) ) 
				fprintf(fout_temp,"%f ",volume->getSTACKS()[r][c]->getEAST()[0]->getReliability(direction(d)));
			else
				fprintf(fout_temp,"%f ",volume->getSTACKS()[r][c]->getSOUTH()[0]->getReliability(direction(d)));
		}
		fprintf(fout_temp,"\n");
	}

	fclose(fout_temp);

	// call python code for performing the global optimization

	std::string python_cmd;
	python_cmd = "python " + python_path + "LQP_HE.py -s=" + fout_name + " -d=" + fin_name + " -noxml=True";
	//printf("---> %s\n",python_cmd.c_str());
	if ( system(python_cmd.c_str()) != 0 ) 
		throw iom::exception(iom::strprintf("in TPAlgoLQP::execute(): cannot execute python command \"%s\" \n\nTry to set the environment variable __LQP_PATH__ to the path of the python script implementing the LQP algorithm", python_cmd.c_str()));

	// read the output file generated by python code

	FILE *fin_temp;

	if ( (fin_temp = fopen(fin_name.c_str(),"r")) == 0 ) {
		throw iom::exception(iom::strprintf("in TPAlgoLQP::execute: cannot open temporary file \"%s\"", fin_name.c_str()));
	}
	
	int *displ = new int[n_vars];

	for ( int d=0; d < 3; d++ ) {

		for ( int j=0; j<n_vars; j++ ) {
			fscanf(fin_temp,"%d",displ+j);
		}

		int pos = 0;
		volume->getSTACKS()[0][0]->setABS(0,direction(d));
		for ( int r=0; r<(rows-1); r++ ) {
			for ( int c=0; c<(cols-1); c++ ) {
				volume->getSTACKS()[r][c+1]->setABS(volume->getSTACKS()[r][c]->getABS(direction(d)) + displ[r*(cols-1)+c],direction(d));
			}
			volume->getSTACKS()[r+1][0]->setABS(pos + displ[(rows * (cols-1)) + r * cols],direction(d));
			pos = volume->getSTACKS()[r+1][0]->getABS(direction(d));
		}
		for ( int c=0; c<(cols-1); c++ ) {
			volume->getSTACKS()[rows-1][c+1]->setABS(volume->getSTACKS()[rows-1][c]->getABS(direction(d)) + displ[(rows-1)*(cols-1)+c],direction(d));
		}
	}

	fclose(fin_temp);

#ifdef _WIN32
	if ( base_path.find('/') != std::string::npos ) {
		std::replace(fin_name.begin(), fin_name.end(), '/', '\\');
		std::replace(fout_name.begin(), fout_name.end(), '/', '\\');
	}
#endif
	iim::delete_file(fin_name.c_str());
	iim::delete_file(fout_name.c_str());
}
