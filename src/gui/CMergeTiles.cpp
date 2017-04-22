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
*       Bria, A., et al., (2012) "Stitching Terabyte-sized 3D Images Acquired in Confocal Ultramicroscopy", Proceedings of the 9th IEEE International Symposium on Biomedical Imaging.
*       Bria, A., Iannello, G., "A Tool for Fast 3D Automatic Stitching of Teravoxel-sized Datasets", submitted on July 2012 to IEEE Transactions on Information Technology in Biomedicine.
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

#include "CMergeTiles.h"
#include "CImportUnstitched.h"
#include <new>
#include <iostream>
#include "StackStitcher.h"
#include "ProgressBar.h"
#include "IM_config.h"
#include "vmBlockVolume.h"
#include "Tiff3DMngr.h"
#include "VolumeConverter.h"

using namespace terastitcher;

CMergeTiles* CMergeTiles::uniqueInstance = NULL;

void CMergeTiles::uninstance()
{
    if(uniqueInstance)
    {
        delete uniqueInstance;
        uniqueInstance = NULL;
    }
}

CMergeTiles::~CMergeTiles()
{
    #ifdef TSP_DEBUG
    printf("TeraStitcher plugin [thread %d] >> CMergeTiles destroyed\n", this->thread()->currentThreadId());
    #endif
}

UnstitchedVolume* CMergeTiles::unstitchedVolume() throw (iim::IOException)
{
	try
	{
		if(!_unstitchedVolume)
			_unstitchedVolume = new UnstitchedVolume(CImportUnstitched::instance()->getVolume());
		return _unstitchedVolume;
	}
	catch (iim::IOException & e)
	{
		throw e;
	}
	catch(...)
	{
		throw iim::IOException("Unhandled exception in CMergeTiles::unstitchedVolume()");
	}
}

//automatically called when current thread is started
void CMergeTiles::run()
{
    #ifdef TSP_DEBUG
    printf("TeraStitcher plugin [thread %d] >> CMergeTiles::run() launched\n", this->thread()->currentThreadId());
    #endif

    try
    {
        //pointer to the <Image4DSimple> object which stores the volume to be shown into Vaa3D
#ifdef VAA3D_TERASTITCHER
        Image4DSimple* img = 0;
#endif

        // check preconditions
        unstitchedVolume();
		if(!pMergeTiles)
			throw iim::IOException("in CMergeTiles::run(): invalid reference to GUI");

#ifdef VAA3D_TERASTITCHER
        V3DPluginCallback* V3D_env = 0;
#endif
        if(pMergeTiles != 0)
        {
			// set libtiff flags
			setLibTIFFcfg(!pMergeTiles->libtiff_uncompressed_checkbox->isChecked(), pMergeTiles->libtiff_bigtiff_checkbox->isChecked());

			// get other info from GUI
			int slice_height =  pMergeTiles->block_height_field->value();
			int slice_width = pMergeTiles->block_width_field->value();
			int slice_depth = pMergeTiles->block_depth_field->value();
			int x0 = pMergeTiles->x0_field->value();
			int x1 = pMergeTiles->x1_field->value()+1;
			int y0 = pMergeTiles->y0_field->value();
			int y1 = pMergeTiles->y1_field->value()+1;
			int z0 = pMergeTiles->z0_field->value();
			int z1 = pMergeTiles->z1_field->value()+1;
			bool restoreSPIM = pMergeTiles->restoreSPIM_cbox->currentIndex() != 0;
			std::string dst_root_dir = pMergeTiles->savedir_field->text().toStdString();
			int restore_direction = pMergeTiles->restoreSPIM_cbox->currentIndex();
			int blending_algo = pMergeTiles->blendingalbo_cbox->currentIndex();
			int img_depth = pMergeTiles->imgdepth_cbox->currentText().section(" ", 0, 0).toInt();
			std::string dst_format = pMergeTiles->vol_format_cbox->currentText().toStdString();
			bool parallel = false;
			bool isotropic = false;
			bool show_progress_bar = true;
			bool timeseries = false;
			bool makeDirs = false;
			bool metaData = false;
			bool halving_method = HALVE_BY_MEAN;
			std::string ch_dir = "";
			std::string mdata_fname = "null";
			std::string outFmt = "RGB";

			// set blending algorithm
			_unstitchedVolume->setBLENDING_ALGO(blending_algo);


			// create volume converter
			VolumeConverter vc;
			vc.setSrcVolume(_unstitchedVolume, outFmt.c_str());
			vc.setSubVolume(y0, y1, x0, x1, z0, z1);


			// make conversion (code taken from teraconverter.cpp)
			if ( dst_format == iim::SIMPLE_RAW_FORMAT )
			{
					vc.generateTilesSimple(dst_root_dir.c_str(),resolutions,
						slice_height,slice_width,halving_method,isotropic,
						show_progress_bar,"raw",img_depth,"",parallel);
			}
			else if ( dst_format == iim::SIMPLE_FORMAT )
				if ( timeseries ) {
					vc.convertTo(dst_root_dir.c_str(),dst_format,img_depth,true,resolutions,
						slice_height,slice_width,slice_depth,halving_method);
				}
				else if ( makeDirs ) {
					vc.createDirectoryHierarchySimple(dst_root_dir.c_str(),resolutions,
						slice_height,slice_width,-1,halving_method,isotropic,
						show_progress_bar,"tif",img_depth,"",parallel);
					// 				vc.createDirectoryHierarchy(dst_root_dir.c_str(),ch_dir,resolutions,
					// 					slice_height,slice_width,-1,halving_method,isotropic,
					// 					show_progress_bar,"tif",img_depth,"",parallel);
				}
				else if ( metaData ) {
					//vc.mdataGenerator(dst_root_dir.c_str(),resolutions,
					//	slice_height,slice_width,-1,halving_method,isotropic,
					//	show_progress_bar,"tif",img_depth,"",parallel);
				}
				else {
					vc.generateTilesSimple(dst_root_dir.c_str(),resolutions,
						slice_height,slice_width,halving_method,isotropic,
						show_progress_bar,"tif",img_depth,"",parallel);
				}
			else if ( dst_format == iim::STACKED_RAW_FORMAT )
				if ( timeseries ) {
					vc.convertTo(dst_root_dir.c_str(),dst_format,img_depth,true,resolutions,
						slice_height,slice_width,slice_depth,halving_method);
				}
				else if ( makeDirs ) {
					vc.createDirectoryHierarchy(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,-1,halving_method,isotropic,
						show_progress_bar,"raw",img_depth,"",parallel);
				}
				else if ( metaData ) {
					vc.mdataGenerator(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,-1,halving_method,isotropic,
						show_progress_bar,"raw",img_depth,"",parallel);
				}
				else {
					vc.generateTiles(dst_root_dir.c_str(),resolutions,
						slice_height,slice_width,halving_method,isotropic,
						show_progress_bar,"raw",img_depth,"",parallel);
				}
			else if ( dst_format == iim::STACKED_FORMAT )
				if ( timeseries ) {
					vc.convertTo(dst_root_dir.c_str(),dst_format,img_depth,true,resolutions,
						slice_height,slice_width,slice_depth,halving_method);
				}
				else if ( makeDirs ) {
					vc.createDirectoryHierarchy(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,-1,halving_method,isotropic,
						show_progress_bar,"tif",img_depth,"",parallel);
				}
				else if ( metaData ) {
					vc.mdataGenerator(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,-1,halving_method,isotropic,
						show_progress_bar,"tif",img_depth,"",parallel);
				}
				else {
					vc.generateTiles(dst_root_dir.c_str(),resolutions,
						slice_height,slice_width,halving_method,isotropic,
						show_progress_bar,"tif",img_depth,"",parallel);
				}
			else if ( dst_format == iim::TILED_FORMAT ) {
				if ( timeseries ) {
					vc.convertTo(dst_root_dir.c_str(),dst_format,img_depth,true,resolutions,
						slice_height,slice_width,slice_depth,halving_method);
				}
				else if ( makeDirs ) {
					vc.createDirectoryHierarchy(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Vaa3DRaw",img_depth,"",parallel);
				}
				else if ( metaData ) {
					vc.mdataGenerator(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Vaa3DRaw",img_depth,"",parallel);
				}
				else {
					vc.generateTilesVaa3DRaw(dst_root_dir.c_str(),resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Vaa3DRaw",img_depth,"",parallel);
				}
			}
			else if ( dst_format == iim::TILED_TIF3D_FORMAT ) {
				if ( timeseries ) {
					vc.convertTo(dst_root_dir.c_str(),dst_format,img_depth,true,resolutions,
						slice_height,slice_width,slice_depth,halving_method);
				}
				else if ( makeDirs ) {
					vc.createDirectoryHierarchy(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Tiff3D",img_depth,"",parallel);
				}
				else if ( metaData ) {
					vc.mdataGenerator(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Tiff3D",img_depth,"",parallel);
				}
				else {
					vc.generateTilesVaa3DRaw(dst_root_dir.c_str(),resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Tiff3D",img_depth,"",parallel);
				}
			}
			else if ( dst_format == iim::TILED_MC_FORMAT )
				if ( timeseries ) {
					vc.convertTo(dst_root_dir.c_str(),dst_format,img_depth,true,resolutions,
						slice_height,slice_width,slice_depth,halving_method);
				}
				else if ( makeDirs ) {
					vc.createDirectoryHierarchy(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Vaa3DRawMC",img_depth,"",parallel);
				}
				else if ( metaData ) {
					vc.mdataGenerator(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Vaa3DRawMC",img_depth,"",parallel);
				}
				else {
					vc.generateTilesVaa3DRawMC(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Vaa3DRaw",img_depth,"",false);
				}
			else if ( dst_format == iim::TILED_MC_TIF3D_FORMAT )
				if ( timeseries ) {
					vc.convertTo(dst_root_dir.c_str(),dst_format,img_depth,true,resolutions,
						slice_height,slice_width,slice_depth,halving_method);
				}
				else if ( makeDirs ) {
					vc.createDirectoryHierarchy(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Tiff3DMC",img_depth,"",parallel);
				}
				else if ( metaData ) {
					vc.mdataGenerator(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Tiff3DMC",img_depth,"",parallel);
				}
				else {
					vc.generateTilesVaa3DRawMC(dst_root_dir.c_str(),ch_dir,resolutions,
						slice_height,slice_width,slice_depth,halving_method,isotropic,
						show_progress_bar,"Tiff3D",img_depth,"",parallel);
				}
			else if ( dst_format == iim::BDV_HDF5_FORMAT )
				vc.generateTilesBDV_HDF5(dst_root_dir.c_str(),resolutions,
				slice_height,slice_width,slice_depth,halving_method,
				show_progress_bar,"Fiji_HDF5",img_depth);
			else if ( dst_format == iim::IMS_HDF5_FORMAT )
				vc.generateTilesIMS_HDF5(dst_root_dir.c_str(),mdata_fname,resolutions,
				slice_height,slice_width,slice_depth,halving_method,
				show_progress_bar,"Fiji_HDF5",img_depth);

            //if a resolution has be selected to be shown in Vaa3D, it is necessary to load each slice and to create a new Image4DSimple object
#ifdef VAA3D_TERASTITCHER
			if(resolution_index_vaa3D != -1)
            {
                //updating progress bar message
                ts::ProgressBar::instance()->start("Loading volume into Vaa3D...");

                //retrieving the directory path where the selected volume is stored
                QString volpath = pMergeTiles->savedir_field->text();
                int height = (stitcher.getV1()-stitcher.getV0())/pow(2.0f,resolution_index_vaa3D);
                int width = (stitcher.getH1()-stitcher.getH0())/pow(2.0f,resolution_index_vaa3D);
                int depth = (stitcher.getD1()-stitcher.getD0())/pow(2.0f,resolution_index_vaa3D);
                volpath.append("/RES(").append(QString::number(height));
                volpath.append("x").append(QString::number(width));
                volpath.append("x").append(QString::number(depth));
                volpath.append(")");
                QDir voldir(volpath);
                QStringList first_level_list = voldir.entryList(QDir::AllDirs | QDir::NoDotAndDotDot, QDir::Name);
                voldir.setPath(voldir.path().append("/").append(first_level_list.first()));
                QStringList second_level_list = voldir.entryList(QDir::AllDirs | QDir::NoDotAndDotDot, QDir::Name);
                voldir.setPath(voldir.path().append("/").append(second_level_list.first()));
                QStringList slices_list = voldir.entryList(QDir::Files);

                // use Vaa3D to load all slices
                std::vector<Image4DSimple*> slices;
                for (int k = 0; k < slices_list.size(); k++)
                    slices.push_back(V3D_env->loadImage(const_cast<char*>(QString(voldir.path().append("/").append(slices_list.at(k).toLocal8Bit().constData())).toStdString().c_str())));

                // allocate image data assuming all slices have the same X, Y, C dimensions and bitdepth
                img = new Image4DSimple();
                img->setFileName(voldir.path().toStdString().c_str());
                V3DLONG slice_dims = slices[0]->getXDim()*slices[0]->getYDim()*slices[0]->getCDim()*slices[0]->getUnitBytes();
                unsigned char* img_data = new iom::uint8[slice_dims * slices_list.size()];

                // copy each loaded slice into the volume
                for (int k = 0; k < slices.size(); k++)
                    for(V3DLONG pc = 0; pc < slice_dims; pc++)
                        img_data[k*slice_dims + pc] = slices[k]->getRawData()[pc];

                // set image data
                img->setData(img_data, slices[0]->getXDim(), slices[0]->getYDim(), slices_list.size(), slices[0]->getCDim(), slices[0]->getDatatype());

                // deallocate data
                for (int k = 0; k < slices.size(); k++)
                    delete slices[k];

                // set image attributes
                img->setRezX(volume->getVXL_H()*pow(2.0f,resolution_index_vaa3D));
                img->setRezY(volume->getVXL_V()*pow(2.0f,resolution_index_vaa3D));
                img->setRezZ(volume->getVXL_D()*pow(2.0f,resolution_index_vaa3D));
                img->setOriginX(stitcher.getMultiresABS_H(resolution_index_vaa3D, stitcher.getH0())/1000.0f);
                img->setOriginY(stitcher.getMultiresABS_V(resolution_index_vaa3D, stitcher.getV0())/1000.0f);
                img->setOriginZ(stitcher.getMultiresABS_D(resolution_index_vaa3D, stitcher.getD0())/1000.0f);
            }
#endif
        }

        //everything went OK
#ifdef VAA3D_TERASTITCHER
        emit sendOperationOutcome(0, img);
#else
		emit sendOperationOutcome(0);
#endif
    }
    catch( iim::IOException& exception)
    {
        /**/ts::warning(strprintf("exception thrown in CMergeTiles::run(): \"%s\"", exception.what()).c_str());
        emit sendOperationOutcome(new iom::exception(exception.what()));
    }
    catch( iom::exception& exception)
    {
        /**/ts::warning(strprintf("exception thrown in CMergeTiles::run(): \"%s\"", exception.what()).c_str());
        emit sendOperationOutcome(new iom::exception(exception.what()));
    }
    catch(const char* error)
    {
        /**/ts::warning(strprintf("exception thrown in CMergeTiles::run(): \"%s\"", error).c_str());
        emit sendOperationOutcome(new iom::exception(error));
    }
    catch(std::bad_alloc& ba)
    {
        /**/ts::warning(strprintf("exception thrown in CMergeTiles::run(): \"%s\"", ba.what()).c_str());
        emit sendOperationOutcome(new iom::exception(ba.what()));
    }
    catch(...)
    {
        /**/ts::warning(strprintf("exception thrown in CMergeTiles::run(): \"%s\"", "Generic error").c_str());
        emit sendOperationOutcome(new iom::exception("Unable to determine error's type"));
    }
}

//reset method
void CMergeTiles::reset()
{
    #ifdef TSP_DEBUG
    printf("TeraStitcher plugin [thread %d] >> CMergeTiles::reset()\n", this->thread()->currentThreadId());
    #endif

    for(int i=0; i<S_MAX_MULTIRES; i++)
        resolutions[i] = i==0;
    pMergeTiles = 0;
	if(_unstitchedVolume)
	{
		delete _unstitchedVolume;
		_unstitchedVolume = 0;
	}
}