#ifndef TIMESERIES_H
#define TIMESERIES_H

#include <vector>
#include "VirtualVolume.h"


class TimeSeries : public iim::VirtualVolume
{
    protected:

        std::vector<VirtualVolume*> frames;     // each time frame corresponds to a complete volumetric image

        // pure virtual method inherithed from abstract class <VirtualVolume>
        void initChannels() ;

    public:

        //CONSTRUCTORS-DESTRUCTOR
        TimeSeries(void){}
        TimeSeries(const char* rootDir, std::string frames_format = "") ;
        virtual ~TimeSeries(void) ;

        // get methods
        std::vector<VirtualVolume*> getFrames() {return frames;}
        std::vector<VirtualVolume*> getActiveFrames() ;
        VirtualVolume* getFrameAt(int t) ;
        iim::uint32 getNFrames(){return static_cast<iim::uint32>(frames.size());}

        // returns a unique ID that identifies the volume format
        std::string getPrintableFormat(){return iim::TIME_SERIES;}

        // added by Alessandro on 2014-02-18: additional info on the reference system (where available)
        float getVXL_1() {return frames.empty() ? 0 : frames[0]->getVXL_1();}
        float getVXL_2() {return frames.empty() ? 0 : frames[0]->getVXL_2();}
        float getVXL_3() {return frames.empty() ? 0 : frames[0]->getVXL_3();}
        iim::axis getAXS_1() {return frames.empty() ? iim::axis_invalid : frames[0]->getAXS_1();}
        iim::axis getAXS_2() {return frames.empty() ? iim::axis_invalid : frames[0]->getAXS_2();}
        iim::axis getAXS_3() {return frames.empty() ? iim::axis_invalid : frames[0]->getAXS_3();}

        // @ADDED by Alessandro on 2016-12-19
        // return true if the given dimension is tiled
        virtual bool isTiled(iim::dimension d) {return false;}
        // return vector of tiles along x-y-z (empty vector if the volume is not tiled)
        virtual std::vector< iim::voi3D<size_t> > tilesXYZ() {return std::vector< iim::voi3D<size_t> >();}

        // @OVERRIDE
        iim::uint32* getActiveChannels(){ if(!frames.empty()) return frames[0]->getActiveChannels(); else return 0;}
        int getNACtiveChannels() { if(!frames.empty()) return frames[0]->getNACtiveChannels(); else return 0;}

		// set active channels (@OVERRIDES VirtualVolume.h by Alessandro on 2014-02-23)
		// WARNING: caller loses ownership of array '_active' 
		void setActiveChannels ( iim::uint32 *_active, int _n_active );


        // pure virtual methods inherithed from abstract class <VirtualVolume>
        iim::real32 *loadSubvolume_to_real32(int V0,int V1, int H0, int H1, int D0, int D1)  ;
        iim::uint8 *loadSubvolume_to_UINT8(int V0=-1,int V1=-1, int H0=-1, int H1=-1, int D0=-1, int D1=-1,
                                                   int *channels=0, int ret_type=iim::DEF_IMG_DEPTH) ;
};

#endif // TIMESERIES_H
