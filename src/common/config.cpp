#include "config.h"

#ifdef WITH_QT
#include <QtGui>
#endif

/*******************************************************************************************************************************
 *   Interfaces, types, parameters and constants   													       *
 *******************************************************************************************************************************/
namespace terastitcher
{
    /*******************
    *    PARAMETERS    *
    ********************
    ---------------------------------------------------------------------------------------------------------------------------*/
	std::string version = num2str<int>(TERASTITCHER_MAJOR) + "." + num2str<int>(TERASTITCHER_MINOR) + "." + num2str<int>(TERASTITCHER_PATCH);
    int DEBUG = NO_DEBUG;                    //debug level
#ifdef WITH_QT
	std::string qtversion = QT_VERSION_STR;
#else	
	std::string qtversion = "";
#endif
    /*-------------------------------------------------------------------------------------------------------------------------*/
}