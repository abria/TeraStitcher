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
* 2019-09-21. Giulio.     @CREATED shared file .inc between bioformats2D and bioformats3D plugins
* 2019-09-13. Giulio.     @ADDED management of endianness: when the data endianness and machine endianness differ bytes must be swapped
* 2018-08-24. Giulio.     @ADDED method 'isInputOnly'
* 2016-05-21. Giulio.     @FIXED the case of non RGB images that on reading return one channel at the time
* 2016-05-21. Giulio.     @CHANGED moved the creation of the reader and the file opening to global variables for effeciency
* 2016-05-11. Giulio.     @CREATED
*/

static JavaVM *jvm;       /* denotes a Java VM */
static JNIEnv *env;       /* pointer to native method interface */
static bool jvm_setup = false;

static jclass ImageReaderClass;
static jobject reader;

static jbyteArray staticByteArray;
//static jbyte *staticJavaPointer;
static bool buffer_setup = false;

static std::string cur_img_path = "";


static
jint setupJVM ( JavaVM **jvm, JNIEnv **env ) {

    JavaVMInitArgs vm_args; /* JDK/JRE 6 VM initialization arguments */

    JavaVMOption* options = new JavaVMOption[1];
    
	// set optional parameters
	std::string bioformats_path;
	char *temp;
	if ( (temp = getenv("__BIOFORMATS_PATH__")) ) {
		bioformats_path = temp;
	}
	else
		throw iom::exception(iom::strprintf("environment variable __BIOFORMATS_PATH__ not defined"), __iom__current__function__);

	if ( bioformats_path[bioformats_path.size()-1] != '/' && bioformats_path[bioformats_path.size()-1] != '\\' ) {
		if ( bioformats_path.find('\\') != std::string::npos )
			bioformats_path += "\\";
		else 
			bioformats_path += "/";
	}
	std::string jar_path = "-Djava.class.path=" + bioformats_path + "bioformats_package.jar";
    options[0].optionString = (char *) jar_path.c_str();

// #ifdef WIN32
//     options[0].optionString = "-Djava.class.path=Z:\\Windows\\bioformat\\bioformats_package.jar";
// #else
//     options[0].optionString = "-Djava.class.path=/Users/iannello/Home/Windows/bioformat/bioformats_package.jar";
// #endif

    vm_args.version = JNI_VERSION_1_6;
    vm_args.nOptions = 1;
    vm_args.options = options;
    vm_args.ignoreUnrecognized = false;

    /* load and initialize a Java VM, return a JNI interface
     * pointer in env */
    jint rc = JNI_CreateJavaVM(jvm, (void **)env, &vm_args);
    delete[] options;

	//jint ver = (*env)->GetVersion();

	return rc;
}

/* function to check endianness of the machine 
 * returns: 
 *  L for little endian machines
 *  B for big endian machines
 *  M for machines that swap 16-bits words
 *  N for unknown endianness
 */
static 
char checkMachineEndian()
{
    char e='N'; //for unknown endianness
    int a[2];

	if ( sizeof(int) == 2 )
    	*((long *)a) = 0x44332211;
    else
    	*a = 0x44332211;
    	
    unsigned char *p = (unsigned char *)a;
    if ((*p==0x11) && (*(p+1)==0x22) && (*(p+2)==0x33) && (*(p+3)==0x44))
        e = 'L';
    else if ((*p==0x44) && (*(p+1)==0x33) && (*(p+2)==0x22) && (*(p+3)==0x11))
        e = 'B';
    else if ((*p==0x22) && (*(p+1)==0x11) && (*(p+2)==0x44) && (*(p+3)==0x33))
        e = 'M';
    else
        e = 'N';

    //printf("[%c] \n", e);
    return e;
}


// insert here your plugin description that will be displayed on the user interface
#ifdef IS_2D_PLUGIN

std::string iomanager::bioformats2D::desc()
{
	return	"******************************************************\n"
			"* bioformats2D v.1.0                                 *\n"
			"******************************************************\n"
			"*                                                    *\n"
			"* 2D image-based I/O plugin that uses the bioformats *\n"
			"* library to read several different formats          *\n"
			"* Implemented methods:                               *\n"
			"* - readMetadata                                     *\n"
			"* - readData                                         *\n"
			"*                                                    *\n"
			"* Tested image format extentions:                    *\n"
			"* - tif                                              *\n"
			"* - JPEG 2000                                        *\n"
			"*                                                    *\n"
			"* Accepted configuration parameters:                 *\n"
			"*  - none                                            *\n"
			"*                                                    *\n"
			"******************************************************\n";
}

#else

std::string iomanager::bioformats3D::desc()
{
	return	"******************************************************\n"
			"* bioformats3D v.1.1                                 *\n"
			"******************************************************\n"
			"*                                                    *\n"
			"* 3D image-based I/O plugin that uses the bioformats *\n"
			"* library to read several different formats          *\n"
			"* Implemented methods:                               *\n"
			"* - readMetadata                                     *\n"
			"* - readData                                         *\n"
			"*                                                    *\n"
			"* Tested image format extentions:                    *\n"
			"* - tif                                              *\n"
			"* - czi (Zeiss File format)                          *\n"
			"*                                                    *\n"
			"* Accepted configuration parameters:                 *\n"
			"*  - none                                            *\n"
			"*                                                    *\n"
			"******************************************************\n";
}

#endif



// Return if channels are interleaved (in case the image has just one channel return value is indefinite)
bool 
#ifdef IS_2D_PLUGIN
	iomanager::bioformats2D::isChansInterleaved( )
{
#else
	iomanager::bioformats3D::isChansInterleaved( )
{
#endif
	return false;
}


// Return if it is an input plugin only
bool 
#ifdef IS_2D_PLUGIN
	iomanager::bioformats2D::isInputOnly( )
{ 
#else
	iomanager::bioformats3D::isInputOnly( ) 
{
#endif
	return true;
}


// read image metadata
#ifdef IS_2D_PLUGIN
void 
	iomanager::bioformats2D::readMetadata(
	std::string img_path,		// (INPUT)	image filepath
	int & img_width,			// (OUTPUT) image width  (in pixels)
	int & img_height,			// (OUTPUT) image height (in pixels)
	int & img_bytes_x_chan,	    // (OUTPUT) number of bytes per channel
	int & img_chans,			// (OUTPUT) number of channels
	const std::string & params)	// (INPUT)	additional parameters <param1=val, param2=val, ...> 
throw (iom::exception)
{
	int img_depth = 1; // default value for a 2D plugin
#else
void 
	iomanager::bioformats3D::readMetadata(
	std::string img_path,		// (INPUT)	image filepath
	int & img_width,			// (OUTPUT) image width  (in pixels)
	int & img_height,			// (OUTPUT) image height (in pixels)
	int & img_depth,			// (OUTPUT) image depth (in pixels)
	int & img_bytes_x_chan,	    // (OUTPUT) number of bytes per channel
	int & img_chans,			// (OUTPUT) number of channels
	const std::string & params)	// (INPUT)	additional parameters <param1=val, param2=val, ...> 
throw (iom::exception)
{
#endif
	/**/iom::debug(iom::LEV3, iom::strprintf("img_path = %s, img_width = %d, img_height = %d, img_depth = %d, img_bytes_x_chan = %d, img_chans = %d, params = \"%s\"",
		img_path.c_str(), img_width, img_height, img_depth,	img_bytes_x_chan, img_chans, params.c_str()).c_str(), __iom__current__function__);

	if ( !jvm_setup ) {		
		if ( setupJVM(&jvm, &env) != JNI_OK ) {
			throw iom::exception(iom::strprintf("Cannot setup Java Virtual Machine"), __iom__current__function__);
		}

		jclass DebugToolsClass = env->FindClass("loci/common/DebugTools");
		jobject dtool = env->NewGlobalRef(env->NewObject(DebugToolsClass, env->GetMethodID(DebugToolsClass, "<init>", "()V")));
		//jboolean res = env->CallBooleanMethod(dtool, env->GetStaticMethodID(DebugToolsClass, "enableIJLogging", "(Z)Z"),false);
		jboolean res = env->CallBooleanMethod(dtool, env->GetStaticMethodID(DebugToolsClass, "enableLogging", "(Ljava/lang/String;)Z"), env->NewStringUTF("OFF"));
		//env->CallVoidMethod(dtool, env->GetMethodID(DebugToolsClass, "close", "()V"));
		env->DeleteGlobalRef(dtool);

		ImageReaderClass = env->FindClass("loci/formats/ImageReader");
		reader = env->NewGlobalRef(env->NewObject(ImageReaderClass, env->GetMethodID(ImageReaderClass, "<init>", "()V")));
	
		jvm_setup = true;
	}

	if ( img_path != cur_img_path ) {
		if ( cur_img_path != "" )
			env->CallVoidMethod(reader, env->GetMethodID(ImageReaderClass, "close", "()V"));
		env->CallVoidMethod(reader, 
			env->GetMethodID(ImageReaderClass, "setId", "(Ljava/lang/String;)V"), env->NewStringUTF(img_path.c_str()));
		cur_img_path = img_path;
	}

	jstring format = (jstring) env->CallObjectMethod(reader, env->GetMethodID(ImageReaderClass, "getFormat", "()Ljava/lang/String;"));
	const char *fdescr = env->GetStringUTFChars(format,NULL);

	int series_no  = -1;
	std::string str_series_no;
	if ( getParamValue(params,"series_no",str_series_no) )
		series_no = atoi(str_series_no.c_str());
	if ( series_no >= 0 ) {
		env->CallVoidMethod(reader, 
			env->GetMethodID(ImageReaderClass, "setSeries", "(I)V"), series_no);
	}

	img_width = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeX", "()I"));
	img_height = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeY", "()I"));
	if ( strcmp(fdescr,"Tagged Image File Format") == 0 ) 
		// multipage TIFFs are interpreted by bioformat as time series instead of a 3D stack
		img_depth = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeT", "()I"));
	else
		img_depth = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeZ", "()I"));
	img_bytes_x_chan = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getBitsPerPixel", "()I"));
	img_bytes_x_chan /= 8;
	img_chans = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeC", "()I"));

	//env->CallVoidMethod(reader, env->GetMethodID(ImageReaderClass, "close", "()V"));
	
	env->ReleaseStringUTFChars(format,fdescr);
	//env->DeleteGlobalRef(reader);
	
// 	printf("---> img_width = %d \n",img_width); 
// 	printf("---> img_height = %d \n",img_height); 
// 	printf("---> img_depth = %d \n",img_depth); 
// 	printf("---> img_bytes_x_chan = %d \n",img_bytes_x_chan); 
// 	printf("---> img_chans = %d \n",img_chans); fflush(stdout);
// 	fflush(stdout);
	//exit(0);

	if ( !buffer_setup ) {
		staticByteArray = (jbyteArray) env->NewGlobalRef(env->NewByteArray(img_width * img_height * img_bytes_x_chan));
	//	staticJavaPointer = env->GetByteArrayElements(staticByteArray, NULL);
		buffer_setup = true;
	}
	
	//if ( jvm->DestroyJavaVM() != JNI_OK ) {
	//	throw iom::exception(iom::strprintf("Cannot destroy Java Virtual Machine"), __iom__current__function__);
	//}
}


// Read 3D image data
#ifdef IS_2D_PLUGIN
unsigned char *						// (OUTPUT) buffer containing the read image data
	iomanager::bioformats2D::readData(
	std::string img_path,			// (INPUT) image filepath
	int & img_width,				// (INPUT/OUTPUT) image width  (in pixels)
	int & img_height,				// (INPUT/OUTPUT) image height (in pixels)
	int & img_bytes_x_chan,			// (INPUT/OUTPUT) number of bytes per channel
	int & img_chans,				// (INPUT/OUTPUT) number of channels to be read
	unsigned char *data,			// (INPUT) image data
	const std::string & params)		// (INPUT) additional parameters <param1=val, param2=val, ...> 
throw (iom::exception) 
{
	int img_depth = 1; // default value for a 2D plugin
	int z0        = 0; // 2D images have only plane 0
	int z1        = 1;
#else
unsigned char *						// (OUTPUT) buffer containing the read image data
	iomanager::bioformats3D::readData(
	std::string img_path,			// (INPUT) image filepath
	int & img_width,				// (INPUT/OUTPUT) image width  (in pixels)
	int & img_height,				// (INPUT/OUTPUT) image height (in pixels)
	int & img_depth,				// (INPUT/OUTPUT) image depth (in pixels)
	int & img_bytes_x_chan,			// (INPUT/OUTPUT) number of bytes per channel
	int & img_chans,				// (INPUT/OUTPUT) number of channels to be read
	unsigned char *data,			// (INPUT) image data
	int z0,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int z1,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	const std::string & params)		// (INPUT) additional parameters <param1=val, param2=val, ...> 
throw (iom::exception) 
{
#endif
	//throw iom::exception(iom::strprintf("not implemented yet"), __iom__current__function__);

	/**/iom::debug(iom::LEV3, iom::strprintf("img_path = %s, img_width = %d, img_height = %d, img_depth = %d, img_bytes_x_chan = %d, img_chans = %d, data = %p, z0 = %d, z1 = %d, params = \"%s\"",
		img_path.c_str(), img_width, img_height, img_depth,	img_bytes_x_chan, img_chans, data, z0, z1, params.c_str()).c_str(), __iom__current__function__);

	if ( !jvm_setup ) {		
		if ( setupJVM(&jvm, &env) != JNI_OK ) {
			throw iom::exception(iom::strprintf("Cannot setup Java Virtual Machine"), __iom__current__function__);
		}

		jclass DebugToolsClass = env->FindClass("loci/common/DebugTools");
		jobject dtool = env->NewGlobalRef(env->NewObject(DebugToolsClass, env->GetMethodID(DebugToolsClass, "<init>", "()V")));
		//jmethodID m = env->GetMethodID(DebugToolsClass, "<init>", "()V");
		//m = env->GetStaticMethodID(DebugToolsClass, "enableIJLogging", "(Z)Z");
		jboolean res = env->CallBooleanMethod(dtool, env->GetStaticMethodID(DebugToolsClass, "enableIJLogging", "(Z)Z"),false);
		//env->CallVoidMethod(dtool, env->GetMethodID(DebugToolsClass, "close", "()V"));
		env->DeleteGlobalRef(dtool);

		ImageReaderClass = env->FindClass("loci/formats/ImageReader");
		reader = env->NewGlobalRef(env->NewObject(ImageReaderClass, env->GetMethodID(ImageReaderClass, "<init>", "()V")));
	
		jvm_setup = true;
	}

	if ( img_path != cur_img_path ) {
		if ( cur_img_path != "" )
			env->CallVoidMethod(reader, env->GetMethodID(ImageReaderClass, "close", "()V"));
		env->CallVoidMethod(reader, 
			env->GetMethodID(ImageReaderClass, "setId", "(Ljava/lang/String;)V"), env->NewStringUTF(img_path.c_str()));
		cur_img_path = img_path;
	}

	int series_no  = -1;
	std::string str_series_no;
	if ( getParamValue(params,"series_no",str_series_no) )
		series_no = atoi(str_series_no.c_str());
	if ( series_no >= 0 ) {
		env->CallVoidMethod(reader, 
			env->GetMethodID(ImageReaderClass, "setSeries", "(I)V"), series_no);
	}

	jboolean isRGB = env->CallBooleanMethod(reader, env->GetMethodID(ImageReaderClass, "isRGB", "()Z"));
	
	unsigned int _width;
	unsigned int _height;
	unsigned int _depth;
	int _bytes_x_chan;
	unsigned int _chans;

	if ( !data ) { // recover the metadata, allocate the buffer and set parameters

		_width = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeX", "()I"));
		_height = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeY", "()I"));
		_depth = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeT", "()I"));
		_bytes_x_chan = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getBitsPerPixel", "()I"));
		_bytes_x_chan /= 8;
		_chans = env->CallIntMethod(reader, env->GetMethodID(ImageReaderClass, "getSizeC", "()I"));

		data = new unsigned char[((sint64)_width) * ((sint64)_height) * ((sint64)_depth) * _bytes_x_chan];
		memset(data,0,sizeof(unsigned char) * ((sint64)_width) * ((sint64)_height) * ((sint64)_depth) * _bytes_x_chan);
		img_width        = _width;
		img_height       = _height;
		img_depth        = _depth;
		img_bytes_x_chan = _bytes_x_chan;
		img_chans        = _chans;
	}

	// set the ROI
	z0 = (z0 < 0) ? 0: z0;
	z1 = (z1 < 0) ? img_depth : z1;

	if ( z0 >= z1 )
		throw iom::exception(iom::strprintf("wrong slice indices (z0 = %d, z1 = %d)",z0, z1), __iom__current__function__);

	// get the image
	//if ( (err_Tiff3Dfmt = readTiff3DFile2Buffer((char *)img_path.c_str(),data,img_width,img_height,z0,z1-1)) != 0 ) {
	//	throw iom::exception(iom::strprintf("(%s) unable to read tiff file %s in page range [%d,%d]",err_Tiff3Dfmt,img_path.c_str(),z0,z1-1), __iom__current__function__);
	//}
	unsigned char *pData = data;                                                               // offset of first channel of current plane
	sint64 sliceSize = ((sint64)img_width) * ((sint64)img_height) * img_bytes_x_chan;          // single channel slice size
	sint64 chanSize = ((sint64)img_width) * ((sint64)img_height) * img_bytes_x_chan * (z1-z0); // single channel buffer size

	// allocates a static buffer - WARNING: buffer size should be checked
	if ( !buffer_setup ) {
		staticByteArray = (jbyteArray) env->NewGlobalRef(env->NewByteArray((isRGB ? img_chans : 1) * sliceSize));
	//	staticJavaPointer = env->GetByteArrayElements(staticByteArray, NULL);
		buffer_setup = true;
	}

	int b_swap = 0; // no swapping by default
	
	// 2019-09-13. Giulio. @ADDED management of endianness when data and machine endianness differ
	jboolean isLittleEndian;
	char endianCodeMachine, endianCodeData;
	//bool isJP2 = (img_path.size() - img_path.find(".jp2")) == strlen(".jp2"); // endianness different is not managed for all formats
	// assume img_bytes_x_chan <= 2
	if ( img_bytes_x_chan == 2 ) { // check endianness
		isLittleEndian = env->CallBooleanMethod(reader, env->GetMethodID(ImageReaderClass, "isLittleEndian", "()Z"));
		endianCodeData = isLittleEndian ? 'L' : 'B';
		endianCodeMachine = checkMachineEndian();
		b_swap = (endianCodeMachine==endianCodeData)?0:1;
	}

	jbyteArray byteArray;
	jbyte *javaPointer;
	jbyte *javaPointerTmp;
	if ( isRGB ) {
		javaPointer = new jbyte[(isRGB ? img_chans : 1) * sliceSize];
		for ( int p=z0; p<z1; p++ ) {
			//byteArray = (jbyteArray) env->CallObjectMethod(reader, env->GetMethodID(ImageReaderClass, "openBytes", "(I)[B"),p);
			//javaPointer = env->GetByteArrayElements(byteArray, NULL);
			byteArray = (jbyteArray) env->CallObjectMethod(reader, env->GetMethodID(ImageReaderClass, "openBytes", "(I[B)[B"),p,staticByteArray);
			env->GetByteArrayRegion(byteArray, 0, (isRGB ? img_chans : 1) * sliceSize, javaPointer);
			unsigned char *ppData = pData; // offset of current channel of current plane in multiple planes buffer
			javaPointerTmp = javaPointer;  // offset of current channel of current plane in single plane buffer
			for ( int c=0; c<img_chans; c++ ) {
				memcpy(ppData, javaPointerTmp, sliceSize);
				ppData += chanSize;          // inter channel offset in multiple planes buffer = single channel buffer size
				javaPointerTmp += sliceSize; // inter channel offset in single plane buffer    = single channel slice size
			}
			//env->ReleaseByteArrayElements(byteArray,javaPointer,0);
			pData += sliceSize;              // inter slice offset in first plane buffer       = single channel slice size
		}
	}
	else {
		javaPointer = new jbyte[sliceSize];
		for ( int p=z0; p<z1; p++ ) {
			unsigned char *ppData = pData; 
			for ( int c=0; c<img_chans; c++ ) { // assumes XYCZT ordering
				//byteArray = (jbyteArray) env->CallObjectMethod(reader, env->GetMethodID(ImageReaderClass, "openBytes", "(I)[B"),(p*img_chans + c));
				//javaPointer = env->GetByteArrayElements(byteArray, NULL);
				//memcpy(ppData, javaPointer, sliceSize);
				//env->ReleaseByteArrayElements(byteArray,javaPointer,JNI_ABORT);
				//byteArray = (jbyteArray) env->CallObjectMethod(reader, env->GetMethodID(ImageReaderClass, "openBytes", "(I[B)[B"),(p*img_chans + c),staticByteArray);
				//memcpy(ppData, staticJavaPointer, sliceSize);
				byteArray = (jbyteArray) env->CallObjectMethod(reader, env->GetMethodID(ImageReaderClass, "openBytes", "(I[B)[B"),(p*img_chans + c),staticByteArray);
				env->GetByteArrayRegion(byteArray, 0, sliceSize, javaPointer);
				memcpy(ppData, javaPointer, sliceSize);
				ppData += chanSize;    // inter channel offset in multiple planes buffer = single channel buffer size
			}
			pData += sliceSize;        // inter slice offset in first plane buffer       = single channel slice size
		}
		delete javaPointer;
	}
	
	if ( b_swap == 1 ) {
		for (int i=0;i<(chanSize * img_chans / 2); i++) {
			iim::swap2bytes((void *)(data+(2*i)));
		}
	}

	//env->CallVoidMethod(reader, env->GetMethodID(ImageReaderClass, "close", "()V"));

	//env->DeleteGlobalRef(reader);

	//if ( jvm->DestroyJavaVM() != JNI_OK ) {
	//	throw iom::exception(iom::strprintf("Cannot destroy Java Virtual Machine"), __iom__current__function__);
	//}

	return data;
}


// Write 3D image data into a single (3D) image file
#ifdef IS_2D_PLUGIN
void 
	iomanager::bioformats2D::writeData(
	std::string img_path,			// (INPUT)	image filepath (it includes the file extension)
	unsigned char * raw_img,		// (INPUT)	image data to be saved into the file
	int img_height,					// (INPUT)	image height (in pixels)
	int img_width,					// (INPUT)	image width (in pixels)
	int img_bytes_x_chan,			// (INPUT)  number of bytes per channel
	int img_chans,					// (INPUT)	number of channels
	int y0,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int y1,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int x0,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int x1,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	const std::string & params)		// (INPUT) additional parameters <param1=val, param2=val, ...> 
throw (iom::exception) 
{
	int img_depth = 1; // default value for a 2D plugin
	int z0        = 0; // 2D images have only plane 0
	int z1        = 1;
#else
void 
	iomanager::bioformats3D::writeData(
	std::string img_path,			// (INPUT)	image filepath (it includes the file extension)
	unsigned char * raw_img,		// (INPUT)	image data to be saved into the file
	int img_height,					// (INPUT)	image height (in pixels)
	int img_width,					// (INPUT)	image width (in pixels)
	int img_depth,					// (INPUT)  image depth (in pixels)
	int img_bytes_x_chan,			// (INPUT)  number of bytes per channel
	int img_chans,					// (INPUT)	number of channels
	int z0,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int z1,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int y0,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int y1,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int x0,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	int x1,							// (INPUT)	region of interest [x0,x1)[y0,y1)[z0,z1) to be set on the image
	const std::string & params)		// (INPUT) additional parameters <param1=val, param2=val, ...> 
throw (iom::exception) 
{
#endif
	throw iom::exception(iom::strprintf("not implemented yet"), __iom__current__function__);

	/**/iom::debug(iom::LEV3, iom::strprintf("img_path = %s, img_width = %d, img_height = %d, img_depth = %d, img_bytes_x_chan = %d, img_chans = %d, z0 = %d, z1 = %d, y0 = %d, y1 = %d, x0 = %d, x1 = %d, params = \"%s\"",
		img_path.c_str(), img_width, img_height, img_depth,	img_bytes_x_chan, img_chans, z0, z1, y0, y1, x0, x1, params.c_str()).c_str(), __iom__current__function__);

}


// Create an empty 3D image 
#ifdef IS_2D_PLUGIN
// method 'create3Dimage' is not defined for 2D plugins
#else
void 
	iomanager::bioformats3D::create3Dimage(
	std::string img_path,			// (INPUT)	image filepath (it includes the file extension)
	int img_height,					// (INPUT)	image height (in pixels)
	int img_width,					// (INPUT)	image width (in pixels)
	int img_depth,					// (INPUT)  image depth (in pixels)
	int img_bytes_x_chan,			// (INPUT)  number of bytes per channel
	int img_chans,					// (INPUT)	number of channels
	const std::string & params)	// (INPUT) additional parameters <param1=val, param2=val, ...> 
throw (iom::exception) 
{
	throw iom::exception(iom::strprintf("not implemented yet"), __iom__current__function__);

	/**/iom::debug(iom::LEV3, iom::strprintf("img_path = %s, img_width = %d, img_height = %d, img_depth = %d, img_bytes_x_chan = %d, img_chans = %d, params = \"%s\"",
		img_path.c_str(), img_width, img_height, img_depth,	img_bytes_x_chan, img_chans, params.c_str()).c_str(), __iom__current__function__);

}
#endif


// Append a single slice at the bottom of a 3D image file
#ifdef IS_2D_PLUGIN
// method 'appendSlice' is not defined for 2D plugins
#else
void 
	iomanager::bioformats3D::appendSlice(
	std::string img_path,			// (INPUT)	image filepath (it includes the file extension)
	unsigned char * raw_img,		// (INPUT)	slice to be saved into the file
	int img_height,					// (INPUT)	slice height (in pixels)
	int img_width,					// (INPUT)	slice width (in pixels)
	int img_bytes_x_chan,			// (INPUT)  number of bytes per channel
	int img_chans,					// (INPUT)	number of channels
	int y0,							// (INPUT)	region of interest [x0,x1)[y0,y1) to be set on the image
	int y1,							// (INPUT)	region of interest [x0,x1)[y0,y1) to be set on the image
	int x0,							// (INPUT)	region of interest [x0,x1)[y0,y1) to be set on the image
	int x1,							// (INPUT)	region of interest [x0,x1)[y0,y1) to be set on the image
	int slice,						// (INPUT)  slice index 
	const std::string & params)		// (INPUT) additional parameters <param1=val, param2=val, ...> 
throw (iom::exception)
{
	throw iom::exception(iom::strprintf("not implemented yet"), __iom__current__function__);

	/**/iom::debug(iom::LEV3, iom::strprintf("img_path = %s, img_width = %d, img_height = %d, img_bytes_x_chan = %d, img_chans = %d, y0 = %d, y1 = %d, x0 = %d, x1 = %d, params = \"%s\"",
		img_path.c_str(), img_width, img_height, img_bytes_x_chan, img_chans, y0, y1, x0, x1, params.c_str()).c_str(), __iom__current__function__);

}
#endif


// read image data (deprecated)
#ifdef IS_2D_PLUGIN
iom::real_t*						// (OUTPUT) a [0.0,1.0]-valued array storing the 3D image in channel->slice->row order
	iomanager::bioformats2D::readData(
	char **files,					// (INPUT)	array of C-strings storing image filenames
	int files_size,					// (INPUT)	size of 'files'
	const char *path /*= 0*/,		// (INPUT)	path to be concatenated to the i-th entry of 'files'
	int first /*= -1*/,				// (INPUT)	selects a range [first, last] of files to be loaded
	int last /*= -1*/,				// (INPUT)	selects a range [first, last] of files to be loaded	
	bool is_sparse /*= false*/,		// (INPUT)	if true, 'files' is a sparse array and null entries should be treated as empty (black) images
	iom::channel chan,				// (INPUT)	channel selection { ALL, R, G, B }. 
	const std::string & params)		// (INPUT)	additional parameters <param1=val, param2=val, ...> 
throw (iom::exception)
{
	throw iom::exception(iom::strprintf("not implemented yet"), __iom__current__function__);
	
	/**/iom::debug(iom::LEV3, iom::strprintf("files_size = %d, path = %s, first = %d, last = %d, is_sparse = %s, chan = %d, params = \"%s\"",
		files_size, path ? path : "null", first, last, is_sparse ? "true" : "false", chan, params.c_str()).c_str(), __iom__current__function__);
}
#else
void
	iomanager::bioformats3D::readData(
	char *finName,					// image filepath
	int XSIZE,						// image width (in pixels)
	int YSIZE,						// image height (in pixels)
	unsigned char *data,			// (OUTPUT) image data
	int first,						// selects a range [first, last] of files to be loaded 
	int last,						// selects a range [first, last] of files to be loaded 
	const std::string & params)		// (INPUT)	additional parameters <param1=val, param2=val, ...> 
throw (iom::exception)
{
	throw iom::exception(iom::strprintf("no more available"), __iom__current__function__);

	/**/iom::debug(iom::LEV3, iom::strprintf("finName = %s, XSIZE = %d, YSIZE = %d, first = %d, last = %d, params = \"%s\"",
		finName, XSIZE, YSIZE, first, last, params.c_str()).c_str(), __iom__current__function__);

}
#endif

// write image data (deprecated)
#ifdef IS_2D_PLUGIN
void 
	iomanager::bioformats2D::writeData(
	std::string img_path,		// (INPUT)	image filepath (it includes the file extension)
	iom::real_t* raw_img,			// (INPUT)	a [0.0,1.0]-valued array storing the 2D image in channel->row order
	int img_height,				// (INPUT)	image height
	int img_width,				// (INPUT)	image width
	int img_chans,				// (INPUT)	number of channels
	int y0,						// (INPUT)	region of interest [x0,x1][y0,y1] to be set on the image
	int y1,						// (INPUT)	region of interest [x0,x1][y0,y1] to be set on the image
	int x0,						// (INPUT)	region of interest [x0,x1][y0,y1] to be set on the image
	int x1,						// (INPUT)	region of interest [x0,x1][y0,y1] to be set on the image
	int bpp,					// (INPUT)	color depth (bits per pixel)
	const std::string & params)	// (INPUT)	additional parameters <param1=val, param2=val, ...> 
throw (iom::exception)
{
	throw iom::exception(iom::strprintf("not implemented yet"), __iom__current__function__);
	
	/**/iom::debug(iom::LEV3, iom::strprintf("img_path = %s, img_height = %d, img_width = %d, y0 = %d, y1 = %d, x0 = %d, x1 = %d, bpp = %d, params = \"%s\"", 
		img_path.c_str(), img_height, img_width, y0, y1, x0, x1, bpp, params.c_str()).c_str(), __iom__current__function__);
}
#endif