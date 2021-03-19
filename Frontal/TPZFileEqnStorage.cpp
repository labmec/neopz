/**
 * @file
 * @brief Contains the implementation of the TPZFileEqnStorage methods.
 */

#include "TPZFileEqnStorage.h"
#include <stdlib.h>

#ifdef WIN32
//#include <dir.h>
#include <io.h>
#endif
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.frontal.tpzfileeqnstorage");
#endif
#include <fstream>
using namespace std;

template<class TVar>
void TPZFileEqnStorage<TVar>::WriteHeaders() {
	/**
	 *Updates fNumBlocks information each time 
	 *WriteHeaders is called
	 */
	fNumBlocks++;
    TPZVec<int64_t> Position(fNumHeaders,0);
	/**
	 *If fCurrentBlock = 0 then a fBlockPos.Push must be called to 
	 *store the first address
	 */
	if(fCurrentBlock==0) {
		int64_t basepos = ftell(fIOStream);
		fBlockPos.Push(basepos);
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout << "basepos = " << basepos;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif	
	}else{
		int64_t tempaddress = ftell(fIOStream);
		fBlockPos.Push(tempaddress);
	}
	
	int64_t firstpos = ftell(fIOStream);
	//if (fCurrentBlock) firstpos = 
	
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Writing the position of the headers, numheaders " << fNumHeaders << " position ";
		for(int i=0; i<fNumHeaders; i++) sout << Position[i] << ' ';
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	/** Writes fNumHeaders positions for the headers */
	fwrite(Position.begin(),sizeof(int64_t),fNumHeaders,fIOStream);

	/** Get starting position of first header */
	int64_t firstaddress = ftell(fIOStream);
	
	/** Writes first position the address of block one */
	fseek(fIOStream,firstpos,SEEK_SET);
	fwrite(&firstaddress,sizeof(int64_t),1,fIOStream);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "At position " << firstpos << " writing the first address " << firstaddress;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	/** Sets fCurBlockPosition to actual address */
	fCurBlockPosition = firstaddress;
	
	/** Return the pointer to the actual position */
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Setting the file position at " << firstaddress;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fseek(fIOStream,firstaddress,SEEK_SET);
	
}

template<class TVar>
void TPZFileEqnStorage<TVar>::Store(int ieq, int jeq, const char *name){
	//Initial tests with C input output files !
	int loop_limit=100;// = jeq-ieq;
	int i;
	FILE *out_file = fopen(name,"wb");
	double number=2.1;
	double val = 0;
	int64_t fPos[5] = {0};
	int64_t firstpos = ftell(out_file);
	fwrite(fPos,sizeof(int64_t),5,out_file);
	double readvec[4][100];
	int iblock =0;
	int64_t sizereturn;
	sizereturn = 0;
	for(iblock=0; iblock<4; iblock++) {
		int64_t currentpos = ftell(out_file);
		fseek(out_file,firstpos+iblock*sizeof(int64_t),SEEK_SET);
		fwrite(&currentpos,sizeof(int64_t),1,out_file);
		fseek(out_file,currentpos,SEEK_SET);
		for(i=0;i<loop_limit;i++){
			val=number*i*(iblock+1);
			fwrite(&val, sizeof(double), 1, out_file);
		}
	}
	fclose(out_file);
	out_file = fopen(name,"rb");
	sizereturn = fread(fPos,sizeof(int64_t),5,out_file);
#ifdef PZDEBUG
	if (sizereturn != 5) DebugStop();
#endif
	for(iblock = 0; iblock<4; iblock++) {
		fseek(out_file,fPos[iblock],SEEK_SET);
		sizereturn = fread(readvec[iblock],sizeof(double),loop_limit,out_file);
#ifdef PZDEBUG
		if (sizereturn != loop_limit) DebugStop();
#endif
	}
}

template<class TVar>
void TPZFileEqnStorage<TVar>::Forward(TPZFMatrix<TVar> &f, DecomposeType dec) const
{
	TPZEqnArray<TVar> REqnArray;
	for(int64_t i=0;i<fBlockPos.NElements();i++) {
		if (fBlockPos[i]) {
			if(fseek(fIOStream,fBlockPos[i],SEEK_SET)){
				cout << "fseek fail on Element " << i << " Position " << fBlockPos[i] << endl;
				cout.flush();
			}
			int64_t position;
			if(!fread(&position,sizeof(int64_t),1,fIOStream)){
				cout << "fread fail on Element " << i << " Position " << position << endl;
				cout << "EOF " << feof(fIOStream) << endl;
				cout << "Error Number " << ferror(fIOStream) << endl;
				cout.flush();
			}
			if(fseek(fIOStream,position,SEEK_SET)){
				cout << "fseek fail on Element " << i << " Position " << position << endl;
				cout.flush();
			}
			
			REqnArray.Read(fIOStream);
			REqnArray.EqnForward(f,dec);
		}
	}
}

template<class TVar>
void TPZFileEqnStorage<TVar>::Backward(TPZFMatrix<TVar> &f, DecomposeType dec) const
{
	TPZEqnArray<TVar> REqnArray;
	int i;
	int64_t sizereturn;
	sizereturn = 0;
	for(i=fBlockPos.NElements()-1;i>=0;i--){
		if (fBlockPos[i]) {
			fseek(fIOStream,fBlockPos[i],SEEK_SET);
			int64_t position;
			sizereturn = fread(&position,sizeof(int64_t),1,fIOStream);
#ifdef PZDEBUG
			if (sizereturn != 1) DebugStop();
#endif
			fseek(fIOStream,position,SEEK_SET);
			REqnArray.Read(fIOStream);
			REqnArray.EqnBackward(f,dec); 
		}
	}
}

template<class TVar>
void TPZFileEqnStorage<TVar>::Reset()
{
	
}

template<class TVar>
void TPZFileEqnStorage<TVar>::Print(const char *name, std::ostream& out) const {
    int i;
	TPZEqnArray<TVar> REqnArray;
	out <<  "Number of entries on File  "<< fBlockPos.NElements() << endl;
	for(i=0;i<fBlockPos.NElements();i++) {
		if (fBlockPos[i]) {
			fseek(fIOStream,fBlockPos[i],SEEK_SET);
			REqnArray.Read(fIOStream);
			REqnArray.Print(name, out);
		}
	}
}

// Redefine TPZFileEqnStorage so it can handle parallel writeing!!!
template<class TVar>
void TPZFileEqnStorage<TVar>::AddEqnArray(TPZEqnArray<TVar> *EqnArray)
{
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "fCurrentBlock "<< fCurrentBlock << " fNumHeaders " << fNumHeaders;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  	
	if(fCurrentBlock%(fNumHeaders-1)==0) {
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            LOGPZ_DEBUG(logger,"writing headers")
        }
#endif
		WriteHeaders();
		//		fBlockPos.Push(fBlockPos[fCurrentBlock-1]+sizeof(int64_t));
		
	}else if(fCurrentBlock!=0){
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
            std::stringstream sout;
            sout << "fCurrentBlock " << fCurrentBlock << " fBlockPos[fCurrentBlock-1] "<< fBlockPos[fCurrentBlock-1] << fBlockPos;
            LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		fBlockPos.Push(fBlockPos[fCurrentBlock-1]+sizeof(int64_t));
	}
	
	EqnArray->Write(fIOStream); 
	
	/**Gets actual position on fIOStream */
	int64_t nextaddress=ftell(fIOStream);

	/** 
	 *Writes this address on next available header block
	 *and sets pointer to its previous position
	 */
	//fseek(fIOStream,fCurBlockPosition+sizeof(int64_t),SEEK_SET);
	fseek(fIOStream,fBlockPos[fCurrentBlock]+sizeof(int64_t),SEEK_SET);
  	fCurBlockPosition=ftell(fIOStream);
  	fwrite(&nextaddress,sizeof(int64_t),1,fIOStream);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "nextaddress " << nextaddress << " fCurBlockPosition " << fCurBlockPosition;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	fseek(fIOStream,nextaddress,SEEK_SET);
	
	fCurrentBlock++;
}

template<class TVar>
TPZFileEqnStorage<TVar>::TPZFileEqnStorage(char option, const std::string & name)
: TPZRegisterClassId(&TPZFileEqnStorage<TVar>::ClassId)
{
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=20;
	fFileName = name;
	if(option=='r'){
		fIOStream = fopen(fFileName.c_str(),"rb"); //open for reading
		/**
		 *Opens binary files and get initial information
		 *use this information for storage requirements
		 */
		int64_t sizereturn;
		sizereturn = 0;
		sizereturn = fread(&fNumHeaders,sizeof(int),1,fIOStream);
#ifdef PZDEBUG
		if (sizereturn != 1) DebugStop();
#endif
		sizereturn = fread(&fNumBlocks,sizeof(int),1,fIOStream);
#ifdef PZDEBUG
		if (sizereturn != 1) DebugStop();
#endif
		ReadBlockPositions();
	}else if(option=='w'){
		fIOStream = fopen(fFileName.c_str(),"wb"); //open for writing
		/**
		 *Writes NumHeaders and NumBlocks information in
		 *the two initial positions on fIOStream
		 */
		int zero = 0;
		fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
		fwrite(&zero,sizeof(int),1,fIOStream);
		//fCurBlockPosition = ftell(fIOStream);
	}
}

template<class TVar>
TPZFileEqnStorage<TVar>::~TPZFileEqnStorage()
{
	remove(fFileName.c_str());
}

template<class TVar>
void TPZFileEqnStorage<TVar>::main()
{
	int Loop_Limit=24;
	//cout << "Loop_Limit <";
	//cin >> Loop_Limit;
	std::string filename;
	filename = "testbinary.txt\0";
	TPZFileEqnStorage FileStoreW('w',filename);
	//	FileStoreW.SetBlockSize(10); 
	TPZEqnArray<TVar> EqnArray;
	
	ifstream input("MatrizInversa.txt");
	TVar aux;
	int i, j;
	TPZFMatrix<TVar> DecMat(Loop_Limit,Loop_Limit);
	for(i=0;i<Loop_Limit;i++){
		for(j=0;j<Loop_Limit;j++){
			input >> aux;
			DecMat(i,j) = aux;
		}
	}

	for(i=0;i<Loop_Limit;i++){
		EqnArray.BeginEquation(i);
		for(j=i;j<Loop_Limit;j++){
			//			aux = DecMat(i,j);
			EqnArray.AddTerm(j,DecMat(i,j));
		}
		EqnArray.EndEquation();
	}
	
	FileStoreW.AddEqnArray(&EqnArray); 
	
	FileStoreW.FinishWriting();
	
	FileStoreW.ReOpen();
	
	ofstream output("testeFileBin.txt");
	
	FileStoreW.Print("Teste",output);
	
	TPZFMatrix<TVar> f(Loop_Limit,1);
	
	for(i=0;i<Loop_Limit;i++) {
		f(i,0) = (TVar)(float)((i+1.)*2.1/23.);
	}
	
	f.Print("Teste"); 
	
	FileStoreW.Forward(f,ECholesky);
	
}

static char filenamestorage[256];

template<class TVar>
TPZFileEqnStorage<TVar>::TPZFileEqnStorage()
: TPZRegisterClassId(&TPZFileEqnStorage<TVar>::ClassId)
{
	strcpy(filenamestorage, "/tmp/binary_frontalXXXXXX");
	int fdtmp = -1;
#ifdef WIN32
	_mktemp(filenamestorage);
#else
	fdtmp = mkstemp(filenamestorage); //returns file description for tmp file
#endif
	
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=20;
	
	fFileName = filenamestorage;
#ifdef WIN32
	fIOStream = fopen(filenamestorage,"wb"); //open for writing
#else
	fIOStream = fdopen(fdtmp,"wb"); //open for writing
#endif
	/**
	 *Writes NumHeaders and NumBlocks information in
	 *the two initial positions on fIOStream
	 */
	int zero = 0;
	fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	fwrite(&zero,sizeof(int),1,fIOStream);
}

template<class TVar>
TPZFileEqnStorage<TVar>::TPZFileEqnStorage(const TPZFileEqnStorage &)
: TPZRegisterClassId(&TPZFileEqnStorage<TVar>::ClassId)
{
	strcpy(filenamestorage, "/tmp/binary_frontalXXXXXX");
#ifdef WIN32
	_mktemp(filenamestorage);
#else
	int fdtmp = -1;
	fdtmp = mkstemp(filenamestorage); //returns file description for tmp file
#endif
	
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=20;
	
	fFileName = filenamestorage;
#ifdef WIN32
	fIOStream = fopen(filenamestorage,"wb"); //open for writing
#else
	fIOStream = fdopen(fdtmp,"wb"); //open for writing
#endif
	/**
	 *Writes NumHeaders and NumBlocks information in
	 *the two initial positions on fIOStream
	 */
	int zero = 0;
	fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	fwrite(&zero,sizeof(int),1,fIOStream);
}

template<class TVar>
void TPZFileEqnStorage<TVar>::Zero()
{
	if(fIOStream) FinishWriting();
	remove(fFileName.c_str());
	
	strcpy(filenamestorage, "/tmp/binary_frontalXXXXXX");
	
	int fdtmp = -1;
#ifdef WIN32
	_mktemp(filenamestorage);
#else
	fdtmp = mkstemp(filenamestorage); //returns file description for tmp file
#endif
	cout << "Temporary file name " << filenamestorage << endl;
	cout.flush();
	
	fBlockPos.Resize(0);
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=20;
	
	fFileName = filenamestorage;
#ifdef WIN32
	fIOStream = fopen(filenamestorage,"wb"); //open for writing
#else
	fIOStream = fdopen(fdtmp,"wb"); //open for writing
#endif
	/**
	 *Writes NumHeaders and NumBlocks information in
	 *the two initial positions on fIOStream
	 */
	int zero = 0;
	fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
	fwrite(&zero,sizeof(int),1,fIOStream);
}

template<class TVar>
void TPZFileEqnStorage<TVar>::ReOpen()
{
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) LOGPZ_DEBUG(logger,"reopening the file")
#endif
	fIOStream = fopen(fFileName.c_str(),"rb"); //open for reading
	/**
	 *Opens binary files and get initial information
	 *use this information for storage requirements
	 */
	int64_t sizereturn;
	sizereturn = 0;
	sizereturn = fread(&fNumHeaders,sizeof(int),1,fIOStream);
#ifdef PZDEBUG
	if (sizereturn != 1) DebugStop();
#endif
	sizereturn = fread(&fNumBlocks,sizeof(int),1,fIOStream);
#ifdef PZDEBUG
	if (sizereturn != 1) DebugStop();
#endif
}

template<class TVar>
void TPZFileEqnStorage<TVar>::OpenGeneric(char option, const char * name)
{
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=11;

	fFileName = name;
	
	if(option=='r'){
		fIOStream = fopen(fFileName.c_str(),"rb"); //open for reading
		/**
		 *Opens binary files and get initial information
		 *use this information for storage requirements
		 */
		int64_t sizereturn;
		sizereturn = 0;
		sizereturn = fread(&fNumHeaders,sizeof(int),1,fIOStream);
#ifdef PZDEBUG
		if (sizereturn != 1) DebugStop();
#endif
		sizereturn = fread(&fNumBlocks,sizeof(int),1,fIOStream);
#ifdef PZDEBUG
		if (sizereturn != 1) DebugStop();
#endif
		ReadBlockPositions();
	}else if(option=='w'){
		fIOStream = fopen(fFileName.c_str(),"wb"); //open for writing
		/**
		 *Writes NumHeaders and NumBlocks information in
		 *the two initial positions on fIOStream
		 */
		int zero = 0;
		fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
		fwrite(&zero,sizeof(int),1,fIOStream);
	}
}

template<class TVar>
void TPZFileEqnStorage<TVar>::FinishWriting()
{
#ifdef PZ_LOG
	if(logger.isDebugEnabled()) LOGPZ_DEBUG(logger,"Closing the binary file")
#endif
	fseek(fIOStream,sizeof(int),SEEK_SET);
	//cout << "Second fseek " << ftell(fIOStream) << endl;
	//fwrite(&fNumBlocks,sizeof(int),1,fIOStream);
	fwrite(&fNumBlocks,sizeof(int),1,fIOStream);
	fclose(fIOStream);
	fIOStream = 0;
}

template<class TVar>
void TPZFileEqnStorage<TVar>::ReadBlockPositions()
{
#ifdef PZ_LOG
	if(logger.isDebugEnabled() ) LOGPZ_DEBUG(logger,"Reading block positions")
#endif
	int aux = fNumBlocks * (fNumHeaders-1);
	fBlockPos.Resize(aux);
	int i, ibl = 0;
	cout << "Reading Block Positions\n";
	cout.flush();
	int64_t sizereturn;
	sizereturn = 0;
	for(i=0;i<fNumBlocks;i++) {
		cout << "*";
		cout.flush();
		if(!(i%20)){
			cout << 100*i/fNumBlocks << "% Read\n";
			cout.flush();
		}
		sizereturn = fread(&fBlockPos[ibl],sizeof(int64_t),fNumHeaders-1,fIOStream);
#ifdef PZDEBUG
		if (sizereturn != fNumHeaders-1) DebugStop();
#endif
		ibl+=fNumHeaders-1;
		int64_t nextpos;
		sizereturn = fread(&nextpos,sizeof(int64_t),fNumHeaders-1,fIOStream);
#ifdef PZDEBUG
		if (sizereturn != fNumHeaders-1) DebugStop();
#endif
		fseek(fIOStream,nextpos,SEEK_SET);
	}
}

template<class TVar>
std::string TPZFileEqnStorage<TVar>::GetStorage() {
    return "File Storage";
}

template class TPZFileEqnStorage<float>;
template class TPZFileEqnStorage<std::complex<float> >;

template class TPZFileEqnStorage<double>;
template class TPZFileEqnStorage<std::complex<double> >;

template class TPZFileEqnStorage<long double>;
template class TPZFileEqnStorage<std::complex<long double> >;
