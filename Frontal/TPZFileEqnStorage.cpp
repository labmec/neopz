/**
 * @file
 * @brief Contains the implementation of the TPZFileEqnStorage methods.
 */
//$Id: TPZFileEqnStorage.cpp,v 1.15 2011-05-11 02:13:57 phil Exp $

#include "TPZFileEqnStorage.h"
#include <stdlib.h>

#ifdef WIN32
#include <dir.h>
#endif
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.frontal.tpzfileeqnstorage"));
#endif
#include <fstream>
using namespace std;

void TPZFileEqnStorage::WriteHeaders(){
	/**
	 *Updates fNumBlocks information each time 
	 *WriteHeaders is called
	 */
	fNumBlocks++;
	int i;
    TPZVec<long int> Position(fNumHeaders,0);
	/**
	 *If fCurrentBlock = 0 then a fBlockPos.Push must be called to 
	 *store the first address
	 */
	if(fCurrentBlock==0){
		long int basepos = ftell(fIOStream);
		fBlockPos.Push(basepos);
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "basepos = " << basepos;
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif	
	}else{
		long int tempaddress = ftell(fIOStream);
		fBlockPos.Push(tempaddress);
	}
	
	long int firstpos = ftell(fIOStream);
	//if (fCurrentBlock) firstpos = 
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Writing the position of the headers, numheaders " << fNumHeaders << " position ";
		for(i=0; i<fNumHeaders; i++) sout << Position[i] << ' ';
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	/**
	 *Writes fNumHeaders positions for the headers
	 */
	fwrite(Position,sizeof(long  int),fNumHeaders,fIOStream);

	/**
	 *Get starting position of first header
	 */
	long int firstaddress = ftell(fIOStream);
	
	/**
	 *Writes first position the address of block one
	 */
	fseek(fIOStream,firstpos,SEEK_SET);
	fwrite(&firstaddress,sizeof(long int),1,fIOStream);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "At position " << firstpos << " writing the first address " << firstaddress;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	/**
	 *Sets fCurBlockPosition to actual address
	 */
	fCurBlockPosition = firstaddress;
	
	/**
	 *Return the pointer to the actual position
	 */
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Setting the file position at " << firstaddress;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	fseek(fIOStream,firstaddress,SEEK_SET);
	
}
/*void TPZFileEqnStorage::SetBlockSize(int bs){
 fBlockSize = bs;
 }*/
void TPZFileEqnStorage::Store(int ieq, int jeq, const char *name){
	//Initial tests with C input output files !
	int loop_limit=100;// = jeq-ieq;
	int i;
	FILE *out_file = fopen(name,"wb");
	//cout << "Loop Limit "; 
	//cin >> loop_limit ;
	//cout << "Block Size";
	//cin >> block_size;
	
	
	
	//	struct _iobuf *temp_i;
	//From MSDN
	/*	
	 
	 char list[30];
	 int  i, numread, numwritten;
	 
	 Open file in text mode: 
	 if( (stream = fopen( "fread.out", "w+t" )) != NULL )
	 {
	 for ( i = 0; i < 25; i++ )
	 list[i] = (char)('z' - i);
	 Write 25 characters to stream 
	 numwritten = fwrite( list, sizeof( char ), 25, stream );
	 fprintf( "Wrote %d items\n", numwritten );
	 fclose( stream );
	 */
	double number=2.1;
	double val = 0;
	long int fPos[5] = {0};
	long int firstpos = ftell(out_file);
	fwrite(fPos,sizeof(long  int),5,out_file);
	double readvec[4][100];
	int iblock =0;
	for(iblock=0; iblock<4; iblock++) {
		long int currentpos = ftell(out_file);
		fseek(out_file,firstpos+iblock*sizeof(long int),SEEK_SET);
		fwrite(&currentpos,sizeof(long int),1,out_file);
		fseek(out_file,currentpos,SEEK_SET);
		for(i=0;i<loop_limit;i++){
			val=number*i*(iblock+1);
			fwrite(&val, sizeof(double), 1, out_file);
		}
	}
	fclose(out_file);
	out_file = fopen(name,"rb");
	fread(fPos,sizeof(long int),5,out_file);
	for(iblock = 0; iblock<4; iblock++) {
		fseek(out_file,fPos[iblock],SEEK_SET);
		fread(readvec[iblock],sizeof(double),loop_limit,out_file);
	}
}

void TPZFileEqnStorage::Forward(TPZFMatrix &f, DecomposeType dec) const
{
	//  cout << "Inside TPZFileEqnStorage::Forward" << endl;
	//  cout << "fBlockPos.NElements() = " << fBlockPos.NElements() << endl;
	
	//if(!fIOStream) SetFileName(fFileName);
	TPZEqnArray REqnArray;
	int i;
	for(i=0;i<fBlockPos.NElements();i++) {
		if (fBlockPos[i]) {
			//     if(!(i%10)) cout << "*";
			//     if(!(i%100)) cout << i << endl;
			
			if(fseek(fIOStream,fBlockPos[i],SEEK_SET)){
				cout << "fseek fail on Element " << i << " Position " << fBlockPos[i] << endl;
				cout.flush();
			}
			long int position;
			if(!fread(&position,sizeof(long int),1,fIOStream)){
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
void TPZFileEqnStorage::Backward(TPZFMatrix &f, DecomposeType dec) const
{
	//  cout << "Inside TPZFileEqnStorage::Backward" << endl;
	//  cout << "fBlockPos.NElements() = " << fBlockPos.NElements() << endl;
	
	TPZEqnArray REqnArray;
	int i;
	for(i=fBlockPos.NElements()-1;i>=0;i--){
		if (fBlockPos[i]) {
			//		     if(!(i%10)) cout << "*";
			//		     if(!(i%100)) cout << i << endl;
			fseek(fIOStream,fBlockPos[i],SEEK_SET);
			long int position;
			fread(&position,sizeof(long int),1,fIOStream);
			fseek(fIOStream,position,SEEK_SET);
			REqnArray.Read(fIOStream);
			REqnArray.EqnBackward(f,dec); 
		}
	}
	
}

void TPZFileEqnStorage::Reset()
{
	
}

void TPZFileEqnStorage::Print(const char *name, std::ostream& out) const {
    int i;
	TPZEqnArray REqnArray;
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

void TPZFileEqnStorage::AddEqnArray(TPZEqnArray *EqnArray)
{
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "fCurrentBlock "<< fCurrentBlock << " fNumHeaders " << fNumHeaders;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  	
	if(fCurrentBlock%(fNumHeaders-1)==0) {
#ifdef LOG4CXX
		LOGPZ_DEBUG(logger,"writing headers")
#endif
		WriteHeaders();
		//		fBlockPos.Push(fBlockPos[fCurrentBlock-1]+sizeof(long int));
		
	}else if(fCurrentBlock!=0){
#ifdef LOG4CXX
		{
            std::stringstream sout;
            sout << "fCurrentBlock " << fCurrentBlock << " fBlockPos[fCurrentBlock-1] "<< fBlockPos[fCurrentBlock-1] << fBlockPos;
            LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		fBlockPos.Push(fBlockPos[fCurrentBlock-1]+sizeof(long int));
	}
	
	EqnArray->Write(fIOStream); 
	
	/**Gets actual position on fIOStream */
	long int nextaddress=ftell(fIOStream);
	
	
	/** 
	 *Writes this address on next available header block
	 *and sets pointer to its previous position
	 */
	//fseek(fIOStream,fCurBlockPosition+sizeof(long int),SEEK_SET);
	fseek(fIOStream,fBlockPos[fCurrentBlock]+sizeof(long int),SEEK_SET);
  	fCurBlockPosition=ftell(fIOStream);
  	fwrite(&nextaddress,sizeof(long int),1,fIOStream);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "nextaddress " << nextaddress << " fCurBlockPosition " << fCurBlockPosition;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	fseek(fIOStream,nextaddress,SEEK_SET);
	//	fBlockPos.Push(fCurBlockPosition);
	
	fCurrentBlock++;
	//	fSubBlockIndex++;
	
	//	cout << "NumBlocks " << fNumBlocks << endl;
	//	cout << "NumHeaders " << fNumHeaders << endl;
}

TPZFileEqnStorage::TPZFileEqnStorage(char option, const std::string & name)
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
		fread(&fNumHeaders,sizeof(int),1,fIOStream);
		fread(&fNumBlocks,sizeof(int),1,fIOStream);
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

TPZFileEqnStorage::~TPZFileEqnStorage()
{
	remove(fFileName.c_str());
}

void TPZFileEqnStorage::main()
{
	int Loop_Limit=24;
	//cout << "Loop_Limit <";
	//cin >> Loop_Limit;
	std::string filename;
	filename = "testbinary.txt\0";
	TPZFileEqnStorage FileStoreW('w',filename);
	//	FileStoreW.SetBlockSize(10); 
	TPZEqnArray EqnArray;
	
	ifstream input("MatrizInversa.txt");
	double aux;
	int i, j;
	TPZFMatrix DecMat(Loop_Limit,Loop_Limit);
	for(i=0;i<Loop_Limit;i++){
		for(j=0;j<Loop_Limit;j++){
			input >> aux;
			DecMat(i,j)=aux;
		}
	}
	
	//DecMat.Print("MatrizInv");
	
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
	
	TPZFMatrix f(Loop_Limit,1);
	
	for(i=0;i<Loop_Limit;i++) {
		f(i,0) = (i+1)*2.1/23;
	}
	
	f.Print("Teste"); 
	
	FileStoreW.Forward(f,ECholesky);
	
}

static char filenamestorage[256];
TPZFileEqnStorage::TPZFileEqnStorage()
{
	strcpy(filenamestorage, "/tmp/binary_frontalXXXXXX");
#ifdef WIN32
	_mktemp(filenamestorage);
#else
	mkstemp(filenamestorage);
#endif
	//     cout << "Temporary file name " << filenamestorage << endl;
	//     cout.flush();
	
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=20;
	//	SetBlockSize(10); 
	//fBlockPos.Resize(fNumHeaders);
	
	
	fFileName = filenamestorage;
	
	fIOStream = fopen(fFileName.c_str(),"wb"); //open for writing
	/**
	 *Writes NumHeaders and NumBlocks information in
	 *the two initial positions on fIOStream
	 */
	int zero = 0;
	fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	fwrite(&zero,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	//fCurBlockPosition = ftell(fIOStream);
	
	
}

TPZFileEqnStorage::TPZFileEqnStorage(const TPZFileEqnStorage &)
{
	strcpy(filenamestorage, "/tmp/binary_frontalXXXXXX");
#ifdef WIN32
	_mktemp(filenamestorage);
#else
	mkstemp(filenamestorage);
#endif
	//     cout << "Temporary file name " << filenamestorage << endl;
	//     cout.flush();
	
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=20;
	//	SetBlockSize(10); 
	//fBlockPos.Resize(fNumHeaders);
	
	
	fFileName = filenamestorage;
	
	fIOStream = fopen(fFileName.c_str(),"wb"); //open for writing
	/**
	 *Writes NumHeaders and NumBlocks information in
	 *the two initial positions on fIOStream
	 */
	int zero = 0;
	fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	fwrite(&zero,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	//fCurBlockPosition = ftell(fIOStream);
	
	
}

void TPZFileEqnStorage::Zero()
{
	if(fIOStream) FinishWriting();
	remove(fFileName.c_str());
	
	strcpy(filenamestorage, "/tmp/binary_frontalXXXXXX");
#ifdef WIN32
	_mktemp(filenamestorage);
#else
	mkstemp(filenamestorage);
#endif
	cout << "Temporary file name " << filenamestorage << endl;
	cout.flush();
	
	fBlockPos.Resize(0);
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=20;
	//	SetBlockSize(10);
	//fBlockPos.Resize(fNumHeaders);
	
	
	fFileName = filenamestorage;
	
	fIOStream = fopen(fFileName.c_str(),"wb"); //open for writing
	/**
	 *Writes NumHeaders and NumBlocks information in
	 *the two initial positions on fIOStream
	 */
	int zero = 0;
	fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	fwrite(&zero,sizeof(int),1,fIOStream);
	//cout << ftell(fIOStream) << endl;
	//fCurBlockPosition = ftell(fIOStream);
	
	
}

void TPZFileEqnStorage::ReOpen()
{
#ifdef LOG4CXX
	LOGPZ_DEBUG(logger,"reopening the file")
#endif
	fIOStream = fopen(fFileName.c_str(),"rb"); //open for reading
	/**
	 *Opens binary files and get initial information
	 *use this information for storage requirements
	 */
	fread(&fNumHeaders,sizeof(int),1,fIOStream);
	fread(&fNumBlocks,sizeof(int),1,fIOStream);
	//	ReadBlockPositions();
	
}

void TPZFileEqnStorage::OpenGeneric(char option, const char * name)
{
	fCurBlockPosition = -1;
	fNumBlocks=0;
	fCurrentBlock=0;
	fNumHeaders=11;
	//	SetBlockSize(10); 
	//fBlockPos.Resize(fNumHeaders);
	
	
	fFileName = name;
	
	if(option=='r'){
		fIOStream = fopen(fFileName.c_str(),"rb"); //open for reading
		/**
		 *Opens binary files and get initial information
		 *use this information for storage requirements
		 */
		fread(&fNumHeaders,sizeof(int),1,fIOStream);
		fread(&fNumBlocks,sizeof(int),1,fIOStream);
		ReadBlockPositions();
	}else if(option=='w'){
		fIOStream = fopen(fFileName.c_str(),"wb"); //open for writing
		/**
		 *Writes NumHeaders and NumBlocks information in
		 *the two initial positions on fIOStream
		 */
		int zero = 0;
		fwrite(&fNumHeaders,sizeof(int),1,fIOStream);
		//cout << ftell(fIOStream) << endl;
		fwrite(&zero,sizeof(int),1,fIOStream);
		//cout << ftell(fIOStream) << endl;
		//fCurBlockPosition = ftell(fIOStream);
	}
}
/*
 //double readvec[4][100];
 int iblock =0;
 for(iblock=0; iblock<4; iblock++) {
 long int currentpos = ftell(fIOStream);
 fseek(fIOStream,firstpos+iblock*sizeof(long int),SEEK_SET);
 fwrite(&currentpos,sizeof(long int),1,fIOStream);
 fseek(fIOStream,currentpos,SEEK_SET);
 for(i=0;i<loop_limit;i++){
 val=number*i*(iblock+1);
 fwrite(&val, sizeof(double), 1, fIOStream);
 }
 }
 fclose(fIOStream);
 fIOStream = fopen(name,"rb");
 fread(fPosition,sizeof(long int),5,fIOStream);
 for(iblock = 0; iblock<4; iblock++) {
 fseek(fIOStream,fPosition[iblock],SEEK_SET);
 fread(readvec[iblock],sizeof(double),loop_limit,fIOStream);
 }*/

void TPZFileEqnStorage::FinishWriting()
{
#ifdef LOG4CXX
	LOGPZ_DEBUG(logger,"Closing the binary file")
#endif
	fseek(fIOStream,sizeof(int),SEEK_SET);
	//cout << "Second fseek " << ftell(fIOStream) << endl;
	//fwrite(&fNumBlocks,sizeof(int),1,fIOStream);
	fwrite(&fNumBlocks,sizeof(int),1,fIOStream);
	fclose(fIOStream);
	fIOStream = 0;
}

void TPZFileEqnStorage::ReadBlockPositions()
{
#ifdef LOG4CXX
	LOGPZ_DEBUG(logger,"Reading block positions")
#endif
	int aux = fNumBlocks * (fNumHeaders-1);
	fBlockPos.Resize(aux);
	int i, ibl = 0;
	cout << "Reading Block Positions\n";
	cout.flush();
	for(i=0;i<fNumBlocks;i++) {
		cout << "*";
		cout.flush();
		if(!(i%20)){
			cout << 100*i/fNumBlocks << "% Read\n";
			cout.flush();
		}
		fread(&fBlockPos[ibl],sizeof(long int),fNumHeaders-1,fIOStream);
		ibl+=fNumHeaders-1;
		long int nextpos;
		fread(&nextpos,sizeof(long int),fNumHeaders-1,fIOStream);
		fseek(fIOStream,nextpos,SEEK_SET);
	}
	
	
}

std::string TPZFileEqnStorage::GetStorage() {return "File Storage";}
