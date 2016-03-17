/**
 *  TestMatLab: File with examples of how to use matlab inside PZ
 *
 *  @author Nathan Shauer
 *  @since Mar 2016
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>

#ifdef USING_MATLAB
#include "engine.h"
#include "mat.h"
#define  BUFSIZE 256
#endif

int mainEngine();
int mainGenerateMat();

// There are two examples. Choose by using the whichMain variable
int main()
{
  int whichMain = 0;
  if (whichMain == 0)
    mainEngine();
  else
    mainGenerateMat();
}

/*
 *	engdemo.c
 *
 *	A simple program to illustrate how to call MATLAB
 *	Engine functions from a C program.
 *
 * Copyright 1984-2011 The MathWorks, Inc.
 * All rights reserved
 */
int mainEngine()
{
  
#ifdef USING_MATLAB
  Engine *ep;
  mxArray *T = NULL, *result = NULL;
  char buffer[BUFSIZE+1];
  double time[10] = { 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 };
  
  /*
   * Call engOpen with a NULL string. This starts a MATLAB process
   * on the current host using the command "matlab".
   */
  if (!(ep = engOpen(""))) {
    fprintf(stderr, "\nCan't start MATLAB engine\n");
    return EXIT_FAILURE;
  }
  
  /*
   * PART I
   *
   * For the first half of this demonstration, send data
   * to MATLAB, analyze the data, and plot the result.
   */
  
  /*
   * Create a variable for the data
   */
  printf("------ PART I - Plotting ------\n");
  
  T = mxCreateDoubleMatrix(1, 10, mxREAL);
  memcpy((void *)mxGetPr(T), (void *)time, sizeof(time));
  /*
   * Place the variable T into the MATLAB workspace
   */
  engPutVariable(ep, "T", T);
  
  /*
   * Evaluate a function of time, distance = (1/2)g.*t.^2
   * (g is the acceleration due to gravity)
   */
  engEvalString(ep, "D = .5.*(-9.8).*T.^2;");
  
  /*
   * Plot the result
   */
  engEvalString(ep, "plot(T,D);");
  engEvalString(ep, "title('Position vs. Time for a falling object');");
  engEvalString(ep, "xlabel('Time (seconds)');");
  engEvalString(ep, "ylabel('Position (meters)');");
  
  /*
   * use fgetc() to pause long enough to be
   * able to see the plot
   */
  printf("Hit return to continue\n\n");
  fgetc(stdin);
  /*
   * We're done for Part I! Free memory, close MATLAB figure.
   */
  printf("Done for Part I.\n\n\n");
  mxDestroyArray(T);
  engEvalString(ep, "close;");
  
  
  /*
   * PART II
   *
   * For the second half of this demonstration, we will request
   * a MATLAB string, which should define a variable X.  MATLAB
   * will evaluate the string and create the variable.  We
   * will then recover the variable, and determine its type.
   */
  
  /*
   * Use engOutputBuffer to capture MATLAB output, so we can
   * echo it back.  Ensure first that the buffer is always NULL
   * terminated.
   */
  printf("------ PART II - Capture Matlab output ------\n");
  buffer[BUFSIZE] = '\0';
  engOutputBuffer(ep, buffer, BUFSIZE);
  while (result == NULL) {
    char str[BUFSIZE+1];
    /*
     * Get a string input from the user
     */
    printf("Enter a MATLAB command to evaluate.  This command should\n");
    printf("create a variable X.  This program will then determine\n");
    printf("what kind of variable you created.\n");
    printf("For example: X = 1:5\n");
    printf(">> ");
    
    fgets(str, BUFSIZE, stdin);
    
    /*
     * Evaluate input with engEvalString
     */
    engEvalString(ep, str);
    
    /*
     * Echo the output from the command.
     */
    printf("%s", buffer);
    
    /*
     * Get result of computation
     */
    printf("\nRetrieving X...\n");
    if ((result = engGetVariable(ep,"X")) == NULL)
      printf("Oops! You didn't create a variable X.\n\n");
    else {
      printf("X is class %s\t\n", mxGetClassName(result));
    }
  }
  
  /*
   * We're done! Free memory, close MATLAB engine and exit.
   */
  printf("\n\nDone for Part II.\n");
  printf("Hit return to continue\n\n");
  fgetc(stdin);
  
  mxDestroyArray(result);
 
  /*
   * PART III - Calculating eigenvalues of a matrix
   *
   */
  printf("------ PART III - Calculate eigenvalues ------\n");
  
  mxArray *myMat = NULL;
  double data[9] = { 3.0, 0.0, 0.0, 0.0, 8.0, 0.0, 0.0, 0.0, 15.0 }; // matrix in columns
  myMat = mxCreateDoubleMatrix(3, 3, mxREAL);
  memcpy((void *)(mxGetPr(myMat)), (void *)data, sizeof(data)); // copy data to myMat

  engPutVariable(ep,"myMat", myMat); // Putting myMat on Engine
  engEvalString(ep,"myMat"); // Printing myMat on Engine
  printf("%s", buffer); // Showing in console what is on the buffer of the Engine
  
  engEvalString(ep,"eigOfMyMat = eig(myMat)"); // Calculating the Eigenvalues
  printf("%s", buffer); // Showing in console what is on the buffer of the Engine
    
  printf("Finished!\n");
  mxDestroyArray(myMat);
  engClose(ep);
  
  printf("\n==> For more examples go to PATH_TO_YOUR_MATLAB/extern/examples\n\n");

#else
  std::cout << "You have to check USING_MATLAB_ENGINE on CMake to run this test!" << std::endl;
#endif
  
  return EXIT_SUCCESS;
}


/*
 * matcreat.cpp - MAT-file creation program
 *
 * See the MATLAB External Interfaces/API Guide for compiling information.
 *
 * Calling syntax:
 *
 *   matcreat
 *
 * Create a MAT-file which can be loaded into MATLAB.
 *
 * This program demonstrates the use of the following functions:
 *
 *  matClose
 *  matGetVariable
 *  matOpen
 *  matPutVariable
 *  matPutVariableAsGlobal
 *
 * Copyright 1984-2007 The MathWorks, Inc.
 */
int mainGenerateMat()
{
  
#ifdef USING_MATLAB
  
  MATFile *pmat;
  mxArray *pa1, *pa2, *pa3;
  std::vector<int> myInts;
  myInts.push_back(1);
  myInts.push_back(2);
  printf("Accessing a STL vector: %d\n", myInts[1]);
  double data[9] = { 1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0 };
  const char *file = "mattest.mat";
  char str[BUFSIZE];
  int status;
  
  printf("Creating file %s...\n\n", file);
  pmat = matOpen(file, "w");
  if (pmat == NULL) {
    printf("Error creating file %s\n", file);
    printf("(Do you have write permission in this directory?)\n");
    return(EXIT_FAILURE);
  }
  
  pa1 = mxCreateDoubleMatrix(3,3,mxREAL);
  if (pa1 == NULL) {
    printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
    printf("Unable to create mxArray.\n");
    return(EXIT_FAILURE);
  }
  
  pa2 = mxCreateDoubleMatrix(3,3,mxREAL);
  if (pa2 == NULL) {
    printf("%s : Out of memory on line %d\n", __FILE__, __LINE__);
    printf("Unable to create mxArray.\n");
    return(EXIT_FAILURE);
  }
  memcpy((void *)(mxGetPr(pa2)), (void *)data, sizeof(data));
  
  pa3 = mxCreateString("MATLAB: the language of technical computing");
  if (pa3 == NULL) {
    printf("%s :  Out of memory on line %d\n", __FILE__, __LINE__);
    printf("Unable to create string mxArray.\n");
    return(EXIT_FAILURE);
  }
  
  status = matPutVariable(pmat, "LocalDouble", pa1);
  if (status != 0) {
    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
    return(EXIT_FAILURE);
  }
  
  status = matPutVariableAsGlobal(pmat, "GlobalDouble", pa2);
  if (status != 0) {
    printf("Error using matPutVariableAsGlobal\n");
    return(EXIT_FAILURE);
  }
  
  status = matPutVariable(pmat, "LocalString", pa3);
  if (status != 0) {
    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
    return(EXIT_FAILURE);
  }
  
  /*
   * Ooops! we need to copy data before writing the array.  (Well,
   * ok, this was really intentional.) This demonstrates that
   * matPutVariable will overwrite an existing array in a MAT-file.
   */
  memcpy((void *)(mxGetPr(pa1)), (void *)data, sizeof(data));
  status = matPutVariable(pmat, "LocalDouble", pa1);
  if (status != 0) {
    printf("%s :  Error using matPutVariable on line %d\n", __FILE__, __LINE__);
    return(EXIT_FAILURE);
  }
  
  /* clean up */
  mxDestroyArray(pa1);
  mxDestroyArray(pa2);
  mxDestroyArray(pa3);
  
  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(EXIT_FAILURE);
  }
  
  /*
   * Re-open file and verify its contents with matGetVariable
   */
  pmat = matOpen(file, "r");
  if (pmat == NULL) {
    printf("Error reopening file %s\n", file);
    return(EXIT_FAILURE);
  }
  
  /*
   * Read in each array we just wrote
   */
  pa1 = matGetVariable(pmat, "LocalDouble");
  if (pa1 == NULL) {
    printf("Error reading existing matrix LocalDouble\n");
    return(EXIT_FAILURE);
  }
  if (mxGetNumberOfDimensions(pa1) != 2) {
    printf("Error saving matrix: result does not have two dimensions\n");
    return(EXIT_FAILURE);
  }
  
  pa2 = matGetVariable(pmat, "GlobalDouble");
  if (pa2 == NULL) {
    printf("Error reading existing matrix GlobalDouble\n");
    return(EXIT_FAILURE);
  }
  if (!(mxIsFromGlobalWS(pa2))) {
    printf("Error saving global matrix: result is not global\n");
    return(EXIT_FAILURE);
  }
  
  pa3 = matGetVariable(pmat, "LocalString");
  if (pa3 == NULL) {
    printf("Error reading existing matrix LocalString\n");
    return(EXIT_FAILURE);
  }
  
  status = mxGetString(pa3, str, sizeof(str));
  if(status != 0) {
    printf("Not enough space. String is truncated.");
    return(EXIT_FAILURE);
  }
  if (strcmp(str, "MATLAB: the language of technical computing")) {
    printf("Error saving string: result has incorrect contents\n");
    return(EXIT_FAILURE);
  }
  
  /* clean up before exit */
  mxDestroyArray(pa1);
  mxDestroyArray(pa2);
  mxDestroyArray(pa3);
  
  if (matClose(pmat) != 0) {
    printf("Error closing file %s\n",file);
    return(EXIT_FAILURE);
  }
  printf("Done\n");

#else
  std::cout << "You have to check USING_MATLAB_ENGINE on CMake to run this test!" << std::endl;
#endif
  
  return(EXIT_SUCCESS);
}
