// This file declares various functions that are useful over the course
// of the semester.  It must be "included" above the main program.
// Version Spring 2016.


// These are some standard C function libraries.
#include <stdlib.h>  // "standard library"
#include <cmath>   // various math functions, such as exp()
#include <stdio.h>   // various "input/output" functions
#include <time.h>    // functions for timing computations

// These functions are found in "Definitions.h".

// - Random number generators:
double MTUniform (unsigned int);
unsigned int Temper (unsigned int);
double LCGUniform (unsigned int);
double MWCUniform (unsigned int);
double MWC2Uniform (unsigned int);
double LCG64Uniform (int unsigned);

// - Histogram functions.
void Histogram (double, double, double, int, int);
void DiscreteHistogram (int, int, int, int);
void NormalHistogram (double, int, int);
void ExponentialHistogram (double, int, int);
void UniformHistogram (double, int, int);

// - Functions related to standard normals.
double PolarNormal (void);
double PsiTS (double);
double Psi (double);
double PsiInv (double);
void   CorrelatedNormals (double, double *);

// - Miscellanceous functions.
void   Pause (void);
void   Exit (void);
double Time (void);
int    Equal (double, double, double);
double **Array (int, int);
unsigned int GetUnsignedInteger (void);
double GetDouble (void);

// - Mortality Functions.
double H (char, double);
double HSlow (char, double);
double F (char, double);
double HConditional (char, double, double);



