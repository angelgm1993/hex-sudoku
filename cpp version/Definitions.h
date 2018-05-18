// This file contains various functions that are useful over the course
//   of the semester.  They must be declared above the main program and this
//   file must be "included" following the main program.
// Version Spring 2016.


// First are some random number generators.

////////////////////////////////////////////////////////////////////////////////
// MERSENNE TWISTER
// By M. Matsumoto and T. Nishimura (1998).
// "Mersenne Twister: a 623-dimensionally equidistributed uniform pseudo-random
//   number generator".
// ACM Transactions of Modeling and Computer Simulation 8(1):3-30.
// Any coding errors introduced are my own (C.D. Howard).

// An unsigned integer is represented by 32 bits in base 2.  The largest unsigned integer is:
// 2^32 - 1 = 4,294,967,295 (base 10);
//          = ffffffff (base 16);
//          = 1111 1111 1111 1111 1111 1111 1111 1111 (base 2).
// The digits in hexadecimal (base 16) are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f.

#include <cstdio>
#include <stdio.h>


double MTUniform (unsigned int seed) {

   static unsigned int X[1248], m[2], initialized = 0, k;
   unsigned int N, Y;

   // Seed the RNG when a new seed is passed or it has not yet been initialized.
   if (seed || !initialized) {
      // If no seed is specified, default is 1.
      X[0] = (seed ? seed : 1);
      // Seed X[1],X[2],...,X[623] with your favorite LCG.
      for (k = 1; k < 624; k++) {
         X[k] = 22695477 * X[k-1] + 1;
      }
      m[0] = 0; m[1] = 0x9908b0df;
      // The counter "k" is now 624.
      initialized = 1;
   }

   // Concatenate the first bit of X[k-624] with the last 31 bits of X[k-623],
   //    "|" is "bit-wise or".
   Y = (X[k-624] & 0x80000000) | (X[k-623] & 0x7fffffff);

   // Now apply the invertible linear transformation A to Y and bit-wise add X[k-227].
   X[k] = ((Y >> 1) ^ m[Y & 1] ^ X[k-227]);

   // Re-load X[0],X[1],...,X[623] as you go.
   X[k-624] = X[k];

   // Apply the tempering function.
   N = Temper (X[k]);

   // Increment the counter; shift vectors when appropriate.
   k ++;
   if (k == 1248) k = 624;

   // Now 0 <= N <= 4,294,967,295; scale it to be on the interval (0,1).
   return ( (N + 0.5) / 4294967296.0 );

}




////////////////////////////////////////////////////////////////////////////////
// Tempering function used by the Mersenne Twister.
unsigned int Temper (unsigned int N) {

   N ^= (N >> 11);
   N ^= (N << 7) & 0x9d2c5680;
   N ^= (N << 15) & 0xefc60000;
   N ^= (N >> 18);

   return N;

}




////////////////////////////////////////////////////////////////////////////////
////// A typical 32 bit Linear Congruential Generator.

double LCGUniform (unsigned int seed) {

   // Static variables retain their values between function calls.
   static unsigned int N=1, a = 22695477, c = 1, initialized = 0;

   // Seed the RNG when a new seed is passed or it has not yet been initialized.
   if (seed || !initialized) {
      // If no seed is specified, default is 1.
      N = (seed ? seed : 1);
      initialized = 1;
   }

   // Generate the next number in the sequence.
   // Here mod 2^32 is automatic since an unsigned integer is represented by 32 bits.
   N = a * N + c;

   // Now 0 <= N <= 4,294,967,295; scale it to be on the interval (0,1).
   return ( (N + 0.5) / 4294967296.0 );

}


////////////////////////////////////////////////////////////////////////////////
////// A Typical Multiply With Carry Generator.

// This function generates random numbers uniformly on [0,1].
// The numbers generated are of the form (n + 0.5) / 2^32, where 0 <= n <= 4,294,967,295 = 2^32 - 1.

// It uses the multiply-with-carry algorithm with m = 2^32, a = 4,294,967,118, and initial c = 1.
// Its period is 9,223,371,654,602,686,463 which is approximately 10^19.

// An unsigned integer is represented by 32 bits in base 2.  The largest unsigned integer is:
// 2^32 - 1 = 4,294,967,295 (base 10);
//          = ffffffff (base 16);
//          = 1111 1111 1111 1111 1111 1111 1111 1111 (base 2).
// The digits in hexadecimal (base 16) are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, a, b, c, d, e, f.


double MWCUniform (unsigned int seed) {

   // Static variables retain their values between function calls.
   static unsigned int a1, a0, n1, n0, c1, c0, initialized=0;

   unsigned int w1, w0, x1, x0, y1, y0, z1, z0, carry, N;

   // This function is called only be MWCUniform and follows below.
   void Split (unsigned int, unsigned int *, unsigned int *);

   // Seed the MWC function when it is passed a non-zero seed or it is not yet
   //    initialized.
   // Represent a, n, and c as with pairs of 16-bit numbers:
   //  a = a1 * 2^16 + a0 (the multiplier)
   //  n = n1 * 2^16 + n0 (initially the seed)
   //  c = c1 * 2^16 + c0 (the increment, initially 1)
   if (seed || !initialized) {
      seed = (seed ? seed : 1);
      Split (4294967118, &a1, &a0);
      Split (seed, &n1, &n0);
      c0 = 1;                                 // first 16 bits of "c"
      c1 = 0;                                 // second 16 bits
      initialized = 1;
   }

   // Generate the next number in the sequence.

   // Compute w1, w0, where a0 * n0 = w1 * 2^16 + w0.
   Split (a0 * n0, &w1, &w0);

   // Compute x1, x0, where a1 * n0 = x1 * 2^16 + x0.
   Split (a1 * n0, &x1, &x0);

   // Compute y1, y0, where a0 * n1 = y1 * 2^16 + y0.
   Split (a0 * n1, &y1, &y0);

   // Compute z1, z0, where a1 * n1 = z1 * 2^16 + z0.
   Split (a1 * n1, &z1, &z0);

   // Compute next n0, n1, c0, c1, where an + c = c1 * 2^48 + c0 * 2^32 + n1 * 2^16 + n0.

   // First compute n0 and whatever is carried.
   Split (w0 + c0, &carry, &n0);

   // Compute n1 and whatever is carried.
   Split (carry + w1 + x0 + y0 + c1, &carry, &n1);

   // Compute c0 and whatever is carried.
   Split (carry + x1 + y1 + z0, &carry, &c0);

   // Compute c1; nothing is carried.
   c1 = z1 + carry;

   // Re-assemble n1 and n0.
   N = (n1 << 16) + n0;

   return ((N + 0.5) / 4294967296.0);

}

// This function computes 16-bit numbers x1 and x0 such that x = x1 * 2^16 + x0.
// Pointers *x1 and *x0 are used because a function can return only one value.
void Split (unsigned int x, unsigned int *x1, unsigned int *x0) {

   // First 16 bits of x.
   *x0 = x & 0xffff;

   // Second 16 bits of x.
   *x1 = x >> 16;

   return;

}




////////////////////////////////////////////////////////////////////////////////
////// 64 bit LCG by Donald Knuth

// This LCG is designed for a 64 bit processor.  Here it is re-written to run
//    on a 32 bit processor.  Since it cannot take full advantage of 64 bit
//    processing this version is time inefficient.
// Here a = 6364136223846793005 = 22609 * 2^48 + 62509 * 2^32 + 19605 * 2^16 + 32557
//  and c = 1442695040888963407 =  5125 * 2^48 + 31614 * 2^32 + 63335 * 2^16 + 33103.

double LCG64Uniform (int unsigned seed) {

   static unsigned int a[] = {32557, 19605, 62509, 22609};
   static unsigned int x[4], y[4], *to, *from, initialized=0;
   unsigned int N, carry, *tmp;
   int i, j, k;

   // Seed the RNG when a new seed is passed or it has not yet been initialized.
   if (seed || !initialized) {
      // If no seed is specified, default is 1.
      seed = (seed ? seed : 1);
      x[1] = seed >> 16;
      x[0] = seed & 0x0000ffff;
      y[0] = y[1] = y[2] = y[3] = x[2] = x[3] = 0;
      to = x;
      from = y;
      initialized=1;
      // This RNG is sensitive to poor initialization.
      LCG64Uniform (0);
      LCG64Uniform (0);
      LCG64Uniform (0);
      LCG64Uniform (0);
   }

   // Toggle the pointers.
   tmp = to;
   to = from;
   from = tmp;

   // Intialize "to" to "c".
   to[0] = 33103; to[1] = 63335; to[2] = 31614; to[3] = 5125;

   // Now multiply a*from and add it to "to".
   for (k = 0; k <= 3; k++) {
      // Work on to[k]:
      for (i = 0; i <= k; i++) {
         // Multiply:
         to[k] += a[i] * from[k-i];
         // Carry any overflow exceeding 16 bits:
         j = k;
         while (1) {
            carry = (to[j] >> 16);
            to[j] = (to[j] & 0x0000ffff);
            j++;
            if (j == 4 || carry == 0) {
               break;
            }
            to[j] += carry;
         } // while
      } // i loop
   } // k loop

   // Re-assemble to[2] and to[3].
   N =  to[2] + (to[3]<<16);

   return ( (N + 0.5) / 4294967296.0 );

}



////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then CONTINUES the program.
void Pause () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to continue program... ");
   fgets (input, 9, stdin);

   return;

}



////////////////////////////////////////////////////////////////////////////////
// This function waits for a user-input "enter", then EXITS the program.
// It prevents the window from closing up before the output can be viewed.
void Exit () {

   char input[100];

   printf ("\n");
   printf ("Hit Enter to exit program... ");
   fgets (input, 9, stdin);

   exit (0);

}


////////////////////////////////////////////////////////////////////////////////
// This function computes elapsed computation time in seconds.
double Time () {

   static clock_t time;
   static int initialized = 0;

   if (!initialized) {
      time = clock ();
      initialized = 1;
   }

   return (clock() - time) / CLOCKS_PER_SEC;

}      


////////////////////////////////////////////////////////////////////////////////
// This function gets a non-negative integer value typed in by the user. ///////
unsigned int GetUnsignedInteger () {

   char input[100];
   unsigned int n;

   fgets (input, 99, stdin);
   sscanf (input, "%d", &n);

   return (n);

}



////////////////////////////////////////////////////////////////////////////////
// This function gets a double precision value typed in by the user. ///////////
double GetDouble () {

   char input[100];
   double x;

   fgets (input, 99, stdin);
   sscanf (input, "%lf", &x);

   return (x);

}




////////////////////////////////////////////////////////////////////////////////
// This function creates a histogram of randomly generated numbers. It is
// typically used with continuous (not discrete)  data.
////////////////////////////////////////////////////////////////////////////////
void   Histogram (double x,          // Item of data to be added to the histogram.
                  double a,          // Lower bound of histogram.
                  double b,          // Upper bound of histogram.
                  int n,             // Number of "buckets".
                  int done)          // 0 if not done, 1 if done.

{

   int k;
   double x0, x1, biggest;
   FILE *fp;

   static int initialized=0;
   static double dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

   }


   // If finished, create the TeX files that generate the plot. This is the
   // final call to the histogram function.
    else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "      From         To  Frequency\n");
      fprintf (fp, "========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f\n", x0, x1, freq[k]);
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      // Generate TeX output files.
      // Scale so that the biggest bucket has 1 (for graphing convenience).

      // First find the largest bucket value.
      biggest = 0;
      for (k = 0; k < n; k++) {
         if (freq[k] > biggest) {
            biggest = freq[k];
         }
      }

      // Now re-scale all the bucket values.
      for (k = 0; k < n; k++) {
         freq[k] /= biggest;
      }

      // Report data to the output file.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1) * dx, freq[k]);
      }
      fclose (fp);

      // Free up the freq[] array.
      free (freq);

      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <%8.3f truein, 2.5 truein>\n", 4.0 / (b-a));
      fprintf (fp, "\\setplotarea x from %8.3f to %8.3f, y from  0 to 1.0\n", a, b);
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from %8.2f to %8.2f by %8.2f\n", a, b, (b-a)/10);
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot %8.3f 0  %8.3f 0  %8.3f 1  %8.3f 1 /\n",
                a-.045*(b-a), a-.03*(b-a), a-.03*(b-a), a-.045*(b-a) );
      fprintf (fp, "\\put {0} [cr] at %8.3f  0\n", a-.055*(b-a));
      fprintf (fp, "\\put {%d} [cr] at %8.3f  1\n", (int) biggest, a-.055*(b-a));
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\put {\\sl Histogram for a Continuous Random Variable} at %8.3f 1.2\n", (a+b)/2.0);
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);


   }

   return;

}



////////////////////////////////////////////////////////////////////////////////
// This function creates a histogram of randomly generated numbers. It is
// typically used with discrete data, like binomially distributed data.
////////////////////////////////////////////////////////////////////////////////
void   DiscreteHistogram (int x,          // Data to be added to the histogram.
                          int a,          // Lower bound of histogram.
                          int b,          // Upper bound of histogram.
                          int done)       // 0 if not done, 1 if done.

{

   static double *freq;
   static int initialized = 0;
   int k;
   double biggest;
   FILE *fp;

   // Allocate space on the first call to this function.
   if (! initialized) {

      freq = (double *) calloc (b-a+2, sizeof (double));

      initialized = 1;

   }

   // If not yet finished, add the current item of data to the histogram.
   if (!done) {

      // Truncate if data is below the lower bound.
      if (x < a) {
         x = a;
      }

      // Truncate if data is above the upper bound.
      if (x > b) {
         x = b;
      }

      // Increment the appropriate bucket.
      freq[x-a] ++;


   }


   // If finished, create the TeX files that generate the plot. This is the
   // final call to the histogram function.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "         X  Frequency\n");
      fprintf (fp, "========== ==========\n");
      for (k = a; k <= b; k++) {
         fprintf (fp, "%10d %10.0f\n", k, freq[k]);
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <%8.3f truein, 2.5 truein>\n", 5.0 / (b-a));
      fprintf (fp, "\\setplotarea x from %d to %d, y from  0 to 1.0\n", a, b);
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from %d to %d by 1\n", a, b);
      fprintf (fp, "/\n");
      fprintf (fp, "\\linethickness=%ftruein\n", 1.0/(b-a));
      fprintf (fp, "\\input HistogramData.txt\n");
      fprintf (fp, "\\put {\\sl Histogram for a Discrete Random Variable} at %8.3f 1.2\n", (a+b)/2.0);
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the biggest bucket has 1 (for graphing convenience).

      // First find the largest bucket.
      biggest = 0;
      for (k = 0; k <= b-a; k++) {
         if (freq[k] > biggest) {
            biggest = freq[k];
         }
      }

      // Now scale all the buckets.
      for (k = 0; k <= b-a; k++) {
         freq[k] /= biggest;
      }

      // Report data to the output file.
      fp = fopen ("HistogramData.txt", "w");
      for (k = 0; k <= b-a; k++) {
         fprintf (fp, "\\putrule from %d 0 to %d %f \n", a+k, a+k, freq[k]);
      }
      fclose (fp);

      // Free up the freq[] array.
      free (freq);

   }

   return;

}


////////////////////////////////////////////////////////////////////////////////
// This histogram function is specially designed for data that should be
// Normal(0,1) in distribution.  The function computes a properly scaled
// standard normal density function for comparison to the data.
////////////////////////////////////////////////////////////////////////////////
void NormalHistogram (double x,          // Data to be added to the histogram.
                      int n,             // Number of "buckets".
                      int done)          // 0 if not done, 1 if done.

{

   int k;
   double z, x0, x1, area;
   FILE *fp;

   static int N=0, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the standard normal.
      a = -5; b = 5;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

   }

   // ...otherwise, when finished, create the TeX output files for viewing.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n",
                       x0, x1, freq[k], N * (Psi(x1) - Psi(x0)));
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in first or last bucket.\n");
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <0.45 truein, 1.5 truein>\n", 5.0 / (b-a));
      fprintf (fp, "\\setplotarea x from -5 to 5, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from -5 to 5 by 1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot \"Normal.txt\"\n");
      fprintf (fp, "\\plot -5.2 0 -5.1 0 -5.1 1  -5.2 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -5.3 0\n");
      fprintf (fp, "\\put {$\\frac{1}{\\sqrt{2\\pi}}$} [cr] at -5.3 1\n");
      fprintf (fp, "\\put {\\sl Standard Normal Histogram} at 0 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the area is sqrt(2pi) (for comparing to the Normal(0,1) density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * dx;
      }
      area /= sqrt(2*3.14159);
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report the histogram to an output file.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1)*dx, freq[k]);
      }
      fclose (fp);

      // Report the appropriately scaled normal density function for comparison.
      fp = fopen ("Normal.txt", "w");
      for (z = -5; z <= 5; z += .01) {
         fprintf (fp, "%8.3f %8.4f\n", z, exp(-z*z/2));
      }
      fclose (fp);

   }

   return;

}


////////////////////////////////////////////////////////////////////////////////
// This function creates a histogram of the variable x and compares it to
// the density function for a mean 1 exponential.
////////////////////////////////////////////////////////////////////////////////
void ExponentialHistogram (double x,          // data to be added to the histogram
                           int n,             // number of "buckets"
                           int done)          // 1 when done, 0 otherwise

{

   int k;
   double z, x0, x1, area;
   FILE *fp;

   static int N=0, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the mean 1 exponential distribution.
      // The cutoff of 5 is arbitrary but it captures most of the action.
      a = 0; b = 5;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data...
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

      // ...otherwise, when finished, create the TeX output files for viewing.
   } else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1)* dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n",
                       x0, x1, freq[k], N * (exp(-x0) - exp(-x1)));
      }
      fprintf (fp, "\n");
      fprintf (fp, "Data outside histogram range is put in last bucket.\n");
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <1 truein, 1.5 truein>\n");
      fprintf (fp, "\\setplotarea x from 0 to 5, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from 0 to 5 by 1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\put {$\\leftarrow$} at 5.1 %8.3f\n", exp(-b) * n / (b-a));
      fprintf (fp, "\\plot \"Exponential.txt\"\n");
      fprintf (fp, "\\plot -.3 0 -.2 0 -.2 1  -.3 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -.35 0\n");
      fprintf (fp, "\\put {1} [cr] at -.35 1\n");
      fprintf (fp, "\\put {\\sl Mean-1 Exponential Histogram} at 2.5 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale frequencies so that the area is 1 (for comparing to the Exp(1) density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * (b-a) / n;
      }
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report histogram data.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n",
                       a + (k + 1) * ((b - a) / n), freq[k]);
      }
      fclose (fp);

      // Report the exp(1) density function.
      fp = fopen ("exponential.txt", "w");
      for (z = 0; z <= b; z += .01) {
         fprintf (fp, "%8.3f %8.4f\n", z, exp(-z));
      }
      fclose (fp);

   }

   return;

}


////////////////////////////////////////////////////////////////////////////////
// This program creates a histogram of the data x, and plots the Uniform[0,1]
// density function for comparison.
////////////////////////////////////////////////////////////////////////////////
void UniformHistogram (double x,         // data to be added to the histogram
                      int n,             // number of "buckets"
                      int done)          // 1 when done, 0 otherwise

{

   int k;
   double x0, x1, area;
   FILE *fp;

   static int N=0, initialized=0;
   static double a, b, dx, *freq;

   // Initialize certain variables on the first call to this function.
   if (!initialized) {

      // Allocate array space for frequencies.
      freq = (double *) calloc (n, sizeof (double));

      // Hardwire upper and lower bounds for the Uniform[0,1].
      a = 0; b = 1;

      // Compute size of each "bucket".
      dx = (b - a) / n;

      // Variables have now been initialized.
      initialized = 1;

   }

   // If adding a new piece of data.
   if (!done) {

      // See if data is below the lower bound.
      if (x <= a) {
         freq[0] ++;
      }

      // See if data is above the upper bound.
      else if (x >= b) {
         freq[n-1] ++;
      }

      // Otherwise find the appropriate bucket.
      else {
         k = (x - a) / dx;
         freq[k] ++;
      }

      // Increment the data counter.
      N ++;

   }

   // When finished, create the output files for viewing.
   else {

      // Report the raw data to a TXT file for viewing without TeX.
      fp = fopen ("HistogramTxtData.txt", "w");
      fprintf (fp, "                                   Expected\n");
      fprintf (fp, "      From         To  Frequency  Frequency\n");
      fprintf (fp, "========== ========== ========== ==========\n");
      for (k = 0; k < n; k++) {
         x0 = a + k * dx;
         x1 = a + (k+1) * dx;
         fprintf (fp, "%10.5f %10.5f %10.0f %10.0f\n", x0, x1, freq[k], N*dx);
      }
      fclose (fp);

      // Create the main TeX file.
      fp = fopen ("Histogram.tex", "w");
      fprintf (fp, "\\input pictex\\magnification=\\magstep1\\nopagenumbers\n");
      fprintf (fp, "\\centerline {\n");
      fprintf (fp, "\\beginpicture\n");
      fprintf (fp, "\\setcoordinatesystem units <4.5 truein, 1.5 truein>\n", 5.0 / (b-a));
      fprintf (fp, "\\setplotarea x from 0 to 1, y from  0 to 1.0\n");
      fprintf (fp, "\\axis bottom\n");
      fprintf (fp, "ticks numbered from 0 to 1 by .1\n");
      fprintf (fp, "/\n");
      fprintf (fp, "\\plot 0 1  1 1 /\n");    // the Uniform[0,1] density function
      fprintf (fp, "\\plot -.06 0 -.05 0 -.05 1  -.06 1 /\n");
      fprintf (fp, "\\put {0} [cr] at -.07 0\n");
      fprintf (fp, "\\put {1} [cr] at -.07 1\n");
      fprintf (fp, "\\put {\\sl Uniform (0,1) Histogram} at 0.5 1.2\n");
      fprintf (fp, "\\sethistograms\n");
      fprintf (fp, "\\plot \"HistogramData.txt\"\n");
      fprintf (fp, "\\endpicture}\\vfill\\end\n");
      fclose (fp);

      // Scale so that the area is 1 (for comparing to the uniform[0,1] density).
      area = 0;
      for (k = 0; k < n; k++) {
         area += freq[k] * dx;
      }
      for (k = 0; k < n; k++) {
         freq[k] /= area;
      }

      // Report the histogram data.
      fp = fopen ("HistogramData.txt", "w");
      fprintf (fp, "%10.5f  %10.5f\n", a, 0.0);
      for (k = 0; k < n; k++) {
         fprintf (fp, "%10.5f  %10.5f\n", a + (k+1)*dx, freq[k]);
      }
      fclose (fp);


   }

   return;

}



////////////////////////////////////////////////////////////////////////////////
// Test to see if a = b up to an allowable error of epsilon.
////////////////////////////////////////////////////////////////////////////////
int Equal (double a, double b, double epsilon) {

   int equal;

   equal = (fabs(a-b) <= epsilon);

   return (equal);

}


////////////////////////////////////////////////////////////////////////////////
// Allocate array space for an n x m array.
////////////////////////////////////////////////////////////////////////////////
double **Array (int n, int m) {

   int i;
   double **A;

   A = (double **) calloc (n+1, sizeof (double *));
   for (i = 0; i <= n; i++) {
      A[i] = (double *) calloc (m+1, sizeof (double));
   }

   return (A);

}


////////////////////////////////////////////////////////////////////////////////
// Normal random number generator (polar method).
////////////////////////////////////////////////////////////////////////////////
double PolarNormal () {

   static int compute_a_new_pair=1;
   double static Y;

   int accepted;
   double U, V, A, B, S, W, X;

   // The flag compute_a_new_pair is 1 if both X and Y  have been used and
   //  a new pair (X,Y) must be computed.  X will then get used first and
   //  compute_a_new_pair is set to 0 until Y is also used.

   // The flag compute_a_new_pair is 0 if (X,Y) have recently been computed and
   //  only X has been used so far.  Y gets used and compute_a_new_pair
   //  is set to 1 so that a new pair will get computed with the next function call.

   if (compute_a_new_pair) {

      // Generate (A,B) uniformly distributed on the unit disk.
      accepted = 0;
      while (!accepted) {

         // Get two fresh uniforms.
         U = MTUniform (0);   V = MTUniform (0);

         // Compute A and B uniformly distributed on [-1,1] x [-1,1].
         A = -1 + 2.0 * U;
         B = -1 + 2.0 * V;

         // Accept if (A,B) lies within the unit disk.
         S = A*A + B*B;
         accepted = (S < 1);

      }

      // Compute X and Y via the polar method.
      W = sqrt ( -2*log(S) / S);
      X = A * W;
      Y = B * W;

      // Flip the flag.
      compute_a_new_pair = 0;

      // Use X.
      return (X);

   } else {

      // Flip the flag.
      compute_a_new_pair = 1;

      // Use Y.
      return (Y);

   }

}


////////////////////////////////////////////////////////////////////////////////
// This function computes Psi(x) via its Taylor series expansion.  With 80
// terms the error will < 1e-10 for -6 <= x <= 6.  For values of x outside this
// range the computation is handled differently, still producing at least 9
// significant digits.
////////////////////////////////////////////////////////////////////////////////
double PsiTS (double x) {

   static double scale = 1.0/sqrt(2*3.14159265358979), r[200];
   static int initialized = 0, N=80;
   double y, sum, term, w;
   int n;

   if (!initialized) {
      for (n = 1; n <= N; n++) {
         r[n] = -(2*n-1.0)/(2*n+1.0)/n/2.0;
      }
      initialized = 1;
   }

   if (x < -6) {
      y = -(scale/x) * exp(-x*x/2);
   }

   else if (x <= 6) {
      sum = term = x;
      w = x*x;
      for (n = 1; n <= N; n++) {
         term *= (w * r[n]);
         sum += term;
      }
      y = 0.5 + scale * sum;
   }

   else {
      y = 1 - (scale/x) * exp(-x*x/2);
   }

   return (y);

}

////////////////////////////////////////////////////////////////////////////////
// Normal (cumulative) distribution function evaluated at "x".
// It runs almost 3 times as fast as the more accurate PsiTS()
// and is accurate to 7 decimal places for all values of x.
////////////////////////////////////////////////////////////////////////////////
double Psi (double x) {

   static double    A =    0.0293868870,
                    B =    0.2161934875,
                    C =    0.6503029656,
                    D =    0.7978845608,
                    E =    0.0594864800,
                    F =    0.4160822924;

   int sign=1;
   double R;

   if (x < 0) {
      sign = -1;
      x *= -1;
   }

   x = fabs(x);
   R = (((A*x+B)*x+C)*x+D)*x / ((E*x+F)*x+1);

   // Return Psi (x)
   return (0.5 + sign * 0.5 * (1.0 - exp (-R)));

}


////////////////////////////////////////////////////////////////////////////////
// This algorithm is due to Peter John Acklam. Any coding errors are my own.
// It is a very good approximation of Psi^(-1).
////////////////////////////////////////////////////////////////////////////////
double PsiInv (double u) {

   static double
    A1 =  -3.969683028665376e+01,
    A2 =   2.209460984245205e+02,
    A3 =  -2.759285104469687e+02,
    A4 =   1.383577518672690e+02,
    A5 =  -3.066479806614716e+01,
    A6 =   2.506628277459239e+00,
    B1 =  -5.447609879822406e+01,
    B2 =   1.615858368580409e+02,
    B3 =  -1.556989798598866e+02,
    B4 =   6.680131188771972e+01,
    B5 =  -1.328068155288572e+01,
    C1 =  -7.784894002430293e-03,
    C2 =  -3.223964580411365e-01,
    C3 =  -2.400758277161838e+00,
    C4 =  -2.549732539343734e+00,
    C5 =   4.374664141464968e+00,
    C6 =   2.938163982698783e+00,
    D1 =   7.784695709041462e-03,
    D2 =   3.224671290700398e-01,
    D3 =   2.445134137142996e+00,
    D4 =   3.754408661907416e+00,
    P0 =  0.02425,
    P1 =  0.97575;

   double N, q, r;

   // Left tail.
   if (u < P0) {
      q = sqrt(-2*log(u));
      N = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   // Middle, general case.
   else if (u <= P1) {
      q = u - 0.5;
      r = q*q;
      N = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
   }

   // Right tail.
   else {
      q = sqrt(-2*log(1.0-u));
      N = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
   }

   return (N);

}




////////////////////////////////////////////////////////////////////////////////
// This function generates two correlated standard normals using the
//    polar method.   They are put into X[1] and X[2].
////////////////////////////////////////////////////////////////////////////////
void CorrelatedNormals (double rho, double *X) {

   // Make sure rho is in the right range.
   if (fabs(rho) > 1) {
      printf ("Rho is too big in CorrelatedNormals.\n");
      Pause ();
   }

   // Generate two independent standard normals.
   X[1] = PolarNormal ();
   X[2] = PolarNormal ();

   // Transform X[2] so that Correlation(X[1],X[2]) = rho.
   X[2] = rho * X[1] + sqrt(1.0 - rho*rho) * X[2];

   return;

}


////////////////////////////////////////////////////////////////////////////////
// The next three functions are mortality functions based on the SSA table for
//   males and females.
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// This function computes F(x) = P[lifespan <= x] from SSA table.
////////////////////////////////////////////////////////////////////////////////
double F(char gender, double t) {

   // Static variables.
   static int initialized=0;
   static double *pmale, *pfemale;

   // Non-static variables.
   int i, i0, NumberOfMaleLives, NumberOfFemaleLives;
   double w0, w1, p0, *p;
   char input[100];
   FILE *fp;

   // Scan in the SSA mortality table for females the first time this function is called.
   if (!initialized) {

      // Allocate array space.
      pmale   = (double *) calloc (120, sizeof (double));
      pfemale = (double *) calloc (120, sizeof (double));

      // Open the SSA life table.
      fp = fopen ("SSALifeTable.txt", "r");

      // Scan in the data (2nd and 3rd columns) and convert to a death probability.
      for (i = 0; i <= 115; i++) {
         fgets (input, 99, fp);
         sscanf (input, "%d %d %d", &i0, &NumberOfMaleLives, &NumberOfFemaleLives);
         pmale[i]   = 1.0 - NumberOfMaleLives / 100000.0;
         pfemale[i] = 1.0 - NumberOfFemaleLives / 100000.0;
      }

      // Close the data file.
      fclose (fp);

      // Now pmale[*] and pfemale[*] have been initialized.
      initialized = 1;

   }

   // Initialization has taken place.  Execute the function call.

   // Point to the correct gender-specific table.
   if (gender == 'f' || gender == 'F') {
      p = pfemale;
   } else {
      p = pmale;
   }

   // No one lives more than 115 years.
   if (t >= 115.0) {
      p0 = 1.0;
   }

   // Otherwise, interpolate between appropriate probabilities.
   else {

      i = (int) t;

      // Now i <= t < i+1, compute the weights.
      w1 = t - i;
      w0 = 1 - w1;

      // Interpolate from the gender specific table as specified above.
      p0 = w0 * p[i] + w1 * p[i+1];

   }

   return (p0);

}


////////////////////////////////////////////////////////////////////////////////
// This function computes the functional inverse of u = F(t).
// It uses interval bisection methodology, which is quite slow.
////////////////////////////////////////////////////////////////////////////////
double HSlow (char gender, double u) {

   int done;
   double t, trial, lower, upper, epsilon=.001;

   lower = 0;
   upper = 115;
   done = 0;

   while (!done) {

      // New trial t.
      trial = (lower + upper) / 2.0;

      // Is t too big???
      if (F(gender, trial) > u) {
         upper = trial;
      }

      // or too small???
      else {
         lower = trial;
      }

      if (upper - lower < 2.0*epsilon) {
         done = 1;
      }

   }

   // The midpoint of [lower,upper] is now within epsilon of the exact t.
   t = (lower + upper) / 2.0;

   return (t);

}


////////////////////////////////////////////////////////////////////////////////
// This is a faster version of Hslow suitable for simulation loops.
// It creates a "look-up" table held in a static array.
////////////////////////////////////////////////////////////////////////////////
double H (char gender, double u){

   // Static variables.
   static int initialized=0;
   static double *Hmale, *Hfemale;

   // Non-static variables.
   double t, u0, w0, w1, *H;
   int i;

   // On the first time through:
   // Create the look-up table for values of u = i / 10000, where 0 <= i <= 10000.
   // Use the slow version HSlow (gender, u) to do this.
   // Create tables for both male and female lives.
   if (!initialized) {
      Hmale   = (double *) calloc (10001, sizeof (double));
      Hfemale = (double *) calloc (10001, sizeof (double));
      for (i = 0; i <= 10000; i++) {
         u0 = i / 10000.0;
         Hmale[i] = HSlow ('m', u0);
         Hfemale[i] = HSlow ('f', u0);
      }
      initialized = 1;
   }

   // Initialization is now complete.  Execute the function call.

   // Point to the correct gender-specific table.
   if (gender == 'f' || gender == 'F') {
      H = Hfemale;
   } else {
      H = Hmale;
   }

   // Determine i such that i/10000 <= u < (i+1)/10000
   i = (int) (10000*u);

   // Will have H[i] <= x < H[i+1]; interpolate.

   // Compute the weights.
   u0 = i / 10000.0;
   w1 = 10000.0 * (u - u0);
   w0 = 1 - w1;

   // Interpolate.
   t = w0 * H[i] + w1 * H[i+1];

   return (t);

}


////////////////////////////////////////////////////////////////////////////////
// This function compute the functional inverse of
//    P[lifespan <= x | lifespan > a].
//    then subtracts off "a" to get the remaining life R.
////////////////////////////////////////////////////////////////////////////////
double HConditional (char gender, double a, double U) {

   double p_a, R;

   p_a = F (gender, a);
   R   = H (gender, p_a + (1.0-p_a) * U) - a;

   return (R);

}


