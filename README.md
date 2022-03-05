# The-Beta-Method-Thesis
The final program for the beta method of iterated exponentials

Coded by James David Nixon, with supplementary code by Sheldon Levenstein and Mike3 of the tetration forum.
Levenstein optimized some protocols I had written, and we use Mike3's graphing protocol.


To begin use of the program one must execute the file in a pari-gp shell.
We can now choose digit precision by taking "\p d" and series precision by "\ps count"; for d,count positive integers.


From here, one must call init(l,b); or init_OFF(l,b) if abs(b) > 1E6. The variable l must have a real part greater than zero; and b cannot equal zero (this produces the trivial beta).

After initialization, the values l = mult and b = base are global variables. It is possible to initialize these objects as polynomials as well; meaning mult = 1+u or base = 1+v are perfectly valid. This protocol internally uses the variable x, so do not use this variable or reassign this variable else where in the program.

This produces 4 function protocols (and a few variations):

beta(z) = exp(base * beta(z-1))/(1+exp(-mult * (z-1)))

tau(z) = log(1+tau(z+1)/beta(z+1))/base - log(1+exp(-mult * z))/base

Sexp_N(z+1) = beta(z+1) + tau(z+1) = exp(base * (beta(z) + tau(z)))

Each of these functions has a period of 2 * Pi * I/mult. They will be accurate to the digit precision, usually; if not, re run the initialization protocol with a higher series precision. Or one can keep the series precision and run init(l,b,COUNT) (resp. init_OFF); where COUNT determines how many iterations to do. One should never have to exceed 1000. There is an optional parameter in tau and Sexp_N which can tell us how many times to do the iteration; the program by default chooses enough for most calculations; sometimes it may be beneficial to induce more iterations.

We can treat these three functions as polynomials as well; so that beta(1+z), tau(1+z), Sexp_N(1+z), are all viable prompts, which produce the polynomial about 1 in z to the desired series precision.

There is also a separate initialization protocol you can run after the first which gives the summative form of the error. This would be:

init_rho(c, {LIMIT = default(seriesprecision)}) which creates up to the series precision a sequence of functions rho, which act as sum(j=1,LIMIT,RHO[j]) = tau(q+CENTER,LIMIT-2). The variable CENTER=c; and this is declared as a global constant. This protocol internally uses the variable q; so it is best warned to not use the variable q elsewhere in the program. 

This spawns the function:

rho(z,j)

which is an evaluation of a polynomial centered about CENTER upto the series precision. The variable j indexes up until how many rho terms were made using LIMIT in init_rho. Or, if there's an overflow it sets the protocol to zero; or if j is outside of the index list, sets the protocol to zero. This satisfies:

rho(z,0) = 0

rho(z,1) = -log(1+exp(-mult * z))/base

sum(j=0,LIMIT-1, rho(z,j)) = tau(z, LIMIT-2)

rho(z,j) = tau(z,j-1) - tau(z,j-2);

And lastly, a protocol:

tau_rho(z) = tau(z,LIMIT-2);

Which is just short hand for the sum.


/*****************************************************************************************************/
The next functions are a Julia set protocol for beta. 

We want to determine the points that are in the weak Julia set, and the points which are not. 
These are points where the orbits escape. A point a belongs to the weak Julia set if, for any D> 0:

limsup n to infinity sup(|z-a|<D, 1/beta(z+n)) = infinity

Sadly, there is no efficient way to calculate this, other than the below expression; which tests if 1/beta(z+depth) grows eggregiously fast in a neighborhood.
This is allowed because where beta(z) is very small, nearby we are arbitrarily close to infinity; whether in further orbits or in the same orbit.


This is the Julia test, which requires a single for-loop; which is accurate enough.

The variable depth must come in the form of a nonzero natural number; it tells us how deep to run the recursion.

The variable rad can be chosen as small as possible; and decides the radius of the circle we sample from; 1E-100 gives a good enough result.

The variable sample must come in the form of a nonzero natural number. It determines how many samples we take from the circle of radius rad about a.


This function takes sample amount of sample points on a disk about a; and tests whether we diverge at depth or not.
The ideal form of this function would be at rad to 0, and sample to infinity; and depth to infinity.
But, nonetheless this is still fairly accurate; and as accurate as possible without unscrupuously shrinking rad and growing the sample; plus, making depth as deep as possible.

Is_Julia(a,{depth = 25},{rad = 1E-100},{sample = 24})


Attached to this function is the Is_Shell_Thron, and Is_Shell_Thron_EXP. This takes a point a and tests upto a depth and a series precision the fixed point (or lack there of) of the exponential exp(log(a) * z) (resp. exp(a * z)); and spits out the constant or a zero for lack of fixed point. This is solely for the purpose of graphing and testing for if a point is in Shell-Thron. 

/*****************************************************************************************************/

The next functions are a tad more involved. We want to find the normalization constant, and this can be difficult.
Expect this normalization process to fail for values like b=0.0001 or b = -1000; or any extreme values.
At this point, it is much better to just run the super exponential that is not normalized.

We effectively run a Newtonian root finder algorithm, but we run it twice; in the neighborhood of z = 0
I have included these functions in the normalization initialization; but they may not initialize properly for b = 0.0001, or  b=-1000.

It can also encounter a good amount of problems with complex b, but it will work; the problem is typically wait times.


init_NORM(l,b)  where l and b must be mult and base or else we fail. Then we have access to:

x_NORM()

such that

Sexp(z) = Sexp_N(z+x_NORM())

in which Sexp(0) = 1; now. Again, this is superfluous and not as well coded as the rest of this program; expect wait times and imprecision. 

/************************************************************************************************************************************

The following is Mike3's readme for his graphing protocol, used to make phase-plot graphs.


/* =============================================================================  
** Color graphing system - mike3 - 20.11.10 04:07
** Hi.
** I thought I'd post the code I use to generate the color graphs from Pari/GP.
** Here it is.
**
** Note: the output is in .PPM format. 
** You'll need something else to convert that to .PNG. (I use GIMP.)
** 
** Also, I might warn you: it takes a LONG time to graph a complex function
** with a significantly complicated calculation procedure, as it must be
** evaluated at every pixel of the graph.
** 
** (updated 12/16/2010 -- program was not all good: 
**   * spurious "func" parameter in MakeGraph and "safetyarg" was missing.)
** ------------------------------------------------------------------------------------------ 
**  
** ============================================================================= */

/* =============================== Code: ==================================================== */
/* Complex function magnitude/phase plotter. */

/* To use:
*     1. Define function to graph as func(z).
*     2. Load this program.
*     3. Execute MakeGraph(width, height, x0, y0, x1, y1, filename) with the parameters given as follows:
*        width, height = width/height of image in pixels
*        x0, y0, x1, y1 = rectangle of complex plane to graph: x0 + y0i in upper-left corner to x1 + y1i in lower-right corner
*        filename = name of file to save as.
* Output is in .PPM format.
*/
