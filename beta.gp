/***************************************************************************

James David Nixon's beta.gp program with supplemetary code by Sheldon Levenstein and Mike3 of the tetration forum.

This code serves as the empirical justification of the beta method thesis on tetration.
***************************************************************************/


/*
This first function is an initialization function. 

The value b = base, is the base of the superexponential;
with the convention that exp(b*z) is the exponential. Therein, our only bad value is b=0.

The value l = mult, is the multiplier of the beta function we choose.
Recall, that 2*Pi*I/mult is the period of the beta function.

count is telling us how deep we want to do the recursion.
By default I have set it to the seriesprecision you choose with \ps; where at least 240 iterations must be done.


Before any of the other functions can be run, we must run init (or we must run the initialization init_OFF).
Which requires a choice of mult and base given by the first and second argument, respectively.

This code is a hybridization between my code and Sheldon's code.
I had originally initialized with the multiplier as a free variable which had its benefits (graphing the multiplier, particularly).
Sheldon's method is significantly faster by requiring us to declare the multiplier; and a quick speed up with subst.
*/


offset = 0; /*this is a global variable that is unused in this protocol unless we call call large values of b (approx abs(b) = 1E6); 
				at which; it prevents a good amount of overflows, and mostly just allows us to still grab a taylor series in the intialization process*/


init(l,b,{count=default(seriesprecision)}) = {
	base = b;
	mult = l;  /*declare these variables to work for all other functions; make them global*/
	
	if(count < 240  && default(seriesprecision) < 240,count=240);  /*This just makes sure we do at least 240 iterations regardless of series depth*/
	
	Cut = 0;  /*this is just a quick cut value, to limit how deep the recursion goes, otherwise we get deep recursion*/

	if(abs(real(mult)) < 1E-3, Cut = offset-1E3, Cut = offset-50 - 1/real(mult));

	beta_taylor = x/exp(mult) + O(x^2);  /*Sheldon is to credit for this polynomial translation*/
	for (n=1,count, 
		beta_taylor=exp(base*subst(beta_taylor,x,x/exp(mult))) / (1 + exp(mult)/x) /*this is a degree faster than my old code; absolutely genius work by Sheldon*/
	); 
	beta_taylor=Pol(beta_taylor); /*we treat these objects as polynomials*/

}

/*
This is a modified initialization process intended for values of b > 1E4 or about; it is slightly slower but more accurate for large values.
*/

init_OFF(l,b,{count=default(seriesprecision)}) = {
	base = b;
	mult = l;  /*declare these variables to work for all other functions; make them global*/
	
	if(count < 240  && default(seriesprecision) < 240,count=240);  /*This just makes sure we do at least 240 iterations regardless of series depth*/
	
		
	Cut = 0;  /*this is just a quick cut value, to limit how deep the recursion goes, otherwise we get deep recursion*/

	if(abs(real(mult)) < 1E-3, Cut = offset-1E3, Cut = offset-50 - 1/real(mult));
	
	/*This variable offsets the nested exponentials so that the iterate will converge regardless how large b is; or how close it is to zero*/
	if(abs(b) > 1E6, offset = log(sqrt(abs(b)))^6);
	if(abs(b) <1E-6, offset = log(sqrt(1/abs(b)))^6);

	beta_taylor = x/exp(mult+offset) + O(x^2)/exp(offset);  /*Sheldon is to credit for this polynomial translation*/
	for (n=1,count, 
		beta_taylor=exp(base*subst(beta_taylor,x,x/exp(mult+offset))) / (1 + exp(mult+offset)/x) /*this is a degree faster than my old code; absolutely genius work by Sheldon*/
	); 
	beta_taylor=Pol(beta_taylor); /*we treat these objects as polynomials*/
}

/*
This truncates a series into the constant coefficient--needed for the recursion
We will use this to generate if conditions.
*/

Const(poly) = {
    my(vars = variables(poly));
    substvec(poly, vars, vector(#vars));
};

/*
A hybridization between mine and Sheldon's code. Sheldon optimized the exponential series defining beta. 
But he didn't account for loss of accuracy in the functional equation.
But otherwise, mine and Sheldon's code are mathematically identical.

This function will satisfy the functional equation:

beta(z+1) = exp(base*beta(z))/(exp(-mult*z) + 1)

It will have singularities at the points mult*(z-j) = (2k+1)*Pi*I for j >= 1 and j,k integers.

It also has a period of 2*Pi*I/mult
*/



beta(z) = {
	
	/*
	The goal of this function: if we're in the radius of convergence of the exponential series (at an accurate level), just substitute the value.
	Otherwise, we pull backwards until we are in the radius of convergence,
	and iterate forwards with the functional equation.
	*/
	
	/*
	The accurate radius of convergence (the radius where we are still sufficiently accurate) is about exp(-50 - 1/real(mult)) for most reasonable values.
	If mult is too small though, this can create recursion errors (Pari doesn't like 1E4 recursive calls)
	You may need to increase depth of iteration/series precision/digit precision for more anomalous values.
	*/
	
	if(real(Const(z)) <= Cut, 
		subst(beta_taylor,x,exp(mult*z)),
		exp(base*beta(z-1))/(1+exp(-mult*(z-1)))
	);
}

/*
This function is specifically for values of base greater than 1E6 or less than 1E-6; and satisfies a different functional equation than beta:
	
exp(base*beta_off(z))/(1+exp(mult*(offset-z))) = beta_off(z+1)

This effectively makes a shift to the argument so that we can still grab the Taylor series. It is related to beta by the relation:

beta_off(z+offset) = beta(z)
*/


beta_off(z) = {
	/*
	The goal of this function: if we're in the radius of convergence of the exponential series (at an accurate level), just substitute the value.
	Otherwise, we pull backwards until we are in the radius of convergence,
	and iterate forwards with the functional equation.
	*/
	
	/*
	The accurate radius of convergence (the radius where we are still sufficiently accurate) is about exp(-50 - 1/real(mult)) for most reasonable values.
	If mult is too small though, this can create recursion errors (Pari doesn't like 1E4 recursive calls)
	You may need to increase depth of iteration/series precision/digit precision for more anomalous values.
	*/

	
	if(real(Const(z)) <= Cut, 
		subst(beta_taylor,x,exp(mult*(z))),
		exp(base*beta_off(z-1))/(1+exp(-mult*(z-offset-1)))
	);
}

if(offset != 0, beta(z) = beta_off(z+offset));  /*if we've initialized init_OFF, instead of init; we'll renormalize beta*/


beta_N(z) = {beta(z-1/mult);} /*This is a rough normalization so that abs(beta_N(0)) < 1; largely unused but helpful*/



/*
This function is the error between tetration and our beta function.

It satisfies the equation beta(z+1) + tau(z+1) = exp(base*(beta(z)+tau(z)))

The value count chooses the depth of iteration. The value 100 is a heuristic value.
At certain times, for certain bases, and certain values, it can be helpful to set count much much larger.

For instance, if we are with base = 1 and mult = 1; and we're looking at z approx -500, we should set count to 600 to be safe.
Technically this value should be infinite, but setting it arbitrary large can increase wait times for bounded bases.
*/

tau(z, {count = 100}) = {
	if(abs(Const(beta(z))) <= 1E6 && count > 0, /*1E6 is chosen because exp(1E6) is about where pari overflows*/
		count--;
		log(1+tau(z+1,count)/beta(z+1))/base - log(1+exp(-mult*z))/base,
		-log(1+exp(-mult*z))/base
	);
}


/*
This is tau but we are using beta_N rather than beta.
This means it gets us closer to the value x_0, where beta_N(x_0) + tau_N(x_0) =1.

A largely unused protocol.
*/

tau_N(z, {count = 100}) = {
	if(abs(Const(beta_N(z))) <= 1E6 && count > 0, /*1E6 is chosen because exp(1E6) is about where pari overflows*/
		count--;
		log(1+tau_N(z+1,count)/beta_N(z+1))/base - log(1+exp(-mult*(z-1/mult)))/base,
		-log(1+exp(-mult*(z-1/mult)))/base
	);
}

/*
This is an alternative initialization; done through rho, that is intended to grab the polynomials describing rho faster.

If you want to analyse tau(center+z,count) = sum(j=1,count+1, rho(z,j)); then this is how its done. 

This will initialize a sequence of polynomials up to LIMIT about CENTER.

This function protocol is very spurious, and unnecessary. If you try to program this protocol using Levenstein's form of the error; we run far too slow.
And if you try and rewrite this protocol to work solely off of polynomials, things like to get too unaccurate too fast; and we diverge rather than converge.
For that run this initialization only if you want to analyse the various rho terms accurately, and be prepared for a wait.
*/

init_rho(c,LIMIT={default(seriesprecision)}) = {
	CENTER = c;
	local(rho,k,j);
	rho=vector(LIMIT,i,0 + O(x^LIMIT));
	global(RHO);

	rho[2] = -log(1+exp(-mult*(x+CENTER)))/base;

	for(k=3,LIMIT, 
		rho[k] = tau(CENTER+x,k-2)-sum(j=1,k-2,rho[j+1]);
		if(rho[k]==0,k=LIMIT+1);
	);
	
	/*initialize the rho function by calling about a point*/
	RHO=rho;
}

rho(z,K) = {
	subst(Pol(RHO[K+1],x),x,z-CENTER); \\shift the root value because PARI indexes start at 1, and I want it to start at 0.
}

tau_rho(z,LIMIT={default(seriesprecision)}) = {
	sum(j=1,LIMIT-1,rho(z,j));
}


/*********************************

This function does pretty much the same thing as tau. 
But instead of quitting at exp(1E6) it quits if there's an overflow error. 
I don't like this function as much. It can be finicky.

*********************************/

tau_Chop(z,{count=100}) ={
	count--;
	if(count > 0,
		iferr(log(1+tau_Chop(z+1,count)/beta(z+1))/base - log(1+exp(-mult*z))/base, 
			Error,
			-log(1+exp(-mult*z))/base, 
		),
		-log(1+exp(-mult*z))/base
	);
}

/*
This is the first super exponential; it is not normalized, so Sexp_N(0) != 1
The next super exponential is the correct normalized one.

This function satisfies Sexp_N(z+1) = exp(base*Sexp_N(z))

Again, for anomalous values we can overiterate by increasing count.
*/

Sexp_N(z, {count = 100}) = {
	beta(z) + tau(z,count);
}



/***********************************************************************
The next function is a Julia set protocol for beta. 

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
*/



Is_Julia(a,{depth = 25},{rad = 1E-100},{sample = 24}) = {
	local(j);
	
	sample = floor(sample);
	a=a-floor(real(a));

	for(j=0, sample-1,
		iferr(1/beta(a+rad*exp(2*Pi*I*j/sample)+depth),Error,return(1E10),);
	);

	return(0);
}

Is_Shell_Thron(a, {depth = 10^(default(realprecision))}, {d = default(realprecision)}) = {


	if(depth > 200, depth = 200);      \\limit the maximum amount of iterations to 200; despite the desired digit precision.
	if(d>3,d=2.5);                    \\make sure to quit once we have ~4 digits. Increasing this bound means you must increase depth;
	                                  \\ and you must be prepared for long wait times.

	local(out,j, norm_test, quit_while);		        \\declare a: recursive variable out; and a counter of the depth of recursion j; 
	out=1; j=0; norm_test = 0; quit_while = 0;			\\a decision whether out is good enough, norm_test; a break variable quit_while.
														
	 
	
	while(j < depth &&  quit_while==0,                   \\while the iteration is not beyond depth; and if we haven't flagged to quit, do a loop.
	
		iferr(out = exp(log(a)*out), Error,out = 0,);    \\assign out = exp(log(a)*out); if this overflows or errors for some reason; out = 0
		
		if(abs(log(out)-log(a)*out) >= 10^(d/2), out = 0);
		
		if(out == 0 || abs(log(out)-log(a)*out) <=10^(-d), \\ if out is zero; we've produced an error; OR if out = exp(log(a)*out) is precise enough; quit the while loop.
			quit_while = 1;
			norm_test = 1;
		);
		
		j++;
	);
	
	if(norm_test == 1,out,0);  \\if we quit the while loop naturally return out; otherwise flag an indiscriminate by returning 0.
}

Is_Shell_Thron_EXP(a, {depth = 10^(default(realprecision))}, {d = default(realprecision)}) = {


	if(depth > 200, depth = 200);      \\limit the maximum amount of iterations to 200; despite the desired digit precision.
	if(d>3,d=2.5);                    \\make sure to quit once we have ~4 digits. Increasing this bound means you must increase depth;
	                                  \\ and you must be prepared for long wait times.

	local(out,j, norm_test, quit_while);		        \\declare a: recursive variable out; and a counter of the depth of recursion j; 
	out=1; j=0; norm_test = 0; quit_while = 0;			\\a decision whether out is good enough, norm_test; a break variable quit_while.
														
	 
	
	while(j < depth &&  quit_while==0,                   \\while the iteration is not beyond depth; and if we haven't flagged to quit, do a loop.
	
		iferr(out = exp(a*out), Error,out = 0,);    \\assign out = exp(a*out); if this overflows or errors for some reason; out = 0
		
		if(abs(log(out)-a*out) >= 10^(d/2), out = 0);
		
		if(out == 0 || abs(log(out)-a*out) <=10^(-d), \\ if out is zero; we've produced an error; OR if out = exp(a*out) is precise enough; quit the while loop.
			quit_while = 1;
			norm_test = 1;
		);
		
		j++;
	);
	
	if(norm_test == 1,out,0);  \\if we quit the while loop naturally return out; otherwise flag an indiscriminate by returning 0.
}
/********************************************************************************************************************/

/*
The next functions are a tad more involved. We want to find the normalization constant, and this can be difficult.
Expect this normalization process to fail for values like b=0.0001 or b = -1000; or any extreme values.
At this point, it is much better to just run the super exponential that is not normalized.

We effectively run a Newtonian root finder algorithm, but we run it twice; in the neighborhood of z = 0
I have included these functions in the normalization initialization; but they may not initialize properly for b = 0.0001, or  b=-1000.

It can also encounter a good amount of problems with complex b, but it will work; the problem is typically wait times.
*/

/********************************************************************************************************************/


/*
This initializes the normalization constant. Choose the same l=mult and b=base as when you run init
I have caught if you don't, and it will not run.

This essentially just finds the closest zero about z=0 and sets this as our normalizer.

This process will take a while, especially for weird values. It is much slower than init.
So if you don't care about normalizing the equation exp(base*Sexp(z)) = Sexp(z+1), ignore this protocol entirely. 
*/

init_Norm(l,b) = {
	if(l == mult && b == base,
		NORMALIZATION_CONSTANT = x_NORM(mult),
		print("You're normalization failed because the initialization multiplier and base did not correspond to your normalization multiplier and base.");
	);
}

/*
This is the first approximation of the normalization constant.
It's essentially a Newtonian root finder.
*/

x_L({y=mult}) = {
    my(out = 0);
    
    /*
    Do a different method depending on if the multiplier is real or not.
    The real algorithm runs much faster, but doesn't work for complex multipliers.
    Also, it makes sure the function is still real valued
    */
    
    /*
    Also, our first guess of the root is 1 - 1/mult
    */
    
	if(imag(y) == 0, 
		out = polrootsreal(real(Pol(Sexp_N(1-1/mult+x),x))), 
		out = polroots(Pol(Sexp_N(1-real(1/mult)+x),x))
	);

	return(out);
}

/*This is the second approximation of the constant*/

x_L2({y=mult}) = {
    my(out = 0);
    
	if(imag(y) == 0, 
		out = vector(#x_L(y),i ,polrootsreal(real(Pol(Sexp_N(1-1/mult+x_L(y)[i]+v),v)))), 
		out = vector(#x_L(y),i ,polroots(Pol(Sexp_N(1-real(1/mult)+x_L(y)[i]+v),v)))
	);

	return(out);
}

/*
This gives the normalization constant; for a good amount of cases, it fails for a good amount of values though.
Again, if you don't care about normalizing your tetration, ignore all of this.
*/

x_NORM({y=mult}) = {
	vecmin(abs(x_L(y)),&i);
	vecmin(abs(x_L2(y)[i]),&j);
	return(2-x_L(y)[i] - x_L2(y)[i][j]);
}



/*
The final super exponential. This works fairly well for a lot of bases and multipliers; but the normalization process may fail.
In such a case, just run the unnormalized tetration.
*/

Sexp(z,{count = 100}) = {
	Sexp_N(z+NORMALIZATION_CONSTANT,count);
}

/***********************************************************************************************************************/

/*
This code is attributed to the user Ember Edison.
Or at least the philosophy.
*/

/***********************************************************************************************************************/

/*
This is what us tetrationers usually refer to as a safe superexponential.
I'm going to call it a chopped superexponential; partially inspired by Ember Edison's code, which doesn't work, unfortunately.

The difference in this super exponential being it assigns 0 to overflows.
Instead of Sexp_N(500) = overflow, Sexp_N(500) = 0.
I've written it without normalization.
*/

Sexp_Chop(z,{count=100}) = {
    iferr(Sexp_N(z,count),Error, 0,errorname(Error) = "e_OVERFLOW");
}

/* HERE ENDS ALL THE CODE I'VE WRITTEN--the rest was procured from mike3 of the tetration forum*/



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


/* Color conversion (HSB to RGB). */

HSB2RGB(hsb) = {
       local(H=hsb[1]);
       local(S=hsb[2]);
       local(B=hsb[3]);
       local(HH);
       local(F);
       local(P);
       local(Q);
       local(T);

       HH = floor(6*H)%6;
       F = (6*H) - floor(6*H);
       P = B*(1 - S);
       Q = B*(1 - (F*S));
       T = B*(1 - (1-F)*S);
       if(B > 1.0, B = 1.0);
       if(B < 0.0, B = 0.0);
       if(P > 1.0, P = 1.0);
       if(P < 0.0, P = 0.0);
       if(Q > 1.0, Q = 1.0);
       if(Q < 0.0, Q = 0.0);
       if(T > 1.0, T = 1.0);
       if(T < 0.0, T = 0.0);

       if(HH == 0, return([B, T, P]));
       if(HH == 1, return([Q, B, P]));
       if(HH == 2, return([P, B, T]));
       if(HH == 3, return([P, Q, B]));
       if(HH == 4, return([T, P, B]));
       if(HH == 5, return([B, P, Q]));
       }

/* Safe argument. */
safetyarg(z) = if(z == 0, 0, arg(z));

/* Make graph. */
MakeGraph(width, height, x0, y0, x1, y1, filename) = {
       xstep = (x1 - x0)/width;
       ystep = (y1 - y0)/height;
       write(filename, "P3");
       write(filename, "# ", filename);
       write(filename, width, " ", height);
       write(filename, "255");

       for(y=0, height-1,
           for(x=0, width-1,
                  xx = x0+(xstep*x);
                  yy = y0+(ystep*y);
               z = xx+yy*I;
               
               /*
               I've added error catching tools on all the following functions,
               This just ensures we aren't overflowing for very large numbers.
               Every potential overflow is set to zero.
               
               I've noticed that even if you feed Mike's program with func, where func catches errors on its own,
               We can still overflow or glitch out. 
               These error catchers are just a second safety protocol layer, 
               to make sure we can still do large graphs carefree of overflow errors.
               
               In short, if the first error catching failed, this will not.
               
               -James
               */
               
               
               funcvalue = iferr(func(z),error,0,);
               mag = iferr(abs(funcvalue),Error,0,);
               phase = iferr(safetyarg(funcvalue),Error,0,);
               H = iferr(phase/(2*Pi),Error,0,);
               S = iferr(1/(1 + 0.3*log(mag + 1)),Error,1,);
               B = iferr(1 - 1/(1.1 + 5*log(mag + 1)),Error,1-1/1.1,);
               
               
               RGB = HSB2RGB([H, S, B]);
                  Red = floor(RGB[1]*255.0);
                  Green = floor(RGB[2]*255.0);
                  Blue = floor(RGB[3]*255.0);
               write1(filename, Red, " ", Green, " ", Blue, "  ");
              );
           write(filename, "");
       );
       print("Done.");
    }
