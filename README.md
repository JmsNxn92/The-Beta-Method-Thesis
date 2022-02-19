# The-Beta-Method-Thesis
The final program for the beta method of iterated exponentials


To begin use of the program one must execute the file in a pari-gp shell.
We can now choose digit precision by taking "\p d" and series precision by "\ps count"; for d,count positive integers.


From here, one must call init(b,l); or init_OFF(b,l) if b > 1E6.

This produces 4 function protocols (and a few variations):

beta(z) = exp(b * beta(z-1))/(1+exp(-l*(z-1)))

tau(z) = log(1+tau(z+1)/beta(z+1))/base - log(1+exp(-l*z))

Sexp_N(1+z) = beta(z+1) + tau(z+1) = exp(b * (beta(z) + tau(z)))


Andc lastly, a protocol:

tau_rho(z+CENTER) = tau(z+CENTER,LIMIT);

which requires an initializatio:

init_rho(CENTER)

This function protocol is less accurate but serves to exposite the mathematics.
