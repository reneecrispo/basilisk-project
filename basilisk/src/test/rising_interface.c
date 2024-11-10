/*
# Rising bubble

A two-dimensional bubble is released in a rectangular box and raises
under the influence of buoyancy. This test case was proposed by
[Hysing et al, 2009](/src/references.bib#hysing2009) (see also [the
FEATFLOW page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble.html)).

We solve the incompressible, variable-density, Navier--Stokes
equations with interfaces and surface tension. We can solve either the
axisymmetric or planar version. We can used standard or "reduced"
gravity. We also test levelset interface tracking and a momentum
formulation. 
*/

#define ADAPT
// #define CASE2

#include "view.h"
#include "vof.h"
#include "draw.h"
#include <math.h>
#include <stdio.h>

#if AXIS
    # include "axi.h" // fixme: does not run with -catch
#endif
#if MOMENTUM
    # include "momentum.h"
#else
    #include "navier-stokes/centered.h"
#if CLSVOF
    # include "two-phase-clsvof.h"
#elif LEVELSET
    # include "two-phase-levelset.h"
#else
    # include "two-phase.h"
#endif
#endif
#if LEVELSET
    # include "integral.h"
#else
    # include "tension.h"
#endif
#if REDUCED
    # include "reduced.h"
#endif

#ifndef LEVEL
    # define LEVEL 10
#endif


/**
The boundary conditions are slip lateral walls (the default) and
no-slip on the right and left walls. */

#if MOMENTUM
    q.t[right] = dirichlet(0);
    q.t[left]  = dirichlet(0);
#else
    u.t[right] = dirichlet(0);
    u.t[left]  = dirichlet(0);
#endif

// Define the adimentional parameters 
double U; // Characteristic velocity
double L; // Characteristic length
double gr = 9.81; // Acceleration due to gravity

// Definition of skirt position and skirt width as global variables
double skirt_position = 0.0; 
double skirt_width = 0.0; 

// Definitin of the variables to calculate the characteristic velocity
double sum = 0.0;
int n_elem_vel = 0;

double t_f = 2.; // final time

int main() {

    size (2 [1]);
    DT = 1. [0, 1];
    init_grid (1 << LEVEL);

    rho1 = 1000., mu1 = 10.; 

    #ifdef CASE2 
        rho2 = 1., mu2 = 0.1;
    #else
        rho2 = 100., mu2 = 1.;  
    #endif 

    #if LEVELSET
        #ifdef CASE2
            const scalar sigma[] = 1.96;
        #else
            const scalar sigma[] = 24.5;
        #endif
        d.sigmaf = sigma;
    #else // !LEVELSET
        #ifdef CASE2
            f.sigma = 1.96;
        #else
            f.sigma = 24.5;
        #endif
    #endif // !LEVELSET

    TOLERANCE = 1e-4 [*];
    #if REDUCED
        G.x = -0.98;
        Z.x = 1.;
    #endif

    run();

    U = sum / n_elem_vel;
    L = 0.5; // diameter = 2 * radius = 2 * 0.25

    double Re = (rho1 * U * L) / mu1; // Reynolds number
    double Bo = (rho1 * gr * sq(L)) / f.sigma; // Bond number

    // Create and open the file to write
    FILE *dati_test = fopen("dati_test.txt", "w");
    if (dati_test == NULL) {
        perror("Error during opening of the file dati_text.txt\n");
        return 1;
    }

    #ifdef CASE2
        fprintf(dati_test, "Executing test case 2: rho1 = %g, mu1 = %g, rho2 = %g, mu2 = %g, sigma = %g, LEVEL = %g\n", rho1, mu1, rho2, mu2, f.sigma, LEVEL);
    #else 
        fprintf(dati_test, "Executing test case 1: rho1 = %g, mu1 = %g, rho2 = %g, mu2 = %g, sigma = %g, LEVEL = %g\n", rho1, mu1, rho2, mu2, f.sigma, LEVEL);
    #endif

    #ifdef ADAPT
        fprintf(dati_test, "Adaptive mesh ON\n\n");
    #endif

    // Output parameters and calculated numbers
    fprintf(dati_test, "Characteristic velocity \t U = %g\n", U);
    fprintf(dati_test, "Characteristic length \t L = %g\n", L);
    fprintf(dati_test, "Reynolds number \t Re = %g\n", Re);
    fprintf(dati_test, "Bond number \t Bo = %g\n", Bo);

    fclose(dati_test);

    return 0;
}

event init (t = 0) {
    /*
    The domain is a rectangle. We only simulate half the bubble. */
    mask (y > 0.5 ? top : none);

    /*
    The bubble is centered on (0.5,0) and has a radius of 0.25. */
    #if LEVELSET
    foreach()
        d[] = sqrt (sq(x - 0.5) + sq(y)) - 0.25;
    #else
    fraction (f, sq(x - 0.5) + sq(y) - sq(0.25));
    #endif
}

/*
We add the acceleration of gravity. */

#if !REDUCED
event acceleration (i++) {
    face vector av = a;
    foreach_face(x)
        av.x[] -= 0.98;
}
#endif

/*
A utility function to check the convergence of the multigrid
solvers. */

void mg_print (mgstats mg) {
    if (mg.i > 0 && mg.resa > 0.)
        printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
            mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
            mg.nrelax);
}

/*
We log the position of the center of mass of the bubble, its velocity
and volume as well as convergence statistics for the multigrid
solvers. */

event logfile (i++) {
    double xb = 0., vb = 0., sb = 0.;
    foreach(reduction(+:xb) reduction(+:vb) reduction(+:sb)) {
        double dv = (1. - f[])*dv();
    #if MOMENTUM
        vb += q.x[]*dv/rho(f[]);
    #else
        vb += u.x[]*dv;
    #endif
        xb += x*dv;
        sb += dv;
    }
    static double sb0 = 0.;
    if (i == 0) {
        printf ("t sb -1 xb vb dt perf.t perf.speed\n");
        sb0 = sb;
    }
    printf ("%g %g %g %g %g %g %g %g ", t, (sb - sb0)/sb0, -1., xb/sb, vb/sb, dt, perf.t, perf.speed);

    FILE *velocity_log = fopen("velocity_log.txt", "a");
    if (velocity_log == NULL) {
        perror("Error during opening of the file velocity_log.txt\n");
        return 1;
    }

    fprintf(velocity_log, "%g %g\n", t, vb/sb);
    if(t >= t_f - t_f * 0.05) {
        sum += vb/sb;
        n_elem_vel++;
    }

    fclose(velocity_log);

    

    #if !MOMENTUM
        mg_print (mgp);
        mg_print (mgpf);
        mg_print (mgu);
    #endif
    putchar ('\n');
    fflush (stdout);
}

/*
At $t=3$ we output the shape of the bubble. */

event interface (t = t_f) {
    output_facets (f, stderr);
}

#ifdef ADAPT
event adapt (i++) {
  double uemax=1e-2;
  adapt_wavelet ({f,u}, (double[]){0.0001,uemax,uemax},LEVEL);
}
#endif

/*
event image (i++) {
  char filename[256];
  sprintf (filename, "out%04d.ppm", i);
  output_ppm (f, file=filename);
}
*/


#ifdef CASE2
event skirt_width_log (i++) {
    double max_y = 0; 
    double max_x = 0;

    foreach() { // for each cell
        if (f[] > 0 && f[] < 1) { // cells at the interface
            if (y > max_y) {
                max_y = y; 
                max_x = x;
            }
        }
    }

    skirt_width = 2 * max_y;
    skirt_position = max_x;

    FILE *skirt_file = fopen("skirt_width_log.txt", "a"); 
    if (skirt_file != NULL) {
        fprintf(skirt_file, "%g %g %g\n", t, skirt_width, skirt_position);
        fclose(skirt_file);
    } else {
        perror("Error during the opening of the file skirt_width_log.txt");
    }
}
#endif


event interface (t = t_f) {
    char filename[100];
    snprintf(filename, sizeof(filename), "interface_%g.dat", t);
    FILE *interface_file = fopen(filename, "a");

    if (interface_file != NULL) {
        output_facets(f, interface_file);
        fclose(interface_file);
    } else {
        perror("Error opening file to save interface\n");
    }
}


/*
event movie (i++) {
    clear();
    
    // orizontal view
    view (tx = -0.5, ty = 0, width = 4800, height = 2400);

    squares ("f", linear = true); 
    draw_vof ("f", lw = 1.5);     

    mirror ({0, 1}) {
        draw_vof ("f", lw = 1.5);
        cells();
    }
    save ("movie.mp4");
}
*/