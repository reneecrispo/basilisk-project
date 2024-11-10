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
#define CASE2

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
    # define LEVEL 7
#endif

#if MOMENTUM
    q.t[right] = dirichlet(0);
    q.t[left]  = dirichlet(0);
#else
    u.t[right] = dirichlet(0);
    u.t[left]  = dirichlet(0);
#endif
   
double rho_bubble1 = 1.0, mu_bubble1 = 0.1;   // first bubble
double rho_bubble2 = 2.0, mu_bubble2 = 0.05;  // second bubble

scalar f1[], f2[]; 

int main() {
    size (2 [1]);  
    DT = 1. [0, 1];
    init_grid (1 << LEVEL);

    #ifdef CASE2
        f.sigma = 1.96;
    #else
        f.sigma = 24.5;
    #endif

    TOLERANCE = 1e-4 [*];
    #if REDUCED
        G.x = -0.98;
        Z.x = 1.;
    #endif

    run();

    return 0;
}

event init (t = 0) {
    mask (y > 0.5 ? top : none);

    #if LEVELSET
        foreach() {
            d[] = sqrt (sq(x - 0.35) + sq(y)) - 0.25; 
            d[] = min(d[], sqrt (sq(x - 0.65) + sq(y)) - 0.25); 
        }
    #else
        scalar f1[], f2[];
        fraction (f1, sq(x - 0.35) + sq(y) - sq(0.25)); 
        fraction (f2, sq(x - 0.65) + sq(y) - sq(0.25)); 
        foreach()
            f[] = max(f1[], f2[]); 
    #endif
}


#if !REDUCED
event acceleration (i++) {
    face vector av = a;
    foreach_face(x)
        av.x[] -= 0.98;
}
#endif

void mg_print (mgstats mg) {
    if (mg.i > 0 && mg.resa > 0.)
        printf ("%d %g %g %g %d ", mg.i, mg.resb, mg.resa,
            mg.resb > 0 ? exp (log (mg.resb/mg.resa)/mg.i) : 0.,
            mg.nrelax);
}


/* At $t=3$ we output the shape of the bubble. */
event interface (t = 2.5) {
    output_facets (f, stderr);
}

#ifdef ADAPT
event adapt (i++) {
    double uemax = 1e-2;
    adapt_wavelet ({f, u}, (double[]){0.0001, uemax, uemax}, LEVEL);
}
#endif

event movie (i++) {
    clear();
    
    // Impostazioni per la visualizzazione orizzontale
    view (tx = -0.5, ty = 0, width = 4800, height = 2400);

    squares ("f1", linear = true, min = 0, max = 1); 
    draw_vof("f1", linear = true);  

    squares ("f2", linear = true, min = 0, max = 1); 
    draw_vof ("f2", lw = 1.5);

    mirror ({0, 1}) {
        draw_vof ("f1", lw = 1.5);
        draw_vof ("f2", lw = 1.5);
        cells();
    }

    save ("movie.mp4");
}



/**
## Results

The final shape of the bubble is compared to that obtained with the
MooNMD Lagrangian solver (see [the FEATFLOW
page](http://www.featflow.de/en/benchmarks/cfdbenchmarking/bubble/bubble_verification.html))
at the highest resolution. We also display the shape of the
axisymmetric version of the test. The axisymmetric bubble moves much
faster.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 1.
set term push
set term @SVG size 640,320
set size ratio -1
set grid
plot [][0:0.4]'c1g3l4s.txt' u 2:($1-0.5) w l t 'MooNMD', \
              'rising/log' u 1:2 w l t 'Basilisk', \
              'rising-levelset' u 1:2 w l t 'Basilisk (levelset)', \
              'rising-clsvof' u 1:2 w l t 'Basilisk (CLSVOF)', \
              'rising-axi' u 1:2 w l t 'Basilisk (axisymmetric)', \
              'rising-axi-clsvof' u 1:2 w l t 'Basilisk (axi + CLSVOF)', \
              'rising-axi-momentum' u 1:2 w l t 'Basilisk (axi + momentum)'
~~~

For test case 2, the mesh in Basilisk is too coarse to accurately
resolve the skirt.

~~~gnuplot Bubble shapes at the final time ($t=3$) for test case 2.
set key bottom left
plot [][0:0.4]'../c2g3l4s.txt' u 2:($1-0.5) w l t 'MooNMD', \
              '../rising2/log' u 1:2 w l t 'Basilisk', \
              '../rising2-levelset/log' u 1:2 w l t 'Basilisk (levelset)', \
              '../rising2-clsvof/log' u 1:2 w l t 'Basilisk (CLSVOF)'
~~~

The agreement for the bubble rise velocity with time is also good.

~~~gnuplot Rise velocity as a function of time for test case 1.
set term pop
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'../c1g3l4.txt' u 1:5 w l t 'MooNMD', \
              'out' u 1:5 w l t 'Basilisk', \
              '../rising-levelset/out' u 1:5 w l t 'Basilisk (levelset)', \
              '../rising-clsvof/out' u 1:5 w l t 'Basilisk (CLSVOF)',     \
              '../rising-axi/out' u 1:5 w l t 'Basilisk (axisymmetric)',  \
              '../rising-axi-clsvof/out' u 1:5 w l t 'Basilisk (axi + CLSVOF)',  \
              '../rising-axi-momentum/out' u 1:5 w l t 'Basilisk (axi + momentum)'
~~~

~~~gnuplot Relative volume difference as a function of time for test case 1.
reset
set grid
set xlabel 'Time'
set ylabel '(vb - vb_0)/vb_0'
set key bottom left
plot [0:3]'out' u 1:2 w l t 'Basilisk', \
          '../rising-levelset/out' u 1:2 w l t 'Basilisk (levelset)',   \
	  '../rising-clsvof/out' u 1:2 w l t 'Basilisk (CLSVOF)',	\
	  '../rising-axi/out' u 1:2 w l t 'Basilisk (axisymmetric)',	\
	  '../rising-axi-clsvof/out' u 1:2 w l t 'Basilisk (axi + CLSVOF)',	\
	  '../rising-axi-momentum/out' u 1:2 w l t 'Basilisk (axi + momentum)'
~~~

~~~gnuplot Rise velocity as a function of time for test case 2.
reset
set grid
set xlabel 'Time'
set key bottom right
plot [0:3][0:]'../c2g3l4.txt' u 1:5 w l t 'MooNMD', \
              '../rising2/out' u 1:5 w l t 'Basilisk', \
              '../rising2-levelset/out' u 1:5 w l t 'Basilisk (levelset)', \
              '../rising2-clsvof/out' u 1:5 w l t 'Basilisk (CLSVOF)'
~~~

~~~gnuplot Relative volume difference as a function of time for test case 2.
reset
set grid
set xlabel 'Time'
set ylabel '(vb - vb_0)/vb_0'
set key top left
plot [0:3]'../rising2/out' u 1:2 w l t 'Basilisk',		       \
          '../rising2-levelset/out' u 1:2 w l t 'Basilisk (levelset)', \
	  '../rising2-clsvof/out' u 1:2 w l t 'Basilisk (CLSVOF)'
~~~
*/