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
    # define LEVEL 6
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

/*
    // Create and open the file to read
    FILE *velocity_log = fopen("velocity_log.txt", "r");
    if (velocity_log == NULL) {
        perror("Error during opening of the file velocity_log.txt\n");
        return 1;
    }

    // Characteristic values
    double *velocities = NULL;
    int count = 0;
    double velocity;

    // Scorri il file e leggi tutte le velocità
    while (fscanf(velocity_log, "%*lf %lf", &velocity) == 1) {
        velocities = realloc(velocities, (count + 1) * sizeof(double));
        velocities[count++] = velocity;
    }
    fclose(velocity_log);

    // Calcola l'indice di inizio dell'ultimo 5% dei dati
    int start_index = count - (int)(0.05 * count);
    double sum = 0.0;
    int num_elements = count - start_index;

    // Somma le velocità nell'ultimo 5%
    for (int i = start_index; i < count; i++) {
        sum += velocities[i];
    }
    
    U = sum / num_elements;

    free(velocities);

    fclose(velocity_log);
*/

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
        fprintf(dati_test, "Executing test case 2: rho1 = %g, mu1 = %g, rho2 = %g, mu2 = %g, sigma = %g\n", rho1, mu1, rho2, mu2, f.sigma);
    #else 
        fprintf(dati_test, "Executing test case 1: rho1 = %g, mu1 = %g, rho2 = %g, mu2 = %g, sigma = %g\n", rho1, mu1, rho2, mu2, f.sigma);
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
    if(t >= 2.3) {
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

event interface (t = 2.5) {
    output_facets (f, stderr);
}

#ifdef ADAPT
event adapt (i++) {
    adapt_wavelet({f,u}, (double[]){5e-4,1e-3,1e-3}, LEVEL);
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
        fprintf(skirt_file, "%g %g %g\n", t, skirt_width, skirt_position); // Tempo, larghezza, posizione
        fclose(skirt_file);
    } else {
        perror("Error during the opening of the file skirt_width_log.txt");
    }
}
#endif


event movie (i++) {
    clear();
    
    // orizontal view
    view (tx = -0.5, ty = 0, width = 4800, height = 2400);

    squares ("f", linear = true); 
    draw_vof ("f", lw = 1.5);     

    #ifdef CASE2
        draw_line(skirt_position, -1.0, 1.0, 0.02, "red");
    #endif

    mirror ({0, 1}) {
        draw_vof ("f", lw = 1.5);
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