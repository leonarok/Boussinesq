#include <stdio.h>

/*
 * Physical quantities
 * -------------------
 * k                      Thermal conductivity     [watt / (meter * kelvin)]
 * rho                    Density                  [kg / meter^3]
 * cp                     Specific heat capacity   [kj / (kg * kelvin)]
 * rho * cp               Volumetric heat capacity [joule / (meter^3 * kelvin)]
 * alpha = k / (rho * cp) Thermal diffusivity      [meter^2 / second]
 *
 * Mercury
 * -------
 * cp = 0.140, rho = 13506, k = 8.69
 * alpha = 8.69 / (0.140 * 13506) = 0.0619 [0.004595841]
 *
 * Copper
 * ------
 * cp = 0.385, rho = 8960, k = 401
 * alpha = 401.0 / (0.385 * 8960) = 0.116 [0.116245369]
 *
 * Tin
 * ---
 * cp = 0.227, k = 67, rho = 7300
 * alpha = 67.0 / (0.227 * 7300) = 0.040 [0.040432080]
 *
 * Aluminum
 * ---------
 * cp = 0.897, rho = 2700, k = 237
 * alpha = 237 / (0.897 * 2700) = 0.098 [0.097857054]
 */

#define MERCURY 0.0619
#define COPPER 0.116
#define TIN 0.040
#define ALUMINIUM 0.098

#define SIZE 256

#define NSTEP 125000
#define CUTOFF 75000

#define BORDER 1

/* Dump rate. 16 is realtime at 25 FPS. */
#define SNAPSHOT 160

/*
 * Indexes the global temperature domain.
 *
 * 1: (int) Buffer row.
 * 2: (int) Buffer column.
 */
#define TEMP(i, j) temperature[(i) * SIZE + (j)]

/* Evaluates valid accesses to local subdomain.
 *
 * 1: (int) Buffer row.
 * 2: (int) Buffer column.
 */
#define BOX(y, x) (                             \
     ((y) >= local_origin[0])                && \
     ((y) < local_origin[0] + local_dims[0]) && \
     ((x) >= local_origin[1])                && \
     ((x) < local_origin[1] + local_dims[1]))   \

/*
 * Indexes the local material domain.
 *
 * 1: (int) Buffer row.
 * 2: (int) Buffer column.
 */
#define LMAT(i, j) local_material[                                    \
        ((i) + BORDER) * (local_dims[1] + 2 * BORDER) + (j) + BORDER] \

/*
 * Indexes the two local temperature subdomains.
 *
 * 1: (int) Step.
 * 2: (int) Buffer row.
 * 3: (int) Buffer column.
 */
#define LTEMP(s, i, j) local_temp[(s) % 2][                       \
    ((i) + BORDER) * (local_dims[1] + 2 * BORDER) + (j) + BORDER] \

/*
 * Solves subdomains by applying the Forward Time Central Space (FTCP) equation.
 *
 * 1: (int) Time step.
 */
void ftcs_solver(int step);

/* Executes Neumann boundary conditions.
 *
 * 1: (int) Time step.
 */
void boundaries(int step);

/*
 * Performs halo exchange.
 *
 * 1: (int) Time step.
 */
void border_exchange(int step);

/*
 * Creates and commits vector types.
 */
void commit_vector_types(void);

/*
 * Applies external temperature to the aluminum bar.
 *
 * 1: (int) Time step.
 */
void external_heat(int step);

/*
 * Initializes the geometry and attributes of the metal blocks.
 */
void configure_geometry(void);

/*
 * Collects local subdomains and dumps the state of the global domain.
 *
 * 1: (int)    Time step.
 * 2: (char *) Output file name.
 */
void collect_area(int step, char *filename);

/*
 * Writes content of an array to file.
 *
 * 1: (FILE *)  File pointer.
 * 2: (float *) Data array.
 */
void write_matrix(FILE *out, float *data);

/*
 * Frees all allocated memory associated with the program.
 */
void free_memory(void);

/*
 * Iterates the global domain NSTEP times.
 */
void iterate(void);
