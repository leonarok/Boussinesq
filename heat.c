#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include <mpi.h>

#include "heat.h"

/* Arrays for the simulation data */
float *temperature;
float *local_temp[2];
float *local_material;

/* Discretization: 5 cm square cells, 2.5 ms time intervals. */
const float h = 5e-2;
const float dt = 2.5e-3;

int size, rank;
int dims[2];
int periods[2] = { 0, 0 };
int coords[2];
int north, south, east, west;
int local_dims[2];
int local_origin[2];

MPI_Comm cart;
MPI_Datatype global_area, local_area;
MPI_Datatype border_row, border_col;

int main(void)
{
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        temperature = (float *) calloc(SIZE * SIZE, sizeof(float));
    }

    MPI_Dims_create(size, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart);
    MPI_Cart_coords(cart, rank, 2, coords);
    MPI_Cart_shift(cart, 0, 1, &north, &south);
    MPI_Cart_shift(cart, 1, 1, &west, &east);

    local_dims[0] = SIZE / dims[0];
    local_dims[1] = SIZE / dims[1];
    local_origin[0] = coords[0] * local_dims[0];
    local_origin[1] = coords[1] * local_dims[1];

    size_t lsize_full =
        (local_dims[0] + 2 * BORDER) * (local_dims[1] + 2 * BORDER);

    local_material = (float *) malloc(lsize_full * sizeof(float));
    local_temp[0] = (float *) malloc(lsize_full * sizeof(float));
    local_temp[1] = (float *) malloc(lsize_full * sizeof(float));

    commit_vector_types();
    configure_geometry();
    iterate();

    free_memory();

    MPI_Finalize();

    return EXIT_SUCCESS;
}

void ftcs_solver(int step)
{
    for (int y = 0; y < local_dims[0]; y++) {
        for (int x = 0; x < local_dims[1]; x++) {
            LTEMP(step + 1, y, x) = LTEMP(step, y, x) + LMAT(y, x) * ((
                LTEMP(step, y - 1, x) + LTEMP(step, y + 1, x) +
                LTEMP(step, y, x - 1) + LTEMP(step, y, x + 1)) -
                4.0 * LTEMP(step, y, x));
        }
    }
}

void boundaries(int step)
{
    int my = local_dims[0] - 1;
    int mx = local_dims[1] - 1;

    /* BEGIN: Calculate inner boundaries. */

    // West boundary exists.
    if (coords[1] == 0) {
        for (int i = 0; i < local_dims[0]; i++) {
            LTEMP(step + 1, i, 0) = LTEMP(step, i, 0) + LMAT(i, 0) * ((
                2 * LTEMP(step, i, 1) + LTEMP(step, i - 1, 0) +
                LTEMP(step, i + 1, 0)) - 4.0 * LTEMP(step, i, 0));
        }
    }

    // East boundary exists.
    if (coords[1] == dims[1] - 1) {
        for (int i = 0; i < local_dims[0]; i++) {
            LTEMP(step + 1, i, mx) = LTEMP(step, i, mx) + LMAT(i, mx) * ((
                2 * LTEMP(step, i, mx - 1) + LTEMP(step, i - 1, mx) +
                LTEMP(step, i + 1, mx)) - 4.0 * LTEMP(step, i, mx));
        }
    }

    // North boundary exists.
    if (coords[0] == 0) {
        for (int i  =0; i < local_dims[1]; i++) {
            LTEMP(step + 1, 0, i) = LTEMP(step, 0, i) + LMAT(0, i) * (
                2 * LTEMP(step, 1, i) + LTEMP(step, 0, i - 1) +
                LTEMP(step, 0, i + 1) - 4.0 * LTEMP(step, 0, i));
        }
    }

    // South boundary exists.
    if (coords[0] == dims[0] - 1) {
        for (int i = 0; i < local_dims[1]; i++) {
            LTEMP(step + 1, my, i) = LTEMP(step, my, i) + LMAT(my, i) * (
                2 * LTEMP(step, my - 1, i) + LTEMP(step, my, i - 1) +
                LTEMP(step, my, i + 1) - 4.0 * LTEMP(step, my, i));
        }
    }

    /* END: Calculate inner boundaries. */

    /* BEGIN: Calculate outer boundaries. */

    if (BOX(0, 0)) {
        LTEMP(step + 1, 0, 0) = LTEMP(step, 0, 0) + LMAT(0, 0) * (
            2 * LTEMP(step, 1, 0) + 2 * LTEMP(step, 0, 1) -
            4.0 * LTEMP(step, 0, 0));
    }

    if (BOX(0, SIZE - 1)) {
        LTEMP(step + 1, 0, mx) = LTEMP(step, 0, mx) + LMAT(0, mx) * (
            2 * LTEMP(step, 1, mx) + 2 * LTEMP(step, 0, mx - 1) -
            4.0 * LTEMP(step, 0, mx));
    }

    if (BOX(SIZE - 1, 0)) {
        LTEMP(step+1,my,0) = LTEMP(step, my, 0) + LMAT(my, 0) * (
            2 * LTEMP(step, my, 1) + 2*LTEMP(step, my - 1, 0) -
            4.0 * LTEMP(step, my, 0));
    }

    if (BOX(SIZE - 1, SIZE - 1)) {
        LTEMP(step + 1, my, mx) = LTEMP(step, my, mx) + LMAT(my, mx) * (
            2 * LTEMP(step, my - 1, mx) + 2 * LTEMP(step, my, mx - 1) -
            4.0 * LTEMP(step, my, mx));
    }

    /* END: Calculate outer boundaries. */
}

void commit_vector_types(void)
{
    MPI_Type_vector(SIZE / dims[0], local_dims[1],
        dims[1] * local_dims[1], MPI_FLOAT, &global_area);

    MPI_Type_vector(SIZE / dims[0], local_dims[1],
        local_dims[1] + 2 * BORDER, MPI_FLOAT, &local_area);

    MPI_Type_vector(BORDER, local_dims[1] + 2 * BORDER,
        local_dims[1] + 2 * BORDER, MPI_FLOAT, &border_row);

    MPI_Type_vector(local_dims[0], BORDER,
        local_dims[1] + 2 * BORDER, MPI_FLOAT, &border_col);

    MPI_Type_commit(&local_area);
    MPI_Type_commit(&global_area);
    MPI_Type_commit(&border_row);
    MPI_Type_commit(&border_col);
}

void border_exchange(int step)
{
    /* Send west. */
    MPI_Sendrecv(
        &LTEMP(step, 0, 0), 1, border_col, west, 0,
        &LTEMP(step, 0, local_dims[1]), 1, border_col, east, 0,
        cart, MPI_STATUS_IGNORE);

    /* Send east. */
    MPI_Sendrecv (
        &LTEMP(step, 0, local_dims[1] - BORDER), 1, border_col, east, 0,
        &LTEMP(step, 0, -BORDER), 1, border_col, west, 0,
        cart, MPI_STATUS_IGNORE);

    /* Send north. */
    MPI_Sendrecv (
        &LTEMP(step, 0, -BORDER), 1, border_row, north, 0,
        &LTEMP(step, local_dims[0], -BORDER), 1, border_row, south, 0,
        cart, MPI_STATUS_IGNORE);

    /* Send south. */
    MPI_Sendrecv (
        &LTEMP(step, local_dims[0] - BORDER, -BORDER), 1, border_row, south, 0,
        &LTEMP(step, -BORDER, -BORDER), 1, border_row, north, 0,
        cart, MPI_STATUS_IGNORE);
}

void iterate(void)
{
    for (int step = 0; step < NSTEP; step++) {
        if (step < CUTOFF) {
            external_heat(step);
        }

        border_exchange(step);
        ftcs_solver(step);
        boundaries(step);

        /* Write snapshot to file. */
        if ((step % SNAPSHOT) == 0) {
            char filename[15];
            sprintf(filename, "data/%.4d.dat", step / SNAPSHOT);
            collect_area(step, filename);
        }
    }
}

void free_memory(void)
{
    if (rank == 0) {
        free (temperature);
    }

    free(local_material);
    free(local_temp[0]);
    free(local_temp[1]);
}

void external_heat(int step)
{
    /* Imposed temperature from outside */
    for (int y = (SIZE / 2) - (SIZE / 16); y <= (SIZE / 2) + (SIZE / 16); y++) {
        for (int x= (SIZE / 4); x <= (3 * SIZE / 4); x++) {
            if (BOX(y, x)) {
                LTEMP(step, y - local_origin[0], x - local_origin[1]) = 100.0;
            }
        }
    }
}

void configure_geometry(void)
{
    /* Fill the pool with mercury. */
    for (int y = 0; y < local_dims[0]; y++) {
        for (int x = 0; x < local_dims[1]; x++) {
            LMAT(y, x) =  MERCURY * dt / (h * h);
            LTEMP(1, y, x) = LTEMP(0, y, x) = 20.0;
        }
    }

    /* Set up the two blocks of copper and tin. */
    for (int y = (SIZE / 8); y < (3 * SIZE / 8); y++) {
        for (int x = (SIZE / 8); x < (SIZE / 2) - (SIZE / 8); x++) {
            if (BOX(y, x)) {
                LMAT(y - local_origin[0], x - local_origin[1]) =
                    COPPER * dt / (h * h);
                LTEMP(0, y - local_origin[0], x - local_origin[1]) = 60.0;
            }
            if (BOX(y, SIZE - x)) {
                LMAT(y - local_origin[0], (SIZE - x) - local_origin[1]) =
                    TIN * dt / (h * h);
                LTEMP(0, y - local_origin[0], (SIZE - x) - local_origin[1]) =
                    60.0;
            }
        }
    }

    /* Set up the block of aluminum in the middle. */
    for (int y = (SIZE / 2) - (SIZE / 16); y <= (SIZE / 2) + (SIZE / 16); y++) {
        for (int x = (SIZE / 4); x <= (3 * SIZE / 4); x++) {
            if (BOX(y, x)) {
                LMAT(y - local_origin[0], x - local_origin[1]) =
                    ALUMINIUM * dt / (h * h);
            }
        }
    }
}

void collect_area(int step, char *filename)
{
    MPI_Request req;
    MPI_Isend(&LTEMP((step % SNAPSHOT), 0, 0), 1, local_area, 0, 0, cart, &req);

    if (rank == 0) {
        int c[2];

        for (int r = 0; r < size; r++) {
            MPI_Cart_coords (cart, r, 2, c);
            MPI_Recv(&TEMP(c[0] * local_dims[0], c[1] * local_dims[1]),
                1, global_area, r, 0, cart, MPI_STATUS_IGNORE);
        }

        FILE *out = fopen(filename, "w");
        write_matrix(out, temperature);
        fclose(out);
        printf("Dump step %d\n", step);
    }

    MPI_Wait(&req, MPI_STATUS_IGNORE);
}

void write_matrix(FILE *out, float *data)
{
    float size = (float) SIZE;
    fwrite(&size, sizeof(float), 1, out);

    for (float x = 0; x < SIZE; x += 1.0) {
        fwrite(&x, sizeof(float), 1, out);
    }

    for (int y = 0; y < SIZE; y++) {
        float len = (float) y;
        fwrite(&len, sizeof(float), 1, out);
        fwrite(&data[y * SIZE], sizeof(float), SIZE, out);
    }
}
