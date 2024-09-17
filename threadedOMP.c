#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <omp.h> // Include OpenMP header
#include "timer.h"

#define SOFTENING 1e-9f

typedef struct {
    float x, y, z;
    float vx, vy, vz;
} Body;

void initializeBodies(Body *bodies, int nBodies) {
    for (int i = 0; i < nBodies; i++) {
        bodies[i].x = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        bodies[i].y = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        bodies[i].z = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
        bodies[i].vx = 0.0f;
        bodies[i].vy = 0.0f;
        bodies[i].vz = 0.0f;
    }
}

void computeForces(Body *bodies, int nBodies, float dt) {
    // Parallelize the outer loop
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < nBodies; i++) {
        float Fx = 0.0f, Fy = 0.0f, Fz = 0.0f;

        for (int j = 0; j < nBodies; j++) {
            if (i != j) { // Avoid self-interaction
                float dx = bodies[j].x - bodies[i].x;
                float dy = bodies[j].y - bodies[i].y;
                float dz = bodies[j].z - bodies[i].z;
                float distSqr = dx * dx + dy * dy + dz * dz + SOFTENING;
                float invDist = 1.0f / sqrtf(distSqr);
                float invDist3 = invDist * invDist * invDist;

                // Accumulate force contributions
                Fx += dx * invDist3;
                Fy += dy * invDist3;
                Fz += dz * invDist3;
            }
        }

        // Update velocities based on accumulated forces
        bodies[i].vx += dt * Fx;
        bodies[i].vy += dt * Fy;
        bodies[i].vz += dt * Fz;
    }
}

void updatePositions(Body *bodies, int nBodies, float dt) {
    // Parallelize the position update loop
    #pragma omp parallel for
    for (int i = 0; i < nBodies; i++) {
        bodies[i].x += bodies[i].vx * dt;
        bodies[i].y += bodies[i].vy * dt;
        bodies[i].z += bodies[i].vz * dt;
    }
}

int main(int argc, char** argv) {
    int nBodies = 30000;
    if (argc > 1) nBodies = atoi(argv[1]);

    const float dt = 0.01f; // Time step
    const int nIters = 10;  // Number of iterations

    // Allocate memory for bodies
    Body *bodies = (Body*)malloc(nBodies * sizeof(Body));
    if (bodies == NULL) {
        printf("Error: Could not allocate memory for bodies.\n");
        return 1;
    }

    srand(time(NULL));
    initializeBodies(bodies, nBodies);

    double totalTime = 0.0;

    // Simulation loop
    for (int iter = 1; iter <= nIters; iter++) {
        StartTimer();

        // Compute forces in parallel
        computeForces(bodies, nBodies, dt);

        // Update positions in parallel
        updatePositions(bodies, nBodies, dt);

        double elapsedTime = GetTimer() / 1000.0; // Time in seconds
        if (iter > 1) { // Skip the first iteration for averaging
            totalTime += elapsedTime;
        }

        printf("Iteration %d: %.3f seconds\n", iter, elapsedTime);
    }

    // Report average time and performance
    double avgTime = totalTime / (double)(nIters - 1);
    double interactionsPerSecond = (nBodies * nBodies) / avgTime * 1e-9;
    printf("Average time per iteration: %.3f seconds\n", avgTime);
    printf("Performance: %.3f Billion Interactions / second\n", interactionsPerSecond);

    free(bodies);
    return 0;
}
