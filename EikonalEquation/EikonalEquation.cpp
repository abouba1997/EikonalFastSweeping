#include "EikonalSolver.h"

double func2d(double x, double y) {
    //return 1.; // RealSolution2d 1
    return 2 * exp(x * x + y * y) * sqrt(x * x + y * y); // Second example RealSolution2d 2

    // From article fast methods
    /*if (x > -0.8 && x < -0.75 && y <= 0.9) {
        return INFINITY;
    }
    return 1.;*/

    /*if ((x > -0.95 && x < -0.9 && y <= 0.9) || (x > -0.75 && x < -0.7 && y > -0.9) || (x > -0.55 && x < -0.5 && y <= 0.9) || (x > -0.35 && x < -0.3 && y > -0.9) || (x > -0.15 && x < -0.1 && y <= 0.9)) {
        return INFINITY;
    }
    return 1.;*/
     
    // From the book
    //return sqrt(pow(sin(PI + x * PI / 2), 2) + pow(sin(PI + y * PI / 2), 2)); // First example
    
    //return 3 * sqrt((x*x + y*y)/pow(9 + x * x + y * y,3)); // Third example
}

double func3d(double x, double y, double z) {
    // First
    return 1.; // Example

    // With walls
    /*if ((x > -0.95 && x < -0.9 && y <= 0.9 && z <= 0.9) || (x > -0.55 && x < -0.5 && y <= 0.9 && z <= 0.9) || (x > -0.15 && x < -0.1 && y <= 0.9 && z <= 0.9)) {
        return INFINITY;
    }
    return 1.;*/

    /*if (x > 0 && y > -0.5 && y < 0.5) {
        return INFINITY;
    }
    return 1.;*/

    // From the book
    //return sqrt(pow(sin(x), 2) + pow(sin(y), 2) + pow(sin(z), 2));
}

double realSolution2d(double x, double y) {
    return -1 + exp(x * x + y * y);
}

double realSolution3d(double x, double y, double z) {
    return sqrt(x * x + y * y + z * z);
}

int main(int argc, char** argv) {

    MPI::Init(argc, argv);
  
    int N = 2000;
    EikonalSolver solver(-1, 1, N, func2d, func3d);

    // 2D - Sequential
    
    //Matrix2dI origins2d = { {N/2, N/2} };
    //solver.setOrigins2d(origins2d);

    // Fast Sweeping Method
    //solver.solveFSM2d(40, "2d_fsm.txt");
    //solver.writeError2D(realSolution2d);
    
    // Fast Marching Method
    //solver.solveFMM2d(10000, "2d_fmm.txt", origins2d);

    // 2D - parallel

    Matrix2dI origins2d = { {N / 2, N / 2} };
    solver.solveFSM2dParallel(40, "parallel_2d.txt", origins2d, realSolution2d, "parallel_2d_error_40.txt");

    // 3D - Sequential
    
    //Matrix2dI origins3d = { {50,50,50} };
    //solver.setOrigins3d(origins3d);

    // Fast Sweeping Method
    //solver.solveFSM3d(5, "3d_fsm.txt");
    //solver.writeRealSolution3D(realSolution3d, "3d_fsm_solution.txt");

    // Fast Marching Method
    //solver.solveFMM3d(2, "3d_fmm.txt", origins3d);

    //3D - Parallel
    /*Matrix2d origins3d = { {50, 50, 50} };
    solver.solveFSM3dParallel(5, "parallel_3d.txt", origins3d);*/

    MPI::Finalize();
    return 0;
}
