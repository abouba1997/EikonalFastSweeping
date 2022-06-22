#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include <numeric>
#include <mpi.h>

#define PI 3.14159265359

using namespace std;

#define Matrix3d vector<vector<vector<double>>>
#define Matrix3dI vector<vector<vector<int>>>
#define Matrix2d vector<vector<double>>
#define Matrix2dI vector<vector<int>>
#define Row vector<double>
#define RowI vector<int>
#define function2d function<double(double, double)>
#define function3d function<double(double, double, double)>

#define min2(a, b) ((a) < (b) ? (a) : (b))
#define min3(a, b, c) (min2(min2((a), (b)), (c)))

class EikonalSolver {
	public:
		EikonalSolver(double, double, int, function2d, function3d);
		~EikonalSolver();
		
		void setOrigins2d(Matrix2dI);
		void setOrigins3d(Matrix2dI);
		
		void solveFSM2d(int, string);
		void solveFSM2dParallel(int, string, Matrix2dI, function2d, string);
		void writeRealSolution2D(function2d, string);
		void writeError2D(function2d);
		void writeError2DParallel(double*,function2d, string);
		void solveFSM3d(int, string);
		void solveFSM3dParallel(int, string, Matrix2d);
		void writeRealSolution3D(function3d, string);

		void solveFMM2d(int, string, Matrix2dI);
		void solveFMM3d(int, string, Matrix2dI);

	private:
		double a, b, h;
		int N;
		Matrix2d space2d;
		Matrix3d space3d;
		Row x, y, z;
		function2d f2d;
		function3d f3d;
		void writeInfile2d(string);
		void writeInfile2dParallel(string, double*);
		void writeInfile3d(string);
		void writeInfile3d_(string);
		void writeInfile3dParallel(string, double*);
		void sweep3d(Matrix3d&, int, int, int);
		void sweep2d(Matrix2d&, int, int);
		void sweep2dParallel (double*, int, int);
		void sweep3dParallel (double*, int, int, int);

		Matrix2dI getIndexesNarrowBand(Matrix2dI);
		Matrix2dI getIndexesNarrowBand(double*);
		Matrix2dI getIndexesNarrowBand3(Matrix2dI);
		Matrix2dI getIndexesNarrowBand3(double*);
		RowI getIndexesMinNarrowBand(Matrix2d, Matrix2dI);
		RowI getIndexesMinNarrowBand3(Matrix3d, Matrix2dI);

		Matrix2dI neighbors = { {1,0},{-1,0},{0,1},{0,-1} };

		Matrix2dI top_left_neighbors = { {1,0},{0,1} };
		Matrix2dI top_right_neighbors = { {1,0},{0,-1} };
		Matrix2dI bottom_left_neighbors = { {-1,0},{0,1} };
		Matrix2dI bottom_right_neighbors = { {-1,0},{0,-1} };

		Matrix2dI top_boundary = { {-1,0},{1,0},{0,1} };
		Matrix2dI bottom_boundary = { {-1,0},{1,0},{0,-1} };
		Matrix2dI left_boundary = { {1,0},{0,1},{0,-1} };
		Matrix2dI right_boundary = { {-1,0},{0,1},{0,-1} };

		Matrix2dI neighbors3d = { {1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1} };
		
		Matrix2d top_far_left = { {1,0,0},{0,1,0},{0,0,-1} };
		Matrix2d top_far_right = { {1,0,0},{0,-1,0},{0,0,-1} };
		Matrix2d top_near_left = { {-1,0,0},{0,1,0},{0,0,-1} };
		Matrix2d top_near_right = { {-1,0,0},{0,-1,0},{0,0,-1} };
		Matrix2d bottom_far_left = { {1,0,0},{0,1,0},{0,0,1} };
		Matrix2d bottom_far_right = { {1,0,0},{0,-1,0},{0,0,1} };
		Matrix2d bottom_near_left = { {-1,0,0},{0,1,0},{0,0,1} };
		Matrix2d bottom_near_right = { {-1,0,0},{0,-1,0},{0,0,1} };

		Matrix2d top_boundary_far = { {1,0,0},{0,1,0},{0,0,-1},{0,-1,0} };
		Matrix2d top_boundary_near = { {-1,0,0},{0,1,0},{0,0,-1},{0,-1,0} };
		Matrix2d top_boundary_left = { {1,0,0},{0,1,0},{-1,0,0},{0,0,-1} };
		Matrix2d top_boundary_right = { {1,0,0},{0,-1,0},{-1,0,0},{0,0,-1} };
		Matrix2d bottom_boundary_far = { {1,0,0},{0,1,0},{0,0,1},{0,-1,0} };
		Matrix2d bottom_boundary_near = { {-1,0,0},{0,1,0},{0,0,1},{0,-1,0} };
		Matrix2d bottom_boundary_left = { {1,0,0},{-1,0,0},{0,0,1},{0,-1,0} };
		Matrix2d bottom_boundary_right = { {1,0,0},{-1,0,0},{0,0,1},{0,1,0} };
		Matrix2d middle_boundary_near_left = { {-1,0,0},{0,1,0},{0,0,1},{0,0,-1} };
		Matrix2d middle_boundary_near_right = { {-1,0,0},{0,-1,0},{0,0,1},{0,0,-1} };
		Matrix2d middle_boundary_far_left = { {1,0,0},{0,1,0},{0,0,1},{0,0,-1} };
		Matrix2d middle_boundary_far_right = { {1,0,0},{0,-1,0},{0,0,1},{0,0,-1} };
};
