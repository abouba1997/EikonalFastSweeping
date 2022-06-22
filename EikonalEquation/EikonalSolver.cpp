#include "EikonalSolver.h"

using namespace std;

EikonalSolver::EikonalSolver(double _a, double _b, int _N, function2d _f2d, function3d _f3d) {
	a = _a; b = _b; f2d = _f2d; f3d = _f3d;
	h = (b - a) / (double)_N;
	N = _N + 1;

	x.resize(N); y.resize(N); z.resize(N);

	int k = 0;
	for (double i = a; i <= b + h; i += h) {
		x[k] = i;
		y[k] = i;
		z[k] = i;
		k++;
	}
}

EikonalSolver::~EikonalSolver() {}

void EikonalSolver::setOrigins2d(Matrix2dI origins) {
	space2d.resize(N);
	for (int i = 0; i < N; i++) {
		space2d[i].resize(N);
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			space2d[i][j] = INFINITY;
		}
	}

	for (int i = 0; i < origins.size(); i++) {
		RowI point = origins[i];
		space2d[point[0]][point[1]] = 0;
	}
}

void EikonalSolver::setOrigins3d(Matrix2dI origins)
{
	space3d.resize(N);
	for (int i = 0; i < N; i++) {
		space3d[i].resize(N);
		for (int j = 0; j < N; j++) {
			space3d[i][j].resize(N);
		}
	}

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				space3d[i][j][k] = INFINITY;
			}
		}
	}

	for (int i = 0; i < origins.size(); i++) {
		vector<int> point = origins[i];
		space3d[point[0]][point[1]][point[2]] = 0;
	}
}

void EikonalSolver::solveFSM2d(int iter, string filename) {
	int iterations = 0;
	Matrix2d sweep_directions = { {1,1},{0,1},{0,0},{1,0} };
	
	double tstart, tfinish;
	tstart = MPI_Wtime();
	
	do {
		for (int sweep = 0; sweep < sweep_directions.size(); sweep++) {
			int iStart = sweep_directions[sweep][0] ? 0 : N - 1;
			int iEnd = sweep_directions[sweep][0] ? N : 0;

			int jStart = sweep_directions[sweep][1] ? 0 : N - 1;
			int jEnd = sweep_directions[sweep][1] ? N: 0;

			for (int i = iStart; i != iEnd; i = (sweep_directions[sweep][0]) ? i + 1 : i - 1) {
				for (int j = jStart; j != jEnd; j = (sweep_directions[sweep][1]) ? j + 1 : j - 1) {
					sweep2d(space2d, i, j);
				}
			}
		}
		
		iterations++;
	} while (iterations != iter);
	
	tfinish = MPI_Wtime() - tstart;

	cout << "Sequential FSM 2D solving" << endl;
	cout << "  FSM - N: \tN = " << N << endl;
	cout << "  FSM - Time: \tT = " << tfinish << endl << endl;
	
	writeInfile2d(filename);
}

void EikonalSolver::solveFSM2dParallel(int iter, string filename, Matrix2dI origins, function2d func, string errorfile)
{
	double* space2dparallel = new double[N * N];
	double* space2dparalleltmp = new double[N * N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			space2dparalleltmp[i * N + j] = INFINITY;
		}
	}

	for (int i = 0; i < origins.size(); i++) {
		RowI point = origins[i];
		space2dparalleltmp[point[0] * N + point[1]] = 0;
	}

	Matrix2d sweep_directions = { {1,1},{0,1},{0,0},{1,0} };

	int size, rank;
	MPI::Status stat;
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	int iterations = 0;

	double tstart, tfinish;
	tstart = MPI_Wtime();

	do{
		for (int my_rank = 0; my_rank < 4; my_rank+=size)
		{
			int iStart = sweep_directions[rank][0] ? 0 : N - 1;
			int iEnd = sweep_directions[rank][0] ? N : 0;

			int jStart = sweep_directions[rank][1] ? 0 : N - 1;
			int jEnd = sweep_directions[rank][1] ? N : 0;

			for (int i = iStart; i != iEnd; i = (sweep_directions[rank][0]) ? i + 1 : i - 1) {
				for (int j = jStart; j != jEnd; j = (sweep_directions[rank][1]) ? j + 1 : j - 1) {
					sweep2dParallel(space2dparalleltmp, i, j);
				}
			}

			// Need to combine computed values
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Allreduce(space2dparalleltmp, space2dparallel, N * N, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

			std::swap(space2dparalleltmp, space2dparallel);
		}
		iterations++;
	} while (iterations != iter);

	// Need to combine computed values
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank == 0) {
		tfinish = MPI_Wtime() - tstart;
		cout << "Parallel FSM 2D solving" << endl;
		cout << "  FSM - N: \tN = " << N << endl;
		cout << "  FSM - Time: \tT = " << tfinish << endl << endl;

		writeInfile2dParallel(filename, space2dparallel);

		writeError2DParallel(space2dparallel, func, errorfile);
	}

	delete[] space2dparallel;
	delete[] space2dparalleltmp;
}

void EikonalSolver::writeRealSolution2D(function2d func, string result)
{
	ofstream file(result);

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			file << x[i] << "\t" << y[j] << "\t" << func(x[i], y[j]) << endl;
		}
		file << endl;
	}
	file.close();
}

void EikonalSolver::writeError2D(function2d func)
{
	double result = 0, sum = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++)
		{
			result = abs(space2d[i][j] - func(x[i], y[j]));
			if (result > sum) {
				sum = result;
			}
		}
	}
	cout << "N = " << N << "\twith norm = " << sum << std::endl;
}

void EikonalSolver::writeError2DParallel(double* space2d, function2d func, string filename)
{
	ofstream file(filename);
 	double result = 0, sum = 0;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			result = abs(space2d[i * N + j] - func(x[i], y[j]));
			//file << result << endl;
			//sum += result * result;
			if (result > sum) {
				sum = result;
			}
		}
		//file << sqrt(sum / (N + 1)) << endl;
		file << sum << endl;
		sum = 0.;
	}
	file.close();

	cout << "N = " << N << "\twith norm = " << sqrt(sum / (N + 1)) << std::endl;
}

void EikonalSolver::solveFSM3d(int iter, string filename)
{
	// Cycle
	int iteration = 0;
	Matrix2d sweep_directions = { {1,1,1},
								  {0,1,0},
								  {0,1,1},
								  {1,1,0},
								  {0,0,0},
								  {1,0,1},
								  {1,0,0},
								  {0,0,1},
								};
	double tstart, tfinish;
	tstart = MPI_Wtime();

	do {
		for (int sweep = 0; sweep < sweep_directions.size(); sweep++) {
			int iStart = sweep_directions[sweep][0] ? 0 : N - 1;
			int iEnd = sweep_directions[sweep][0] ? N : 0;

			int jStart = sweep_directions[sweep][1] ? 0 : N - 1;
			int jEnd = sweep_directions[sweep][1] ? N : 0;

			int kStart = sweep_directions[sweep][2] ? 0 : N - 1;
			int kEnd = sweep_directions[sweep][2] ? N : 0;

			for (int i = iStart; i != iEnd; i = (sweep_directions[sweep][0]) ? i + 1 : i - 1)
			{
				for (int j = jStart; j != jEnd; j = (sweep_directions[sweep][1]) ? j + 1 : j - 1)
				{
					for (int k = kStart; k != kEnd; k = (sweep_directions[sweep][2]) ? k + 1 : k - 1)
					{
						sweep3d(space3d, i, j, k);
					}
				}
			}
			
		}

		iteration++;
	} while (iteration != iter);
	
	tfinish = MPI_Wtime() - tstart;
	cout << "Sequential FSM 3D solving" << endl;
	cout << "  FSM - N: \tN = " << N << endl;
	cout << "  FSM - Time: \tT = " << tfinish << endl << endl;

	writeInfile3d(filename);
}

void EikonalSolver::solveFSM3dParallel(int iter, string filename, Matrix2d origins)
{
	double* space3dparallel = new double[N * N * N];
	double* space3dparalleltmp = new double[N * N * N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				space3dparalleltmp[i * N * N + j * N + k] = INFINITY;
			}
		}
	}

	for (int i = 0; i < origins.size(); i++) {
		Row point = origins[i];
		space3dparalleltmp[(int)point[0] * N  * N + (int)point[1] * N + (int)point[2]] = 0;
	}

	Matrix2d sweep_directions = { {1,1,1},
								  {0,1,0},
								  {0,1,1},
								  {1,1,0},
								  {0,0,0},
								  {1,0,1},
								  {1,0,0},
								  {0,0,1},
	};

	int size, rank;
	double tstart, tfinish;
	MPI::Status stat;
	size = MPI::COMM_WORLD.Get_size();
	rank = MPI::COMM_WORLD.Get_rank();

	int iterations = 0;

	MPI_Barrier(MPI_COMM_WORLD);
	tstart = MPI_Wtime();

	do {
		for (int my_rank = rank; my_rank < 8; my_rank += size)
		{
			int iStart = sweep_directions[my_rank][0] ? 0 : N - 1;
			int iEnd = sweep_directions[my_rank][0] ? N : 0;

			int jStart = sweep_directions[my_rank][1] ? 0 : N - 1;
			int jEnd = sweep_directions[my_rank][1] ? N : 0;


			int kStart = sweep_directions[my_rank][2] ? 0 : N - 1;
			int kEnd = sweep_directions[my_rank][2] ? N : 0;

			for (int i = iStart; i != iEnd; i = (sweep_directions[my_rank][0]) ? i + 1 : i - 1) {
				for (int j = jStart; j != jEnd; j = (sweep_directions[my_rank][1]) ? j + 1 : j - 1) {
					for (int k = kStart; k != kEnd; k = (sweep_directions[my_rank][2]) ? k + 1 : k - 1)
					{
						sweep3dParallel(space3dparalleltmp, i, j, k);
					}
				}
			}

		}
		
		// Need to combine computed values
		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Allreduce(space3dparalleltmp, space3dparallel, N * N * N, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

		std::swap(space3dparalleltmp, space3dparallel);

		iterations++;
	} while (iterations != iter);
	MPI_Barrier(MPI_COMM_WORLD);
	tfinish = MPI_Wtime() - tstart;

	if (rank == 0) {
		cout << "Parallel FSM 3D solving" << endl;
		cout << "  FSM - Iteration: Iter = " << iterations << endl;
		cout << "  FSM - Time: \tT = " << tfinish << endl << endl;
		writeInfile3dParallel(filename, space3dparallel);
	}

	delete[] space3dparallel;
	delete[] space3dparalleltmp;
}

void EikonalSolver::writeRealSolution3D(function3d func, string result)
{
	ofstream file(result);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				file << x[i] << "\t" << y[j] << "\t" << z[k] << "\t" << func(x[i], y[j], z[k]) << endl;
			}
			file << endl;
		}
	}
	file.close();
}

void EikonalSolver::solveFMM2d(int iters, string filename, Matrix2dI origins)
{
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			space2d[i][j] = INFINITY;
		}
	}

	// Beginning
	double* far_away = new double[N * N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			far_away[i * N + j] = 2;
		}
	}

	for (int i = 0; i < origins.size(); i++) {
		RowI point = origins[i];
		space2d[point[0]][point[1]] = 0;
		far_away[point[0] * N + point[1]] = 0;
	}

	// Get Neighbors of Origins points
	Matrix2dI narrowBandIndexOfOrigins = getIndexesNarrowBand(origins);

	// Label Neighbors of Origins points as narrowBand
	for (int i = 0; i < narrowBandIndexOfOrigins.size(); i++) {
		RowI row = narrowBandIndexOfOrigins[i];
		far_away[row[0] * N + row[1]] = 1;

		sweep2d(space2d, row[0], row[1]);
	}

	int iter = 1;
	double tstart, tfinish;
	tstart = MPI_Wtime();
	while (true) {
		// Get all Indexes of narrowBand
		Matrix2dI indexesOnesNarrowBand = getIndexesNarrowBand(far_away);
		// When not finished all element inside narrowband
		if (!indexesOnesNarrowBand.empty()) {
			// Get the minimum element index
			RowI minIndex = getIndexesMinNarrowBand(space2d, indexesOnesNarrowBand);
			// Add the minimum to the frozen elements
			far_away[minIndex[0] * N + minIndex[1]] = 0;

			// Get the neighbors of the minvalue
			Matrix2dI newNeighborsIndexes;
			for (int j = 0; j < neighbors.size(); j++) {
				RowI ind;
				ind.push_back(minIndex[0] + neighbors[j][0]);
				ind.push_back(minIndex[1] + neighbors[j][1]);
				newNeighborsIndexes.push_back(ind);
			}

			// Put 1 where 2 in far_away (label far away as narrowBand)
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				RowI current = newNeighborsIndexes[i];
				if (far_away[current[0] * N + current[1]] == 2) {
					far_away[current[0] * N + current[1]] = 1;
				}
			}
			// We got new narrowBand... must get all indexes and recalculate distance field value for new added...

			// Calculate U tentative of that
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				RowI currentIndex = newNeighborsIndexes[i];
				sweep2d(space2d, currentIndex[0], currentIndex[1]);
			}
		}

		if (iter == iters) {
			break;
		}
		iter++;
	}
	
	tfinish = MPI_Wtime() - tstart;
	cout << "Sequential FMM - Time: \tt = " << tfinish << endl;
	cout << "  FMM - N: \tN = " << N << endl;
	cout << "  FMM solving" << endl;

	writeInfile2d(filename);

	delete[] far_away;
}

void EikonalSolver::solveFMM3d(int iters, string filename, Matrix2dI origins)
{
	// Beginning
	double* far_away = new double[N * N * N];

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				far_away[i * N * N + j * N + k] = 2;
			}
		}
	}

	for (int i = 0; i < origins.size(); i++) {
		RowI point = origins[i];
		far_away[point[0] * N * N + point[1] * N + point[2]] = 0;
	}


	// Debugging
	/*for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				cout << far_away[i * N * N + j * N + k] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}*/

	// Get Neighbors of Origins points
	Matrix2dI narrowBandIndexOfOrigins = getIndexesNarrowBand3(origins);

	// Label Neighbors of Origins points as narrowBand
	for (int i = 0; i < narrowBandIndexOfOrigins.size(); i++) {
		RowI row = narrowBandIndexOfOrigins[i];
		far_away[row[0] * N * N + row[1] * N + row[2]] = 1;

		sweep3d(space3d, row[0], row[1], row[2]);
	}

	int iter = 1;
	double tstart, tfinish;
	tstart = MPI_Wtime();
	while (true) {
		// Get all Indexes of narrowBand
		Matrix2dI indexesOnesNarrowBand = getIndexesNarrowBand3(far_away);
		// When not finished all element inside narrowband
		if (!indexesOnesNarrowBand.empty()) {
			// Get the minimum element index
			RowI minIndex = getIndexesMinNarrowBand3(space3d, indexesOnesNarrowBand);
			// Add the minimum to the frozen elements
			far_away[minIndex[0] * N * N + minIndex[1] * N + minIndex[2]] = 0;

			// Get the neighbors of the minvalue
			Matrix2dI newNeighborsIndexes;

			for (int j = 0; j < neighbors3d.size(); j++) {
				RowI ind;
				ind.push_back(minIndex[0] + neighbors3d[j][0]);
				ind.push_back(minIndex[1] + neighbors3d[j][1]);
				ind.push_back(minIndex[2] + neighbors3d[j][2]);
				newNeighborsIndexes.push_back(ind);
			}

			// Put 1 where 2 in far_away (label far away as narrowBand)
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				RowI current = newNeighborsIndexes[i];
				if (far_away[current[0] * N * N + current[1] * N + current[2]] == 2) {
					far_away[current[0] * N * N + current[1] * N + current[2]] = 1;
				}
			}
			// We got new narrowBand... must get all indexes and recalculate distance field value for new added...

			// Calculate U tentative of that
			for (int i = 0; i < newNeighborsIndexes.size(); i++) {
				RowI currentIndex = newNeighborsIndexes[i];
				sweep3d(space3d, currentIndex[0], currentIndex[1], currentIndex[2]);
			}
		}

		if (iter == iters) {
			break;
		}
		iter++;
		break;
	}

	tfinish = MPI_Wtime() - tstart;
	cout << "Sequential FMM 3D solving" << endl;
	cout << "Sequential FMM - Time: \tt = " << tfinish << endl;

	//writeInfile3d(filename);

	delete[] far_away;
}

void EikonalSolver::writeInfile2d(string filename) {
	ofstream file(filename);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			file << x[i] << "\t" << y[j] << "\t" << space2d[i][j] << std::endl;
		}
		file << "\n";
	}
	file.close();
}

void EikonalSolver::writeInfile2dParallel(string filename, double* space2dparallel)
{
	ofstream file(filename);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			file << x[i] << "\t" << y[j] << "\t" << space2dparallel[i * N + j] << std::endl;
		}
		file << "\n";
	}
	file.close();
}

void EikonalSolver::writeInfile3d(string filename)
{
	ofstream file(filename);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				file << x[i] << "\t" << y[j] << "\t" << z[k] << "\t" << space3d[i][j][k] << std::endl;
			}
			file << "\n";
		}
		file << "\n";
	}
	file.close();
}

void EikonalSolver::writeInfile3d_(string filename)
{
	ofstream file(filename);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				file << space3d[i][j][k] << "\t";
			}
			file << "\n";
		}
		file << "\n";
	}
	file.close();
}

void EikonalSolver::writeInfile3dParallel(string filename, double* space3dparallel)
{
	ofstream file(filename);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				file << x[i] << "\t" << y[j] << "\t" << z[k] << "\t" << space3dparallel[i * N * N + j * N + k] << std::endl;
			}
			file << "\n";
		}
		file << "\n";
	}
	file.close();

	std::cout << "Result writting successfully!!!" << std::endl;
}

void EikonalSolver::sweep3d(Matrix3d& space3d, int i, int j, int k)
{
	double txmin, tymin, tzmin, t;
	double b, c;

	if (!isnan(space3d[i][j][k])) {
		// by i
		if (i == 0) {
			tymin = space3d[i + 1][j][k];
		}
		else if (i == N - 1) {
			tymin = space3d[i - 1][j][k];
		}
		else {
			tymin = min(space3d[i - 1][j][k], space3d[i + 1][j][k]);
		}

		// by j
		if (j == 0) {
			txmin = space3d[i][j + 1][k];
		}
		else if (j == N - 1) {
			txmin = space3d[i][j - 1][k];
		}
		else {
			txmin = min(space3d[i][j - 1][k], space3d[i][j + 1][k]);
		}

		// by k
		if (k == 0) {
			tzmin = space3d[i][j][k + 1];
		}
		else if (k == N - 1) {
			tzmin = space3d[i][j][k - 1];
		}
		else {
			tzmin = min(space3d[i][j][k - 1], space3d[i][j][k + 1]);
		}

		//result
		double l = txmin * txmin + tymin * tymin + tzmin * tzmin;
		b = pow((txmin + tymin + tzmin), 2);
		c = 3 * (l - h * h * f3d(x[i], y[j], z[k]));

		if ((b - c) > 0) {
			t = ((txmin + tymin + tzmin) + sqrt(b - c)) / 3.;
		}
		else {
			t = min3(txmin, tymin, tzmin) + f3d(x[i], y[j], z[k]) * h;
		}

		// set space data
		space3d[i][j][k] = min(space3d[i][j][k], t);
	}
}

void EikonalSolver::sweep2d(Matrix2d& space2d, int i, int j)
{
	double txmin, tymin, t, eps = 0;
	if (!isnan(space2d[i][j])) {
		// by i
		if (i == 0) {
			tymin = space2d[i + 1][j];
		}
		else if (i == N - 1) {
			tymin = space2d[i - 1][j];
		}
		else {
			tymin = min(space2d[i - 1][j], space2d[i + 1][j]);
		}

		// by j
		if (j == 0) {
			txmin = space2d[i][j + 1];
		}
		else if (j == N - 1) {
			txmin = space2d[i][j - 1];
		}
		else {
			txmin = min(space2d[i][j - 1], space2d[i][j + 1]);
		}

		//result
		if (abs(txmin - tymin) >= f2d(x[i], y[j]) * h) {
			t = min(txmin, tymin) + f2d(x[i], y[j]) * h;
		}
		else {
			double sqrtP = 2 * pow(f2d(x[i], y[j]), 2) * h * h - pow(txmin - tymin, 2);
			t = (txmin + tymin + sqrt(sqrtP)) / 2.;
		}
		// set space data
		space2d[i][j] = min(space2d[i][j], t);
	}
}

void EikonalSolver::sweep2dParallel(double* space2d, int i, int j)
{
	double txmin, tymin, t;
	if (!isnan(space2d[i * N + j])) {
		// by i
		if (i == 0) {
			tymin = space2d[(i + 1) * N + j];
		}
		else if (i == N - 1) {
			tymin = space2d[(i - 1) * N + j];
		}
		else {
			tymin = min(space2d[(i - 1) * N + j], space2d[(i + 1) * N + j]);
		}

		// by j
		if (j == 0) {
			txmin = space2d[i * N + (j + 1)];
		}
		else if (j == N - 1) {
			txmin = space2d[i * N + (j - 1)];
		}
		else {
			txmin = min(space2d[i * N + (j - 1)], space2d[i * N + (j + 1)]);
		}

		//result
		if (abs(txmin - tymin) >= f2d(x[i], y[j]) * h) {
			t = min(txmin, tymin) + f2d(x[i], y[j]) * h;
		}
		else {
			double sqrtP = 2 * pow(f2d(x[i], y[j]), 2) * h * h - pow(txmin - tymin, 2);
			t = (txmin + tymin + sqrt(sqrtP)) / 2.;
		}
		// set space data
		space2d[i * N + j] = min(space2d[i * N + j], t);
	}
}

void EikonalSolver::sweep3dParallel(double* space3d, int i, int j, int k)
{
	double txmin, tymin, tzmin, t;
	double b, c;

	if (!isnan(space3d[i * N * N + j * N + k])) {
		// by i
		if (i == 0) {
			tymin = space3d[(i + 1) * N * N + j * N + k];
		}
		else if (i == N - 1) {
			tymin = space3d[(i - 1) * N * N + j * N + k];
		}
		else {
			tymin = min(space3d[(i - 1) * N * N + j * N + k], space3d[(i + 1) * N * N + j * N + k]);
		}

		// by j
		if (j == 0) {
			txmin = space3d[i * N * N + (j + 1) * N + k];
		}
		else if (j == N - 1) {
			txmin = space3d[i * N * N + (j - 1) * N + k];
		}
		else {
			txmin = min(space3d[i * N * N + (j - 1) * N + k], space3d[i * N * N + (j + 1) * N + k]);
		}

		// by k
		if (k == 0) {
			tzmin = space3d[i * N * N + j * N + k + 1];
		}
		else if (k == N - 1) {
			tzmin = space3d[i * N * N + j * N + k - 1];
		}
		else {
			tzmin = min(space3d[i * N * N + j * N + k - 1], space3d[i * N * N + j * N + k + 1]);
		}

		//result
		double l = txmin * txmin + tymin * tymin + tzmin * tzmin;
		b = pow((txmin + tymin + tzmin), 2);
		c = 3 * (l - h * h * f3d(x[i], y[j], z[k]));

		if ((b - c) > 0) {
			t = ((txmin + tymin + tzmin) + sqrt(b - c)) / 3.;
		}
		else {
			t = min3(txmin, tymin, tzmin) + f3d(x[i], y[j], z[k]) * h;
		}

		// set space data
		space3d[i * N * N + j * N + k] = min(space3d[i * N * N + j * N + k], t);
	}
}

Matrix2dI EikonalSolver::getIndexesNarrowBand(Matrix2dI origins)
{
	vector<vector<int>> narrowBandIndexOfOrigins;

	for (int i = 0; i < origins.size(); i++) {
		if (origins[i][0] == 0 && origins[i][1] == 0) {
			for (int k = 0; k < top_left_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_left_neighbors[k][0]);
				ind.push_back(origins[i][1] + top_left_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1) {
			for (int k = 0; k < top_right_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_right_neighbors[k][0]);
				ind.push_back(origins[i][1] + top_right_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0) {
			for (int k = 0; k < bottom_left_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_left_neighbors[k][0]);
				ind.push_back(origins[i][1] + bottom_left_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1) {
			for (int k = 0; k < bottom_right_neighbors.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_right_neighbors[k][0]);
				ind.push_back(origins[i][1] + bottom_right_neighbors[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0) {
			for (int k = 0; k < top_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary[k][0]);
				ind.push_back(origins[i][1] + top_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1) {
			for (int k = 0; k < bottom_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == 0) {
			for (int k = 0; k < left_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + left_boundary[k][0]);
				ind.push_back(origins[i][1] + left_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == N - 1) {
			for (int k = 0; k < right_boundary.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + right_boundary[k][0]);
				ind.push_back(origins[i][1] + right_boundary[k][0]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else {
			for (int j = 0; j < neighbors.size(); j++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + neighbors[j][0]);
				ind.push_back(origins[i][1] + neighbors[j][1]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
	}

	return narrowBandIndexOfOrigins;
}

Matrix2dI EikonalSolver::getIndexesNarrowBand3(Matrix2dI origins)
{
	Matrix2dI narrowBandIndexOfOrigins;

	for (int i = 0; i < origins.size(); i++) {
		if (origins[i][0] == 0 && origins[i][1] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_far_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_far_left[k][0]);
				ind.push_back(origins[i][1] + bottom_far_left[k][1]);
				ind.push_back(origins[i][2] + bottom_far_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_near_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_near_left[k][0]);
				ind.push_back(origins[i][1] + bottom_near_left[k][1]);
				ind.push_back(origins[i][2] + bottom_near_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_far_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_far_right[k][0]);
				ind.push_back(origins[i][1] + bottom_far_right[k][1]);
				ind.push_back(origins[i][2] + bottom_far_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_near_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_near_right[k][0]);
				ind.push_back(origins[i][1] + bottom_near_right[k][1]);
				ind.push_back(origins[i][2] + bottom_near_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_far_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_far_left[k][0]);
				ind.push_back(origins[i][1] + top_far_left[k][1]);
				ind.push_back(origins[i][2] + top_far_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_near_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_near_left[k][0]);
				ind.push_back(origins[i][1] + top_near_left[k][1]);
				ind.push_back(origins[i][2] + top_near_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_far_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_far_right[k][0]);
				ind.push_back(origins[i][1] + top_far_right[k][1]);
				ind.push_back(origins[i][2] + top_far_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_near_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_near_right[k][0]);
				ind.push_back(origins[i][1] + top_near_right[k][1]);
				ind.push_back(origins[i][2] + top_near_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		} //top
		else if (origins[i][0] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_far.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_far[k][0]);
				ind.push_back(origins[i][1] + top_boundary_far[k][1]);
				ind.push_back(origins[i][2] + top_boundary_far[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_near.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_near[k][0]);
				ind.push_back(origins[i][1] + top_boundary_near[k][1]);
				ind.push_back(origins[i][2] + top_boundary_near[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == N - 1 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_left[k][0]);
				ind.push_back(origins[i][1] + top_boundary_left[k][1]);
				ind.push_back(origins[i][2] + top_boundary_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == 0 && origins[i][2] == N - 1) {
			for (int k = 0; k < top_boundary_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + top_boundary_right[k][0]);
				ind.push_back(origins[i][1] + top_boundary_right[k][1]);
				ind.push_back(origins[i][2] + top_boundary_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}// Bottom
		else if (origins[i][0] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_far.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_far[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_far[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_far[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_near.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_near[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_near[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_near[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == 0 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_left[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_left[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][1] == N - 1 && origins[i][2] == 0) {
			for (int k = 0; k < bottom_boundary_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + bottom_boundary_right[k][0]);
				ind.push_back(origins[i][1] + bottom_boundary_right[k][1]);
				ind.push_back(origins[i][2] + bottom_boundary_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		} // Middle
		else if (origins[i][0] == 0 && origins[i][1] == 0) {
			for (int k = 0; k < middle_boundary_far_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_far_left[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_far_left[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_far_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == 0 && origins[i][1] == N - 1) {
			for (int k = 0; k < middle_boundary_far_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_far_right[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_far_right[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_far_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == 0) {
			for (int k = 0; k < middle_boundary_near_left.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_near_left[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_near_left[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_near_left[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
		else if (origins[i][0] == N - 1 && origins[i][1] == N - 1) {
			for (int k = 0; k < middle_boundary_near_right.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + middle_boundary_near_right[k][0]);
				ind.push_back(origins[i][1] + middle_boundary_near_right[k][1]);
				ind.push_back(origins[i][2] + middle_boundary_near_right[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}// other
		else {
			for (int k = 0; k < neighbors3d.size(); k++) {
				vector<int> ind;
				ind.push_back(origins[i][0] + neighbors3d[k][0]);
				ind.push_back(origins[i][1] + neighbors3d[k][1]);
				ind.push_back(origins[i][2] + neighbors3d[k][2]);
				narrowBandIndexOfOrigins.push_back(ind);
			}
		}
	}

	return narrowBandIndexOfOrigins;
}

Matrix2dI EikonalSolver::getIndexesNarrowBand(double* far)
{
	vector<vector<int>> indexes;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (far[i * N + j] == 1) {
				vector<int> ind;
				ind.push_back(i);
				ind.push_back(j);
				indexes.push_back(ind);
				/*if (i != 0 && i != N - 1 && j != 0 && j != N - 1) {
					vector<int> ind;
					ind.push_back(i);
					ind.push_back(j);
					indexes.push_back(ind);
				}
				else {
					if (i == 0 && j == 0) {
						vector<int> ind;
						ind.push_back(i);
						ind.push_back(j);
						indexes.push_back(ind);
					}
				}*/
			}
		}
	}

	return indexes;
}

Matrix2dI EikonalSolver::getIndexesNarrowBand3(double* far)
{
	Matrix2dI indexes;

	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				if (far[i * N * N + j * N + k] == 1) {
					vector<int> ind;
					ind.push_back(i);
					ind.push_back(j);
					ind.push_back(k);
					indexes.push_back(ind);

					/*if (i != 0 && i != N - 1 && j != 0 && j != N - 1 && k != 0 && k != N - 1) {
						vector<int> ind;
						ind.push_back(i);
						ind.push_back(j);
						ind.push_back(k);
						indexes.push_back(ind);
					}
					else {
						if (i == 0 && j == 0 && k == 0) {
							vector<int> ind;
							ind.push_back(i);
							ind.push_back(j);
							ind.push_back(k);
							indexes.push_back(ind);
						}
					}*/
				}
			}
		}
	}

	return indexes;
}

RowI EikonalSolver::getIndexesMinNarrowBand(Matrix2d T, Matrix2dI neigh)
{
	RowI minIndex = neigh[0];
	double minT = T[minIndex[0]][minIndex[1]];

	for (int i = 1; i < neigh.size(); i++) {
		RowI currentIndex = neigh[i];

		double currentMin = T[currentIndex[0]][currentIndex[1]];

		if (currentMin < minT) {
			minT = currentMin;
			minIndex = currentIndex;
		}
	}

	return minIndex;
}

RowI EikonalSolver::getIndexesMinNarrowBand3(Matrix3d T, Matrix2dI neigh)
{
	RowI minIndex = neigh[0];
	double minT = T[minIndex[0]][minIndex[1]][minIndex[2]];

	for (int i = 1; i < neigh.size(); i++) {
		RowI currentIndex = neigh[i];

		double currentMin = T[currentIndex[0]][currentIndex[1]][currentIndex[2]];

		if (currentMin < minT) {
			minT = currentMin;
			minIndex = currentIndex;
		}
	}

	return minIndex;
}
