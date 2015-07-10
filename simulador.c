#include <stdio.h>
#include <string.h>
#include <math.h>

typedef struct
{
	double x;
	double y;
} vector2_t;

typedef struct
{
	double k;          ///< Rigidity coefficient
	double len;        ///< Relaxed length
	vector2_t anchor;  ///< Fixed anchor point
	int vertex;        ///< Attached vertex
} spring_t;            ///< Spring properties

typedef struct
{
	double m;     ///< Mass
	vector2_t v;  ///< Current Velocity
	vector2_t p;  ///< Current Position
} vertex_t;       ///< Vertex properties / state

/* Simulation parameters */
#define EDO_STEP 0.1
#define STEP_COUNT 200
//#define GRAVITY 9.8
//#define VERTEX_COUNT 4
//#define BAR_COUNT 4
//#define SPRING_COUNT 2
#define GRAVITY 9.8
#define VERTEX_COUNT 2
#define BAR_COUNT 1
#define SPRING_COUNT 1

const spring_t springs[SPRING_COUNT] = {
//	{.k = 0.1, .len = 1, .anchor = {.x = 0, .y = 0}, .vertex = 0},
//	{.k = 0.2, .len = 1, .anchor = {.x = 5, .y = 0}, .vertex = 1}
	{.k = 1.8, .len = 3, .anchor = {.x = 0, .y = 0}, .vertex = 0}
};

vertex_t vertexes[VERTEX_COUNT] = {
//	{.m = 0.10, .v = {.x = 0, .y = 0}, .p = {.x = 1, .y = 1}},
//	{.m = 0.05, .v = {.x = 0, .y = 0}, .p = {.x = 3, .y = 1}},
//	{.m = 0.05, .v = {.x = 0, .y = 0}, .p = {.x = 1, .y = 2}},
//	{.m = 0.02, .v = {.x = 0, .y = 0}, .p = {.x = 3, .y = 2}}
	{.m = 0.5, .v = {.x = 0, .y = 0}, .p = {.x = 0, .y = 2}},
	{.m = 0.5, .v = {.x = 0, .y = 0}, .p = {.x = 1, .y = 2}}
};

// Vertex pairs (i, i+1) ; i even
/*
const int bars[2 * BAR_COUNT] = {0, 1,
                                 1, 3,
                                 2, 3,
                                 0, 2};
*/

const int bars[2 * BAR_COUNT] = {0, 1};

/* Restriction forces calculation matrixes */

// Derivative of Jacobian restriction matrix
double dJ[BAR_COUNT * 2 * VERTEX_COUNT];

// Lambda (Lagrange coefs) calculations:  JWJT * Lambda = -JWf -dJv
double JT[BAR_COUNT * 2 * VERTEX_COUNT];
double JW[BAR_COUNT * 2 * VERTEX_COUNT];  // = - J * W
double JWJT[4 * BAR_COUNT * BAR_COUNT];   // = - J * W * J(transposed) = JW * JT

// Restriction forces
double fc[2 * VERTEX_COUNT];

// Vertex forces
static double forces[2 * VERTEX_COUNT];

/**
 * \brief Calculates a linear combination of 2 vectors
 * \param ca A coefficient
 * \param A Vector A
 * \param cb B coefficient
 * \param B Vector B
 * \param result Combined vector
 * \param N Element count
 */
void LinearCombination(double ca, double * A, double cb, double * B, double * result, int N)
{
	for (int i=0; i<N; i++)
		result[i] = ca * A[i] + cb * B[i];
}

/**
 * \brief Solve a linear system where a is an upper triangular matrix
 * \param n Variable count
 * \param a Coefficient matrix
 * \param b Independent vector
 * \param x Solution vector
 */
void BackSubstitution(int n, double *a, double *b, double *x)
{
	double s;
	//x[n-1] = a[(n-1)*(n-1)];
	for (int i=n-1; i>=0; i--)
	{
		s=0;
		for (int j=i+1; j<n; j++)
			s += a[i*n + j] * x[j];

		x[i] = (b[i] - s) / a[i*n + i];
	}
}

/**
 * \brief Gauss Elimination: Make a into an upper triangular matrix
 * \param n Variable count
 * \param a Coefficient matrix
 * \param b Independent vector
 */
void GaussElimination(int n, double *a, double *b)
{
	double f;
	int biggest;
	for (int j=0; j<n-1; j++)
	{
		biggest = j;
		for (int z=j+1; z<n; z++)
			if (fabs(a[z*n + j]) > fabs(a[biggest*n +j]))
				biggest = z;

		if (biggest != j)
		{
			// Swap Lines
			double tmp;
			tmp = b[j];
			b[j] = b[biggest];
			b[biggest] = tmp;

			for (int y=j; y<n; y++)
			{
				tmp = a[j*n + y];
				a[j*n + y] = a[biggest*n + y];
				a[biggest*n + y] = tmp;
			}
		}

		for (int i=j+1; i<n; i++)
		{
			f = a[i*n + j] / a[j*n + j];
			for (int k=j; k<n; k++)
				a[i*n + k] -= f * a[j*n + k];

			b[i] -= f * b[j];
		}
	}
}

/**
 * \brief Multiply 2 matrixes A * B
 * \param m A line count
 * \param n A col count (B line count)
 * \param l B col count
 * \param A A matrix
 * \param B B matrix
 * \param C result matrix
 */
void MatrixMultiply(int m, int n, int l, double *A, double *B, double *C)
{
	for (int a=0; a<m; a++)        // line
		for (int b=0; b<l; b++)    // column
		{
			C[a*l + b] = 0;

			for (int i=0; i<n; i++)
				C[a*l + b] += A[a*n + i] * B[i*l + b];
		}
}

/**
 * \brief Transposes a m x n matrix
 * \param m Lines
 * \param n Columns
 * \param A Matrix
 * \param X Transposed result
 */
void MatrixTranspose(int m, int n, double *A, double *X)
{
	for (int x=0; x<n; x++)
		for (int y=0; y<m; y++)
			X[x*m + y] = A[y*n + x];
}

/**
 * \brief Calculates the dot product of 2 vectors
 * \param A Vector A
 * \param B Vector B
 * \return Dot product
 */
double DotProduct(vector2_t A, vector2_t B)
{
	double result;
	result = A.x * B.x + A.y * B.y;
	return result;
}

/**
 * \brief Calculates values for the simulation matrixes JW, dJ, JT, JWJT
 */
void CalcMatrixes(void)
{
	// Jacobian restriction matrix
	double J[BAR_COUNT * 2 * VERTEX_COUNT];

	// Weight macro. Finds vertex weight. x is the vertex index * 2 (even) or
	// vertex index * 2 + 1 (odd)
	#define W(x) (vertexes[(int)((x) / 2)].m)
  
	// First vertex (index) in a bar
	#define vi(x) (bars[2 * (x)])

	// Second vertex (index) in a bar
	#define vj(x) (bars[2 * (x) + 1])

	// Calc J, JW and dJ
	memset(J, 0, sizeof(J));
	memset(dJ, 0, sizeof(dJ));

	for (int k=0; k<BAR_COUNT; k++)
	{
		vector2_t l, dl;
		double ml, mdl;    // l and dl modulus
		vector2_t dJinc;   // dJ element pair increment value
		double mlsq;       // ml squared
		double ldl;        // l . dl dot product
		l.x  = vertexes[vi(k)].p.x - vertexes[vj(k)].p.x;
		l.y  = vertexes[vi(k)].p.y - vertexes[vj(k)].p.y;
		dl.x = vertexes[vi(k)].v.x - vertexes[vj(k)].v.x;
		dl.y = vertexes[vi(k)].v.y - vertexes[vj(k)].v.y;

		ml = sqrt(l.x * l.x + l.y * l.y);
		mdl = sqrt(dl.x * dl.x + dl.y * dl.y);
		
		J[k * 2 * VERTEX_COUNT + 2 * vi(k)]     += l.x / ml;
		J[k * 2 * VERTEX_COUNT + 2 * vi(k) + 1] += l.y / ml;
		J[k * 2 * VERTEX_COUNT + 2 * vj(k)]     -= l.x / ml;
		J[k * 2 * VERTEX_COUNT + 2 * vj(k) + 1] -= l.y / ml;

		mlsq = ml * ml;
		ldl = DotProduct(l, dl);
		dJinc.x = (mlsq * dl.x - ldl * l.x) / (mlsq * ml);
		dJinc.y = (mlsq * dl.y - ldl * l.y) / (mlsq * ml);

		dJ[k * 2 * VERTEX_COUNT + 2 * vi(k)]     += dJinc.x;
		dJ[k * 2 * VERTEX_COUNT + 2 * vi(k) + 1] += dJinc.y;
		dJ[k * 2 * VERTEX_COUNT + 2 * vj(k)]     -= dJinc.x;
		dJ[k * 2 * VERTEX_COUNT + 2 * vj(k) + 1] -= dJinc.y;
	}

	for (int i=0; i<2*VERTEX_COUNT; i++)   // Column
		for (int k=0; k<BAR_COUNT; k++)    // Line
			JW[k*2*VERTEX_COUNT + i] = J[k*2*VERTEX_COUNT + i] / W(i);

	MatrixTranspose(BAR_COUNT, 2 * VERTEX_COUNT, J, JT);

	// Calc JWJT
	MatrixMultiply(BAR_COUNT, 2 * VERTEX_COUNT, BAR_COUNT, JW, JT, JWJT);
}

void CalcForces(void)
{
	double *f = forces;
	memset(forces, 0, sizeof(forces));

	// Springs
	for (int i=0; i<SPRING_COUNT; i++)
	{
		int vertex = springs[i].vertex;
		vector2_t sd; // spring direction
		double d;     // spring ends distance
		sd.x = vertexes[vertex].p.x - springs[i].anchor.x;
		sd.y = vertexes[vertex].p.y - springs[i].anchor.y;
		d = sqrt(sd.x * sd.x + sd.y * sd.y);
		
		f[2*vertex]     = - springs[i].k * (d - springs[i].len) * sd.x / d;
		f[2*vertex + 1] = - springs[i].k * (d - springs[i].len) * sd.y / d;
	}

	// Gravity
	for (int i=0; i<VERTEX_COUNT; i++)
		f[2*i + 1] -= GRAVITY * vertexes[i].m;

	// Calc lambda
	static double lambda[BAR_COUNT];
	static double b[BAR_COUNT];

	// b = - JW * f - dJ * v
	for (int i=0; i<BAR_COUNT; i++)
	{
		b[i] = 0;
		for (int k=0; k<2*VERTEX_COUNT; k+=2)
		{
			b[i] -= JW[i*2*VERTEX_COUNT + k]   * f[k]   + dJ[i*2*VERTEX_COUNT + k]   * vertexes[(int)(k/2)].v.x;
			b[i] -= JW[i*2*VERTEX_COUNT + k+1] * f[k+1] + dJ[i*2*VERTEX_COUNT + k+1] * vertexes[(int)(k/2)].v.y;
		}
	}

	// Solve system (change this to LU decomposition)
	static double tmp[BAR_COUNT * BAR_COUNT];
	memcpy(tmp, JWJT, sizeof(tmp));
	GaussElimination(BAR_COUNT, tmp, b);
	BackSubstitution(BAR_COUNT, tmp, b, lambda);
	
	// Calc restriction forces
	// fc = JT * lambda
	MatrixMultiply(2*VERTEX_COUNT, BAR_COUNT, 1, JT, lambda, fc);

	for (int i=0; i<2*VERTEX_COUNT; i++)
		forces[i] += fc[i];
}

/**
 * \brief First order acceleration calculation for all vertexes (will recalculate forces)
 * \param v0 Current velocities
 * \param result Will contain the calculated accelerations
 */
void CalcAccelerationsStep(double * v0, double * result)
{
	for (int i=0; i<VERTEX_COUNT; i++)
	{
		vertexes[i].v.x = v0[2*i];
		vertexes[i].v.y = v0[2*i +1];
	}

	CalcMatrixes();
	CalcForces();

	for (int i=0; i<VERTEX_COUNT; i++)
	{
		result[2*i]    = forces[2*i]    / vertexes[i].m;
		result[2*i +1] = forces[2*i +1] / vertexes[i].m;
	}
}

/**
 * \brief 4th order Runge-Kutta velocity calculation for all vertexes
 * \param dt Integration time interval
 * \param result Calculated velocities vector
 */
void CalcVelocities(double dt, double * result)
{
	static double  v0[2*VERTEX_COUNT];
	static double  s1[2*VERTEX_COUNT];
	static double  s2[2*VERTEX_COUNT];
	static double  s3[2*VERTEX_COUNT];
	static double  s4[2*VERTEX_COUNT];
	static double tmp[2*VERTEX_COUNT];

	for (int i=0; i<VERTEX_COUNT; i++)
	{
		v0[2*i]    = vertexes[i].v.x;
		v0[2*i +1] = vertexes[i].v.y;
	}

	CalcAccelerationsStep(v0, s1);
	LinearCombination(dt, s1, 0, s1, s1, 2*VERTEX_COUNT);

	LinearCombination(1, v0, 0.5, s1, tmp, 2*VERTEX_COUNT);
	CalcAccelerationsStep(tmp, s2);
	LinearCombination(dt, s2, 0, s2, s2, 2*VERTEX_COUNT);
	
	LinearCombination(1, v0, 0.5, s2, tmp, 2*VERTEX_COUNT);
	CalcAccelerationsStep(tmp, s3);
	LinearCombination(dt, s3, 0, s3, s3, 2*VERTEX_COUNT);
	
	LinearCombination(1, v0, 1, s3, tmp, 2*VERTEX_COUNT);
	CalcAccelerationsStep(tmp, s4);
//	LinearCombination(dt, s4, 0, s4, s4, 2*VERTEX_COUNT);  <-- Done in the sum below

	for (int i=0; i<2*VERTEX_COUNT; i++)
		result[i] = v0[i] + (s1[i] + 2.0 * s2[i] + 2.0 * s3[i] + dt * s4[i])/6.0;
}

/**
 * \brief First order velocity calculation for all vertexes (will recalculate velocities and forces)
 * \param dt Integration interval
 * \param p0 Current positions
 * \param result Will contain the calculated velocities
 */
void CalcVelocitiesStep(double dt, double * p0, double * result)
{
	for (int i=0; i<VERTEX_COUNT; i++)
	{
		vertexes[i].p.x = p0[2*i];
		vertexes[i].p.y = p0[2*i +1];
	}

	CalcVelocities(dt, result);
}

/**
 * \brief Calculates the next position for each vertex
 * \param dt Integration interval
 */
void CalcPositions(double dt)
{
	static double  p0[2*VERTEX_COUNT];
	static double  s1[2*VERTEX_COUNT];
	static double  s2[2*VERTEX_COUNT];
	static double  s3[2*VERTEX_COUNT];
	static double  s4[2*VERTEX_COUNT];
	static double tmp[2*VERTEX_COUNT];

	for (int i=0; i<VERTEX_COUNT; i++)
	{
		p0[2*i]    = vertexes[i].p.x;
		p0[2*i +1] = vertexes[i].p.y;
	}

	CalcVelocitiesStep(dt, p0, s1);
	LinearCombination(dt, s1, 0, s1, s1, 2*VERTEX_COUNT);

	LinearCombination(1, p0, 0.5, s1, tmp, 2*VERTEX_COUNT);
	CalcVelocitiesStep(dt, tmp, s2);
	LinearCombination(dt, s2, 0, s2, s2, 2*VERTEX_COUNT);
	
	LinearCombination(1, p0, 0.5, s2, tmp, 2*VERTEX_COUNT);
	CalcVelocitiesStep(dt, tmp, s3);
	LinearCombination(dt, s3, 0, s3, s3, 2*VERTEX_COUNT);
	
	LinearCombination(1, p0, 1, s3, tmp, 2*VERTEX_COUNT);
	CalcVelocitiesStep(dt, tmp, s4);
//	LinearCombination(dt, s4, 0, s4, s4, 2*VERTEX_COUNT);  <-- Done in the sum below

	for (int i=0; i<VERTEX_COUNT; i++)
	{
		vertexes[i].p.x = p0[2*i]    + (s1[2*i]    + 2.0 * s2[2*i]    + 2.0 * s3[2*i]    + dt * s4[2*i])/6.0;
		vertexes[i].p.y = p0[2*i +1] + (s1[2*i +1] + 2.0 * s2[2*i +1] + 2.0 * s3[2*i +1] + dt * s4[2*i +1])/6.0;
	}
}

void maple()
{
	char * colors[] = {"red", "blue", "green", "black", "brown", "grey", "maroon", "yellow", "orange", "purple", "purple", NULL};

	printf("with(plots):\n");
	for (int i=0; i<STEP_COUNT; i++)
	{
		CalcPositions(EDO_STEP);

		double d = sqrt((vertexes[0].p.x - vertexes[1].p.x)*(vertexes[0].p.x - vertexes[1].p.x) + (vertexes[0].p.y - vertexes[1].p.y)*(vertexes[0].p.y - vertexes[1].p.y));
		printf("#%0.5lf\n", d);

		printf("%c:=<", 'a'+i);
		for (int v=0; v<VERTEX_COUNT; v++)
		{
			printf("<%0.5lf,%0.5lf>", vertexes[v].p.x, vertexes[v].p.y);
			if (v != VERTEX_COUNT-1)
				printf("|");
		}

		printf(">;\n");
	}

	for (int v=0; v<VERTEX_COUNT; v++)
	{
		printf("p%d:=[", v);
		for (int i=0; i<STEP_COUNT; i++)
		{
			printf("[%c[1,%d],%c[2,%d]]", 'a'+i, v+1, 'a'+i, v+1);
			//printf("[%d,%c[2,%d]]", i, 'a'+i, v+1);
			if (i != 15)
				printf(",");
		}
		printf("];\n");
	}

	for (int v=0; v<VERTEX_COUNT; v++)
		if (colors[v] == NULL)
			break;
		else
			printf("g%d:=listplot(p%d,color=%s);\n", v, v, colors[v]);


	printf("display(");
	for (int v=0; v<VERTEX_COUNT; v++)
		if (v != VERTEX_COUNT-1)
			printf("g%d,", v);
		else
			printf("g%d);\n", v);
}

void matlab(void)
{
	for (int v=0; v<VERTEX_COUNT; v++)
		printf("p%d=zeros(2, %d);\n", v, STEP_COUNT);

	for (int i=0; i<STEP_COUNT; i++)
	{
		CalcPositions(EDO_STEP);
		for (int v=0; v<VERTEX_COUNT; v++)
			printf("p%d(1:2, %d)=[%0.5lf %0.5lf];\n", v, i+1, vertexes[v].p.x, vertexes[v].p.y);
	}
}

int main()
{
	matlab();
	return 0;
}
