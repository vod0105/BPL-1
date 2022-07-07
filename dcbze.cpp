ï»¿
#define TESTMODE
#include "BezierCurveReconstruction.h"
#include <malloc.h>
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <string>
typedef Point2* BezierCurve;
#define MAXPOINTS	1000		/* The most points you can have */
using namespace std;

void DrawBezierCurve(int n, BezierCurve bcurve)
{
	int i;
	printf("\nKet qua 4 diem tao duong cong: \n");
	for (i = 0; i <= n; i++) {
		printf("%0.2f, %0.2f\n", bcurve[i].x, bcurve[i].y);
	}
	ofstream output;
	output.open("D:/DUT/DoAnCoSo/newbezier/newbezier/output.txt", std::fstream::out| std::ofstream::app);
	output<< "Day la 4 diem tao dung nen duong cong:\n";
		for (i = 0; i <= n; i++)
			output << "x="<< bcurve[i].x << ";y=" <<  bcurve[i].y<< endl;
     output.close();
}

/*
 *  B0, B1, B2, B3 :
 *	Bezier multipliers
 */
static double B0(double	u)
{
	double tmp = 1.0 - u;
	return (tmp * tmp * tmp);
}


static double B1(double	u)
{
	double tmp = 1.0 - u;
	return (3 * u * (tmp * tmp));
}

static double B2(double	u)
{
	double tmp = 1.0 - u;
	return (3 * u * u * tmp);
}

static double B3(double	u)
{
	return (u * u * u);
}

static Vector2 V2AddII(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x + b.x;  c.y = a.y + b.y;
	return (c);
}
static Vector2 V2ScaleIII(Vector2 v, double s) {
	Vector2 result;
	result.x = v.x * s; result.y = v.y * s;
	return (result);
}

static Vector2 V2SubII(Vector2 a, Vector2 b) {
	Vector2 c;
	c.x = a.x - b.x; c.y = a.y - b.y;
	return (c);
}

/* return the distance between two points */
double V2DistanceBetween2Points(Point2* a, Point2* b) {
	double dx = a->x - b->x;
	double dy = a->y - b->y;
	return(sqrt((dx * dx) + (dy * dy)));
}

/* negates the input vector and returns it */
Vector2* V2Negate(Vector2* v) {
	v->x = -v->x;  v->y = -v->y;
	return(v);
}

double V2SquaredLength(Vector2* a) {
	return((a->x * a->x) + (a->y * a->y));
}


/* returns length of input vector */
double V2Length(Vector2* a) {
	return(sqrt(V2SquaredLength(a)));
}

/* normalizes the input vector and returns it */
Vector2* V2Normalize(Vector2* v) {
	double len = V2Length(v);
	if (len != 0.0) { v->x /= len;  v->y /= len; }
	return(v);
}

/* scales the input vector to the new length and returns it */
Vector2* V2Scale(Vector2* v, double newlen) {
	double len = V2Length(v);
	if (len != 0.0) { v->x *= newlen / len;   v->y *= newlen / len; }
	return(v);
}

/* return vector sum c = a+b */
Vector2* V2Add(Vector2* a, Vector2* b, Vector2* c) {
	c->x = a->x + b->x;  c->y = a->y + b->y;
	return(c);
}

/* return vector difference c = a-b */
Vector2* V2Sub(Vector2* a, Vector2* b, Vector2* c) {
	c->x = a->x - b->x;  c->y = a->y - b->y;
	return(c);
}

/* return the dot product of vectors a and b */
double V2Dot(Vector2* a, Vector2* b) {
	return((a->x * b->x) + (a->y * b->y));
}




/*
 *  ChordLengthParameterize :
 *  Assign parameter values to digitized points
 *  using relative distances between points.
 */
static double* ChordLengthParameterize(Point2* d, int first, int last)
{
	int     i;
	double* u;         /*  Parameterization        */

	u = (double*)malloc((unsigned)(last - first + 1) * sizeof(double));  //uuuuuuuuu
	

	u[0] = 0.0;
	for (i = first + 1; i <= last; i++) {
		u[i - first] = u[i - first - 1] +
			V2DistanceBetween2Points(&d[i], &d[i - 1]);
	}

	for (i = first + 1; i <= last; i++) {
		u[i - first] = u[i - first] / u[last - first];
	}

	return(u);
}



/*
 *  GenerateBezier :
 *  Use least-squares method to find Bezier control points for region.
 *
 */
static BezierCurve  GenerateBezier(
	Point2* d,			/*  Array of digitized points	*/
	int		first, int last,		/*  Indices defining region	*/
	double* uPrime,		/*  Parameter values for region */
	Vector2	tHat1, Vector2 tHat2)	/*  Unit tangents at endpoints	*/
{
	int 	i;
	Vector2 	A[MAXPOINTS][2];	/* Precomputed rhs for eqn	*/
	int 	nPts;			/* Number of pts in sub-curve */
	double 	C[2][2];			/* Matrix C		*/
	double 	X[2];			/* Matrix X			*/
	double 	det_C0_C1,		/* Determinants of matrices	*/
		det_C0_X,
		det_X_C1;
	double 	alpha_l,		/* Alpha values, left and right	*/
		alpha_r;
	Vector2 	tmp;			/* Utility variable		*/
	BezierCurve	bezCurve;	/* RETURN bezier curve ctl pts	*/

	bezCurve = (Point2*)malloc(4 * sizeof(Point2));
	nPts = last - first + 1;


	/* Compute the A's	*/
	for (i = 0; i < nPts; i++) {
		Vector2		v1, v2;
		v1 = tHat1;
		v2 = tHat2;
		V2Scale(&v1, B1(uPrime[i]));
		V2Scale(&v2, B2(uPrime[i]));
		A[i][0] = v1;
		A[i][1] = v2;
	}

	/* Create the C and X matrices	*/
	C[0][0] = 0.0;
	C[0][1] = 0.0;
	C[1][0] = 0.0;
	C[1][1] = 0.0;
	X[0] = 0.0;
	X[1] = 0.0;

	for (i = 0; i < nPts; i++) {
		C[0][0] += V2Dot(&A[i][0], &A[i][0]);
		C[0][1] += V2Dot(&A[i][0], &A[i][1]);
		/*					C[1][0] += V2Dot(&A[i][0], &A[i][1]);*/
		C[1][0] = C[0][1];
		C[1][1] += V2Dot(&A[i][1], &A[i][1]);

		tmp = V2SubII(d[first + i],
			V2AddII(
				V2ScaleIII(d[first], B0(uPrime[i])),
				V2AddII(
					V2ScaleIII(d[first], B1(uPrime[i])),
					V2AddII(
						V2ScaleIII(d[last], B2(uPrime[i])),
						V2ScaleIII(d[last], B3(uPrime[i]))))));


		X[0] += V2Dot(&A[i][0], &tmp);
		X[1] += V2Dot(&A[i][1], &tmp);
	}

	/* Compute the determinants of C and X	*/
	det_C0_C1 = C[0][0] * C[1][1] - C[1][0] * C[0][1];
	det_C0_X = C[0][0] * X[1] - C[0][1] * X[0];
	det_X_C1 = X[0] * C[1][1] - X[1] * C[0][1];

	/* Finally, derive alpha values	*/
	if (det_C0_C1 == 0.0) {
		det_C0_C1 = (C[0][0] * C[1][1]) * 10e-12;
	}
	alpha_l = det_X_C1 / det_C0_C1;
	alpha_r = det_C0_X / det_C0_C1;


	/*  If alpha negative, use the Wu/Barsky heuristic (see text) */
	/* (if alpha is 0, you get coincident control points that lead to
	 * divide by zero in any subsequent NewtonRaphsonRootFind() call. */
	if (alpha_l < 1.0e-6 || alpha_r < 1.0e-6) {
		double	dist = V2DistanceBetween2Points(&d[last], &d[first]) /
			3.0;

		bezCurve[0] = d[first];
		bezCurve[3] = d[last];
		V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
		V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
		return (bezCurve);
	}

	/*  First and last control points of the Bezier curve are */
	/*  positioned exactly at the first and last data points */
	/*  Control points 1 and 2 are positioned an alpha distance out */
	/*  on the tangent vectors, left and right, respectively */
	bezCurve[0] = d[first];
	bezCurve[3] = d[last];
	V2Add(&bezCurve[0], V2Scale(&tHat1, alpha_l), &bezCurve[1]);
	V2Add(&bezCurve[3], V2Scale(&tHat2, alpha_r), &bezCurve[2]);
	return (bezCurve);
}

/*
 *  Bezier :
 *  	Evaluate a Bezier curve at a particular parameter value
 *
 */
static Point2 BezierII(
	int		degree,		/* The degree of the bezier curve	*/
	Point2* V,		/* Array of control points		*/
	double 	t)		/* Parametric value to find point for	*/
{
	int 	i, j;
	Point2 	Q;	        /* Point on curve at parameter t	*/
	Point2* Vtemp;		/* Local copy of control points		*/

	/* Copy array	*/
	Vtemp = (Point2*)malloc((unsigned)((degree + 1)
		* sizeof(Point2)));
	for (i = 0; i <= degree; i++) {
		Vtemp[i] = V[i];
	}

	/* Triangle computation	*/
	for (i = 1; i <= degree; i++) {
		for (j = 0; j <= degree - i; j++) {
			Vtemp[j].x = (1.0 - t) * Vtemp[j].x + t * Vtemp[j + 1].x;
			Vtemp[j].y = (1.0 - t) * Vtemp[j].y + t * Vtemp[j + 1].y;
		}
	}

	Q = Vtemp[0];
	free((void*)Vtemp);
	return Q;
}



/*
 *  ComputeMaxError :
 *	Find the maximum squared distance of digitized points
 *	to fitted curve.
*/
static double ComputeMaxError(
	Point2* d,			/*  Array of digitized points	*/
	int		first, int last,		/*  Indices defining region	*/
	BezierCurve	bezCurve,		/*  Fitted Bezier curve		*/
	double* u,			/*  Parameterization of points	*/
	int* splitPoint)		/*  Point of maximum error	*/
{
	int		i;
	double	maxDist;		/*  Maximum error		*/
	double	dist;		/*  Current error		*/
	Point2	P;			/*  Point on curve		*/
	Vector2	v;			/*  Vector from point to curve	*/

	*splitPoint = (last - first + 1) / 2;
	maxDist = 0.0;
	for (i = first + 1; i < last; i++) {
		P = BezierII(3, bezCurve, u[i - first]);
		v = V2SubII(P, d[i]);
		dist = V2SquaredLength(&v);
		if (dist >= maxDist) {
			maxDist = dist;
			*splitPoint = i;
		}
	}
	return (maxDist);
}


/*
 *  NewtonRaphsonRootFind : cÃ´ng cá»¥ newton-raph  cho viá»‡c tÃ¬m kiáº¿m cÄƒn sá»‘
 *	Use Newton-Raphson iteration to find better root. : sá»­ dá»¥ng cÃ´ng cá»¥ Ä‘Ã³ láº·p Ä‘i láº·p láº·p láº¡i
	Ä‘á»ƒ tÃ¬m cÄƒn sá»‘ tá»‘t hÆ¡n
 */
static double NewtonRaphsonRootFind(
	BezierCurve	Q,			/*  Current fitted curve: Ä‘Æ°á»ng cong thÃ­ch há»£p hiá»‡n táº¡i	*/
	Point2 		P,		/*  Digitized point	: Ä‘iá»ƒm sá»‘ 2	*/
	double 		u)		/*  Parameter value for "P" : giÃ¡ trá»‹ tham sá»‘ cho P	*/
{
	double 		numerator, denominator;
	Point2 		Q1[3], Q2[2];	/*  Q' and Q'':	Q vÃ  Q'		*/
	Point2		Q_u, Q1_u, Q2_u; /*u evaluated at: u Ä‘Æ°á»£c xac Ä‘á»‹nh táº¡i Ä‘iá»ƒm Q, Q', & Q''	*/
	double 		uPrime;		/*  Improved u			*/
	int 		i;

	/* Compute Q(u)	: tÃ­nh toÃ¡n Q(u)*/
	Q_u = BezierII(3, Q, u);

	/* Generate control vertices for Q' : táº¡o cÃ¡c Ä‘á»‰nh Ä‘iá»u khiá»ƒn cho Q'	*/
	for (i = 0; i <= 2; i++) {
		Q1[i].x = (Q[i + 1].x - Q[i].x) * 3.0;
		Q1[i].y = (Q[i + 1].y - Q[i].y) * 3.0;
	}

	/* Generate control vertices for Q'' :táº¡o cÃ¡c Ä‘á»‰nh Ä‘iá»u khiá»ƒn cho Q'' */
	for (i = 0; i <= 1; i++) {
		Q2[i].x = (Q1[i + 1].x - Q1[i].x) * 2.0;
		Q2[i].y = (Q1[i + 1].y - Q1[i].y) * 2.0;
	}

	/* Compute Q'(u) and Q''(u)	: tÃ­nh toÃ¡n Q'(U) vÃ  U'' (u)*/
	Q1_u = BezierII(2, Q1, u);
	Q2_u = BezierII(1, Q2, u);

	/* Compute f(u)/f'(u) : tÃ­nh toÃ¡n F(u) vÃ  / f(U)' */
	numerator = (Q_u.x - P.x) * (Q1_u.x) + (Q_u.y - P.y) * (Q1_u.y);
	denominator = (Q1_u.x) * (Q1_u.x) + (Q1_u.y) * (Q1_u.y) +
		(Q_u.x - P.x) * (Q2_u.x) + (Q_u.y - P.y) * (Q2_u.y);

	/* u = u - f(u)/f'(u) */
	uPrime = u - (numerator / denominator);
	return (uPrime);
}



/*
 *  Reparameterize: xÃ¡c Ä‘á»‹nh láº¡i tham sá»‘
 *	Given set of points and their parameterization, try to find: táº­p há»£p cÃ¡c Ä‘iá»ƒm vÃ  tham sá»‘ cá»§a chÃºng, cá»‘ gáº¯ng tÃ¬m má»™t tham sá»‘ chÃ­nh xÃ¡c hÆ¡n
 *   a better parameterization.
 *
 */
static double* Reparameterize(
	Point2* d,			/*  Array of digitized points : máº£ng cÃ¡c Ä‘iá»ƒm sá»‘ hoÃ¡	*/
	int		first, int last,		/*  Indices defining region	: chá»‰ sá»‘ xÃ¡c Ä‘á»‹nh vÃ¹ng*/
	double* u,			/*  Current parameter values : giÃ¡ trá»‹ tham sá»‘ hiá»‡n táº¡i	*/
	BezierCurve	bezCurve)	/*  Current fitted curve : Ä‘Æ°á»ng cong thÃ­ch há»£p hiá»‡n táº¡i	*/
{
	int 	nPts = last - first + 1;
	int 	i;
	double* uPrime;		/*  New parameter values: giÃ¡ trá»‹ tham sÃ³ má»›i 	*/

	uPrime = (double*)malloc(nPts * sizeof(double));
	for (i = first; i <= last; i++) {
		uPrime[i - first] = NewtonRaphsonRootFind(bezCurve, d[i], u[i -
			first]);
	}
	return (uPrime);
}






/*
 * ComputeLeftTangent, ComputeRightTangent, ComputeCenterTangent : tÃ­nh toÃ¡n tiáº¿p tuyáº¿n trÃ¡i, tÃ­nh toÃ¡n tiáº¿p tuyáº¿n pháº£i
 tÃ­nh toÃ¡n tiáº¿p tuyáº¿n giá»¯a
 *Approximate unit tangents at endpoints and "center" of digitized curve
 CÃ¡c tiáº¿p tuyáº¿n Ä‘Æ¡n vá»‹ gáº§n Ä‘Ãºng táº¡i cÃ¡c Ä‘iá»ƒm cuá»‘i vÃ  "tÃ¢m" cá»§a Ä‘Æ°á»ng cong sá»‘ hÃ³a
 */
static Vector2 ComputeLeftTangent(
	Point2* d,		/*  Digitized points : Ä‘iá»ƒm tham sá»‘ hoÃ¡*/
	int		end)		/*  Index to "left" end of region :Chá»‰ má»¥c Ä‘áº¿n "bÃªn trÃ¡i" vÃ¹ng cuá»‘i */
{
	Vector2	tHat1;
	tHat1 = V2SubII(d[end + 1], d[end]);
	tHat1 = *V2Normalize(&tHat1);
	return tHat1;
}

static Vector2 ComputeRightTangent(
	Point2* d,		/*  Digitized points	: Ä‘iá»ƒm sá»‘ hoÃ¡	*/
	int		end)		/*  Index to "right" end of region : chá»‰ má»¥c Ä‘áº¿n "bÃªn pháº£i"  vÃ¹ng cuá»‘i */
{
	Vector2	tHat2;
	tHat2 = V2SubII(d[end - 1], d[end]);
	tHat2 = *V2Normalize(&tHat2);
	return tHat2;
}
static Vector2 ComputeCenterTangent(
	Point2* d,		/*  Digitized points	: Ä‘iá»ƒm sá»‘ hoÃ¡		*/
	int		center)		/*  Index to point inside region	Chá»‰ má»¥c Ä‘á»ƒ chá»‰ trong khu vá»±c */
{
	Vector2	V1, V2, tHatCenter;

	V1 = V2SubII(d[center - 1], d[center]);
	V2 = V2SubII(d[center], d[center + 1]);
	tHatCenter.x = (V1.x + V2.x) / 2.0;
	tHatCenter.y = (V1.y + V2.y) / 2.0;
	tHatCenter = *V2Normalize(&tHatCenter);
	return tHatCenter;
}
/*
 *  FitCubic :
 *  	Fit a Bezier curve to a (sub)set of digitized points :Khá»›p Ä‘Æ°á»ng cong Bezier vá»›i táº­p há»£p (phá»¥) cÃ¡c Ä‘iá»ƒm Ä‘Æ°á»£c sá»‘ hÃ³a
 */
static void FitCubic(
	Point2* d,			/*  Array of digitized points */
	int		first, int last,	/* Indices of first and last pts in region */
	Vector2	tHat1, Vector2 tHat2,	/* Unit tangent vectors at endpoints */
	double	error)		/*  User-defined error squared	   */
{
	BezierCurve	bezCurve; /*Control points of fitted Bezier curve*/
	double* u;		/*  Parameter values for point :GiÃ¡ trá»‹ tham sá»‘ cho Ä‘iá»ƒm */
	double* uPrime;	/*  Improved parameter values : Cáº£i thiá»‡n giÃ¡ trá»‹ tham sá»‘ */
	double	maxError;	/*  Maximum fitting error : Lá»—i sáº¯p Ä‘Äƒt tá»‘i Ä‘a	 */
	int		splitPoint;	/*  Point to split point set at	: Äiá»ƒm Ä‘á»ƒ chia Ä‘iá»ƒm Ä‘Æ°á»£c Ä‘áº·t táº¡i */
	int		nPts;		/*  Number of points in subset : Sá»‘ Ä‘iá»ƒm trong táº­p há»£p con */
	double	iterationError; /*Error below which you try iterating :Lá»—i bÃªn dÆ°á»›i mÃ  báº¡n thá»­ láº·p láº¡i */
	int		maxIterations = 4; /*  Max times to try iterating :Sá»‘ láº§n tá»‘i Ä‘a Ä‘á»ƒ thá»­ láº·p */
	Vector2	tHatCenter;   	/* Unit tangent vector at splitPoint :ÄÆ¡n vá»‹ vectÆ¡ tiáº¿p tuyáº¿n táº¡i Ä‘iá»ƒm phÃ¢n chia*/
	int		i;

	iterationError = error * error;
	nPts = last - first + 1;

	/*  Use heuristic if region only has two points in it :Sá»­ dá»¥ng heuristic náº¿u khu vá»±c chá»‰ cÃ³ hai Ä‘iá»ƒm trong Ä‘Ã³ */
	if (nPts == 2) {
		double dist = V2DistanceBetween2Points(&d[last], &d[first]) / 3.0;

		bezCurve = (Point2*)malloc(4 * sizeof(Point2));
		bezCurve[0] = d[first];
		bezCurve[3] = d[last];
		V2Add(&bezCurve[0], V2Scale(&tHat1, dist), &bezCurve[1]);
		V2Add(&bezCurve[3], V2Scale(&tHat2, dist), &bezCurve[2]);
		//	DrawBezierCurve(3, bezCurve);
		free((void*)bezCurve);
		return;
	}

	/*  Parameterize points, and attempt to fit curve :Tham sá»‘ hÃ³a cÃ¡c Ä‘iá»ƒm vÃ  cá»‘ gáº¯ng khá»›p Ä‘Æ°á»ng cong */
	u = ChordLengthParameterize(d, first, last);
	bezCurve = GenerateBezier(d, first, last, u, tHat1, tHat2);

	/*  Find max deviation of points to fitted curve :TÃ¬m Ä‘á»™ lá»‡ch tá»‘i Ä‘a cá»§a cÃ¡c Ä‘iá»ƒm Ä‘á»‘i vá»›i Ä‘Æ°á»ng cong Ä‘Æ°á»£c trang bá»‹ */
	maxError = ComputeMaxError(d, first, last, bezCurve, u, &splitPoint);
	if (maxError < error) {
		DrawBezierCurve(3, bezCurve);
		free((void*)u);
		free((void*)bezCurve);
		return;
	}


	/*  If error not too large, try some reparameterization :Náº¿u lá»—i khÃ´ng quÃ¡ lá»›n, hÃ£y thá»­ má»™t sá»‘ tham sá»‘ láº¡i */
	/*  and iteration  :vÃ  láº·p Ä‘i láº·p láº¡i*/
	if (maxError < iterationError) {
		for (i = 0; i < maxIterations; i++) {
			uPrime = Reparameterize(d, first, last, u, bezCurve);
			bezCurve = GenerateBezier(d, first, last, uPrime, tHat1, tHat2);
			maxError = ComputeMaxError(d, first, last,
				bezCurve, uPrime, &splitPoint);
			if (maxError < error) {
				DrawBezierCurve(3, bezCurve);
				free((void*)u);
				free((void*)bezCurve);
				return;
			}
			free((void*)u);
			u = uPrime;
		}
	}

	/* Fitting failed -- split at max error point and fit recursively: sáº¯p Ä‘áº·t tháº¥t báº¡i - phÃ¢n chia táº¡i Ä‘iá»ƒm lá»—i tá»‘i Ä‘a vÃ  khá»›p Ä‘á»‡ quy*/
	free((void*)u);
	free((void*)bezCurve);
	tHatCenter = ComputeCenterTangent(d, splitPoint);
	FitCubic(d, first, splitPoint, tHat1, tHatCenter, error);
	V2Negate(&tHatCenter);
	FitCubic(d, splitPoint, last, tHatCenter, tHat2, error);
}
/*
 *  FitCurve :
 *  	Fit a Bezier curve to a set of digitized points :Khá»›p Ä‘Æ°á»ng cong Bezier vá»›i táº­p há»£p cÃ¡c Ä‘iá»ƒm Ä‘Æ°á»£c sá»‘ hÃ³a
 */
void FitCurve(
	Point2 d[],			/*  Array of digitized points : máº£ng cá»§a Ä‘iá»ƒm 	*/
	int		nPts,		/*  Number of digitized points	: Sá»‘ Ä‘iá»ƒm Ä‘Æ°á»£c sá»‘ hÃ³a*/
	double	error)		/*  User-defined error squared	: BÃ¬nh phÆ°Æ¡ng lá»—i do ngÆ°á»i dÃ¹ng xÃ¡c Ä‘á»‹nh*/
{
	Vector2	tHat1, tHat2;	/*  Unit tangent vectors at endpoints :CÃ¡c vectÆ¡ tiáº¿p tuyáº¿n Ä‘Æ¡n vá»‹ táº¡i cÃ¡c Ä‘iá»ƒm cuá»‘i */
		tHat1 = ComputeLeftTangent(d, 0);
		tHat2 = ComputeRightTangent(d, nPts - 1);
		FitCubic(d, 0, nPts - 1, tHat1, tHat2, error);

}
#ifdef TESTMODE
/*
 *  main:
 *	Example of how to use the curve-fitting code.  Given an array
 *   of points and a tolerance (squared error between points and
 *	fitted curve), the algorithm will generate a piecewise
 *	cubic Bezier representation that approximates the points.
 VÃ­ dá»¥ vá» cÃ¡ch sá»­ dá»¥ng mÃ£ khá»›p Ä‘Æ°á»ng cong. ÄÆ°a ra má»™t máº£ng cÃ¡c Ä‘iá»ƒm vÃ  dung sai (sai sá»‘ bÃ¬nh phÆ°Æ¡ng giá»¯a
cÃ¡c Ä‘iá»ƒm vÃ  Ä‘Æ°á»ng cong Ä‘Æ°á»£c trang bá»‹), thuáº­t toÃ¡n sáº½ táº¡o ra má»™t biá»ƒu diá»…n Bezier khá»‘i láº­p phÆ°Æ¡ng gáº§n Ä‘Ãºng vá»›i cÃ¡c Ä‘iá»ƒm.
 *	When a cubic is generated, the routine "DrawBezierCurve"
 *	is called, which outputs the Bezier curve just created
 *	(arguments are the degree and the control points, respectively).
 *	Users will have to implement this function themselves
 *   ascii output, etc.
 Khi má»™t khá»‘i Ä‘Æ°á»£c táº¡o, thÆ°á»ng trÃ¬nh "DrawBezierCurve" Ä‘Æ°á»£c gá»i, sáº½ táº¡o ra Ä‘Æ°á»ng cong Bezier vá»«a táº¡o
 (cÃ¡c Ä‘á»‘i sá»‘ láº§n lÆ°á»£t lÃ  Ä‘á»™ vÃ  cÃ¡c Ä‘iá»ƒm kiá»ƒm soÃ¡t). NgÆ°á»i dÃ¹ng sáº½ pháº£i tá»± thá»±c hiá»‡n chá»©c nÄƒng nÃ y Ä‘áº§u ra ascii, v.v.
 *
 */
int x;
int numlines = 0;
int count= 0;
int check() {
	ifstream file;
	char c;
	int numchars = 0;
	int num = 0;
	file.open("D:/DUT/DoAnCoSo/newbezier/newbezier/input.txt");
	//file >> x;
	file.get(c);
	while (file) {
		while (file && c != '\n') {
			numchars = numchars + 1;
			num = num + 1;
			file.get(c);
		}
			numlines = numlines + 1;
			file.get(c);
	}
	numlines = numlines - 1;
	file.close();
	return numlines;
}
void input(Point2Struct  d[]) {
	ifstream file;
	char c;
	int numchars = 0;
	file.open("D:/DUT/DoAnCoSo/newbezier/newbezier/input.txt");
	//for (string line; getline(file, line); ){
		for (int i = 0; i < numlines; i++) {
				file >> d[i].x;
				::count++;
				file>> d[i].y;
				::count++;
		}
	//}
	file.close();
}
void display(Point2Struct  d[]) {
	/*cout << "Tap diem thuoc duong cong bezier: ";
	for (int i = 0; i < numlines; i++) {
		d[i].show();
	}*/
	cout << "\n";
	ofstream output;
	output.open("output.txt");
	/*output << "\n Nhom 18N15-N12";
	output << "\n Thanh vien : -Truong Cong Hien Nhan";
	output << "\n              -Tang Ba Hong Phuc";
	output << "\n Lop : 18TCLC_DT3";*/
	output << "\n De tai cua nhom laï¿½: Noi suy duong cong Bezier";
	output << "\n Tap diem thuoc duong cong benzier: " << endl;
	for (int i = 0; i < numlines; i++) {
		output << "x=" << d[i].x << ";y=" << d[i].y << endl;
	}
		output.close();
}
int main()
{
  
    Point2Struct  d[MAXPOINTS];/*    {	  Digitized points  : táº­p há»£p Ä‘iá»ƒm*/
	int choose = 0; 
	double	error = 4.0;/*  Squared error : lá»—i bÃ¬nh phÆ°Æ¡ng*/
	cout << "\n-----------------Noi suy duong cong bezier-------------";
	cout << "\n enter 1. lay du lieu tu file input co san";
	cout << "\n enter 2. nhap du lieu tu ban phim";
	cout << "\n enter 0. thoat" << endl;
	cin >> choose;
		switch (choose)
		{
		case 0:
			break;
		case 1:
			check();
			input(d);
			display(d);
			FitCurve(d, ::numlines, error);		/*  Fit the Bezier curves :khá»›p vá»›i Ä‘Æ°á»ng cong bezier */
			break;
		case 2:
			int n;
			cout << "\n Nhap so diem thuoc duong cong bezier n=";
			cin >> n;
			for (int i = 0; i < n; i++) {
				cout << "\nx= "; cin >> d[i].x;
				cout << "y= "; cin >> d[i].y;
			}
			FitCurve(d, n, error);
			break;
		}
		getchar();
	return 0;
}
#endif						 /* TESTMODE */