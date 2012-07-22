// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#ifndef POINTS_H
#define POINTS_H

#include <math.h>

// ----------------------------------------------------------------
#define PERIODIC_BOUNDARY_CONDITIONS
//#undef PERIODIC_BOUNDARY_CONDITIONS

// ================================================================
// points_t data structure: contains lattice points, wormhole point, and a
// doubly linked list of cycle-information structure.  Each cycle-info
// structure contains the cycle length and a pointer to one site in the cycle.
// Each site contains a pointer to the cycle-info structure for the cycle it is
// in.
//
//  pcycinfo_list_head                                       pcycinfo_list_tail
//             pprev            pprev                    pprev
//   O  <-----------  O  <-----------  O  .... O  <-----------  O
//      -----------> ^|  ----------->             ----------->
//      pnext        ||  pnext                    pnext
//                   ||
//                   || psite
//         pcycinfo  ||
//                   |v
//               k=nz o----o----o----o----o----o
//                   /    /    /    /    /    /|
//             k=4  o----o----o----o----o----o o      o <== wormhole point w
//                 /    /    /    /    /    /|/|            (non-spatial)
//           k=3  o----o----o----o----o----o o o
//               /    /    /    /    /    /|/| |
//         k=2  o----o----o----o----o----o o o o
//             /    /    /    /    /    /|/| | |
//        k=1 o----o----o----o----o----o o o o o
//           /    /    /    /    /    /|/| | | |
// j=ny k=0 o----o----o----o----o----o o x======> Stored at each site x:
//          |    |    |    |    |    |/|/|/|/|/   * c:       Coords. (i,j,k)
// j=4      o----o----o----o----o----o o o o o    * pfwd:    Pointer to pi(x)
//          |    |    |    |    |    |/|/|/|/     * pbwd:    Pointer to piinv(x)
// j=3      o----o----o----o----o----o o o o      * fwd_dsq: |x-pi(x)|^2
//          |    |    |    |    |    |/|/|/       * fwd_d:   |x-pi(x)|
// j=2      o----o----o----o----o----o o o        * pcycinfo: Pointer to cycinfo
//          |    |    |    |    |    |/|/         * mark:    scratch space
// j=1      o----o----o----o----o----o o
//          |    |    |    |    |    |/
// j=0      o----o----o----o----o----o
//
//          i=0  i=1  i=2  i=3  ...  i=nx
// ================================================================

// ----------------------------------------------------------------
// Integer lattice coordinates.
typedef struct _coord_t {
	int    selfi;
	int    selfj;
	int    selfk;
} coord_t;

// Integer lattice dimensions.
typedef struct _dims_t {
	int    nx;
	int    ny;
	int    nz;
} dims_t;

// Point data structure.
typedef struct _point_t {
	coord_t  c;

	struct  _point_t* pfwd;
	struct  _point_t* pbwd;
	int      fwd_dsq;
	double   fwd_d;
	struct  _cycinfo_t* pcycinfo; // Ptr to cycle info for this site's cycle.
	int      mark;                // Scratch space for various routines.
} point_t;

// Cycle-information data structure.
typedef struct _cycinfo_t {
	struct _cycinfo_t* pprev;  // Backward pointer for doubly linked list.
	struct _cycinfo_t* pnext;  // Forward  pointer for doubly linked list.
	point_t* psite;  // Points to one of the lattice sites in the cycle.
	int      cyclen; // Number of lattice sites in the cycle.
} cycinfo_t;

// Points data structure.
typedef struct _points_t {
	dims_t  dims;
	int     N;     // nx * ny * nz

	// d-dimensional array of lattice sites.
	point_t*** lattice;
	point_t    wormhole;

	// Doubly linked list of cycle-info structs.
	cycinfo_t*   pcycinfo_list_head;
	cycinfo_t*   pcycinfo_list_tail;
} points_t;

// ----------------------------------------------------------------
// Does not populate the cycinfo pointers.
points_t* get_cubic_lattice_points(int L, int d);
void free_points (points_t* ppoints);
void print_point(char* pre, point_t* px, char* post);
void fprint_point(FILE* fp, char* pre, point_t* px, char* post);

// Populates the cycinfo pointers.
void set_up_cycinfo_list(points_t* ppoints);
// For rvrcm's use only:
void recompute_cycinfo_list(points_t* ppoints);

void free_cycinfo_list(points_t* ppoints);
void cycinfo_list_insert(points_t* ppoints, cycinfo_t* pcycinfo);
void cycinfo_list_remove(points_t* ppoints, cycinfo_t* pcycinfo);
void print_cycinfo_list(points_t* ppoints);
void sanity_check_cycinfo_list(points_t* ppoints);
void find_dxy_dyx_xoy(points_t* ppoints, point_t* x, point_t* y,
	int* pdxy, int* pdyx, int* px_y_same_cycle);
// This should be called after the proposed Metropolis changed is accepted, in
// the sense of modifying pfwd and pbwd pointers at the sites x and y.
void update_cycinfo(points_t* ppoints, point_t* x, point_t* y,
	int dxy, int dyx, int x_y_previously_in_same_cycle);

// ----------------------------------------------------------------
// Zero-centered floating-point modulus.
// E.g. if L = 10: -5 <= m <= 5;
//      if L = 11: -5 <= m <= 5.
// See the dissertation for the explanation.
// (Short answer: this is needed for symmetry of the dot product.)
static inline int zero_centered_mod(int x, int L)
{
	int half_L = L / 2;
	while (x >  half_L)
		x -= L;
	while (x < -half_L)
		x += L;
	return x;
}

static inline int get_distance_squared(
	int xi, int xj, int xk,
	int yi, int yj, int yk,
	int Lx, int Ly, int Lz)
{
	int dsq;
	int tmp;

#ifdef PERIODIC_BOUNDARY_CONDITIONS
	tmp = zero_centered_mod(xi - yi, Lx);
	dsq  = tmp*tmp;

	tmp = zero_centered_mod(xj - yj, Ly);
	dsq += tmp*tmp;

	tmp = zero_centered_mod(xk - yk, Lz);
	dsq += tmp*tmp;
#else
	tmp = xi - yi;
	dsq  = tmp*tmp;

	tmp = xj - yj;
	dsq += tmp*tmp;

	tmp = xk - yk;
	dsq += tmp*tmp;
#endif

	return dsq;
}

// ----------------------------------------------------------------
static inline int get_distance_squaredp(
	coord_t* px, coord_t* py, dims_t* pdims)
{
	return get_distance_squared(
		px->selfi, px->selfj, px->selfk,
		py->selfi, py->selfj, py->selfk,
		pdims->nx, pdims->ny, pdims->nz);
}

// ----------------------------------------------------------------
static inline int get_Zd_distance_squared(
	int xi, int xj, int xk,
	int yi, int yj, int yk)
{
	int dsq;
	int tmp;

	tmp  = xi - yi;
	dsq  = tmp*tmp;

	tmp  = xj - yj;
	dsq += tmp*tmp;

	tmp  = xk - yk;
	dsq += tmp*tmp;

	return dsq;
}

// ----------------------------------------------------------------
static inline int get_Zd_distance_squaredp(
	coord_t* px, coord_t* py)
{
	return get_Zd_distance_squared(
		px->selfi, px->selfj, px->selfk,
		py->selfi, py->selfj, py->selfk);
}

// ----------------------------------------------------------------
static inline double get_distancep(
	coord_t* px, coord_t* py, dims_t* pdims)
{
	return sqrt(get_distance_squaredp(px, py, pdims));
}

// ----------------------------------------------------------------
static inline double get_dot(coord_t* pu, coord_t* pv)
{
	return
		pu->selfi * pv->selfi +
		pu->selfj * pv->selfj +
		pu->selfk * pv->selfk;
}

// ----------------------------------------------------------------
static inline double get_normsq(coord_t* pu)
{
	return
		pu->selfi * pu->selfi +
		pu->selfj * pu->selfj +
		pu->selfk * pu->selfk;
}

// ----------------------------------------------------------------
static inline double get_norm(coord_t* pu)
{
	return sqrt(
		pu->selfi * pu->selfi +
		pu->selfj * pu->selfj +
		pu->selfk * pu->selfk);
}

// ----------------------------------------------------------------
static inline coord_t get_coord_diff(coord_t* pu, coord_t* pv, dims_t* pdims)
{
	coord_t w;
#ifdef PERIODIC_BOUNDARY_CONDITIONS
	w.selfi = zero_centered_mod(pu->selfi - pv->selfi, pdims->nx);
	w.selfj = zero_centered_mod(pu->selfj - pv->selfj, pdims->ny);
	w.selfk = zero_centered_mod(pu->selfk - pv->selfk, pdims->nz);
#else
	w.selfi = pu->selfi - pv->selfi;
	w.selfj = pu->selfj - pv->selfj;
	w.selfk = pu->selfk - pv->selfk;
#endif
	return w;
}

// ----------------------------------------------------------------
static inline coord_t get_coord_sum(coord_t* pu, coord_t* pv, dims_t* pdims)
{
	coord_t w;
#ifdef PERIODIC_BOUNDARY_CONDITIONS
	w.selfi = zero_centered_mod(pu->selfi + pv->selfi, pdims->nx);
	w.selfj = zero_centered_mod(pu->selfj + pv->selfj, pdims->ny);
	w.selfk = zero_centered_mod(pu->selfk + pv->selfk, pdims->nz);
#else
	w.selfi = pu->selfi + pv->selfi;
	w.selfj = pu->selfj + pv->selfj;
	w.selfk = pu->selfk + pv->selfk;
#endif
	return w;
}

// ----------------------------------------------------------------
double get_jump_dot(point_t* px, point_t* py, points_t* ppoints);
// Returns 0 if the angle is undefined, i.e. if |u| or |v| is zero.
// Else, returns 1 and sets *pcos_angle = cos(theta).
int get_jump_cos_angle(point_t* px, point_t* py, points_t* ppoints,
	double* pcos_angle);

double get_mean_jump_length(points_t* ppoints, double* pmaxjumplen);

#endif // POINTS_H
