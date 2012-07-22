// ================================================================
// Copyright (c) 2009 John Kerl.
// kerl.john.r@gmail.com
// Conditions for use of this file are expressed in the file ./README.txt.
//================================================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "points.h"
#include "pmt.h"
#include "util.h"
#include "wormidx.h"
#include "rcmrand.h"
#include "checks.h"

static void check_d(char * func_name, int d);

// ----------------------------------------------------------------
static void populate_points(points_t* ppoints, int nx, int ny, int nz, int d)
{
	int i, j, k;
	point_t* px;
	point_t* pw = &ppoints->wormhole;

	if (nx < 1) {
		fprintf(stderr,
			"get_cubic_lattice_points:  nx must be >= 1; got %d.\n", nx);
		exit(1);
	}
	if (ny < 1) {
		fprintf(stderr,
			"get_cubic_lattice_points:  ny must be >= 1; got %d.\n", ny);
		exit(1);
	}
	if (nz < 1) {
		fprintf(stderr,
			"get_cubic_lattice_points:  nz must be >= 1; got %d.\n", nz);
		exit(1);
	}

	ppoints->dims.nx = nx;
	ppoints->dims.ny = ny;
	ppoints->dims.nz = nz;
	ppoints->N = nx * ny * nz;

	ppoints->lattice = (point_t ***)malloc_or_die(nx * sizeof(point_t**));
	for (i = 0; i < nx; i++) {
		ppoints->lattice[i] = (point_t **)malloc_or_die(ny * sizeof(point_t*));
		for (j = 0; j < ny; j++) {
			ppoints->lattice[i][j]=(point_t*)malloc_or_die(nz*sizeof(point_t));
			for (k = 0; k < nz; k++) {
				px = &ppoints->lattice[i][j][k];
				px->c.selfi  = i;
				px->c.selfj  = j;
				px->c.selfk  = k;
				px->pfwd     = px;
				px->pbwd     = px;
				px->fwd_dsq  = 0;
				px->fwd_d    = 0.0;
				px->pcycinfo = 0;
				px->mark     = 0;
			}
		}
	}

	pw->c.selfi  = WORMIDX;
	pw->c.selfj  = WORMIDX;
	pw->c.selfk  = WORMIDX;
	pw->pfwd     = pw;
	pw->pbwd     = pw;
	pw->fwd_dsq  = 0;
	pw->fwd_d    = 0.0;
	pw->pcycinfo = 0;
	pw->mark     = 0;

	// These will be populated by set_up_cycinfo_list().
	ppoints->pcycinfo_list_head = 0;
	ppoints->pcycinfo_list_tail = 0;
}

// ----------------------------------------------------------------
points_t* get_cubic_lattice_points(int L, int d)
{
	points_t* ppoints = (points_t*)malloc_or_die(sizeof(points_t));

	check_d("get_cubic_lattice_points", d);
	if (d == 1) populate_points(ppoints, L, 1, 1, d);
	if (d == 2) populate_points(ppoints, L, L, 1, d);
	if (d == 3) populate_points(ppoints, L, L, L, d);

	return ppoints;
}

// ----------------------------------------------------------------
void free_points(points_t* ppoints)
{
	int i, j;

	free_cycinfo_list(ppoints);

	for (i = 0; i < ppoints->dims.nx; i++) {
		for (j = 0; j < ppoints->dims.ny; j++) {
			free(ppoints->lattice[i][j]);
		}
		free(ppoints->lattice[i]);
	}
	free(ppoints->lattice);

	ppoints->dims.nx = 0;
	ppoints->dims.ny = 0;
	ppoints->dims.nz = 0;
	ppoints->N       = 0;
	ppoints->lattice = 0;
	free(ppoints);
}

// ----------------------------------------------------------------
void free_cycinfo_list(points_t* ppoints)
{
	cycinfo_t* pcycinfo = ppoints->pcycinfo_list_head;

	while (pcycinfo != 0) {
		cycinfo_t* ptemp = pcycinfo;
		pcycinfo = pcycinfo->pnext;
		free(ptemp);
	}

	ppoints->pcycinfo_list_head = 0;
	ppoints->pcycinfo_list_tail = 0;
}

// ----------------------------------------------------------------
void print_point(char* pre, point_t* px, char* post)
{
	fprint_point(stdout, pre, px, post);
}

// ----------------------------------------------------------------
void fprint_point(FILE* fp, char* pre, point_t* px, char* post)
{
	fprintf(fp, "%s", pre);
	if (px->c.selfi == WORMIDX)
		fprintf(fp, "*,*,*");
	else
		fprintf(fp, "%d,%d,%d", px->c.selfi, px->c.selfj, px->c.selfk);
	fprintf(fp, "%s", post);
}

// ----------------------------------------------------------------
static void check_d(char * func_name, int d)
{
	if ((d < 1) || (d > 3)) {
		fprintf(stderr, "%s:  d out of range 1-3 (got %d).\n", func_name, d);
		exit(1);
	}
}

// ----------------------------------------------------------------
// u.v = |u| |v| cos(theta)
// Returns 0 if cos(theta) is undefined, i.e. if |u| or |v| is zero.
// Else, returns 1 and sets *pcos_angle = cos(theta).
int get_jump_cos_angle(point_t* px, point_t* py, points_t* ppoints,
	double* pcos_angle)
{
	point_t* pix = px->pfwd;
	point_t* piy = py->pfwd;
	coord_t u = get_coord_diff(&pix->c, &px->c, &ppoints->dims);
	coord_t v = get_coord_diff(&piy->c, &py->c, &ppoints->dims);
	double normu = get_norm(&u);
	double normv = get_norm(&v);
	double udotv = get_dot(&u, &v);

	*pcos_angle = 0.0;
	if (normu == 0.0)
		return 0;
	if (normv == 0.0)
		return 0;
	*pcos_angle = udotv / (normu*normv);
	return 1;
}

// ----------------------------------------------------------------
double get_jump_dot(point_t* px, point_t* py, points_t* ppoints)
{
	point_t* pix = px->pfwd;
	point_t* piy = py->pfwd;
	coord_t u = get_coord_diff(&pix->c, &px->c, &ppoints->dims);
	coord_t v = get_coord_diff(&piy->c, &py->c, &ppoints->dims);
	return get_dot(&u, &v);
}

// ----------------------------------------------------------------
double get_mean_jump_length(points_t* ppoints, double* pmaxjumplen)
{
	double mean_jump_length = 0.0;
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int i, j, k;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				double jumplen = ppoints->lattice[i][j][k].fwd_d;
				mean_jump_length += jumplen;
				if (jumplen > *pmaxjumplen)
					*pmaxjumplen = jumplen;
			}
		}
	}

	return mean_jump_length / ppoints->N;
}

// ----------------------------------------------------------------
// Allocates a list of cycinfo_t structs, one for each permutation cycle.
// Sets the point_t pcycinfo_list_head to the head of that list.
// Sets each lattice site's pcycinfo pointer.

//   head                                                      tail
//    O <---------> O <---------> O <---------> ... <---------> O
//    |             |             |                             |
//    |             |             |                             |
//    |             |             |                             |
//    v             v             |                             |
//    o-->o         o             v                             v
//  /      \       / \            o                       o<--->o
// o        o     /   \.
// ^        |    o<----o
// |        v
// o        o
//   \    /
//    o<-o
//
// An 8-cycle    a 3-cycle    a 1-cycle         ...     a 2-cycle
//
//

void set_up_cycinfo_list(points_t* ppoints)
{
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int i, j, k;
	point_t* px     = 0;
	point_t* pw     = &ppoints->wormhole;
	point_t* pcurr  = 0;
	point_t* pnext  = 0;
	cycinfo_t* pcycinfo = 0;

	if (ppoints->pcycinfo_list_head != 0) {
		fprintf(stderr,
		"set_up_cycinfo_list:  list has already been populated.  Exiting.\n");
		exit(1);
	}

	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			for (k = 0; k < nz; k++)
				ppoints->lattice[i][j][k].mark = 0;
	pw->mark = 0;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				px = &ppoints->lattice[i][j][k];
				if (px->mark)
					continue; // The site has already been visited.

				// Allocate a cycle structure.
				pcycinfo = (cycinfo_t*)malloc_or_die(sizeof(cycinfo_t));

				// The cycle-info structure points one lattice site in the
				// cycle, and all sites in the cycle point to the cycle-info
				// structure.
				pcycinfo->psite  = px;
				pcycinfo->pprev  = 0;
				pcycinfo->pnext  = 0;
				pcycinfo->cyclen = 0;

				pcurr = px;
				pnext = 0;

				while (1) {
					pcycinfo->cyclen++;
					pcurr->pcycinfo = pcycinfo;
					pcurr->mark = 1;
					pnext = pcurr->pfwd;
					pcurr = pnext;
					if (pnext == px)
						break;
				}

				// Insert this cycle-info structure into the list of
				// cycle-info structures.
				cycinfo_list_insert(ppoints, pcycinfo);
			}
		}
	}

	if (!pw->mark) {
		// The wormhole point has not already been visited.

		// Allocate a cycle structure.
		pcycinfo = (cycinfo_t*)malloc_or_die(sizeof(cycinfo_t));

		// The cycle-info structure points one lattice site in the
		// cycle, and all sites in the cycle point to the cycle-info
		// structure.
		pcycinfo->psite  = pw;
		pcycinfo->pprev  = 0;
		pcycinfo->pnext  = 0;
		pcycinfo->cyclen = 0;

		pcurr = pw;
		pnext = 0;

		while (1) {
			pcycinfo->cyclen++;
			pcurr->pcycinfo = pcycinfo;
			pcurr->mark = 1;
			pnext = pcurr->pfwd;
			pcurr = pnext;
			if (pnext == pw)
				break;
		}

		// Insert this cycle-info structure into the list of
		// cycle-info structures.
		cycinfo_list_insert(ppoints, pcycinfo);
	}
}

// ----------------------------------------------------------------
void recompute_cycinfo_list(points_t* ppoints)
{
	free_cycinfo_list(ppoints);
	set_up_cycinfo_list(ppoints);
}

// ----------------------------------------------------------------
void cycinfo_list_insert(points_t* ppoints, cycinfo_t* pcycinfo)
{
	if (ppoints->pcycinfo_list_tail == 0) {
		// First insertion into an empty list.
		ppoints->pcycinfo_list_head = pcycinfo;
		ppoints->pcycinfo_list_tail = pcycinfo;
		pcycinfo->pprev = 0;
		pcycinfo->pnext = 0;
	}
	else {
		// Append to the end of an existing list.
		ppoints->pcycinfo_list_tail->pnext = pcycinfo;
		pcycinfo->pprev = ppoints->pcycinfo_list_tail;
		pcycinfo->pnext = 0; // End of list.
		ppoints->pcycinfo_list_tail = pcycinfo;
	}
}

// ----------------------------------------------------------------
void cycinfo_list_remove(points_t* ppoints, cycinfo_t* pcycinfo)
{
	if (ppoints->pcycinfo_list_head == 0) {
		// The list is already empty and we are asked to remove something
		// else.  This is certainly an error.
		fprintf(stderr,
			"cycinfo_list_remove:  list is already empty.  Exiting.\n");
		exit(1);
	}

	if (pcycinfo->pprev == 0) {
		// Removal from head of list
		if (pcycinfo->pnext == 0) {
			// Removal of the only item from the list.  This would be a coding
			// error somewhere:  There is always at least one cycle present in
			// a permutation on N sites.
			fprintf(stderr,
				"cycinfo_list_remove:  cycle list will be empty.  Exiting.\n");
			exit(1);
		}

		// Before:  X <--> B <--> C <--> ...
		//          ^
		//          |
		//          head
		//
		// After:          B <--> C <--> ...
		//                 ^
		//                 |
		//                 head

		pcycinfo->pnext->pprev = 0; // New list head
		ppoints->pcycinfo_list_head = pcycinfo->pnext;
	}
	else if (pcycinfo->pnext == 0) {
		// Removal from tail of list

		// Before:  ... <--> A <--> B <--> X
		//                                 ^
		//                                 |
		//                                 tail
		//
		// After:   ... <--> A <--> B
		//                          ^
		//                          |
		//                          tail

		pcycinfo->pprev->pnext = 0; // New list tail
		ppoints->pcycinfo_list_tail = pcycinfo->pprev;
	}
	else {
		// Removal from middle of list

		// Before:  ... <--> A <--> X <--> B <--> ...
		//
		// After:   ... <--> A <---------> B <--> ...
		pcycinfo->pprev->pnext = pcycinfo->pnext;
		pcycinfo->pnext->pprev = pcycinfo->pprev;

	}

	pcycinfo->cyclen = -1;
	pcycinfo->psite  =  0;
	pcycinfo->pprev  =  0;
	pcycinfo->pnext  =  0;
	free(pcycinfo);
}

// ----------------------------------------------------------------
void print_cycinfo_list(points_t* ppoints)
{
	cycinfo_t* pci = 0;
	point_t* px, * pw;
	int num_cycles = 0;
	int i, j, k;

	printf("\n");
	printf("CYCLE-INFORMATION LIST DUMP:\n");

	printf("Lattice sites:\n");
	for (i = 0; i < ppoints->dims.nx; i++) {
		for (j = 0; j < ppoints->dims.ny; j++) {
			for (k = 0; k < ppoints->dims.nz; k++) {
				px = &ppoints->lattice[i][j][k];
				print_point("  (", px, ") ");
				printf("bwd 0x%p at 0x%p fwd 0x%p "
					"pci 0x%p ccl %d\n",
					px->pbwd, px, px->pfwd,
					px->pcycinfo, px->pcycinfo->cyclen);
			}
		}
	}
	pw = &ppoints->wormhole;
	print_point("  (", pw, ") ");
	printf("bwd 0x%p at 0x%p fwd 0x%p pci 0x%p ccl %d\n",
			pw->pbwd, pw, pw->pfwd,
			pw->pcycinfo, pw->pcycinfo->cyclen);

	if (ppoints->pcycinfo_list_head == 0) {
		printf("  (empty)\n");
		return;
	}
	printf("List head 0x%p tail 0x%p\n",
		ppoints->pcycinfo_list_head,
		ppoints->pcycinfo_list_tail);
	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		num_cycles++;
		px = pci->psite;
		printf("  bwd 0x%p at 0x%p fwd 0x%p len %d ste 0x%p ",
			pci->pprev, pci, pci->pnext, pci->cyclen, px);
		print_point("(", px, ")\n");
	}
	printf("Number of cycles: %d.\n", num_cycles);
	printf("\n");
}

// ----------------------------------------------------------------
// The cycinfo sanity check consists of the following tests:

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void check_cycinfo_partition(points_t* ppoints)
{
	int nx = ppoints->dims.nx;
	int ny = ppoints->dims.ny;
	int nz = ppoints->dims.nz;
	int i, j, k;
	point_t*   px  = 0, * pnext = 0;
	point_t*   pw  = &ppoints->wormhole;
	cycinfo_t* pci = 0;

	// Mark all lattice sites as unvisited.
	for (i = 0; i < nx; i++)
		for (j = 0; j < ny; j++)
			for (k = 0; k < nz; k++)
				ppoints->lattice[i][j][k].mark = 0;
	pw->mark = 0;

	// Follow all cycles in the cycle-info structure.
	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		px = pci->psite;
		px->mark = 1;
		pnext = px->pfwd;
		while (pnext != 0) {
			if (pnext == px)
				break;
			pnext->mark++;
			pnext = pnext->pfwd;
		}
	}

	// Make sure all lattice sites have been visited exactly once.
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			for (k = 0; k < nz; k++) {
				px = &ppoints->lattice[i][j][k];
				if (px->mark != 1) {
					fprintf(stderr, "check_cycinfo_partition: site %d,%d,%d "
						"# visits = %d, not 1.  Exiting.\n", i, j, k, px->mark);
					exit(1);
				}
			}
		}
	}

	if (pw->mark != 1) {
		fprintf(stderr, "check_cycinfo_partition: site *,*,* "
			"# visits = %d, not 1.  Exiting.\n", pw->mark);
		exit(1);
	}
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void check_cycinfo_cyclens(points_t* ppoints)
{
	point_t*   px  = 0, * pnext = 0;
	cycinfo_t* pci = 0;

	// Follow all cycles in the cycle-info structure.
	// Compute actual cycle length.
	// Check against the cycinfo cycle length.
	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		int actual_cyclen = 1;
		px = pci->psite;
		pnext = px->pfwd;
		while (pnext != 0) {
			if (pnext == px)
				break;
			actual_cyclen++;
			pnext = pnext->pfwd;
		}
		if (actual_cyclen != pci->cyclen) {
			fprintf(stderr,
				"check_cycinfo_cyclens:  actual cyclen %d at site ",
				actual_cyclen);
			fprint_point(stderr, "(", px, ")");
			fprintf(stderr, " != cached cyclen %d.\n",
				pci->cyclen);
			exit(1);
		}
	}
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Check to see if all sites in the cycle point to the same cycinfo struct.
void check_lattice_cycinfo_pointers(points_t* ppoints)
{
	point_t*   px  = 0, * pnext = 0;
	cycinfo_t* pci = 0;

	for (pci = ppoints->pcycinfo_list_head; pci != 0; pci = pci->pnext) {
		px = pci->psite;
		if (px->pcycinfo != pci) {
			fprintf(stderr,
				"check_lattice_cycinfo_pointers: pointer mismatch at site");
			fprint_point(stderr, "(", px, ").\n");
			exit(1);
		}
		pnext = px->pfwd;
		while (pnext != 0) {
			if (pnext == px)
				break;
			if (px->pcycinfo != pci) {
				fprintf(stderr,
					"check_lattice_cycinfo_pointers: pointer mismatch at site ");
				fprint_point(stderr, "(", px, ").\n");
				exit(1);
			}
			pnext = pnext->pfwd;
		}
	}
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void sanity_check_cycinfo_list(points_t* ppoints)
{
#ifdef CHECK_CYCINFO // in checks.h
	check_cycinfo_partition(ppoints);
	check_cycinfo_cyclens(ppoints);
	check_lattice_cycinfo_pointers(ppoints);
#endif
}

// ----------------------------------------------------------------
//int find_dxy_slowly(point_t* x, point_t* y)
//{
//	int dxy = 0;
//	point_t* z = x;
//	while (z != y) {
//		z = z->pfwd;
//		dxy++;
//	}
//	return dxy;
//}

// ----------------------------------------------------------------
void find_dxy_dyx_xoy(points_t* ppoints, point_t* x, point_t* y,
	int* pdxy, int* pdyx, int* px_y_same_cycle)
{
	if (x->pcycinfo == y->pcycinfo) {
		//int slow_dxy = find_dxy_slowly(x, y);
		//int slow_dyx = find_dxy_slowly(y, x);
		int d = 0;
		point_t* zfwd = x;
		point_t* zbwd = x;
		*px_y_same_cycle = 1;

		// Empirically, in these simulations it is found that almost always y
		// is a few hops from x or x is a few hops from y.  So, search forward
		// and backward simultaneously to find the smallest.  Use dxy + dyx =
		// ellx = elly to find the other.

		while (1) {
			if (zfwd == y) {
				// We found dxy.  dyx must be xcyclen - dxy.
				*pdxy = d;
				*pdyx = x->pcycinfo->cyclen - d;
				break;
			}
			if (zbwd == y) {
				// We found dyx.  dxy must be xcyclen - dyx.
				*pdyx = d;
				*pdxy = x->pcycinfo->cyclen - d;
				break;
			}
			zfwd = zfwd->pfwd;
			zbwd = zbwd->pbwd;
			d++;
		}

		//if ((slow_dxy != *pdxy) || (slow_dyx != *pdyx)) {
		//	fprintf(stderr, "dxy b0rk!\n");
		//	exit(1);
		//}
	}
	else {
		*px_y_same_cycle = 0;
		*pdxy = -1;
		*pdyx = -1;
	}
}

// ----------------------------------------------------------------
// update_cycinfo
//
// This routine is called by the Metropolis-update logic immediately after
// making a permutation swap.  The permutation pointers for x and y reflect the
// new configuration; x->pcycinfo->cyclen, y->pcycinfo->cyclen, dxy, dyx, and
// x<-->y all reflect the old configuration.  Our task is to update the cycinfo
// structs.
//
// We are simply applying the change which was proposed in the
// Metropolis-update logic, and evaluated in get_Delta_rell in energy.c.  See
// get_Delta_rell for more information.
//
// In the figures, O denotes a cycinfo struct; x and y denote the points x and
// y; o denotes other sites.
//
// ----------------------------------------------------------------
// SPLIT
//      BEFORE                                      AFTER
// ...    O  ...                        ...    O <---------> O  ...
//        |                                   |              |
//        |                                   |              |
//        |                                   |              |
//        v                                   v              v
//        x------>y                           x              y
//       /       /                           / \.
//      /      /                            /   \.
//     o<----o                             o<----o
// ... a 4-cycle  ...                  ... a 3-cycle    a 1-cycle ...
//
// ----------------------------------------------------------------
// MERGE
//            BEFORE                        AFTER
// ...    O <---------> O  ...         ...    O  ...
//        |             |                     |
//        |             |                     |
//        |             |                     |
//        v             v                     v
//        x             y                     x------>y
//       / \                                 /       /
//      /   \                               /      /
//     o<----o                             o<----o
// ... a 3-cycle    a 1-cycle ...      ... a 4-cycle  ...
//
// ----------------------------------------------------------------
// Again, note that this routine is called *after* the permutation pointers
// themselves have been changed, but the cached cycle information does not yet
// reflect that change.  (It is the task of this routine to make that update.)

void update_cycinfo(points_t* ppoints, point_t* x, point_t* y,
	int dxy, int dyx, int x_y_previously_in_same_cycle)
{
	point_t* z = 0;

	if (x_y_previously_in_same_cycle) { // SPLIT CASE
		// This is a split: x and y were previously in the same cycle; now they
		// are not.
		int new_x_cyclen = dyx;
		int new_y_cyclen = dxy;
		cycinfo_t* pnewcycinfo = 0;

		// To minimize the number of computations, find the shorter new cycle.
		// Suppose x's new cycle is longer than y's.  (If not, swap them
		// to make this so.)
		if (new_x_cyclen < new_y_cyclen) {
			point_t* ptemp;
			int itemp;
			ptemp = x; x = y; y = ptemp;
			itemp=new_x_cyclen; new_x_cyclen=new_y_cyclen; new_y_cyclen=itemp;
		}

		// All points z in y's new cycle must now point to a new cycinfo
		// struct.  The cycle has been split, so follow from the new pi(y)
		// (which was pi(x) before the merge) around to and including y.
		//
		// Also, the cycinfo structs for both new cycles need to have their
		// cyclens updated.
		pnewcycinfo = (cycinfo_t*)malloc_or_die(sizeof(cycinfo_t));
		pnewcycinfo->pprev  = 0;
		pnewcycinfo->pnext  = 0;
		pnewcycinfo->psite  = y;
		pnewcycinfo->cyclen = new_y_cyclen;

		y->pcycinfo = pnewcycinfo;
		z = y->pfwd;
		while (z != y) {
			z->pcycinfo = pnewcycinfo;
			z = z->pfwd;
			if (z == z->pfwd) {
				printf("b1rk!\n");
				exit(1);
			}
		}
		x->pcycinfo->cyclen = new_x_cyclen;
		x->pcycinfo->psite  = x;

		// Then, the cycinfo struct for the y cycle must be added to
		// the cycinfo list.
		cycinfo_list_insert(ppoints, pnewcycinfo);

	}
	else { // MERGE CASE
		// This is a merge: x and y were not in the same cycle; now they are.
		// New x & y cyclen = old x cyclen + old y cyclen.
		cycinfo_t* poldcycinfo = 0;

		// To minimize the number of computations, find the shorter old cycle.
		// Suppose x's old cycle is longer than y's.  (If not, swap them
		// to make this so.)
		if (x->pcycinfo->cyclen < y->pcycinfo->cyclen) {
			point_t* ptemp;
			int itemp;
			ptemp = x; x = y; y = ptemp;
			itemp = dxy; dxy = dyx; dyx = itemp;
		}
		poldcycinfo = y->pcycinfo;

		// The cycinfo struct for the merged cycle needs to have its cyclen
		// updated.
		x->pcycinfo->cyclen += y->pcycinfo->cyclen;

		// All points z in y's old cycle must have their cycinfo structs
		// now point to x's cycinfo struct.  The two cycles have been
		// merged, though, so follow from the new pi(x) (which was pi(y)
		// before the merge) around to and including y.

		z = x->pfwd;
		while (z != y) {
			z->pcycinfo = x->pcycinfo;
			z = z->pfwd;
			if (z == z->pfwd) {
				printf("b2rk!\n");
				exit(1);
			}
		}
		y->pcycinfo = x->pcycinfo;

		// Then, the cycinfo struct for the old y cycle must be removed from
		// the cycinfo list and freed.
		cycinfo_list_remove(ppoints, poldcycinfo);

	}
}
