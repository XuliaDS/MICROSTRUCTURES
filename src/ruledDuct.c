/*
 *      EGADS: Electronic Geometry Aircraft Design System
 *
 *             EGADS Tessellation using wv with Quad Tessellation
 *
 *      Copyright 2011-2018, Massachusetts Institute of Technology
 *      Licensed under The GNU Lesser General Public License, version 2.1
 *      See http://www.opensource.org/licenses/lgpl-2.1.php
 *
 */

#include <math.h>
#include <string.h>
#include <time.h>
#ifdef  WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <winsock2.h>
#endif
#include "egads.h"

#define DEBUG

/* IRIT includes */
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/user_lib.h"

extern int EG_spline2dAppx( ego context,     int    endc,   /*@null@*/ const double *uknot,
			    /*@null@*/ const double *vknot,  /*@null@*/ const int    *vdata,
			    /*@null@*/ const double *wesT,   /*@null@*/ const double *easT,
			    /*@null@*/ const double *south,  /*@null@*/       double *snor,
			    /*@null@*/ const double *north,  /*@null@*/       double *nnor, int imax, int jmax,
			    const double *xyz,  double tol,   ego *esurf );

static int makeIRITsrf(ego eobj, CagdSrfStruct **surf);

static int
EG_getCommonEdge ( int ne1, ego *srf1Edges, int ne2, ego *srf2Edges, int *idx ) {
  int stat, i, j, per;
  double ts[2], xyz[18];
  *idx = -1;
  for ( i = 0 ; i < ne1; i++ ) {
      EG_getRange( srf1Edges[i], ts, &per );
      ts[0] = ( ts[0] + ts[1] ) * 0.5;
      stat = EG_evaluate ( srf1Edges[i], &ts[0], xyz);
      for ( j = 0; j < ne2; j++ ) {
	  stat =  EG_inTopology ( srf2Edges[j], xyz );
	  if ( stat == EGADS_SUCCESS) {
	      *idx = i;
	      break;
	  }
      }
      if ( *idx != -1 ) return EGADS_SUCCESS;
  }
  return EGADS_GEOMERR;
}

//makes them unclamped ( multiplicity at 0,1 = 1) because EG_2dApproximate already does the clamping
// nK corresponds to the knot sequence + clamping ( num Knots + 2 * splines_degreeree)
static int
setKnotSeq ( int n, double data[], int nK, int splines_degree, double *knots ) {
  int k = 0, j = 0, i = 0, nCP;
  double *t = NULL, arc = 0, diam = 0.0, alpha, d;
  nCP = nK - ( splines_degree + 1 );
  t   = (double * ) EG_alloc ( n * sizeof ( double ));
  if ( t == NULL ) return EGADS_MALLOC;

  for (k = 1; k < n; k++) {
      arc = 0.0;
      for(j = 0; j < 3; ++j)
	arc += pow ( ( data[3 * k + j ] - data[ 3 * ( k - 1 ) + j ]), 2);
      diam += pow ( arc, 0.5 );
  }
  t[0] = 0.0; t[n - 1] = 1.0;
  for ( k = 1; k < n - 1; k++) {
      arc = 0.0;
      for ( j = 0; j < 3; ++j)
	arc += pow ( ( data[3 * k + j ] - data[ 3 * ( k - 1 ) + j ]), 2);
      t[k] = t[k-1] + pow ( arc, 0.5 )/diam;
  }

  d = (double)n/(double)( nCP - splines_degree );
  k = 0;
  knots[ k++ ] = t[0];
  for ( j = 1; j < nCP - splines_degree; ++j) {
      i            = floor((double) j * d );
      alpha        =       (double) j * d - (double) i ;
      knots[ k++ ] = ( 1.0 - alpha ) * t[i-1] + alpha * t[i];
  }
  knots[ k++ ] = t[ n - 1 ];
  EG_free (t ) ;
  if ( k != nK - 2 * splines_degree ) {
      printf(" I messed up unclamping functino\n");
      return EGADS_INDEXERR;
  }
  return EGADS_SUCCESS;
}


/* We are assuming that curves form a loop: uv[0] -> uv[1] -> uv[2] -> uv[3] 
 * and are oriented accordingly. E.g., position uv[0][i] is aligned (vertically)
 * with position uv[3][n - i] (n = total points in u dir)
 */
static void
EG_setGridTFI ( int n[], double **uv, double *grid ){
  int    i, j, k, d ;
  double  yj, xi;
  if ( n[0] != n[2] || n[1] != n[3] ) {
      printf(" need same number of points in opposite directions\n");
      return;
  }
  for (k = j = 0; j < n[1]; j++) {
      yj = ((double) (j+1)) / ((double) n[1]);
      grid[2 * k     ] = uv[3][2* ( n[1] - 1 - j )    ];
      grid[2 * k + 1 ] = uv[3][2* ( n[1] - 1 - j ) + 1];
      k++;
      for (i = 1; i < n[0] - 1; i++) {
	  if ( j == 0 ) {
	      grid[2 * k    ] = uv[0][2*i    ];
	      grid[2 * k + 1] = uv[0][2*i + 1];
	      k++;
	  }
	  else if  ( j == n[1] - 1 ) {
	      grid[2 * k    ] = uv[2][2 * ( n[0] - 1 - i )    ];
	      grid[2 * k + 1] = uv[2][2 * ( n[0] - 1 - i ) + 1];
	      k++;
	  }
	  else {
	      xi          = ((double) (i+1)) / ((double) n[0]);
	      for (d = 0 ; d < 2; d++) {
		  grid[2 * k + d] = (1.0 - yj )  * uv[0][2 * i  + d  ] +
		      (      yj ) * uv[2][2      * ( n[0] - 1 - i ) + d  ] +
		      (1.0 - xi ) * uv[1][2      *              j   + d  ] +
		      (      xi ) * uv[3][2      * ( n[1] - 1 - j ) + d  ] -
		      (1.0 - xi ) * ( 1.0 - yj ) * uv[1][d] -
		      (      xi ) * (       yj ) * uv[3][d] -
		      (      xi ) * ( 1.0 - yj ) * uv[0][d] -
		      (1.0 - xi ) * (       yj ) * uv[2][d];
	      }
	      k++;
	  }
      }
      grid[2 * k     ] = uv[1][2* j    ];
      grid[2 * k + 1 ] = uv[1][2* j + 1];
      k++;
  }
}

int main(int argc, char *argv[])
{
  int        i, j, k, l, ibody, nbody, stat, oclass, mtype,  iface, nface, nloop, iedge, nedge, periodic, dims[4];
  int        *senses, **esenses = NULL, ce[2][4], idsF[4], idsO[4], nT[3], uvwID[3], splines_degree, n, e1, e2, e3, e4;
  double     data[18], trange[6][2 * 4], arc[ 2 ][ 4], **XYZ = NULL,  t[4],  **knots = NULL, *surf_grid = NULL, *uvTFI = NULL;
  double     dt, tEval[8], **tuvs = NULL, paux[4][500];
  const char *OCCrev;
  char       fname[100], fnItd[100], *ErrStr;
  ego        model, geom, *bodies , *dum = NULL, eref,  context, *faces ,
      *surf = NULL, *loops , **edges = NULL, *trimmed_surf = NULL;
#ifdef DEBUG
  FILE *fil;
#endif
  CagdSrfStruct **IRITsurf = NULL;
  TrivTVStruct   *TVMap    = NULL;
  /*  Creates a grid around four edges and test the TFI routine */
  if ( argc != 2 && argc != 5 ) {
      printf("\n Usage: vApproxTrim filename\n\n");
      return 1;
  }
  nT[0] = 10; nT[1] = 10; nT[2] = 5;
  if ( argc == 5 ) {
      nT[0] = atoi ( argv[2] );
      nT[1] = atoi ( argv[3] );
      nT[2] = atoi ( argv[4] );
  }
  splines_degree = 3;
  printf( " ATTENTION !!!!!!!!!!  ONLY ACCEPTING CUBIC SPLINES \n");
  /* look at EGADS revision */
  EG_revision(&i, &j, &OCCrev);
  printf("\n Using EGADS %2d.%02d with %s\n\n", i, j, OCCrev);

  /* initialize */
  printf(" EG_open           = %d\n", EG_open(&context));
  stat = 0;
  if(argv[1])  stat = EG_loadModel(context, 0, argv[1], &model);
  else return EGADS_EMPTY;
  printf(" EG_loadModel      = %d\n", stat);

  /* get all bodies */
  stat = EG_getTopology(model, &geom, &oclass, &mtype, NULL, &nbody,
			&bodies, &senses);
  if (stat != EGADS_SUCCESS) {
      printf(" EG_getTopology = %d\n", stat);
      return 1;
  }
  printf(" EG_getTopology:     nBodies = %d\n", nbody);
  /* fill our structure a body at at time */
  for (ibody = 0; ibody < nbody; ibody++) {
      mtype = 0;
      stat  = EG_getTopology(bodies[ibody], &geom, &oclass,
			     &mtype, NULL, &j, &dum, &senses);
      if ( stat != EGADS_SUCCESS) goto cleanup;
      stat  = EG_getBodyTopos(bodies[ibody], NULL, FACE, &nface, &faces);
      if ( stat != EGADS_SUCCESS || nface != 6 ) {
	  stat = EGADS_INDEXERR;
	  goto cleanup;
      }
      edges         = (ego ** ) EG_alloc ( 6 * sizeof ( ego * ) ) ;
      esenses       = (int ** ) EG_alloc ( 6 * sizeof ( int * ) ) ;
      surf          = (ego *  ) EG_alloc ( 6 * sizeof ( ego   ) ) ;
      trimmed_surf  = (ego *  ) EG_alloc ( 6 * sizeof ( ego   ) ) ;
      if ( edges == NULL || surf == NULL || trimmed_surf == NULL || esenses == NULL ) {
	  stat = EGADS_MALLOC;
	  goto cleanup;
      }
      for ( iface = 0 ; iface < nface; iface++ ) {
	  stat = EG_getTopology(faces[iface], &surf[iface], &oclass, &mtype,
				data, &nloop, &loops, &senses);
	  //stat = EG_convertToBSpline(faces[iface], &faces[iface]);
	  if ( stat != EGADS_SUCCESS) goto cleanup;
	  if (nloop != 1) {
	      printf(" Face %d has more than one Loop\n", iface+1);
	      stat = EGADS_TOPOERR;
	      goto cleanup;
	  }
	  stat = EG_getTopology(loops[0], &eref, &oclass, &mtype,
				data, &nedge, &edges[iface], &esenses[iface]);
	  if ( stat != EGADS_SUCCESS ) {
	      printf(" EG_getTopology for face %d  = %d !!\n", iface + 1, stat ) ;
	      goto cleanup;
	  }
	  if (nedge != 4) {
	      printf(" Face %d is not bounded by 4 Edges\n", iface+1);
	      stat = EGADS_TOPOERR;
	      goto cleanup;
	  }
	  for ( iedge = 0 ; iedge < 4; iedge++ ) {
	      stat = EG_getRange(edges[iface][iedge], &trange[iface][2 * iedge], &periodic);
	      if ( stat != EGADS_SUCCESS) {
		  printf(" EG_getRange = %d\n", stat);
		  goto cleanup;
	      }
	      if ( esenses[iface][iedge] == -1 ) {
		  t      [0                     ] = trange [ iface ][2 * iedge + 1];
		  trange [ iface ][2 * iedge + 1] = trange [ iface ][2 * iedge    ];
		  trange [ iface ][2 * iedge    ] = t      [0                     ];
	      }
	  }
      }
      /* FACE 0 = TOP   ; FACE 5 = BOTTOM
         FACE 1 = RIGHT ; FACE 3 = LEFT
         FACE 2 = BACK  ; FACE 4 = FRONT */
      idsF[0] = 0; idsF[1] = 5; // u direction uses faces 0 5 (idsF) vs 2 4 (idsO)
      idsF[2] = 0; idsF[3] = 5; // v direction uses faces 0 5 (idsF) vs 1 3 (idsO)

      idsO[0] = 2; idsO[1] = 4;  // front + back -> u dir
      idsO[2] = 1; idsO[3] = 3;  // right + left -> v dir
      /* We are going to set knot sequences for u, v: find the largest side in each direction */
      for ( i = 0 ; i < 2; ++i ) { // 0 = u, 1 = v
	  /* face 1 */
	  stat = EGADS_SUCCESS;
	  for ( k = 0 ; k < 2; ++k ) {
	      stat += EG_getCommonEdge ( 4, edges [ idsF[ 2 * i     ] ],
					 4, edges [ idsO[ 2 * i + k ] ],   &ce[i][k]);
	      stat += EG_arcLength     (    edges [ idsF[ 2 * i]    ] [     ce[i][k]    ],
					    trange[ idsF[ 2 * i]    ] [ 2 * ce[i][k]    ],
					    trange[ idsF[ 2 * i]    ] [ 2 * ce[i][k] + 1], &arc[i][k] );
	  }
	  /* face 2 : opposite to face 1: back or left */
	  for ( k = 0 ; k < 2; ++k ) {
	      stat += EG_getCommonEdge ( 4, edges [ idsF[ 2 * i + 1 ] ],
					 4, edges [ idsO[ 2 * i + k ] ],    &ce[i][k + 2]);
	      stat += EG_arcLength     (    edges [ idsF[ 2 * i + 1 ] ][     ce[i][k + 2]    ],
					    trange[ idsF[ 2 * i + 1 ] ][ 2 * ce[i][k + 2]    ],
					    trange[ idsF[ 2 * i + 1 ] ][ 2 * ce[i][k + 2] + 1], &arc[i][k + 2] );
	  }
	  if ( stat != EGADS_SUCCESS  ) {
	      printf(" errr, I am assuming that faces %d and %d are common to %d and %d \n", idsF[2 * i], idsF[2 * i + 1], idsO[ 2 * i], idsO[2 * i + 1] );
	      stat = EGADS_INDEXERR; goto cleanup;
	  }
#ifdef DEBUG
	  printf(" -----------    FACE   %d   DIR  %d WILL USE EDGES %d  and %d  \n", idsF[2 * i ]    , i,  ce[i][0], ce[i][1]);
	  printf(" -----------    FACE   %d   DIR  %d WILL USE EDGES %d  and %d  \n", idsF[2 * i + 1 ], i,  ce[i][2], ce[i][3]);
#endif
      }
      uvwID[0] = 0; uvwID[1] = 0; uvwID[2] = 0;
      for ( i = 1 ; i < 4; i++ )  {
	  if ( arc[0][i] > arc[0][uvwID[0] ] ) uvwID[0] = i;
	  if ( arc[1][i] > arc[1][uvwID[1] ] ) uvwID[1] = i;
	  if ( arc[2][i] > arc[2][uvwID[2] ] ) uvwID[2] = i;
      }
      XYZ = ( double **) EG_alloc ( 3 * sizeof ( double *) ); //u, v, w directions
      if ( XYZ == NULL ) {
	  stat = EGADS_MALLOC;
	  goto cleanup;
      }
      // Now sample data along chosen edges and construct knot sequence
      /* sample curves for building a suitable knot sequence*/
      for ( j = 0 ; j < 3; j++) {
	  n      = 10 * nT[j];
	  dt     = 1.0 / (double)( n - 1 ) ;
	  XYZ[j] = ( double * ) EG_alloc ( 3 * n * sizeof ( double ) );
	  if ( XYZ[j] == NULL ) {
	      stat = EGADS_MALLOC;
	      goto cleanup;
	  }
	  if ( uvwID[j] < 2 ) k = 0;
	  else                k = 1;
	  tEval[0] = trange [ idsF[ 2 * j + k ] ] [ 2 * ce[j][ uvwID[j] ]    ];
	  tEval[1] = trange [ idsF[ 2 * j + k ] ] [ 2 * ce[j][ uvwID[j] ] + 1];
	  for ( i = 0 ; i < n; i++ ) {
	      // check the directino of the t0 t1
	      t[0] = tEval[0] * ( 1.0 - (double) i * dt ) + tEval[1] * (double) i * dt ;
	      EG_evaluate( edges[ idsF[ 2 * j + k ] ] [     ce[j][ uvwID[j] ] ], &t[0],data);
	      XYZ[j][ 3 * i ] = data[0]; XYZ[j][ 3 * i + 1 ] = data[1]; XYZ[j][ 3 * i + 2 ] = data[2];
	  }
      }
      knots  = (double **) EG_alloc ( 3 * sizeof ( double * ) ) ;
      if ( knots == NULL ) {
	  stat = EGADS_MALLOC;
	  goto cleanup;
      }
      for ( i = 0 ; i < 3; i++ ) {
	  knots[i] = ( double* ) EG_alloc ( (nT[i] + 1)  * sizeof ( double ) );
	  if ( knots[i] == NULL ) {
	      stat = EGADS_MALLOC;
	      goto cleanup;
	  }
	  stat  = setKnotSeq ( n, XYZ[i], nT[i] - 1 + 2 * ( splines_degree + 1 ), splines_degree, knots[i]);  // num knots (3d argument) "pretend they are clamped and ignore the multiplicity" /
	  if ( stat != EGADS_SUCCESS ) {
	      printf(" setKnotSeq STAT %d !!!\n", stat ) ;
	      goto cleanup;
	  }

#ifdef DEBUG
	  printf("=====================\n");
	  printf(" KNOT SEQUENCE FOR %d DIR \n", i);
	  for ( j = 0 ; j <= nT[i]; j++ ) {
	      printf(" K[%d] = %f\t", j, knots[i][j]);
	      if ( j % 5 == 0  ) printf("\n");
	  }
	  printf("\n");
#endif
      }
#ifdef DEBUG
      printf("=====================\n");
#endif
      // Approximate Bspline Surface u x v and then produce ruled 
      /* TFI FILLING NOT DONE */
      /* u x v = top   and bottom 0,1 */
      /* u x w = front and back   2,3 */
      /* v x w = left  and right  4,5 */
      idsF[0] = 0; idsF[1] = 5; idsF[2] = 2; idsF[3] = 4;
      idsO[0] = 1; idsO[1] = 1; idsO[2] = 0; idsO[3] = 0;
      /* vector containig the grid info along the four edges */
      tuvs    = (double **) EG_alloc (4 * sizeof ( double *));
      if ( tuvs == NULL) {
	  stat = EGADS_MALLOC;
	  goto cleanup;
      }
      for ( i = 0 ; i < 4; i++ ) tuvs[i] = NULL;
      /* Using Face 2 = front as leading */ 
      for ( k = 0; k < 6; k = k+5 ) {
	  stat = EG_getCommonEdge (  4, edges [ k ], 4, edges [ 2 ], &e1 );
	  if ( stat != EGADS_SUCCESS ) {
	      printf(" ERROR EG_getCommonEdge DIR %d FACES %d  %d have nobody in common \n", j, idsF[k], 2);
	      goto cleanup;
	  }
	  e2 = ( e1 + 1 ) % 4;
	  e3 = ( e1 + 2 ) % 4;
	  e4 = ( e1 + 3 ) % 4;
	  // set min max values
	  tEval[0] = trange [ k ] [ 2 * e1    ];
	  tEval[1] = trange [ k ] [ 2 * e1 + 1];
	  tEval[2] = trange [ k ] [ 2 * e2    ];
	  tEval[3] = trange [ k ] [ 2 * e2 + 1];
	  //
	  tEval[4] = trange [ k ] [ 2 * e3    ];
	  tEval[5] = trange [ k ] [ 2 * e3 + 1];
	  tEval[6] = trange [ k ] [ 2 * e4    ];
	  tEval[7] = trange [ k ] [ 2 * e4 + 1];
	  //
	  for ( i = 0 ; i < 4; i++ ) {
	      if ( i %2 == 0 ) tuvs[i] = EG_reall ( tuvs[i], (nT[0] + 1) * 2 * sizeof (double ));
	      else             tuvs[i] = EG_reall ( tuvs[i], (nT[1] + 1) * 2 * sizeof (double ));
	      if ( tuvs[i] == NULL ) {
		  stat = EGADS_MALLOC;
		  goto cleanup;
	      }
	  }
	  // get corresponding UVs from t in edge
#ifdef DEBUG
	  printf(" RANGE FOR EDGE %d = [ %f   %f ] \n", e1, tEval[0], tEval[1]);
	  printf(" RANGE FOR EDGE %d = [ %f   %f ] \n", e2, tEval[2], tEval[3]);
	  printf(" RANGE FOR EDGE %d = [ %f   %f ] \n", e3, tEval[4], tEval[5]);
	  printf(" RANGE FOR EDGE %d = [ %f   %f ] \n", e4, tEval[6], tEval[7]);
#endif
	  for ( i = 0 ; i <= nT[0] ; i++ ) { // e0
	      t[0] = tEval[0] * ( 1.0 - knots[ 0 ][i] ) + tEval[1] * knots[ 0 ][i];  // ->
	      stat = EG_getEdgeUV ( faces  [ k ], edges  [ k ][ e1 ], esenses[ k ][ e1 ], 
				    t[0], &tuvs[0][2 * i  ]  );
	      t[0] = tEval[5] * ( 1.0 - knots[ 0 ][i] ) + tEval[4] * knots[ 0 ][i]; // <-
	      stat = EG_getEdgeUV ( faces  [ k ], edges  [ k ][ e3 ], esenses[ k ][ e3 ],
				    t[0], &tuvs[2][2 * ( nT[0] - i ) ]  );
	  }
	  for ( i = 0 ; i <= nT[1] ; i++ ) { //v
	      t[0] = tEval[2] * ( 1.0 - knots[ 1 ][i] ) + tEval[3] * knots[ 1 ][i]; // ->
	      stat = EG_getEdgeUV ( faces  [ k ], edges  [ k ][ e2 ], esenses[ k ][ e2 ], 
				    t[0], &tuvs[1][2 * i ]  );
	      t[0] = tEval[7] * ( 1.0 - knots[ 1 ][i] ) + tEval[6] * knots[ 1 ][i]; // <-
	      stat = EG_getEdgeUV ( faces  [ k ], edges  [ k ][ e4 ],
				    esenses[ k ][ e4 ], t[0], &tuvs[3][2 * ( nT[1] - i) ]  );
	  }
	  if ( k == 5 ) { // reverse all for ruling 
	      for ( i = 0 ; i < 4; i++ ) {
		  if ( i % 2 == 0 ) n = nT[0] + 1;
		  else              n = nT[1] + 1;
		  for ( l = 0 ; l < 2*n; l++ ) paux[i][l] = tuvs[i][l];
	      }
	      for ( i = 0 ; i < 4; i++ ) {
		  if      ( i == 1 ) j = 3;
		  else if ( i == 3 ) j = 1; 
		  else               j = i;
		  if ( i % 2 == 0 )  n = nT[0] + 1;
		  else               n = nT[1] + 1;
		  for ( l = 0 ; l < n; l++ ) {
		      tuvs[j][2 * l    ] = paux[i][2 * (n - 1 - l)     ] ;
		      tuvs[j][2 * l + 1] = paux[i][2 * (n - 1 - l)  + 1] ;
		  }
	      }
	  }
	  EG_free ( uvTFI );
	  EG_free ( surf_grid );
	  dims[0]   = nT[0] + 1;
	  dims[1]   = nT[1] + 1;
	  dims[2]   = nT[0] + 1;
	  dims[3]   = nT[1] + 1;
	  uvTFI     = ( double * ) EG_alloc ( 2 * dims[0] * dims[1] * sizeof ( double ) ) ;
	  surf_grid = ( double * ) EG_alloc ( 3 * dims[0] * dims[1] * sizeof ( double ) ) ;
	  if ( uvTFI == NULL || surf_grid == NULL) {
	      stat = EGADS_MALLOC;
	      goto cleanup;
	  }
	  /*TFI Creates a uv grid using the four edges and the discretization along each edge */
	  EG_setGridTFI ( dims, tuvs, uvTFI );
	  /* now get physical coordinates */
	  for ( l = 0 ; l <= nT[0] ; l++ ) {
	  for ( i = 0 ; i <= nT[1] ; i++ ) {
	       stat = EG_evaluate ( faces[ k ], &uvTFI[2 * ( ( nT[0] + 1 ) * i + l )  ], data );
		  if ( stat != EGADS_SUCCESS ) {
		      printf(" PROBLEM EVALUATING %f  %f  \n", uvTFI[2 * ( ( nT[0] + 1 ) * i + l )  ], uvTFI[2 * ( ( nT[0] + 1 ) * i + l )  + 1] );
		      goto cleanup;
		  }
		  /*surf_grid[ 3 * ( ( nT[0] + 1 ) * i + l )     ] = data[0];
		  		  surf_grid[ 3 * ( ( nT[0] + 1 ) * i + l )  + 1] = data[1];
		  		  surf_grid[ 3 * ( ( nT[0] + 1 ) * i + l )  + 2] = data[2];*/
		  surf_grid[ 3 * ( ( nT[1] + 1 ) * l + i) ] = data[0];
		  surf_grid[ 3 * ( ( nT[1] + 1 ) * l + i) + 1] = data[1];
		  surf_grid[ 3 * ( ( nT[1] + 1 ) * l + i) + 2] = data[2];
	      }
	  }
#ifdef DEBUG
	  snprintf(fname,100, "TFIdata_%d.txt", k);
	  fil = fopen (fname, "w");
	  if ( fil == NULL ) {
	      printf(" I COULDN'T OPEN FILE = %p name %s\n", fil,  fname);
	      stat = EGADS_MALLOC;
	      goto cleanup;
	  }

	  for ( i = 0 ; i <= nT[1] ; i++ ) {
	      for ( l = 0 ; l <= nT[0] ; l++ ) {
		  fprintf(fil, "%lf %lf %lf %lf \t\t %lf\t%lf\t%lf\n",knots[0][l], knots[1][i],
			  uvTFI[2 * ( ( nT[0] + 1 ) * i + l )  ],uvTFI[2 * ( ( nT[0] + 1 ) * i + l ) + 1 ],
			  surf_grid[ 3 * ( ( nT[0] + 1 ) * i + l )    ],
			  surf_grid[ 3 * ( ( nT[0] + 1 ) * i + l ) + 1],
			  surf_grid[ 3 * ( ( nT[0] + 1 ) * i + l ) + 2] );
	      }
	      fprintf(fil,"\n\n");
	  }
	  fclose (fil);
#endif
	  stat =  EG_spline2dAppx ( context, 1 , knots[1], knots[0], 
				    NULL, NULL, NULL, NULL, NULL, NULL, NULL,
				    nT[1] + 1, nT[0] + 1,  surf_grid, 1.E-05,
				    &trimmed_surf[ k ] );
	  if ( stat != EGADS_SUCCESS ) {
	      printf(" EG_spline2dAppx %d !!!! \n", stat );
	      goto cleanup;
	  }
      }
      /* Save a trimmed egas file */
      /* SAVE SURFACES FOR IRIT *.itd */
      /* ********************** */
      IRITsurf = EG_alloc ( 2 * sizeof(CagdSrfStruct));
      if (IRITsurf == NULL ) {
	  stat =  EGADS_MALLOC;
	  goto cleanup;
      }
      stat  = makeIRITsrf(trimmed_surf[0], &IRITsurf[0]);
      stat  = makeIRITsrf(trimmed_surf[5], &IRITsurf[1]);
      sprintf ( fname, "surf%d"    ,1 );
      sprintf ( fnItd, "surf%d.itd",1 );
      BspSrfWriteToFile(IRITsurf[0], fnItd, 0, fname, &ErrStr);
      sprintf ( fname, "surf%d"    ,2);
      sprintf ( fnItd, "surf%d.itd",2);
      BspSrfWriteToFile(IRITsurf[1], fnItd, 0, fname, &ErrStr);
      if (ErrStr != NULL) {
	  EG_free (IRITsurf);
	  printf(" PROBLEM IN  surf1 ErrStr: %s\n", ErrStr);
	  stat = EGADS_GEOMERR;
	  goto cleanup;
      }
      TVMap = TrivRuledTV(IRITsurf[0], IRITsurf[1], 2, nT[2]);
      sprintf ( fname, "wingduct_T%dx%dx%d",nT[0], nT[1], nT[2] ) ;
      sprintf ( fnItd, "wingduct_T%dx%dx%d.itd",nT[0], nT[1], nT[2] ) ;

      stat  = TrivBspTVWriteToFile(TVMap, fnItd, 0, fname,  &ErrStr);
      printf("=========================================================\n");
      printf("=========================================================\n");
      printf( " Designed a duct for %d x %d x %d tiles \n", nT[0], nT[1], nT[2] ) ;
      printf("=========================================================\n");
      
      EG_free (IRITsurf);
  }
  cleanup:
  EG_free  (surf_grid);
  EG_free  (uvTFI    );
  if ( knots ) {
      for ( i = 0 ; i < 3; i++)
	EG_free  (knots[i] );
  }
  if ( XYZ) {
      for ( i = 0 ; i < 3; i++)
	EG_free  (XYZ[i]   );
  }
  EG_free (knots );
  EG_free (XYZ);
  if ( tuvs ) {
      for ( i = 0 ; i < 4; i++ ) EG_free (tuvs[i]);
  }
  EG_free (tuvs);
  EG_free (trimmed_surf);
  EG_free (edges);
  EG_free (esenses);
  EG_free (faces);
  EG_free (surf);
  EG_close (context);
  return 0;

}



/*
 ************************************************************************
 *                                                                      *
 *   makeIRITsrf -- convert EGADS BSPLINE into IRIT Surface              *
 *                                                                      *
 ************************************************************************
 */

static int
makeIRITsrf(ego           eobj,         /* (in)  Egads surface of Face */
	    CagdSrfStruct **Srf)        /* (in)  IRIT  surface (or NULL) */
{
  int     status = EGADS_SUCCESS;

  int     i, j, k, indx, oclass, mtype, *ivec=NULL, nloop, nedge, nnode, nchild, *senses, *senses2;
  double  *rvec=NULL, data[4], xyzsw[4], xyzse[4], xyznw[4], xyzne[4];
  ego     topRef, eprev, enext, eref, esurf, *eloops, *eedges, *enodes, *echilds;



  /* --------------------------------------------------------------- */

  /* default return */
  *Srf = NULL;

  /* get type of eobj */
  status = EG_getInfo(eobj, &oclass, &mtype, &topRef, &eprev, &enext);
  if ( status != EGADS_SUCCESS ) {
      printf( " EG_getInfo stat = %d\n", status);
      return status;
  }

  /* eobj is a Surface */
  if (oclass == SURFACE) {
      esurf = eobj;

      /* make sure esurf is a BSpline (and not a NURB and not periodic) */
      status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
      if ( status != EGADS_SUCCESS ) {
	  printf( " EG_getGeometry stat = %d\n", status);
	  return status;
      }
      if (mtype != BSPLINE) {
	  printf(" makeIRITsrf: esurf NOT a BSpline = %d\n", mtype);
	  status = EGADS_GEOMERR;
	  goto cleanup;
      } else if (ivec[0] != 0) {
	  printf(" makeIRITsrf: BSpline flags = %d\n", ivec[0]);
	  status = EGADS_GEOMERR;
	  goto cleanup;
      }
      /* start a new IRIT BSpline surface */
      *Srf = BspSrfNew(ivec[2], ivec[5], ivec[1]+1, ivec[4]+1, CAGD_PT_E3_TYPE);

      /* Populate the control points in Srf */
      k = ivec[3] + ivec[6];
      for (j = 0; j < ivec[5]; j++) {
	  for (i = 0; i < ivec[2]; i++) {
	      indx = CAGD_MESH_UV(*Srf, i, j);
	      (*Srf)->Points[1][indx] = rvec[k++];   /* X coefs */
	      (*Srf)->Points[2][indx] = rvec[k++];   /* Y coefs */
	      (*Srf)->Points[3][indx] = rvec[k++];   /* Z coefs */
	  }
      }

      /* Copy the knot vectors into Srf */
      for (i = 0; i < ivec[3]; i++) {
	  (*Srf)->UKnotVector[i] = rvec[i];
      }
      for (i = 0; i < ivec[6]; i++) {
	  (*Srf)->VKnotVector[i] = rvec[i+ivec[3]];
      }

      /* eobj is a Face */
  } else if (oclass == FACE) {
      status = EG_getTopology(eobj, &eref, &oclass, &mtype, data,
			      &nloop, &eloops, &senses);
      if ( status != EGADS_SUCCESS ) {
	  printf( " EG_getTopology stat = %d\n", status);
	  return status;
      }



      status = EG_getTopology(eloops[0], &eref, &oclass, &mtype, data,
			      &nedge, &eedges, &senses);
      if ( status != EGADS_SUCCESS ) {
	  printf( " EG_getTopology stat = %d\n", status);
	  return status;
      }

      /* get the 4 corners as the beginning/end of each Edge */
      if (nedge == 4) {
	  status = EG_getTopology(eedges[0], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[0] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyzsw,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyzsw,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  status = EG_getTopology(eedges[1], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[1] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyzse,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyzse,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  status = EG_getTopology(eedges[2], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[2] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyzne,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyzne,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  status = EG_getTopology(eedges[3], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[3] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyznw,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyznw,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  /* make an untrimmed Bspline Surface */
	  *Srf = BspSrfNew(2, 2, 2, 2, CAGD_PT_E3_TYPE);

	  /* knot vectors */
	  (*Srf)->UKnotVector[0] = 0;
	  (*Srf)->UKnotVector[1] = 0;
	  (*Srf)->UKnotVector[2] = 1;
	  (*Srf)->UKnotVector[3] = 1;

	  (*Srf)->VKnotVector[0] = 0;
	  (*Srf)->VKnotVector[1] = 0;
	  (*Srf)->VKnotVector[2] = 1;
	  (*Srf)->VKnotVector[3] = 1;

	  /* control points */
	  indx = CAGD_MESH_UV(*Srf, 0, 0);
	  (*Srf)->Points[1][indx] = xyzsw[0];
	  (*Srf)->Points[2][indx] = xyzsw[1];
	  (*Srf)->Points[3][indx] = xyzsw[2];

	  indx = CAGD_MESH_UV(*Srf, 1, 0);
	  (*Srf)->Points[1][indx] = xyzse[0];
	  (*Srf)->Points[2][indx] = xyzse[1];
	  (*Srf)->Points[3][indx] = xyzse[2];

	  indx = CAGD_MESH_UV(*Srf, 0, 1);
	  (*Srf)->Points[1][indx] = xyznw[0];
	  (*Srf)->Points[2][indx] = xyznw[1];
	  (*Srf)->Points[3][indx] = xyznw[2];

	  indx = CAGD_MESH_UV(*Srf, 1, 1);
	  (*Srf)->Points[1][indx] = xyzne[0];
	  (*Srf)->Points[2][indx] = xyzne[1];
	  (*Srf)->Points[3][indx] = xyzne[2];

	  /* if Face does not have 4 Nodes but is a planar
           Face (hopefully with only a notch), so that we can
           convert Face to BSpline and then use it directly */
      } else {
	  status = EG_convertToBSpline(eobj, &esurf);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_convertToSPline stat = %d\n", status);
	      return status;
	  }

	  status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getGeometry stat = %d\n", status);
	      return status;
	  }

	  /* start a new IRIT BSpline surface */
	  *Srf = BspSrfNew(ivec[2], ivec[5], ivec[1]+1, ivec[4]+1, CAGD_PT_E3_TYPE);

	  /* Populate the control points in Srf */
	  k = ivec[3] + ivec[6];
	  for (j = 0; j < ivec[5]; j++) {
	      for (i = 0; i < ivec[2]; i++) {
		  indx = CAGD_MESH_UV(*Srf, i, j);
		  (*Srf)->Points[1][indx] = rvec[k++];   /* X coefs */
		  (*Srf)->Points[2][indx] = rvec[k++];   /* Y coefs */
		  (*Srf)->Points[3][indx] = rvec[k++];   /* Z coefs */
	      }
	  }

	  /* Copy the knot vectors into Srf */
	  for (i = 0; i < ivec[3]; i++) {
	      (*Srf)->UKnotVector[i] = rvec[i];
	  }
	  for (i = 0; i < ivec[6]; i++) {
	      (*Srf)->VKnotVector[i] = rvec[i+ivec[3]];
	  }
      }

      /* eobj is an unknown type */
  } else {
      printf(" makeIRITsrf: eobj is neither Surface nor Face (oclass=%d)\n", oclass);
      status = EGADS_GEOMERR;
      goto cleanup;
  }

  cleanup:
  if (ivec != NULL) EG_free(ivec);
  if (rvec != NULL) EG_free(rvec);

  return status;
}







