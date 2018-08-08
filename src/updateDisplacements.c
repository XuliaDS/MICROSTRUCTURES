
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"


char tiledName[100], ductName[100];
#define DEBUG
#define EPS04  1.e-4

typedef struct TiledDuct{
  TrivTVStruct **trivars, *duct;
  int nTV;
  int trivsPerTile[3];
  int tilesPerKnot[3];
} TiledDuct;

// ASSUMING DUCT TRIVAR LIVES IN [0,1]^3
CagdBBoxStruct getTileBBox ( TrivTVStruct *TVAR, int id, int tilesPerKnot[] ) {
  int ijk[3], m, n, d, res1, knotsPerDir[3], totalTiles[3];
  CagdRType knotCoords[2], tileCoords[2];
  CagdBBoxStruct BBox;
  id--;
  knotsPerDir[0] = TVAR -> ULength - TVAR -> UOrder + 1 ;
  knotsPerDir[1] = TVAR -> VLength - TVAR -> VOrder + 1 ;
  knotsPerDir[2] = TVAR -> WLength - TVAR -> WOrder + 1 ;
  totalTiles [0] = tilesPerKnot[0] * knotsPerDir[0];
  totalTiles [1] = tilesPerKnot[1] * knotsPerDir[1];
  totalTiles [2] = tilesPerKnot[2] * knotsPerDir[2];
  ijk        [2] = id / ( totalTiles [0] * totalTiles [1]);
  res1           = id % ( totalTiles [0] * totalTiles [1]);
  ijk        [1] = res1 / totalTiles [0] ;
  ijk        [0] = res1 % totalTiles [0] ;
  printf(" TILE COORDS %d %d %d\n",ijk[0], ijk[1], ijk[2] );
  // Get tile and knot coordinates:
  for ( d = 0 ; d < 3; d++ ) {
      m = ijk[d] / tilesPerKnot[d];
      n = ijk[d] % tilesPerKnot[d];
      tileCoords[0] = (double)  m      / (double) knotsPerDir[d]; // min
      tileCoords[1] = (double) (m + 1) / (double) knotsPerDir[d]; // max
      knotCoords[0] = (double)  n      / (double) tilesPerKnot[d]; // min
      knotCoords[1] = (double) (n + 1) / (double) tilesPerKnot[d]; // max
      BBox.Min  [d] = (1.0 - knotCoords[0] ) * tileCoords[0] + knotCoords[0] * tileCoords[1];
      BBox.Max  [d] = (1.0 - knotCoords[1] ) * tileCoords[0] + knotCoords[1] * tileCoords[1];
#ifdef DEBUGG
      printf(" DIR  %d : TILE %d of knot %d \n", d, n, m ) ;
      printf(" TILE COORDS %f  %f   KNOT COORDS %f  %f \n",
	     tileCoords[0], tileCoords[1], knotCoords[0], knotCoords[1] );
      printf(" DUCT COORDINATES %f  %f  \n", BBox.Min[d], BBox.Max[d] );
#endif
  }
  return BBox;
}


void getDuctUVW ( double uvw[], int ID, TiledDuct *TV, double *uvwDuct ) ;

void updatePressure ( char *ductFile, char *tileFile, char *dataFile,  int tilesPerKnot[]) {
  const char *name , *name2;
  char *retstr, *lastdot;
  char  s[100],  outName[100];
  int i, j, nP[2], **trivIDs = NULL;
  FILE *fil, *fout;
  CagdRType **pts = NULL, uvw[3];
  TiledDuct  *TD  = NULL;
  TD   = ( TiledDuct * ) malloc ( sizeof ( TiledDuct ) );
  if ( TD == NULL ) {
      fprintf(stderr," MALLOC!! unable to allocate memo in updatePressure\n");
      return ;
  }
  name2                 = tileFile;
  printf(" LOAD TIlES %s\n",name2);
  IPObjectStruct *TObj = IPGetDataFiles(&name2, 1, TRUE, TRUE);
  if ( TObj == NULL ) {
      fprintf(stderr," NULL OBJECT. Failed to load tiled duct \n");
      return ;
  }
  if (IP_IS_OLST_OBJ(TObj)) {
      IPObjectStruct *PTmp;
      for (i = 0; (PTmp = IPListObjectGet(TObj, i++)) != NULL; )
	TD -> trivars = ( TrivTVStruct ** ) malloc ( i * sizeof ( TrivTVStruct * )  );
      if ( TD -> trivars == NULL ) {
	  fprintf(stderr," MALLOC!! unable to allocate memo in updatePressure\n");
	  return ;
      }
      for (i = 0; (PTmp = IPListObjectGet(TObj, i++)) != NULL; ) {
	  if (IP_IS_TRIVAR_OBJ(PTmp))  TD -> trivars[i] = PTmp -> U.Trivars;
      }
  }
  else {
      fprintf(stderr, "ERROR. I AM ASSUMING A OBJECT LIST !\n");
      return ;
  }
  name                 = ductFile;
  printf(" LOAD DUCT %s\n", ductFile);
  IPObjectStruct *duct = IPGetDataFiles(&name, 1, TRUE, TRUE);
  if ( duct == NULL ) {
      printf(" NULL OBJECT. Failed to duct \n");
      return;
  }
  TD -> duct            = duct -> U.Trivars;
  TD -> tilesPerKnot[0] = tilesPerKnot[0];
  TD -> tilesPerKnot[1] = tilesPerKnot[1];
  TD -> tilesPerKnot[2] = tilesPerKnot[2];
  // create pressure file name
  if ((retstr = malloc (strlen (dataFile) + 1)) == NULL)   return;
  strcpy (retstr, dataFile);
  lastdot = strrchr (retstr, '.');
  if (lastdot != NULL)
    *lastdot = '\0';
  snprintf(outName, 100,"%s_pressures.dat", retstr);
  printf(" Pressure file name  %s \n", outName);
  //
  trivIDs = ( int ** )      malloc ( 2 * sizeof ( int       *) );
  pts     = ( CagdRType **) malloc ( 2 * sizeof ( CagdRType *) );
  if ( trivIDs == NULL || pts == NULL ) {
      fprintf(stderr," MALLOC !! failed to allocate MEMO updatePressure\n");
      return ;
  }
  // SURFACE 1:
  fil = fopen (dataFile,"r");
  if ( fil == NULL ) {
      fprintf(stderr," FATAL ERROR: Unable to open file %s\n", dataFile);
      return ;
  }
  for ( j = 0 ; j < 2; j++ ) {
      fscanf(fil, "%d", &nP[j] );
      pts    [j] = (CagdRType * ) malloc ( 3 * nP[j] * sizeof ( CagdRType ) ) ;
      trivIDs[j] = (int       * ) malloc (     nP[j] * sizeof ( int       ) ) ;
      if ( pts[j] == NULL || trivIDs[j] == NULL ) {
	  fprintf(stderr, " MALLOC PROBLEM IN updatePressure function!! ABORT\n");
	  fclose ( fil );
	  IPFreeObject ( TObj    );
	  free         ( trivIDs );
	  free         ( pts     );
	  return ;
      }
      for ( i = 0 ; i < nP[j]; i++ ) {
	  fscanf ( fil, "%lf", &pts    [j][3 * i     ]);
	  fscanf ( fil, "%lf", &pts    [j][3 * i + 1 ]);
	  fscanf ( fil, "%lf", &pts    [j][3 * i + 2 ]);
	  fscanf ( fil, "%d",  &trivIDs[j][    i     ]);
      }
  }
  fout = fopen( outName, "w" );
  for ( j = 0 ; j < 2; j++ ) {
      fprintf(fout,"%d\n", nP[j] );
      printf( " TOTAL POINTS  %d  \n", nP[j] );
      for ( i = 0 ; i < nP[j] ; i++ ) {
	  printf(" TRIV %d UVW %f  %f  %f  \n", trivIDs[j][i], pts[j][ 3 * i],
		 pts[j][ 3 * i + 1], pts[j][ 3 * i + 2]);
	  getDuctUVW ( &pts[j][3 * i], trivIDs[j][i], TD, uvw );
	  //EG_evaluatePressure ( );
	  double press = 0.001;
	  fprintf(fout,"%lf %lf %lf %lf %d\n", pts[j][3*i], pts[j][3*i + 1], pts[j][3*i + 2], press, trivIDs[j][i] );
      }
  }
  fclose ( fil  );
  fclose ( fout );
  return;
}



void getDuctUVW ( double uvw[], int ID, TiledDuct *TD, double *uvwDuct ) {
  int i, j, tileID;
  int tileIDbounds[3], ductIDbounds[3];
  CagdRType *P3 = NULL, *PP3 = NULL, uvw0[3];
  CagdBBoxStruct BBox;
  CagdSrfStruct  **BndSrf = NULL;
  void *handler;
  BndSrf = ( CagdSrfStruct ** ) malloc (6 * sizeof ( CagdSrfStruct * ));
  if ( BndSrf == NULL ) {
      fprintf(stderr," MALLOC AT BndSrf !!!\n");
      return ;
  }
  for ( i = 0 ; i < 6; i++ )  {
      BndSrf[i] = (CagdSrfStruct * ) malloc ( sizeof ( CagdSrfStruct  ));
      BndSrf[i] = IRIT_ZAP_MEM ( BndSrf[i], sizeof ( CagdSrfStruct  ));
  }
  P3  = IritMalloc(CAGD_MAX_PT_SIZE * sizeof(CagdRType));
  PP3 = IritMalloc(CAGD_MAX_PT_SIZE * sizeof(CagdRType));
  // Get six faces
  printf(" GET DUCT 6 SURFACES \n");
  TrivBndrySrfsFromTVToData( TD -> duct, FALSE, BndSrf );
  // get tile ID
  tileID = AttrGetIntAttrib( TD -> trivars[ID-1] ->Attr, "MSGeomID");
  printf(" TRIVAR %d BELONGS TO TILE %d \n", ID, tileID);
  if ( tileID < 0 ) {
      fprintf(stderr, "FATAL ERROR:  I can't find a tile that holds trivariate %d \n", ID);
      exit(1);
  }
  BBox = getTileBBox (TD -> duct, tileID, TD -> tilesPerKnot);
  printf("\n\n TILE %d -> UVW BBOX %f %f %f MAX %f %f %f \n", tileID, BBox.Min[0], BBox.Min[1], BBox.Min[2],
	 BBox.Max[0], BBox.Max[1], BBox.Max[2] );
  // check that indeed we have the right tile:
  tileIDbounds[0] = AttrGetIntAttrib( TD -> trivars[ID -1 ] -> Attr, "MSLclTVBndryU" );
  tileIDbounds[1] = AttrGetIntAttrib( TD -> trivars[ID -1 ] -> Attr, "MSLclTVBndryV" );
  tileIDbounds[2] = AttrGetIntAttrib( TD -> trivars[ID -1 ] -> Attr, "MSLclTVBndryW" );
  ductIDbounds[0] = AttrGetIntAttrib( TD -> trivars[ID -1 ] -> Attr, "MSDfrmTVBndryU");
  ductIDbounds[1] = AttrGetIntAttrib( TD -> trivars[ID -1 ] -> Attr, "MSDfrmTVBndryV");
  ductIDbounds[2] = AttrGetIntAttrib( TD -> trivars[ID -1 ] -> Attr, "MSDfrmTVBndryW");
  for ( i = 0 ; i < 3; i++) {
      if ( tileIDbounds[i] < 0 ) continue;
      printf("\n\n TILE ID BOUND %d ( SURFACE %d ) \n", tileIDbounds[i], ductIDbounds[i]);
      printf(" LOOKING FOR %f  %f  %f  from TRIV (%d) in DUCT \n", uvw[0], uvw[1], uvw[2], ID );
      // get uv0s using tiles centre wrt trivariate
      for ( j = 0 ; j < 3; ++j)
	uvw0[j] = 0.5 * ( BBox.Min[j] + BBox.Max[j] );
      if ( tileIDbounds[i] % 2 == 0 ) uvw0[i] = BBox.Min[i];
      else                            uvw0[i] = BBox.Max[i];
      TRIV_TV_EVAL_E3(TD -> trivars[ID -1]  , uvw [0], uvw [1], uvw [2], P3);
      TRIV_TV_EVAL_E3(TD -> duct, uvw0[0], uvw0[1], uvw0[2], PP3); // unused ATM
      handler = MvarDistSrfPointPrep( BndSrf[ ductIDbounds[i] ] );
      uvwDuct = MvarDistSrfPoint    ( BndSrf[ ductIDbounds[i] ], handler, P3, TRUE, 0.001, 1.e-10, NULL );
      MvarDistSrfPointFree(handler);
#ifdef DEBUG
      printf(" POINT COORDINATES IN TRIVARIATE: %f  %f  %f \n ", P3[0], P3[1], P3[2]);
      printf(" Initial guess using tile: %f  %f  %f  -> %f   %f   %f (UNUSED)\n", uvw0[0], uvw0[1], uvw0[2], PP3[0], PP3[1], PP3[2]);
      // inv evaluate: get global uvw
      printf(" USING VAR %d  MIN MAX %d -> SURFACE %d \n", i, tileIDbounds[i], ductIDbounds[i] ) ;
      TRIV_TV_EVAL_E3(TD ->duct, uvwDuct[0], uvwDuct[1], uvwDuct[2], PP3);
      printf(" I have found DUCT Coords  %f  %f   (%f ~ 0 ) -> %f  %f  %f\n\n\n", uvwDuct[0], uvwDuct[1], uvwDuct[2],
	     PP3[0], PP3[1], PP3[2]);
#endif
      for ( i = 0 ; i < 6 ; i++ )
	free( BndSrf[i]);
      free ( BndSrf);
      IritFree ( PP3);
      IritFree (P3);
      return ;
  }
  for ( i = 0 ; i < 6 ; i++ )
    free( BndSrf[i]);
  free ( BndSrf);
  IritFree ( PP3);
  IritFree (P3);
}

int main(int argc, char **argv)
{
  if ( argc < 3 ) {
      printf(" I NEED A TILE FILE + DISPLACEMENT FILE \n");
      return -1;
  }
  // LOAD DISPLACEMENT TILES


  DispTV = TrivTVAdd(GeomTV, PropTV);
  DispSrf = CagdSrfAdd(GeomSrf, propSrf);
CagdSrfStruct *CagdMergeSrfSrf(const CagdSrfStruct *CSrf1,

                        const CagdSrfStruct *CSrf2,

                        CagdSrfDirType Dir,

                        TRUE, FALSE);


  return 0;
}








