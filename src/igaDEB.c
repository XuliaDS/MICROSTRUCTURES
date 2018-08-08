
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"
#include <time.h>
char tiledName[100], ductName[100];
#define DEBUG
#define EPS04  1.e-4
#define NUMERIC_TOL 1e-6
double MAXERR;

FILE *filSurf;
char buf[100];

int POINTCOUNT, TOTPTS;
typedef struct TiledDuct{
  TrivTVStruct **trivars, *duct;
  IPObjectStruct **tiles;
  int nTV;
  int trivsPerTile[3];
  int tilesPerKnot[3];
  int *IDTrivTile;
} TiledDuct;


int getTileIDFromTriv ( IPObjectStruct *MS, int trivID, int trivsPerTile ) {
  int i, j, offset;
  IPObjectStruct *tile;
  TrivTVStruct *TV;
  if ( trivsPerTile != 0 ) {
      offset = trivID % ( trivsPerTile + 1) ; // ids have +1 bias
      tile   = IPListObjectGet (MS, offset ) ;
      if ( AttrGetObjectIntAttrib(tile,"MSTileID") != offset + 1) {
	  fprintf(stderr, " I don't understand. I am grabbing the wrong tile....\n");
	  return -1;
      }
      // check that indeed we have the right tile:
      for ( TV= tile -> U.Trivars; TV !=NULL; TV = TV -> Pnext ) {
	  i = AttrGetIntAttrib( TV -> Attr, "MSGeomID");
	  if ( i == trivID ) return offset ;
      }
      fprintf(stderr, " I shouldn't have made it this far. messed up with the tile id\n");
      return -1;
  }
  else {
      for ( offset = 0; offset < IPListObjectLength(MS); offset++) {
	  tile  = IPListObjectGet (MS, offset ) ;
	  // check that indeed we have the right tile:
	  for ( TV = tile -> U.Trivars; TV !=NULL; TV = TV -> Pnext ) {
	      i = AttrGetIntAttrib( TV -> Attr, "MSGeomID");
	      if ( i == trivID ) {
		  printf(" YEP  %d = %d in tile %d \n", i, trivID, offset + 1);
		  return offset;
	      }
	  }
      }
      fprintf(stderr, " I shouldn't have made it this far. messed up with the tile id\n");
      return -1;
  }
}



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


void getDuctUVW ( double uvw[], int ID, int faceID, TiledDuct *TV, double *uvwDuct ) ;

void updatePressure ( char *ductFile, char *tileFile, char *dataFile,  int tilesPerKnot[]) {
  const char *name , *name2;
  char *retstr, *lastdot;
  char  s[100],  outName[100];
  int i, j, sum, sum2,nP[2], **trivIDs = NULL;
  FILE *fil, *fout;
  CagdRType **pts = NULL, uvw[3];
  TiledDuct  *TD  = NULL;
  TrivTVStruct *TVs;
  CagdSrfStruct  **BndSrf = NULL;
  BndSrf = ( CagdSrfStruct ** ) malloc (6 * sizeof ( CagdSrfStruct * ));
  if ( BndSrf == NULL ) {
      fprintf(stderr," MALLOC AT BndSrf !!!\n");
      return ;
  }
  for ( i = 0 ; i < 6; i++ )  {
      BndSrf[i] = (CagdSrfStruct * ) malloc ( sizeof ( CagdSrfStruct  ));
      BndSrf[i] = IRIT_ZAP_MEM ( BndSrf[i], sizeof ( CagdSrfStruct  ));
  }
  TD   = ( TiledDuct * ) malloc ( sizeof ( TiledDuct ) );
  if ( TD == NULL ) {
      fprintf(stderr," MALLOC!! unable to allocate memo in updatePressure\n");
      return ;
  }
  name2                 = tileFile;
  printf(" LOAD TIlES %s\n",name2);
  IPSetFlattenObjects(FALSE);
  IPObjectStruct *TObj = IPGetDataFiles(&name2, 1, TRUE, TRUE);
  if ( TObj == NULL ) {
      fprintf(stderr," NULL OBJECT. Failed to load tiled duct \n");
      return ;
  }
  if (IP_IS_OLST_OBJ(TObj)) {
      TD -> tiles = ( IPObjectStruct ** ) malloc ( IPListObjectLength(TObj) * sizeof ( IPObjectStruct * )  );
      if ( TD -> tiles == NULL ) {
	  fprintf(stderr," MALLOC!! unable to allocate memo in updatePressure\n");
	  return ;
      }
      for ( sum = i = 0; i <  IPListObjectLength(TObj); i++ ) {
	  TD -> tiles[i] = IPListObjectGet (TObj, i ) ;
	  for ( TVs = TD -> tiles[i]->U.Trivars; TVs !=NULL; TVs = TVs -> Pnext ) {
	      sum++;
	  }
      }
      printf(" TOTAL TRIVARS %d \n", sum);
      TD -> IDTrivTile = ( int * ) malloc ( sum * sizeof ( int ) );
      if (  TD -> IDTrivTile == NULL ) {
	  fprintf(stderr, "MALLOC!!! IDTrivTile\n");
	  return;
      }
      sum2 = 0;
      for ( i = 0; i <  IPListObjectLength(TObj); i++ ) {
	  int tID = AttrGetObjectIntAttrib(TD -> tiles[i] ,"MSTileID");
	  sum = 0;
	  printf(" TILE %d\n", tID);
	  for ( TVs = TD -> tiles[i]->U.Trivars; TVs !=NULL; TVs = TVs-> Pnext ) {
	      sum++;
	      j = AttrGetIntAttrib( TVs -> Attr, "MSGeomID");
	      TD -> IDTrivTile [ j - 1 ] = tID;
	      printf(" TRIV %d IN TILE %d \n", j, tID ) ;
	  }
	  sum2 +=sum;
	  printf(" TILE HAS %d trivars ACC %d\n", sum, sum2 );

      }
      printf(" TOT TRIVARS %d \n", sum2);
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

  // SURFACE 1:
  fil = fopen (dataFile,"r");
  if ( fil == NULL ) {
      fprintf(stderr," FATAL ERROR: Unable to open file %s\n", dataFile);
      return ;
  }

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
  int faceID;
  char *ErrStr;
  printf(" GET DUCT 6 SURFACES \n");
  TrivBndrySrfsFromTVToData( TD -> duct, FALSE, BndSrf );
  for ( j = 0 ; j < 2; j++) {
      if ( j == 0 ) faceID = 4;
      else          faceID = 5;
      fprintf(fout,"%d\n", nP[j] );
      printf( " TOTAL POINTS  %d  \n", nP[j] );
      MAXERR = 0.0;
      snprintf(buf, 100, "BoundingSurf_%d.itd", faceID);
      CagdSrfWriteToFile(BndSrf[faceID], buf, 0, buf,  &ErrStr);
#ifdef DEBUG
      printf ( " SAMPLE SURFACES: TOP AND BOTTOM \n");
      int ii, jj, N = 100;
      CagdRType *P3 = NULL;
      P3  = IritMalloc(CAGD_MAX_PT_SIZE * sizeof(CagdRType));
	  FILE *fs;
      char buff[100];
      double uv[2], Umin, Umax, Vmin, Vmax;
      CagdSrfDomain ( BndSrf[faceID], &Umin, &Umax, &Vmin, &Vmax );
      snprintf(buff,100,"ITDsurface%d.dat", faceID );
      fs   = fopen ( buff, "w" );
      for ( ii = 0 ; ii <= N; ii++ ) {
	  for ( jj = 0 ; jj <= N; jj++ ) {
	      uv[0] = (double ) jj * (Umax - Umin)/ (double)N;
	      uv[1] = (double ) ii * (Vmax - Vmin) / (double)N;
	      CAGD_SRF_EVAL_E3(BndSrf[faceID], uv[0], uv[1], P3 ) ;
	      fprintf(fs, "%lf %lf %lf\n", P3[0], P3[1], P3[2] );
	  }
      }
      fclose ( fs ) ;
#endif
      snprintf(buf, 100, "SURFPTS_%d", faceID);
      filSurf = fopen ( buf, "w" ) ;
      if ( filSurf == NULL ) {
	  fprintf(stderr, "FAILED TO OPEN FILE %s. ABORT\n", buf ) ;
	  return;
      }
      for ( i = 0 ; i < nP[j] ; i++ ) {
	  POINTCOUNT = i + 1;
	  TOTPTS     = nP[j];
	  printf(" CALLING FOR FACE %d \n", faceID );
	  printf(" TRIV %d UVW %f  %f  %f  \n", trivIDs[j][i], pts[j][ 3 * i],
		 pts[j][ 3 * i + 1], pts[j][ 3 * i + 2]);
	  getDuctUVW ( &pts[j][3 * i], trivIDs[j][i], faceID, TD, uvw );
	  //EG_evaluatePressure ( );
	  double press = 0.001;
	  fprintf(fout,"%lf %lf %lf %lf %d\n", pts[j][3*i], pts[j][3*i + 1], pts[j][3*i + 2], press, trivIDs[j][i] );
      }
      fclose ( filSurf ) ;
      fprintf(stderr," MAX ERROR IN INV PROJECTION FOR FACE %d = %.15e\n",faceID, MAXERR);
  }
  fclose ( fil  );
  fclose ( fout );
  return;
}



void getDuctUVW ( double uvw[], int ID, int faceID, TiledDuct *TD, double *uvwDuct ) {
  int i, j, tileID;
  int tileIDbounds[3], ductIDbounds[3];
  clock_t starT, endT;
  double task1, task2, dist;
  CagdRType *P3 = NULL, *PP3 = NULL, uvw0[3];
  CagdBBoxStruct BBox;
  CagdSrfStruct  **BndSrf = NULL;
  TrivTVStruct *TVs;
  MvarPtStruct *seed = MvarPtNew(2), *outSeed = MvarPtNew(2);
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

  tileID = TD -> IDTrivTile[ID -1 ];
  printf(" TRIVAR %d BELONGS TO TILE %d \n", ID, tileID);
  if ( tileID < 0 ) {
      fprintf(stderr, "FATAL ERROR:  I can't find a tile that holds trivariate %d \n", ID);
      exit(1);
  }
  BBox = getTileBBox (TD -> duct, tileID, TD -> tilesPerKnot);
  printf("\n\n TILE %d -> UVW BBOX %f %f %f MAX %f %f %f \n", tileID, BBox.Min[0], BBox.Min[1], BBox.Min[2],
	 BBox.Max[0], BBox.Max[1], BBox.Max[2] );
  // check that indeed we have the right tile:
  j = - 1;
  for ( TVs = TD -> tiles[tileID -1]->U.Trivars; TVs !=NULL; TVs = TVs-> Pnext ) {
      if ( AttrGetIntAttrib( TVs -> Attr, "MSGeomID") == ID ) {
	  j = 1 ; break;
      }
  }
  if ( j == -1 ) {
      fprintf(stderr," FATAL ERROR : TRIVAR %d doesn't belong to tile %d \n", ID, tileID);
      return;
  }
  tileIDbounds[0] = AttrGetIntAttrib( TVs -> Attr, "MSLclTVBndryU" );
  tileIDbounds[1] = AttrGetIntAttrib( TVs -> Attr, "MSLclTVBndryV" );
  tileIDbounds[2] = AttrGetIntAttrib( TVs -> Attr, "MSLclTVBndryW" );
  ductIDbounds[0] = AttrGetIntAttrib( TVs -> Attr, "MSDfrmTVBndryU");
  ductIDbounds[1] = AttrGetIntAttrib( TVs -> Attr, "MSDfrmTVBndryV");
  ductIDbounds[2] = AttrGetIntAttrib( TVs -> Attr, "MSDfrmTVBndryW");
  for ( i = 0 ; i < 3; i++) {
      if ( ductIDbounds[i] != faceID ) continue;
      // get uv0s using tiles centre wrt trivariate
      for ( j = 0 ; j < 3; ++j)
	uvw0[j] = 0.5 * ( BBox.Min[j] + BBox.Max[j] );
      /* ONLY FOR BOTTOM AND TOP SURFACES. WON't WORK IN GENERAL. I AM ASSUMING THAT FACE 4 IS BOTTOM SO
       * INITIAL GUESS USES TILE MIN VALUE AND FACE 5 IS TOP SO INITIAL GUESSS USES TILE MAX
       */
      if ( faceID == 4 ) uvw0[i] = BBox.Min[i];
      else               uvw0[i] = BBox.Max[i];
      TRIV_TV_EVAL_E3( TVs   , uvw [0], uvw [1], uvw [2], P3);
      seed -> Pt[0] = uvw0[0];
      seed -> Pt[1] = uvw0[1];
      handler       = MvarDistSrfPointPrep( BndSrf[ ductIDbounds[i] ] );
      // USING SEED
      starT = clock();
      outSeed       = MvarLclDistSrfPoint    ( BndSrf[ ductIDbounds[i] ], handler, P3, seed, 0.01, NUMERIC_TOL );
      endT = clock();
      task1 = ( endT - starT ) / CLOCKS_PER_SEC;
      //
      uvwDuct[0]    = outSeed -> Pt[0];
      uvwDuct[1]    = outSeed -> Pt[1];
      uvwDuct[2]    = uvw0[2];
#ifdef DEBUG
      printf(" Time Ellapsed  WITH SEED %f s\n", task1 );
      printf(" POINT COORDINATES IN TRIVARIATE: %f  %f  %f \n ", P3[0], P3[1], P3[2]);
      printf(" Initial guess using tile: %f  %f   -> %f   %f   %f (UNUSED)\n", uvw0[0], uvw0[1],  PP3[0], PP3[1], PP3[2]);
      // inv evaluate: get global uvw
      printf(" USING VAR %d  MIN MAX %d -> SURFACE %d \n", i, tileIDbounds[i], ductIDbounds[i] ) ;
      CAGD_SRF_EVAL_E3( BndSrf[ ductIDbounds[i] ], uvwDuct[0], uvwDuct[1], PP3);
      printf(" I have found DUCT Coords  %f  %f  -> %f  %f  %f\n\n\n", uvwDuct[0], uvwDuct[1],
      	     PP3[0], PP3[1], PP3[2]);
      dist = 0.0;
      if ( filSurf ) {
     	  fprintf(filSurf, " %.16f  %.16f   %.16f   %.16f   %.16f   %.16f \n", P3[0], P3[1], P3[2], PP3[0], PP3[1], PP3[2] ) ;
           }
            for ( j = 0 ; j < 3; j++) dist += ( PP3[j] - P3[j] ) * ( PP3[j] - P3[j] );
            dist = sqrt ( dist ) ;
            if ( dist > MAXERR) MAXERR = dist;
            if ( dist > NUMERIC_TOL ) {
      	  fprintf(stderr, "\n\n SURFACE EVALUATION LOCAL AT POINT %d / %d APPROXIMATION FAILED !\n", POINTCOUNT, TOTPTS ) ;
      	  fprintf( stderr, " TRIV %d  IN TILE %d  FACE %d  POINT FROM TRIVAR (%.16f, %.16f, %.16f ) != (%.16f, %.16f, %.16f )  DIST %.16e  TASK %.16f\n",
      		   ID, tileID, faceID, P3[0],  P3[1], P3[2], PP3[0], PP3[1], PP3[2], dist, task2);
      	    }
      TRIV_TV_EVAL_E3(TD ->duct, uvwDuct[0], uvwDuct[1], uvwDuct[2], PP3);
      printf(" I have found DUCT Coords  %f  %f   (%f ~ 0 ) -> %f  %f  %f\n\n\n", uvwDuct[0], uvwDuct[1], uvwDuct[2],
	     PP3[0], PP3[1], PP3[2]);
     dist = 0.0;
      for ( j = 0 ; j < 3; j++) dist += ( PP3[j] - P3[j] ) * ( PP3[j] - P3[j] );
      dist = sqrt ( dist ) ;
      if ( dist > MAXERR) MAXERR = dist;
      if ( dist > NUMERIC_TOL ) {
	  fprintf(stderr, "\n\n TRIVAR EVALUATION LOCAL AT POINT %d / %d APPROXIMATION FAILED !\n", POINTCOUNT, TOTPTS ) ;
	  fprintf( stderr, " TRIV %d  IN TILE %d  FACE %d  POINT FROM TRIVAR (%.16f, %.16f, %.16f ) != (%.16f, %.16f, %.16f )  DIST %.16e  TASK %.16f\n",
		   ID, tileID, faceID, P3[0],  P3[1], P3[2], PP3[0], PP3[1], PP3[2], dist, task2);
	  printf( " AT POINT %d / %d APPROXIMATION FAILED !\n", POINTCOUNT, TOTPTS ) ;
	  	  printf(  " TRIV %d  IN TILE %d  FACE %d  POINT FROM TRIVAR (%.16f, %.16f, %.16f ) != (%.16f, %.16f, %.16f )  DIST %.16e  TASK %.16f\n",
	  		   ID, tileID, faceID, P3[0],  P3[1], P3[2], PP3[0], PP3[1], PP3[2], dist, task2);
      }
#endif
      starT = clock();
      uvwDuct       = MvarDistSrfPoint    ( BndSrf[ ductIDbounds[i] ], handler, P3,  NULL, TRUE, 0.001, NUMERIC_TOL, NULL );
      endT = clock();
      task2 = ( endT - starT ) / CLOCKS_PER_SEC;
#ifdef DEBUG
      printf(" POINT COORDINATES IN TRIVARIATE: %f  %f  %f \n ", P3[0], P3[1], P3[2]);
      printf(" WITHOUT SEEDING \n");
      printf(" USING VAR %d  MIN MAX %d -> SURFACE %d \n", i, tileIDbounds[i], ductIDbounds[i] ) ;
      CAGD_SRF_EVAL_E3( BndSrf[ ductIDbounds[i] ], uvwDuct[0], uvwDuct[1], PP3);
      printf(" I have found DUCT Coords  %f  %f  -> %f  %f  %f\n\n\n", uvwDuct[0], uvwDuct[1],
	     PP3[0], PP3[1], PP3[2]);
      dist = 0.0;
      for ( j = 0 ; j < 3; j++) dist += ( PP3[j] - P3[j] ) * ( PP3[j] - P3[j] );
      dist = sqrt ( dist ) ;
      if ( dist > NUMERIC_TOL ) {
	  fprintf(stderr, " NO SEED AT POINT %d / %d APPROXIMATION FAILED !\n", POINTCOUNT, TOTPTS ) ;
	  fprintf( stderr, " TRIV %d  IN TILE %d  FACE %d  POINT FROM TRIVAR (%.16f, %.16f, %.16f ) != (%.16f, %.16f, %.16f )  DIST %.16e  TASK %.16f\n",
		   ID, tileID, faceID, P3[0],  P3[1], P3[2], PP3[0], PP3[1], PP3[2], dist, task1);
	  //exit (1);
      }
#endif

      MvarDistSrfPointFree(handler);

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
  int tilesPerKnot[3];
  if ( argc < 7 ) {
      printf(" I NEED A DUCT FILE + TILE FILE (*.itd) + Data Points FILE + TILES PER DIR \n");
      return -1;
  }
  tilesPerKnot[0] = atoi ( argv[4] ) ;
  tilesPerKnot[1] = atoi ( argv[5] ) ;
  tilesPerKnot[2] = atoi ( argv[6] ) ;
  updatePressure ( argv[1], argv[2], argv[3], tilesPerKnot );
  return 0;
}





