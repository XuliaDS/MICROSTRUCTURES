
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"
#include <nlopt.h>

int CIRC  = 0;
int BOUNDSKIN = 0 ;
char tiledName[100], ductName[100];
#define DEBUG
#define EPS04  1.e-4

typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
  TrivTVStruct *DefMap;
  TrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
  CagdRType fMax[6], boundingBox[6];
} UserMicroLocalDataStruct;


typedef struct UserTotalMicroStruct {
  UserMicroRegularParamStruct *MSRP;
  UserMicroParamStruct MSP;
  UserMicroLocalDataStruct LDS;
  int nTperKnot[3], nTperDir[3];
} UserTotalMicroStruct ;

static CagdRType setMaxVal(  TrivTVStruct  *TV, CagdRType *uvwMinMax ) {
  CagdRType *P3, uvw[3], max = 0.0, aux;
  CagdPType E3;
  int i, j, k, nT = 10;
  for ( i = 0 ; i <= nT ; i++ ) {
      uvw[0] = uvwMinMax[0] + (uvwMinMax[1] - uvwMinMax[0] ) * (double) i/(double) nT;
      for ( j = 0 ; j <= nT ; j++ ) {
	  uvw[1] = uvwMinMax[2] + (uvwMinMax[3] - uvwMinMax[2] ) * (double) j/(double) nT;
	  for ( k = 0 ; k <= nT ; k++ ) {
	      uvw[2] = uvwMinMax[4] + (uvwMinMax[5] - uvwMinMax[4] ) * (double) k/(double) nT;
	      P3     = TrivTVEvalMalloc(TV, uvw[0], uvw[1], uvw[2] );
	      CagdCoerceToE3(E3, &P3 , -1, TV -> PType);
	      aux    = sqrt ( E3[0] * E3[0] + E3[1] * E3[1] + E3[2] * E3[2]);
	      if ( aux > max ) max = aux;
	  }
      }
  }
  return max;
}

/*****************************************************************************
 * DESCRIPTION:                                                               *
 *   Converts a hierarchy of geometries into a linear list of geometries.     *
 *                                                                            *
 *                                                                            *
 * PARAMETERS:                                                                *
 *   MS:     Input hierarchy.                                                 *
 *   LinMS:  The resulting equivalent linear geometry.                        *
 *                                                                            *
 * RETURN VALUE:                                                              *
 *   void                                                                     *
 *****************************************************************************/
static void FlattenHierarchy2LinearList(const IPObjectStruct *MS,
					IPObjectStruct *LinMS)
{
  int i;
  switch (MS -> ObjType) {
    case IP_OBJ_LIST_OBJ:
      for (i = 0; i < IPListObjectLength(MS); i++)
	FlattenHierarchy2LinearList(IPListObjectGet(MS, i), LinMS);
      break;
    default:
      IPListObjectAppend(LinMS, IPCopyObject(NULL, MS, TRUE));
      break;
  }
}

static double partialDerMag ( TrivTVStruct *field, CagdRType u,  CagdRType v,  CagdRType w ) {
  CagdRType *P3, f;
  CagdPType E3;
  P3 = TrivTVEvalMalloc( field, u, v, w ) ;
  CagdCoerceToE3(E3, &P3 , -1, field -> PType);
  f = sqrt ( E3[0] * E3[0] + E3[1] * E3[1] + E3[2] * E3[2] ) ;
  return f;
}

static void getFieldConditions( UserMicroPreProcessTileCBStruct *CBData ,  UserMicroTileBndryPrmStruct *params )
{
  int       i, j, isBound = 0;
  CagdRType norm[4], centre[3], w, f[6], uvw_min[3], uvw_max[3], boundingBox[6] ;
  UserMicroLocalDataStruct
  *TVfield = (UserMicroLocalDataStruct *) CBData -> CBFuncData;
  /* Find coordinates in the deformation map: GLOBAL */
  for ( i = 0 ; i < 3; i++ ) {
      uvw_min[i] = (1.0 - CBData -> TileIdxsMin[i] ) * CBData -> TileIdxsMinOrig[i] +
	  CBData -> TileIdxsMin[i] * CBData -> TileIdxsMaxOrig[i];
      uvw_max[i] = (1.0 - CBData -> TileIdxsMax[i] ) * CBData -> TileIdxsMinOrig[i] +
	  CBData -> TileIdxsMax[i] * CBData -> TileIdxsMaxOrig[i];
      centre [i] = 0.5 * ( uvw_min[i] + uvw_max[i]) ;
  }

#ifdef DEBUGG
  printf(" Tiles Coords %f  %f  %f  and %f  %f  %f\n",
	 CBData -> TileIdxsMin[0],  CBData -> TileIdxsMin[1],  CBData -> TileIdxsMin[2],
	 CBData -> TileIdxsMax[0],  CBData -> TileIdxsMax[1],  CBData -> TileIdxsMax[2]);
  printf(" ORI Tiles Coords %f  %f  %f  and %f  %f  %f\n",
	 CBData -> TileIdxsMinOrig[0],  CBData -> TileIdxsMinOrig[1],  CBData -> TileIdxsMinOrig[2],
	 CBData -> TileIdxsMaxOrig[0],  CBData -> TileIdxsMaxOrig[1],  CBData -> TileIdxsMaxOrig[2]);
  printf(" GLOBAL COORDINATES: %f  %f  %f  and  %f  %f  %f\n",
	 uvw_min[0], uvw_min[1], uvw_min[2], uvw_max[0], uvw_max[1], uvw_max[2] );
#endif

  for (i = 0 ; i < 6; i++ ) {
      for ( j = 0 ; j < 4; ++j )
	params[i].Bndry[j] = 0.0;
      params[i].Circular   = CIRC;
  }
  f[0] = partialDerMag ( TVfield -> DuDefMap, uvw_min[0], centre[1] , centre[2] );
  f[1] = partialDerMag ( TVfield -> DuDefMap, uvw_max[0], centre[1] , centre[2] );
  f[2] = partialDerMag ( TVfield -> DvDefMap, centre[0] , uvw_min[1], centre[2] );
  f[3] = partialDerMag ( TVfield -> DvDefMap, centre[0] , uvw_max[1], centre[2] );
  f[4] = partialDerMag ( TVfield -> DwDefMap, centre[0] , centre[1] , uvw_min[2] ) ;
  f[5] = partialDerMag ( TVfield -> DwDefMap, centre[0] , centre[1] , uvw_max[2] ) ;

#ifdef DEBUGG
  printf(" ------------- FIELD DERIVATIVES ------------\n");
  for ( i = 0 ; i < 3; i ++) {
      printf(" i = %d  ->  param =  %f   Df %f \n", 2*i    , uvw_min[i], f[2*i    ] );
      printf(" i = %d  ->  param =  %f   Df %f \n", 2*i + 1, uvw_max[i], f[2*i + 1] );
  }
  printf("---------------------------------------------\n");
#endif

  /* Use Trivariate derivatives dT/du, dT/dV, dT/dW to measure the field "strength"
   MAX VARIATION IN [0 1]  MAX "LOAD" VARYING EACH FACE
   STUDY THE DOMAIN BOUNDARY AND DETERMINE THE WALL THICKNESS */

  for ( i = 0 ; i < 2; i++ ) {
      if ( i == 0 ) w = uvw_min[2];
      else          w = uvw_max[2];
      if  ( fabs ( w - TVfield -> boundingBox[4] ) <  EPS04 ||
	  fabs ( w - TVfield -> boundingBox[5] ) <  EPS04  ) {
#ifdef DEBUGG
	  printf(" w %f   bound box %f  %f \n", w , TVfield -> boundingBox[4] , TVfield -> boundingBox[5] );
#endif

	  norm[0] = partialDerMag ( TVfield -> DwDefMap, uvw_min[0], uvw_min[1], w );
	  norm[1] = partialDerMag ( TVfield -> DwDefMap, uvw_max[0], uvw_min[1], w );
	  norm[2] = partialDerMag ( TVfield -> DwDefMap, uvw_min[0], uvw_max[1], w );
	  norm[3] = partialDerMag ( TVfield -> DwDefMap, uvw_max[0], uvw_max[1], w );
	  /* If you set boundaries to value, it won't print attributes. Test it commenting the for loop below */
	  if ( BOUNDSKIN == 1 ) {
	      for ( j = 0 ; j < 4; ++j )
		params[4 + i].Bndry[j] = 0.1 ;//norm[j] / TVfield -> fMax[i] * 0.25;
	  }
	  isBound = 1;
      }
  }
  for ( j = i = 0 ; i < 6; i++ ) {
      // for wing tiles, hollowed tiles are very slow to visualize. Uncomment in the future
      //if ( isBound )
      //      params[i].InnerRadius = 0.0;
      //else
      if ( i %2 == 0 && i > 0 ) j++;
      params[i].InnerRadius = 0.0;//0.03 *  (  j  + 1)* f[i]/TVfield ->fMax[i] ;
      params[i].OuterRadius = 0.1 * f[i]/TVfield ->fMax[i] ;
  }

}


static IPObjectStruct *PreProcessTile( IPObjectStruct *Tile,  UserMicroPreProcessTileCBStruct *CBData )
{
  int i ;
  CagdPType  DfDwMin, DfDwMax, DfDw0, DfDw1;
  UserMicroTileBndryPrmStruct *UVWparams = NULL;
  // BOUND BOX [0,1]^3 //

#ifdef DEBUGG
  printf("\n==========================  Tile[%d,%d,%d]  ===========\n",
	 CBData -> TileIdxs[0],  CBData -> TileIdxs[1],  CBData -> TileIdxs[2]);
  printf(" Locally from (%0.2f, %0.2f %0.2f) to (%0.2f, %0.2f, %.2f)\t || globally from (%0.2f, %0.2f %0.2f) to (%0.2f, %0.2f, %.2f)\n",
	 CBData -> TileIdxsMin[0],
	 CBData -> TileIdxsMin[1],
	 CBData -> TileIdxsMin[2],
	 CBData -> TileIdxsMax[0],
	 CBData -> TileIdxsMax[1],
	 CBData -> TileIdxsMax[2],
	 CBData -> TileIdxsMinOrig[0],
	 CBData -> TileIdxsMinOrig[1],
	 CBData -> TileIdxsMinOrig[2],
	 CBData -> TileIdxsMaxOrig[0],
	 CBData -> TileIdxsMaxOrig[1],
	 CBData -> TileIdxsMaxOrig[2]);
#endif

  IPFreeObject(Tile);                       /* Free the old (dummy) tile. */
  UVWparams = (UserMicroTileBndryPrmStruct*) malloc ( 6 * sizeof(UserMicroTileBndryPrmStruct) );
  IRIT_ZAP_MEM(UVWparams, 6 * sizeof(UserMicroTileBndryPrmStruct) );
  if ( UVWparams == NULL ) {
      fprintf(stderr, "Failed to allocate memo\n");
      return Tile;
  }
  getFieldConditions (CBData, UVWparams);

#ifdef DEBUGG
  for ( i = 0 ; i < 6; i ++ ) {
      printf(" Parameters  %d : CIRCULAR %d \t"
	  "BDRY %lf %lf %lf %lf RADII OUTER  %lf  INNER  %lf \n",
	  i, UVWparams[i].Circular,
	  UVWparams[i].Bndry[0],    UVWparams[i].Bndry[1],
	  UVWparams[i].Bndry[2],    UVWparams[i].Bndry[3],
	  UVWparams[i].OuterRadius, UVWparams[i].InnerRadius );
  }
  printf("==========================================================================\n");
#endif
  Tile = UserMicro3DCrossTile(&UVWparams[0], &UVWparams[1], &UVWparams[2],
			      &UVWparams[3], &UVWparams[4], &UVWparams[5] );

  free (UVWparams);
  Tile = GMTransformObject(Tile, CBData -> Mat);
  return Tile;
}


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


void getDuctUVW ( double uvw[], int ID,   TrivTVStruct *TVduct, int tilesPerKnot[], IPObjectStruct *MS, double *uvwDuct ) {
  int i, j, tileID;
  int tileIDbounds[3], ductIDbounds[3];
  CagdRType *P3 = NULL, *PP3 = NULL, uvw0[3];
  IPObjectStruct *tile;
  TrivTVStruct *TV;
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
  TrivBndrySrfsFromTVToData( TVduct, FALSE, BndSrf );
  // get tile ID
  tileID = getTileIDFromTriv ( MS, ID, 0 );
  if ( tileID < 0 ) {
      fprintf(stderr, "FATAL ERROR:  I can't find a tile that holds trivariate %d \n", ID);
      exit(1);
  }
  tile = IPListObjectGet (MS, tileID ) ;
  BBox = getTileBBox (TVduct, tileID + 1, tilesPerKnot);
  printf("\n\n TILE %d -> UVW BBOX %f %f %f MAX %f %f %f \n", tileID, BBox.Min[0], BBox.Min[1], BBox.Min[2],
	 BBox.Max[0], BBox.Max[1], BBox.Max[2] );
  // check that indeed we have the right tile:
  for ( TV = tile -> U.Trivars; TV !=NULL; TV = TV -> Pnext ) {
      if (  AttrGetIntAttrib( TV -> Attr, "MSGeomID") == ID )  break;
  }
  tileIDbounds[0] = AttrGetIntAttrib( TV -> Attr, "MSDfrmTVBndryU");
  tileIDbounds[1] = AttrGetIntAttrib( TV -> Attr, "MSDfrmTVBndryV");
  tileIDbounds[2] = AttrGetIntAttrib( TV -> Attr, "MSDfrmTVBndryW");
  ductIDbounds[0] = AttrGetIntAttrib( TV -> Attr, "MSLclTVBndryU" );
  ductIDbounds[1] = AttrGetIntAttrib( TV -> Attr, "MSLclTVBndryV" );
  ductIDbounds[2] = AttrGetIntAttrib( TV -> Attr, "MSLclTVBndryW" );
  for ( i = 0 ; i < 3; i++) {
      if ( tileIDbounds[i] < 0 ) continue;
      printf("\n\n TILE ID BOUND %d ( SURFACE %d ) \n", tileIDbounds[i], ductIDbounds[i]);
      printf(" LOOKING FOR %f  %f  %f  from TRIV (%d) in DUCT \n", uvw[0], uvw[1], uvw[2], ID );
      // get uv0s using tiles centre wrt trivariate
      for ( j = 0 ; j < 3; ++j)
	uvw0[j] = 0.5 * ( BBox.Min[j] + BBox.Max[j] );
      if ( tileIDbounds[i] % 2 == 0 ) uvw0[i] = BBox.Min[i];
      else                            uvw0[i] = BBox.Max[i];
      TRIV_TV_EVAL_E3(TV    , uvw [0], uvw [1], uvw [2], P3);
      TRIV_TV_EVAL_E3(TVduct, uvw0[0], uvw0[1], uvw0[2], PP3); // unused ATM
      handler = MvarDistSrfPointPrep( BndSrf[ ductIDbounds[i] ] );
      uvwDuct = MvarDistSrfPoint    ( BndSrf[ ductIDbounds[i] ], handler, P3, TRUE, 0.0001, 1.e-10, NULL );
      MvarDistSrfPointFree(handler);
#ifdef DEBUG
      printf(" POINT COORDINATES IN TRIVARIATE: %f  %f  %f \n ", P3[0], P3[1], P3[2]);
      printf(" Initial guess using tile: %f  %f  %f  -> %f   %f   %f (UNUSED)\n", uvw0[0], uvw0[1], uvw0[2], PP3[0], PP3[1], PP3[2]);
      // inv evaluate: get global uvw
      printf(" USING VAR %d  MIN MAX %d -> SURFACE %d \n", i, tileIDbounds[i], ductIDbounds[i] ) ;
      TRIV_TV_EVAL_E3(TVduct, uvwDuct[0], uvwDuct[1], uvwDuct[2], PP3);
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





void GenerateMicroStructure(int tilesPerKnot[3] ) {
  char *ErrStr;
  const char *name = ductName;
  int i, j, ErrLine, Handler, IDTriv, IDbounds[3], IDboundsLoc[3];
  IPObjectStruct *pTmp;
  TrivTVStruct   *TVout;
  CagdRType  KeepStackStep = 1,  TrackAllocID = 0;
  IPObjectStruct *MS, *LinMS;
  MvarMVStruct *DeformMV;
  TrivTVStruct *TVMap;
  UserMicroParamStruct MSParam;
  UserMicroRegularParamStruct *MSRegularParam;
  UserMicroLocalDataStruct LclData;
  IPObjectStruct *TObj = IPGetDataFiles(&name, 1, TRUE, TRUE);
  if ( TObj == NULL ) {
      printf(" NULL OBJECT. Failed to load duct \n");
      return;
  }
  TVMap    = TObj -> U.Trivars;
  DeformMV = MvarCnvrtTVToMV(TVMap);

  /* Create the structure to be passed to the call back function. */
  LclData.DefMap   = TVMap;
  LclData.DuDefMap = TrivTVDerive(TVMap, TRIV_CONST_U_DIR);
  LclData.DvDefMap = TrivTVDerive(TVMap, TRIV_CONST_V_DIR);
  LclData.DwDefMap = TrivTVDerive(TVMap, TRIV_CONST_W_DIR);
  TrivTVDomain(LclData.DefMap, &LclData.boundingBox[0], &LclData.boundingBox[1],
	       &LclData.boundingBox[2], &LclData.boundingBox[3],
	       &LclData.boundingBox[4], &LclData.boundingBox[5]);

#ifdef DEBUGG
  printf(" BOUNDING BOX  u [%f x %f] v [%f x %f] w [%f x %f] \n",
	 LclData.boundingBox[0],LclData.boundingBox[1],LclData.boundingBox[2],
	 LclData.boundingBox[3],LclData.boundingBox[4],LclData.boundingBox[5]);
#endif
  LclData.fMax[0] = setMaxVal ( LclData.DuDefMap , LclData.boundingBox);
  LclData.fMax[1] = LclData.fMax[0];
  LclData.fMax[2] = setMaxVal ( LclData.DvDefMap , LclData.boundingBox);
  LclData.fMax[3] = LclData.fMax[2];
  LclData.fMax[4] = setMaxVal ( LclData.DwDefMap , LclData.boundingBox);
  LclData.fMax[5] = LclData.fMax[4];

  IRIT_ZAP_MEM(&MSParam, sizeof(UserMicroParamStruct));
  MSParam.TilingType = USER_MICRO_TILE_REGULAR;

  MSParam.ShellThickness           = 0.1;           /* In tile [0, 1]^3 space/size. */
  MSParam.DeformMV                 = DeformMV;
  MSRegularParam                   = &MSParam.U.RegularParam;
  MSRegularParam -> Tile           = NULL;/* Tile is synthesized by call back func. */
  MSRegularParam -> TilingStepMode = TRUE;
  MSRegularParam -> MaxPolyEdgeLen = 0.1;
  MSParam.ApproxLowOrder           = 4 ;
  MSParam. ShellCapBits            = USER_MICRO_BIT_CAP_ALL;
  for (i = 0; i < 3; ++i) {
      MSRegularParam -> TilingSteps[i]    =
	  (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
      MSRegularParam -> TilingSteps[i][0] = 1;
  }
  MSRegularParam -> TilingSteps[0][1] = tilesPerKnot[0];
  MSRegularParam -> TilingSteps[1][1] = tilesPerKnot[1];
  MSRegularParam -> TilingSteps[2][1] = tilesPerKnot[2];


  /* Call back function - will be called for each tile in the grid just   */
  /* before it is mapped through the deformation function, with the tile  */
  /* (that can be modified) and call back data.                           */
  MSRegularParam -> PreProcessCBFunc = PreProcessTile;
  MSRegularParam -> CBFuncData       = &LclData;         /* The call back data. */

  MS = UserMicroStructComposition(&MSParam); /* Construct microstructure. */

  LinMS = IPGenLISTObject(NULL);     /* Create an empty list to populate. */
  FlattenHierarchy2LinearList(MS, LinMS);  /* Convert into a linear list. */
  IPFreeObject(MS);
  MS = LinMS;
  fprintf(stderr, "%d surfaces created\n", IPListObjectLength(MS));
  Handler = IPOpenDataFile(tiledName, FALSE, 1);
  IPPutObjectToHandler(Handler, MS);
  /* Pretend to have a uvw parameters list. Use boundary trivars centres as evaluation.
   * We are goint to obtain global uvw (duct coordinates)
   */
  for ( i = 0; i <  IPListObjectLength(MS); i++ ) {
      pTmp = IPListObjectGet (MS, i ) ;
      int ID = AttrGetObjectIntAttrib(pTmp,"MSTileID");
      int k = 0 , ii;
      for ( TVout = pTmp -> U.Trivars; TVout !=NULL; TVout = TVout -> Pnext ) {
	  IDTriv         = AttrGetIntAttrib( TVout -> Attr, "MSGeomID"      );
	  IDbounds   [0] = AttrGetIntAttrib( TVout -> Attr, "MSDfrmTVBndryU");
	  IDbounds   [1] = AttrGetIntAttrib( TVout -> Attr, "MSDfrmTVBndryV");
	  IDbounds   [2] = AttrGetIntAttrib( TVout -> Attr, "MSDfrmTVBndryW");
	  IDboundsLoc[0] = AttrGetIntAttrib( TVout -> Attr, "MSLclTVBndryU" );
	  IDboundsLoc[1] = AttrGetIntAttrib( TVout -> Attr, "MSLclTVBndryV" );
	  IDboundsLoc[2] = AttrGetIntAttrib( TVout -> Attr, "MSLclTVBndryW" );
	  double uvw[3], uvwDuct[3];
	  for ( j = 0 ; j < 3; ++j ) {
	      if ( IDbounds[j] >= 0 ) {
		  printf("----------- BOUND %d -> (At surf %d )\n", IDbounds[j], IDboundsLoc[j]);
		  uvw[0] = 0.5 ; uvw[1] = 0.5 ; uvw[2] = 0.5 ;
		  if ( IDbounds[j] %2 == 0 ) uvw[j] = 0.0;
		  else uvw[j] = 1.0;
		  getDuctUVW( uvw, IDTriv, TVMap, tilesPerKnot, MS, uvwDuct );
		  printf("-------------------------------------\n");
	      }
	  }
      }
  }
  MvarMVFree(DeformMV);
  UserMicroTileFree(MSRegularParam -> Tile);
  for (i = 0; i < 3; ++i)
    IritFree(MSRegularParam -> TilingSteps[i]);

  IPFreeObject(MS);
  IPCloseStream(Handler, TRUE);
  /* Free the structure for the call back function. */
  TrivTVFree(TVMap);
  TrivTVFree(LclData.DuDefMap);
  TrivTVFree(LclData.DvDefMap);
  TrivTVFree(LclData.DwDefMap);


}




int main(int argc, char **argv)
{
  int nT[3];
  if ( argc < 6 ) {
      printf(" I NEED A DOMAIN FILE (*.itd) + number of tiles per dir + circ (1) or squared (0) \n");
      return -1;
  }
  snprintf( ductName, 100,"%s", argv[1]);
  snprintf(tiledName, 100,"tiled_%s", ductName);
  nT[0] = atoi ( argv[2] ) ;
  nT[1] = atoi ( argv[3] ) ;
  nT[2] = atoi ( argv[4] ) ;
  CIRC  = atoi ( argv[5] );
  BOUNDSKIN = 1;
  GenerateMicroStructure(nT);
  return 0;
}








