
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
      printf(" HOL    TOTAL TILES ?  %d \n", IPListObjectLength(MS));
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

#ifdef DEBUG
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

#ifdef DEBUG
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
#ifdef DEBUG
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

#ifdef DEBUG
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
  UVWparams = malloc ( 6 * sizeof(UserMicroTileBndryPrmStruct) );
  IRIT_ZAP_MEM(UVWparams, 6 * sizeof(UserMicroTileBndryPrmStruct));
  if ( UVWparams == NULL ) {
      fprintf(stderr, "Failed to allocate memo\n");
      return Tile;
  }
  getFieldConditions (CBData, UVWparams);

#ifdef DEBUG
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
  Tile = GMTransformObject(Tile, CBData -> Mat);

  return Tile;
}



CagdBBoxStruct getTileBBox ( TrivTVStruct *TVAR, int id, int tilesPerKnot[3] ) {
  int ijk[3], m, n, d, res1, knotsPerDir[3],totalTiles[3];
  CagdRType knotCoords[2], tileCoords[2];
  CagdBBoxStruct BBox;
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

int ULength, VLength, WLength;/* Mesh size in tri-variate tensor product.*/
int UVPlane;	  /* Should equal ULength * VLength for fast access. */
int UOrder, VOrder, WOrder;      /* Order in trivariate (B-spline only). */
CagdBType UPeriodic, VPeriodic, WPeriodic;   /* Valid only for B-spline. */
CagdRType *Points[CAGD_MAX_PT_SIZE];     /* Pointer on each axis vector. */
CagdRType *UKnotVector, *VKnotVector, *WKnotVector;

void GenerateMicroStructure(int tilesPerKnot[3] ) {
  char *ErrStr;
  const char *name = ductName;
  int i, ErrLine, Handler;
  CagdRType
  KeepStackStep = 1,
  TrackAllocID = 0;
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

#ifdef DEBUG
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
  IPObjectStruct *pTmp;
  TrivTVStruct   *TVout;
  CagdBBoxStruct BBox, TrivBBox;
  CagdSrfStruct  **BndrySrf, **TrivarBnds;
  int IDbounds[3], IDboundsLoc[3], j, IDTriv;
  TrivBndrySrfsFromTVToData( TVMap, FALSE, BndrySrf);
  for (i = 0; i < 6; i++) {
      CagdSrfBBox(BndrySrf[i], &BBox);
      printf(" SURFACE %d BOUND BOX : %f %f  %f  %f  %f  %f \n",
	     i, BBox.Min[0], BBox.Max[0], BBox.Min[1], BBox.Max[1], BBox.Min[2], BBox.Max[2] ) ;
  }
  Handler = IPOpenDataFile(tiledName, FALSE, 1);
    IPPutObjectToHandler(Handler, MS);

  for ( i = 0; i <  IPListObjectLength(MS); i++ ) {
      pTmp = IPListObjectGet (MS, i ) ;
      int ID = AttrGetObjectIntAttrib(pTmp,"MSIndexID");
      BBox = getTileBBox ( TVMap, i, tilesPerKnot);
      printf(" TILE %d ->  BOUND BOX : %f %f  %f  %f  %f  %f \n", i,
	     BBox.Min[0], BBox.Max[0], BBox.Min[1], BBox.Max[1], BBox.Min[2], BBox.Max[2] ) ;
      int k = 0 , ii;
      for ( TVout = pTmp -> U.Trivars; TVout !=NULL; TVout = TVout -> Pnext ) {
	  TrivBndrySrfsFromTVToData( TVout, FALSE, TrivarBnds );
	  CagdSrfBBox(TrivarBnds[i], &TrivBBox);
	        printf(" TRIVAR BBOX: %f %f  %f  %f  %f  %f \n",
	  	     TrivBBox.Min[0], TrivBBox.Max[0], TrivBBox.Min[1], TrivBBox.Max[1],
		     TrivBBox.Min[2], TrivBBox.Max[2] ) ;

	  printf(" NEW TRIVAR %d \n",k++);
	  int IDTriv      = AttrGetIntAttrib( TVout -> Attr, "MSGeomID");
	  IDbounds   [0] = AttrGetIntAttrib( TVout -> Attr, "MSDfrmTVBndryU");
	  IDbounds   [1] = AttrGetIntAttrib( TVout -> Attr, "MSDfrmTVBndryV");
	  IDbounds   [2] = AttrGetIntAttrib( TVout -> Attr, "MSDfrmTVBndryW");
	  IDboundsLoc[0] = AttrGetIntAttrib( TVout -> Attr, "MSLclTVBndryU");
	  IDboundsLoc[1] = AttrGetIntAttrib( TVout -> Attr, "MSLclTVBndryV");
	  IDboundsLoc[2] = AttrGetIntAttrib( TVout -> Attr, "MSLclTVBndryW");
	  for ( j = 0 ; j < 3; ++j )
	    if ( IDbounds[j] >= 0 ) printf(" BOUND %d -> TVB %d LOCAL %d\n", j, IDbounds[j], IDboundsLoc[j]);

	}
      //TRIV_TV_EVAL_E3(TV, uvw[0], uvw[1], uvw[2], P3 );
      //CAGD_SRF_EVAL_E3 ( SURF, u, v , P3);


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








