#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"
//#include <nlopt.h>

int CIRC  = 0;
char tiledName[100], ductName[100];
//#define DEBUG
#define EPS04  1.e-4

typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
  TrivTVStruct *DefMap;
  TrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
  CagdRType fMax[6], boundingBox[6];
  CagdRType *RoutUVW, *RinUVW, *wBdry;
} UserMicroLocalDataStruct;


typedef struct UserTotalMicroStruct {
  UserMicroRegularParamStruct *MSRP;
  UserMicroParamStruct MSP;
  UserMicroLocalDataStruct LDS;
  int nTperKnot[3], nTperDir[3];
  IPObjectStruct *MS;
} UserTotalMicroStruct ;


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



static IPObjectStruct *PreProcessTile( IPObjectStruct *Tile,  UserMicroPreProcessTileCBStruct *CBData ) {
  int i, i, j, isBound = 0,tileID = 0, wID = 0  ;
  CagdPType  DfDwMin, DfDwMax, DfDw0, DfDw1;
  UserMicroTileBndryPrmStruct *UVWparams = NULL;
  CagdRType  centre[3], w, uvw_min[3], uvw_max[3];
  UserMicroLocalDataStruct
  *TVfield = (UserMicroLocalDataStruct *) CBData -> CBFuncData;
  IPFreeObject(Tile);                       /* Free the old (dummy) tile. */
  IRIT_ZAP_MEM(UVWparams, 6 * sizeof(UserMicroTileBndryPrmStruct));


  if ( UVWparams == NULL ) {
      fprintf(stderr, "Failed to allocate memo\n");
      return Tile;
  }
  for ( i = 0 ; i < 6; i++ ) {
      UVWparams[i].OuterRadius = CBData -> RoutUVW[ tileID + i];
  }

  /* Find coordinates in the deformation map: GLOBAL */
  for ( i = 0 ; i < 3; i++ ) {
      uvw_min[i] = (1.0 - CBData -> TileIdxsMin[i] ) * CBData -> TileIdxsMinOrig[i] +
	  CBData -> TileIdxsMin[i] * CBData -> TileIdxsMaxOrig[i];
      uvw_max[i] = (1.0 - CBData -> TileIdxsMax[i] ) * CBData -> TileIdxsMinOrig[i] +
	  CBData -> TileIdxsMax[i] * CBData -> TileIdxsMaxOrig[i];
      centre [i] = 0.5 * ( uvw_min[i] + uvw_max[i]) ;
  }
  for ( i = 0 ; i < 2; i++ ) {
      if ( i == 0 ) w = uvw_min[2];
      else          w = uvw_max[2];
      if  ( fabs ( w - TVfield -> boundingBox[4] ) <  EPS04 ||
	  fabs ( w - TVfield -> boundingBox[5] ) <  EPS04  ) {
	  printf(" w %f   bound box %f  %f \n", w , TVfield -> boundingBox[4] , TVfield -> boundingBox[5] );

	  for ( j = 0 ; j < 4; ++j )
	    UVWparams[4 + i].Bndry[j] =CBData -> wBdry[tileID + j];
      }
  }



  // BOUND BOX [0,1]^3 //
#ifdef DEBUG
  printf("\n==========================  Tile[%d,%d,%d] ([ %d, %d, %d] )  ===========\n",
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
  Tile = UserMicro3DCrossTile(&UVWparams[0], &UVWparams[1], &UVWparams[2],
			      &UVWparams[3], &UVWparams[4], &UVWparams[5] );
  Tile = GMTransformObject(Tile, CBData -> Mat);

  return Tile;
}



void setTypeParams ( CagdRType *sigma, int *nT, CagdRType *radius ) {
  int i, j, k, m, sum;
  for (  i = 0 ; i < nT[0]; i++) {
      for ( j = 0 ; j < nT[1]; j++) {
	  for ( k = 0 ; k < nT[2]; k++) {
	      sum = k * ( nT[0] * nT[1] ) + j * nT[1] + i;
	      for ( m = 0 ; m < 3; m++ )
		radius[ sum + 2 * m ] =  sigma[0];
	      if ( k > 0 && k < nT[2] -1 )
		radius[ sum + 5 ] =  radius [ ( k - 1)  * ( nT[0] * nT[1] ) +   j        * nT[1] + i        + 4 ];
	      if ( j > 0 && j < nT[1] -1 )
		radius[ sum + 3 ] =  radius [  k        * ( nT[0] * nT[1] ) + ( j - 1  ) * nT[1] + i        + 2 ];
	      if ( i > 0 && i < nT[0] -1 )
		radius[ sum + 1 ] =  radius [  k        * ( nT[0] * nT[1] ) +   j        * nT[1] + ( i - 1)  ];
	  }
      }
  }
}

double optimizeTiles (  unsigned n, double *Ruvw, /*@unused@*/ grad, void *inPtr  ) {
  int i, nT;
  CagdRType *sigma;
  UserTotalMicroStruct *tileInfo;
  tileInfo = (UserTotalMicroStruct*) inPtr;
  tileInfo -> MS = UserMicroStructComposition(&tileInfo -> MSP); /* Construct microstructure. */
  IPObjectStruct *LinMS;
  LinMS = IPGenLISTObject(NULL);     /* Create an empty list to populate. */
  FlattenHierarchy2LinearList(tileInfo -> MS, LinMS);  /* Convert into a linear list. */
  IPFreeObject(tileInfo -> MS);
  tileInfo -> MS = LinMS;
  nT             = IPListObjectLength(tileInfo -> MS);
  fprintf(stderr, "%d surfaces created\n", IPListObjectLength(tileInfo -> MS));

  double error;
  // EVALUATE FIELD
  IRIT_ZAP_MEM(&sigma,  6 * nT * sizeof ( CagdRType ) ) ;

  for ( i = 0 ; i< IPListObjectLength(tileInfo -> MS); i++ ) {
      //getIGAresultForTile( i )
      printf(" Calling IGA to re-set radius \n");
      setTilesParams ( sigma, nT, Ruvw ) ;
  }

  // Get Estimated Ruvw

  return error;
}




void GenerateMicroStructure(int nu, int nv, int nw) {
  char *ErrStr;
  const char *name = ductName;
  int i, ErrLine, Handler, NLOPTvars;
  CagdRType  KeepStackStep = 1,  TrackAllocID = 0, *epsABS = NULL, *lowerBounds = NULL, *upperBounds = NULL, fOpti, *Ruvw = NULL;
  IPObjectStruct *LinMS;
  MvarMVStruct *DeformMV;
  UserTotalMicroStruct *tileInfo ;
  nlopt_opt opt;
  nlopt_result res ;
  IPObjectStruct *TObj = IPGetDataFiles(&name, 1, TRUE, TRUE);
  if ( TObj == NULL ) {
      printf(" NULL OBJECT. Failed to load duct \n");
      return;
  }
  tileInfo -> LDS.DefMap    = TObj -> U.Trivars;
  DeformMV = MvarCnvrtTVToMV( tileInfo -> LDS.DefMap );

  /* Create the structure to be passed to the call back function. */
  tileInfo -> nTperKnot[0] = nu;
  tileInfo -> nTperKnot[1] = nv;
  tileInfo -> nTperKnot[2] = nw;
  tileInfo -> nTperDir [0] = tileInfo -> LDS.DefMap -> ULength * nu;
  tileInfo -> nTperDir [1] = tileInfo -> LDS.DefMap -> VLength * nv;
  tileInfo -> nTperDir [2] = tileInfo -> LDS.DefMap -> WLength * nw;
  for ( i = 0; i < 3; i++)
    printf(" DIR %d TILES PER KNOT %d TOT KNOTS %d  --> TILES %d \n", i, tileInfo -> nTperKnot[i], tileInfo -> nTperDir[i]);
  tileInfo -> LDS.DuDefMap = TrivTVDerive(tileInfo -> LDS.DefMap, TRIV_CONST_U_DIR);
  tileInfo -> LDS.DvDefMap = TrivTVDerive(tileInfo -> LDS.DefMap, TRIV_CONST_V_DIR);
  tileInfo -> LDS.DwDefMap = TrivTVDerive(tileInfo -> LDS.DefMap, TRIV_CONST_W_DIR);
  TrivTVDomain(   tileInfo -> LDS.DefMap,
		  &tileInfo -> LDS.boundingBox[0], &tileInfo -> LDS.boundingBox[1],
		  &tileInfo -> LDS.boundingBox[2], &tileInfo -> LDS.boundingBox[3],
		  &tileInfo -> LDS.boundingBox[4], &tileInfo -> LDS.boundingBox[5]);

#ifdef DEBUG
  printf(" BOUNDING BOX  u [%f x %f] v [%f x %f] w [%f x %f] \n",
	 tileInfo -> boundingBox[0], tileInfo -> boundingBox[1],
	 tileInfo -> boundingBox[2], tileInfo -> boundingBox[3],
	 tileInfo -> boundingBox[4], tileInfo -> boundingBox[5]);
#endif

  tileInfo -> LDS.fMax[0] = setMaxVal ( tileInfo -> LDS.DuDefMap ,
					tileInfo -> LDS.boundingBox);
  tileInfo -> LDS.fMax[1] = tileInfo -> LDS.fMax[0];
  tileInfo -> LDS.fMax[2] = setMaxVal ( tileInfo -> LDS.DvDefMap ,
					tileInfo -> LDS.boundingBox);
  tileInfo -> LDS.fMax[3] = tileInfo -> LDS.fMax[2];
  tileInfo -> LDS.fMax[4] = setMaxVal ( tileInfo -> LDS.DwDefMap ,
					tileInfo -> LDS.boundingBox);
  tileInfo -> LDS.fMax[5] = tileInfo -> LDS.fMax[4];

  IRIT_ZAP_MEM(&tileInfo -> MSP, sizeof(UserMicroParamStruct));
  tileInfo -> MSP.TilingType = USER_MICRO_TILE_REGULAR;

  tileInfo -> MSP.ShellThickness     = 0.1;           /* In tile [0, 1]^3 space/size. */
  tileInfo -> MSP.DeformMV           = DeformMV;
  tileInfo -> MSRP                   = &tileInfo -> MSP.U.RegularParam;
  tileInfo -> MSRP -> Tile           = NULL;/* Tile is synthesized by call back func. */
  tileInfo -> MSRP -> TilingStepMode = TRUE;
  tileInfo -> MSRP -> MaxPolyEdgeLen = 0.1;
  tileInfo -> MSRP -> ApproxLowOrder = 4 ;
  tileInfo -> MSP. ShellCapBits      = USER_MICRO_BIT_CAP_ALL;
  for (i = 0; i < 3; ++i) {
      tileInfo -> MSRP -> TilingSteps[i]    =
	  (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
      tileInfo -> MSRP -> TilingSteps[i][0] = tileInfo -> nTperKnot[i];
  }

  /* Call back function - will be called for each tile in the grid just   */
  /* before it is mapped through the deformation function, with the tile  */
  /* (that can be modified) and call back data.                           */
  NLOPTvars                            = nu + nv + nw ;
  tileInfo -> MSRP -> PreProcessCBFunc = PreProcessTile;
  tileInfo -> MSRP -> CBFuncData       = &tileInfo -> LDS;         /* The call back data. */
  tileInfo -> LDS  -> RoutUVW          = IRIT_ZAP_MEM(Ruvw, 6 * 4 * NLOPTvars * sizeof (CagdRType));
  tileInfo -> LDS  -> RInUVW           = IRIT_ZAP_MEM(Ruvw, 6 * 4 * NLOPTvars * sizeof (CagdRType));
  lowerBounds                          = IRIT_ZAP_MEM(Ruvw, 2     * NLOPTvars * sizeof (CagdRType));
  upperBounds                          = IRIT_ZAP_MEM(Ruvw, 2     * NLOPTvars * sizeof (CagdRType));
  epsABS                               = IRIT_ZAP_MEM(Ruvw, 2     * NLOPTvars * sizeof (CagdRType));
  if (lowerBounds == NULL || upperBounds == NULL || Ruvw == NULL || epsABS == NULL )  {
      fprintf(stderr, " MALLOC PROBLEM IN GenerateMicroStructure\n");
      return;
  }
  opt    = nlopt_create(NLOPT_LN_COBYLA, 2 * NLOPTvars );
  for ( i = 0 ; i < 2 * NLOPTvars; ++i) epsABS[i] = EPS04;
  i = nlopt_set_lower_bounds (opt, lowerBounds);
  i = nlopt_set_upper_bounds (opt, upperBounds);

  fOpti = optimizeTiles ( NLOPTvars, Ruvw, NULL, tileInfo ) ;

  LinMS = IPGenLISTObject(NULL);     /* Create an empty list to populate. */
  FlattenHierarchy2LinearList(tileInfo -> MS, LinMS);  /* Convert into a linear list. */
  IPFreeObject(tileInfo -> MS);
  tileInfo -> MS = LinMS;
  fprintf(stderr, "%d surfaces created\n", IPListObjectLength(tileInfo -> MS));

  MvarMVFree(DeformMV);
  UserMicroTileFree(tileInfo -> MSRP -> Tile);
  for (i = 0; i < 3; ++i)
    IritFree(tileInfo -> MSRP -> TilingSteps[i]);

  Handler = IPOpenDataFile(tiledName, FALSE, 1);
  IPPutObjectToHandler(Handler, tileInfo -> MS);
  IPFreeObject(tileInfo -> MS);
  IPCloseStream(Handler, TRUE);

  /* Free the structure for the call back function. */
  TrivTVFree(TVMap);
  TrivTVFree(tileInfo -> LDS.DuDefMap);
  TrivTVFree(tileInfo -> LDS.DvDefMap);
  TrivTVFree(tileInfo -> LDS.DwDefMap);

}









int main(int argc, char **argv)
{
  int nu, nv, nw;
  if ( argc < 6 ) {
      printf(" I NEED A DOMAIN FILE (*.itd) + number of tiles per dir + circ (1) or squared (0) \n");
      return -1;
  }
  snprintf( ductName, 100,"%s", argv[1]);
  snprintf(tiledName, 100,"tiled_%s", ductName);
  nu   = atoi ( argv[2] ) ;
  nv   = atoi ( argv[3] ) ;
  nw   = atoi ( argv[4] ) ;
  CIRC = atoi ( argv[5] );
  GenerateMicroStructure(nu, nv, nw);
  return 0;
}
