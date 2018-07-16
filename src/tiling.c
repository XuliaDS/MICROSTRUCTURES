#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

int CIRC  = 0;
const char *tileFile;
char defMap[100], msFile[100];


typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
  TrivTVStruct *DefMap;
  TrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
  CagdRType fMax[6];
} UserMicroLocalDataStruct;

static CagdRType setMaxVal(  TrivTVStruct  *TV ) {
  CagdRType *P3, uvw[3], max = 0.0, aux;
  CagdPType E3;
  int i, j, k, nT = 10;
  for ( i = 0 ; i <= nT ; i++ ) {
      uvw[0] = (double) i/(double) nT;
      for ( j = 0 ; j <= nT ; j++ ) {
	  uvw[1] = (double) j/(double) nT;
	  for ( k = 0 ; k <= nT ; k++ ) {
	      uvw[2] = (double) k/(double) nT;
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
  CagdRType norm[4], centre[3], w, f[6], EPS = 1.e-4;
  UserMicroLocalDataStruct
  *TVfield = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

  for ( i = 0 ; i < 3; i++ )
    centre[i] = 0.5 * ( CBData -> TileIdxsMax[i] + CBData -> TileIdxsMin[i] ) ;
  for (i = 0 ; i < 6; i++ ) {
      for ( j = 0 ; j < 4; ++j )
	params[i].Bndry[j] = 0.0;
      params[i].Circular   = CIRC;
  }
  f[0] = partialDerMag ( TVfield -> DuDefMap, CBData -> TileIdxsMin [0], centre[1]                , centre[2] );
  f[1] = partialDerMag ( TVfield -> DuDefMap, CBData -> TileIdxsMax [0], centre[1]                , centre[2] );
  f[2] = partialDerMag ( TVfield -> DvDefMap, centre[0]                , CBData -> TileIdxsMin [1], centre[2] );
  f[3] = partialDerMag ( TVfield -> DvDefMap, centre[0]                , CBData -> TileIdxsMax [1], centre[2] );
  f[4] = partialDerMag ( TVfield -> DwDefMap, centre[0]                , centre[1]                , CBData -> TileIdxsMin[2] );
  f[5] = partialDerMag ( TVfield -> DwDefMap, centre[0]                , centre[1]                , CBData -> TileIdxsMax[2] );
  printf(" ------------- FIELD DERIVATIVES ------------\n");
  for ( i = 0 ; i < 3; i ++) {
      printf(" i = %d  ->  param =  %f   Df %f \n", 2*i    , CBData -> TileIdxsMin[i], f[2*i    ] );
      printf(" i = %d  ->  param =  %f   Df %f \n", 2*i + 1, CBData -> TileIdxsMax[i], f[2*i + 1] );
  }
  printf("---------------------------------------------\n");
  // Use Trivariate derivatives dT/du, dT/dV, dT/dW to measure the field "strength"
  // MAX VARIATION IN [0 1]
  printf( " MAX FIELD Du %f  Dv %f  Dw %f \n", TVfield ->fMax[0], TVfield ->fMax[2], TVfield ->fMax[4]);
  // MAX "LOAD" VARYING EACH FACE
  //  STUDY THE DOMAIN BOUNDARY AND DETERMINE THE WALL THICKNESS

  for ( i = 0 ; i < 2; i++ ) {
      if ( i == 0 ) w = CBData -> TileIdxsMinOrig[2];
      else          w = CBData -> TileIdxsMaxOrig[2];
      printf(" w = %f  \n", w );
      if      ( w <= EPS || w >= 1.0 - EPS ) {
	  printf(" IS BOUNDARY \n");
	  norm[0] = partialDerMag ( TVfield -> DwDefMap, CBData -> TileIdxsMin[0], CBData -> TileIdxsMin[1], w );
	  norm[1] = partialDerMag ( TVfield -> DwDefMap, CBData -> TileIdxsMax[0], CBData -> TileIdxsMin[1], w );
	  norm[2] = partialDerMag ( TVfield -> DwDefMap, CBData -> TileIdxsMin[0], CBData -> TileIdxsMax[1], w );
	  norm[3] = partialDerMag ( TVfield -> DwDefMap, CBData -> TileIdxsMax[0], CBData -> TileIdxsMax[1], w );
	  for ( j = 0 ; j < 4; ++j )
	    params[4 + i].Bndry[j] = norm[j] / TVfield ->fMax[i] * 0.25;
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
      params[i].OuterRadius = 0.1;//0.05 * ( 2*j + 1 ) * f[i]/TVfield ->fMax[i] ;
  }

}


static IPObjectStruct *PreProcessTile( IPObjectStruct *Tile,  UserMicroPreProcessTileCBStruct *CBData )
{
  int i ;
  CagdPType  DfDwMin, DfDwMax, DfDw0, DfDw1;
  UserMicroTileBndryPrmStruct *UVWparams = NULL;
  // BOUND BOX [0,1]^3 //

  fprintf(stderr, "Tile[%d,%d,%d] locally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n\t   globally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n",
	  CBData -> TileIdxs[0],
	  CBData -> TileIdxs[1],
	  CBData -> TileIdxs[2],
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

  IPFreeObject(Tile);                       /* Free the old (dummy) tile. */
  UVWparams = malloc ( 6 * sizeof(UserMicroTileBndryPrmStruct) );
  if ( UVWparams == NULL ) {
      fprintf(stderr, "Failed to allocate memo\n");
      return Tile;
  }
  getFieldConditions (CBData, UVWparams);
  for ( i = 0 ; i < 6; i ++ ) {
      printf(" TILE CONDS %d : CIRCULAR %d \n"
	  "BDRY %lf %lf %lf %lf RADII OUTER  %lf  INNER  %lf \n",
	  i, UVWparams[i].Circular,
	  UVWparams[i].Bndry[0],    UVWparams[i].Bndry[1],
	  UVWparams[i].Bndry[2],    UVWparams[i].Bndry[3],
	  UVWparams[i].OuterRadius, UVWparams[i].InnerRadius );
  }
  Tile = UserMicro3DCrossTile(&UVWparams[0], &UVWparams[1], &UVWparams[2],
			      &UVWparams[3], &UVWparams[4], &UVWparams[5] );
  Tile = GMTransformObject(Tile, CBData -> Mat);

  return Tile;
}

static int GenerateMicroStructure(int nu, int nv, int nw) {
  char *ErrStr;
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

  if ((TVMap = TrivTVReadFromFile(defMap,
				  &ErrStr, &ErrLine)) == NULL) {
      fprintf(stderr, "Failed to load deformation map\n");
      printf(" ERRSTR %d  ERRLINE %d  defmAP %s TVMPA %p \n", ErrLine, Handler, defMap, TVMap);
      return -1;
  }
  DeformMV = MvarCnvrtTVToMV(TVMap);

  /* Create the structure to be passed to the call back function. */
  LclData.DefMap   = TVMap;
  LclData.DuDefMap = TrivTVDerive(TVMap, TRIV_CONST_U_DIR);
  LclData.DvDefMap = TrivTVDerive(TVMap, TRIV_CONST_V_DIR);
  LclData.DwDefMap = TrivTVDerive(TVMap, TRIV_CONST_W_DIR);

  LclData.fMax[0] = setMaxVal ( LclData.DuDefMap );
  LclData.fMax[1] = LclData.fMax[0];
  LclData.fMax[2] = setMaxVal ( LclData.DvDefMap );
  LclData.fMax[3] = LclData.fMax[2];
  LclData.fMax[4] = setMaxVal ( LclData.DwDefMap );
  LclData.fMax[5] = LclData.fMax[4];

  IRIT_ZAP_MEM(&MSParam, sizeof(UserMicroParamStruct));
  MSParam.TilingType = USER_MICRO_TILE_REGULAR;

  MSParam.ShellThickness           = 0.1;           /* In tile [0, 1]^3 space/size. */
  MSParam.DeformMV                 = DeformMV;
  MSRegularParam                   = &MSParam.U.RegularParam;
  MSRegularParam -> Tile           = NULL;/* Tile is synthesized by call back func. */
  MSRegularParam -> TilingStepMode = TRUE;
  MSRegularParam -> MaxPolyEdgeLen = 0.1;

  for (i = 0; i < 3; ++i) {
      MSRegularParam -> TilingSteps[i]    =
	  (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
      MSRegularParam -> TilingSteps[i][0] = 1;
  }
  MSRegularParam -> TilingSteps[0][1] = nu;
  MSRegularParam -> TilingSteps[1][1] = nv;
  MSRegularParam -> TilingSteps[2][1] = nw;


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

  MvarMVFree(DeformMV);
  UserMicroTileFree(MSRegularParam -> Tile);
  for (i = 0; i < 3; ++i)
    IritFree(MSRegularParam -> TilingSteps[i]);

  Handler = IPOpenDataFile(msFile, FALSE, 1);
  IPPutObjectToHandler(Handler, MS);
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
  int nu, nv, nw;
  if ( argc < 6 ) {
      printf(" I NEED A DOMAIN FILE (*.itd) + number of tiles per dir + circ (1) or squared (0) \n");
      return -1;
  }
  snprintf(defMap,100,"%s",argv[1]);
  snprintf(msFile,100,"tiled");
  nu   = atoi ( argv[2] ) ;
  nv   = atoi ( argv[3] ) ;
  nw   = atoi ( argv[4] ) ;
  CIRC = atoi ( argv[5] );
  GenerateMicroStructure(nu, nv, nw);
  return 0;
}
