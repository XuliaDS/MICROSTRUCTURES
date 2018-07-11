#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"


const char *tileFile;
char defMap[100], msFile[100];


typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
  TrivTVStruct *DefMap;
  TrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
} UserMicroLocalDataStruct;


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

static void getFieldConditions( UserMicroLocalDataStruct *TVfield,  CagdRType tMin[], CagdRType tMax[],  UserMicroTileBndryPrmStruct *params )
{
  // params = 6D structure: params[0] = face_u0, params[1] = face_u_1; params[2] = face_v0; params[3] = face_v1, etc.
  int       i, j, isBound = 0;
  CagdRType norm[4], centre[3], w, f[6];
  for (i = 0 ; i < 6; i++ ) {
      for ( j = 0 ; j < 4; ++j )
	params[i].Bndry[j] = 0.0;
      params[i].Circular   = FALSE;
  }
  for ( i = 0 ; i < 2; i++ ) {
      if ( i == 0 ) w = tMin[2];
      else          w = tMax[2];
      if      ( w <= 0.0 || w >= 1.0 ) {
	  for ( j = 0 ; j < 4; ++j )
	    params[4 + i].Bndry[j] = 0.1;
      }
  }
  for ( j = i = 0 ; i < 6; i++ ) {
      params[i].InnerRadius = 0.0;
      if ( i %2 == 0 && i > 0 ) ++j;
      if ( i %2 == 0 ) params[i].OuterRadius   = 0.1 + 0.1 * tMin[j];
      else             params[i].OuterRadius   = 0.1 + 0.1 * tMax[j];
  }
}


static IPObjectStruct *PreProcessTile( IPObjectStruct *Tile,  UserMicroPreProcessTileCBStruct *CBData )
{
  int i ;
  CagdPType  DfDwMin, DfDwMax, DfDw0, DfDw1;
  UserMicroLocalDataStruct
  *LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;
  UserMicroTileBndryPrmStruct *UVWparams = NULL;
  // BOUND BOX [0,1]^3 //

  printf("\n\n\n\n Tile[%d,%d,%d] locally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n\t   "
      "globally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n ",
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
  getFieldConditions (LclData, CBData -> TileIdxsMin,CBData -> TileIdxsMax, UVWparams);
  for ( i = 0 ; i < 6; i ++ ) {
      printf(" TILE CONDS %d : CIRCULAR %d \t RADII OUTER  %lf  INNER  %lf \t"
	  "BDRY %lf %lf %lf %lf \n", i,
	  UVWparams[i].Circular, UVWparams[i].OuterRadius,  UVWparams[i].InnerRadius,
	  UVWparams[i].Bndry[0], UVWparams[i].Bndry[1],
	  UVWparams[i].Bndry[2], UVWparams[i].Bndry[3]
      );
  }
  Tile = UserMicro3DCrossTile(&UVWparams[0], &UVWparams[1], &UVWparams[2],
			      &UVWparams[3], &UVWparams[4], &UVWparams[5] );
  Tile = GMTransformObject(Tile, CBData -> Mat);

  return Tile;
}

static int GenerateMicroStructure() {
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
  MSRegularParam -> TilingSteps[0][1] = 10;
  MSRegularParam -> TilingSteps[1][1] = 8;
  MSRegularParam -> TilingSteps[2][1] = 5;



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
  if ( argc < 1 ) {
      printf(" I NEED A DOMAIN FILE (*.itd) \n");
  }
  snprintf(defMap,100,"%s",argv[1]);
  tileFile = argv[2];
  printf (" DEF MAP %s\n", defMap);
  printf( " TILE %s \n",tileFile);
  snprintf(msFile,100,"tiled");
  // call tiles
  GenerateMicroStructure();
  return 0;
}
