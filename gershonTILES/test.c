#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
    TrivTVStruct *DefMap;
    TrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
} UserMicroLocalDataStruct;


static IPObjectStruct *PreProcessTile(IPObjectStruct *Tile,
				      UserMicroPreProcessTileCBStruct *CBData)
{
    static IPObjectStruct
	*ReadTile = NULL;
    const char *tileFile = "/temp/tile.itd";
    CagdRType  *R, UMin, UMax, VMin, VMax, WMin, WMax;
    CagdPType  DfDw;
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

    TrivTVDomain(LclData -> DefMap, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);
    R = TrivTVEvalToDyn(LclData -> DwDefMap,
		   CBData -> TileIdxsMin[0],
		   CBData -> TileIdxsMin[1],
		   CBData -> TileIdxsMin[2]);
    CagdCoerceToE3(DfDw, &R, -1, LclData -> DwDefMap -> PType);
    IritFree(R);

    fprintf(stderr, "Tile[%d,%d,%d] locally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n\t   globally from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f)\n\t   DfDw = (%.4f %.4f %.4f)\n",
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
	    CBData -> TileIdxsMaxOrig[2],
	    DfDw[0], DfDw[1], DfDw[2]);

    IPFreeObject(Tile);                       /* Free the old (dummy) tile. */
    
    IPSetFlattenObjects(FALSE);
    if (ReadTile == NULL) {
	if ((ReadTile = IPGetDataFiles(&tileFile, 1, 1, 1)) == NULL) {
	    fprintf(stderr, "Failed to load the tile\n");
	    exit(1);
	}
    }
    Tile = GMTransformObject(ReadTile, CBData -> Mat);

    return Tile;
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

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Create a micro structure with a varying-in-size tiling example.          *
*                                                                            *
* PARAMETERS:                                                                *
*   None                                                                     *
*                                                                            *
* RETURN VALUE:                                                              *
*   void                                                                     *
*****************************************************************************/
static void GenerateMicroStructure(void)
{
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

    if ((TVMap = TrivTVReadFromFile("/temp/DefMap.itd",
				     &ErrStr, &ErrLine)) == NULL) {
	 fprintf(stderr, "Failed to load deformation map\n");
	 exit(1);
    }
    DeformMV = MvarCnvrtTVToMV(TVMap);

    /* Create the structure to be passed to the call back function. */
    LclData.DefMap = TVMap;
    LclData.DuDefMap = TrivTVDerive(TVMap, TRIV_CONST_U_DIR);
    LclData.DvDefMap = TrivTVDerive(TVMap, TRIV_CONST_V_DIR);
    LclData.DwDefMap = TrivTVDerive(TVMap, TRIV_CONST_W_DIR);

    IRIT_ZAP_MEM(&MSParam, sizeof(UserMicroParamStruct));
    MSParam.TilingType = USER_MICRO_TILE_REGULAR;

    /* Sets boundary end conditions on the geometry - cap the tiles in five */
    /* boundaries and do a solid shell base at the (WMin) bottom.	    */
    /* No capping/shelling for now.
    MSParam.ShellCapBits = USER_MICRO_BIT_CAP_UMIN |
                           USER_MICRO_BIT_CAP_UMAX |
                           USER_MICRO_BIT_CAP_VMIN |
                           USER_MICRO_BIT_CAP_VMAX |
                           USER_MICRO_BIT_SHELL_WMIN |
                           USER_MICRO_BIT_CAP_WMAX;
    */
    MSParam.ShellThickness = 0.1;           /* In tile [0, 1]^3 space/size. */
    MSParam.DeformMV = DeformMV;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;/* Tile is synthesized by call back func. */
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i] = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i][0] = 1;
	MSRegularParam -> TilingSteps[i][1] = 2;
    }

    /* Call back function - will be called for each tile in the grid just   */
    /* before it is mapped through the deformation function, with the tile  */
    /* (that can be modified) and call back data.                           */
    MSRegularParam -> PreProcessCBFunc = PreProcessTile;
    MSRegularParam -> CBFuncData = &LclData;         /* The call back data. */
    
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

    Handler = IPOpenDataFile("/temp/MS.itd", FALSE, 1); 
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
    GenerateMicroStructure();
    
    exit(0);
}
