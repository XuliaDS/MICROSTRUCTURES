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


static IPObjectStruct *PreProcessTile(IPObjectStruct *Tile,
				      UserMicroPreProcessTileCBStruct *CBData)
{
    static IPObjectStruct
	*ReadTile = NULL;
    CagdRType  *R, UMin, UMax, VMin, VMax, WMin, WMax, Du, Dv, Dw;
    CagdPType  DfDw;
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

    TrivTVDomain(LclData -> DefMap, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);
    R = TrivTVEval(LclData -> DwDefMap,
		   CBData -> UVWMin[0], CBData -> UVWMin[1], CBData -> UVWMin[2]);
    CagdCoerceToE3(DfDw, &R, -1, LclData -> DwDefMap -> PType);
    Du = UMax - UMin;
    Dv = VMax - VMin;
    Dw = WMax - WMin;

    fprintf(stderr, "Tile[%d,%d,%d] from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f), DfDw = (%.4f %.4f %.4f)\n",
	    CBData -> UVWIndices[0],
	    CBData -> UVWIndices[1],
	    CBData -> UVWIndices[2],
	    CBData -> UVWMin[0], CBData -> UVWMin[1], CBData -> UVWMin[2],
	    CBData -> UVWMax[0], CBData -> UVWMax[1], CBData -> UVWMax[2],
	    DfDw[0], DfDw[1], DfDw[2]);

    IPFreeObject(Tile);                       /* Free the old (dummy) tile. */
    
    IPSetFlattenObjects(FALSE);
    if (ReadTile == NULL)
	ReadTile = IPGetDataFiles(&tileFile, 1, 1, 1);
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
    int ErrLine, Handler;    
    CagdRType
        KeepStackStep = 1,
        TrackAllocID = 0;
    IPObjectStruct *MS, *LinMS;
    MvarMVStruct *DeformMV;
    TrivTVStruct *TVMap;
    UserMicroParamStruct MSParam;
    UserMicroRegularParamStruct *MSRegularParam;
    UserMicroLocalDataStruct LclData;

    TVMap = TrivTVReadFromFile("DefMap.itd", &ErrStr, &ErrLine);
    DeformMV = MvarTVToMV(TVMap);

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
    MSRegularParam -> Tile =            /* Dummy tile - will not be used... */
        UserMicroTileNew(IPGenCRVObject(BspCrvCreateUnitPCircle()));
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> TilingSteps[0] = 
	MSRegularParam -> TilingSteps[1] = 
	    MSRegularParam -> TilingSteps[2] = 5;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

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
   GenerateMicroStructure();
   if ( argc != 3 ) {
      printf(" I NEED A DOMAIN FILE (*.itd) AND A TILE FILE (*.itd) \n");
    }
    snprintf(defMap,100,"%s",argv[1]);
    tileFile = argv[2];
    printf (" DEF MAP %s\n", defMap);
    exit(0);
}
