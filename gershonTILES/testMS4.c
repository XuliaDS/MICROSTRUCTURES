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
    const char *tileFile = "/temp/tile.itd";
    CagdRType  *R, UMin, UMax, VMin, VMax, WMin, WMax, Du, Dv, Dw;
    CagdPType  DfDw;
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

    TrivTVDomain(LclData -> DefMap, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);
    R = TrivTVEval(LclData -> DwDefMap,
		   CBData -> TileIdxsMin[0],
		   CBData -> TileIdxsMin[1],
		   CBData -> TileIdxsMin[2]);
    CagdCoerceToE3(DfDw, &R, -1, LclData -> DwDefMap -> PType);
    Du = UMax - UMin;
    Dv = VMax - VMin;
    Dw = WMax - WMin;

    fprintf(stderr, "Tile[%d,%d,%d] from (%.3f, %.3f %.3f) to (%.3f, %.3f, %.4f), DfDw = (%.4f %.4f %.4f)\n",
	    CBData -> TileIdxs[0],
	    CBData -> TileIdxs[1],
	    CBData -> TileIdxs[2],
	    CBData -> TileIdxsMin[0],
	    CBData -> TileIdxsMin[1],
	    CBData -> TileIdxsMin[2],
	    CBData -> TileIdxsMax[0],
	    CBData -> TileIdxsMax[1],
	    CBData -> TileIdxsMax[2],
	    DfDw[0], DfDw[1], DfDw[2]);

    IPSetFlattenObjects(FALSE);
    Tile = IPGetDataFiles(&tileFile, 1, 1, 1);
    Tile = GMTransformObjectInPlace(Tile, CBData -> Mat);

    return Tile;
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
    IPObjectStruct *MS;
    MvarMVStruct *DeformMV;
    TrivTVStruct *TVMap;
    UserMicroParamStruct MSParam;
    UserMicroRegularParamStruct *MSRegularParam;
    UserMicroLocalDataStruct LclData;

    TVMap = TrivTVReadFromFile("/temp/DefMap.itd", &ErrStr, &ErrLine);
    DeformMV = MvarCnvrtTVToMV(TVMap);

    /* Create the structure to be passed to the call back function. */
    LclData.DefMap = TVMap;
    LclData.DuDefMap = TrivTVDerive(TVMap, TRIV_CONST_U_DIR);
    LclData.DvDefMap = TrivTVDerive(TVMap, TRIV_CONST_V_DIR);
    LclData.DwDefMap = TrivTVDerive(TVMap, TRIV_CONST_W_DIR);

    IRIT_ZAP_MEM(&MSParam, sizeof(UserMicroParamStruct));
    MSParam.TilingType = USER_MICRO_TILE_REGULAR;

    /* Sets boundary end conditions on the geometry - cap the tiles in all  */
    /* boundaries and as a side effect color trivar tiles on boundaries,    */
    /* so one can set boundary conditions (i.e. toward analysis).           */
    MSParam.ShellCapBits = USER_MICRO_BIT_CAP_UMIN |
                           USER_MICRO_BIT_CAP_UMAX |
                           USER_MICRO_BIT_CAP_VMIN |
                           USER_MICRO_BIT_CAP_VMAX |
                           USER_MICRO_BIT_CAP_WMIN |
			   USER_MICRO_BIT_CAP_WMAX;
    MSParam.ShellThickness = 0.1;           /* In tile [0, 1]^3 space/size. */
    MSParam.DeformMV = DeformMV;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;       /* Tile is synthesized on the fly. */
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

    MvarMVFree(DeformMV);
    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i]);
    UserMicroTileFree(MSRegularParam -> Tile);

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
