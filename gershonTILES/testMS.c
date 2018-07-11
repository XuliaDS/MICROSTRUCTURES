/*****************************************************************************
* Microstructures test example.						     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, Dec 2017    *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

const char *TileFileName;

char defMap[100], msFile[100];
/* Create a tile from a list of surfaces.		*/
UserMicroTileStruct *UserMicroCreateTileFrmSrfs(CagdSrfStruct *SrfList)
{
    IPObjectStruct *IPObj,
	*IPListObj = IPGenLISTObject(NULL);
    UserMicroTileStruct *Tile;
    CagdSrfStruct *Srf, *SrfNext;

    for (Srf = SrfList; Srf != NULL; Srf = SrfNext) {
	IPObj = IPGenSRFObject(Srf);
	IPListObjectAppend(IPListObj, IPObj);
	SrfNext = Srf -> Pnext;
	Srf -> Pnext = NULL;
    }

    Tile = UserMicroTileNew(IPListObj);
    return Tile;
}

/* Create a tile from a list of trimmed surfaces.	*/
UserMicroTileStruct *UserMicroCreateTileFrmTrimSrfs(TrimSrfStruct *TrmSrfList)
{
    IPObjectStruct *IPObj,
	*IPListObj = IPGenLISTObject(NULL);
    UserMicroTileStruct *Tile;
    TrimSrfStruct *TrmSrf, *TrmSrfNext;

    for (TrmSrf = TrmSrfList; TrmSrf != NULL; TrmSrf = TrmSrfNext) {
	IPObj = IPGenTRIMSRFObject(TrmSrf);
	IPListObjectAppend(IPListObj, IPObj);
	TrmSrfNext = TrmSrf -> Pnext;
	TrmSrf -> Pnext = NULL;
    }

    Tile = UserMicroTileNew(IPListObj);

    return Tile;
}

/* Generate a list of two surfaces.		*/
CagdSrfStruct *GenSrfList(void) 
{
    CagdSrfStruct *Srf, *Srf2;
    CagdVType Center;

    Center[0] = Center[1] = 0.5;
    Center[2] = 0;

    Srf = CagdPrimCylinderSrf(Center, 0.5, 1, FALSE, CAGD_PRIM_CAPS_NONE);
    Srf2 = CagdPrimCylinderSrf(Center, 0.25, 1, FALSE, CAGD_PRIM_CAPS_NONE);
    Srf -> Pnext = Srf2;

    return Srf;
}

/* Generate a list of two trimmed-surfaces.	*/
TrimSrfStruct *GenTrimSrfList(void)
{
    TrimSrfStruct *TrmSrf, *TrmSrf2;
    TrimCrvStruct *TrCrv;
    TrimCrvSegStruct *TrCrvSeg;
    CagdSrfStruct *Srf;
    CagdVType Center;
    CagdPtStruct CircCenter;
    CagdCrvStruct *TCrv;

    Center[0] = Center[1] = 0.5;
    Center[2] = 0;

    IRIT_PT_COPY(CircCenter.Pt, Center);

    Srf = CagdPrimCylinderSrf(Center, 0.5, 1, FALSE, CAGD_PRIM_CAPS_NONE);
    TCrv = BspCrvCreatePCircle(&CircCenter, 0.45);    
    TrCrvSeg = TrimCrvSegNew(TCrv, NULL);
    TrCrv = TrimCrvNew(TrCrvSeg);
    TrmSrf = TrimSrfNew(Srf, TrCrv, FALSE);


    Srf = CagdPrimCylinderSrf(Center, 0.25, 1, FALSE, CAGD_PRIM_CAPS_NONE);
    TCrv = BspCrvCreatePCircle(&CircCenter, 0.25);    
    TrCrvSeg = TrimCrvSegNew(TCrv, NULL);
    TrCrv = TrimCrvNew(TrCrvSeg);
    TrmSrf2 = TrimSrfNew(Srf, TrCrv, FALSE);

    TrmSrf -> Pnext = TrmSrf2;

    return TrmSrf;
}

/* A tiling example.			*/
void TestFunction(CagdBType TrimFlag)
{
    char *ErrStr;
    int i, ErrLine, Handler;    
    CagdRType
        Start         = 1,
        KeepStackStep = 1,
        TrackAllocID  = 0;
    CagdSrfStruct *SrfList;
    TrimSrfStruct *TrimSrfList;
    IPObjectStruct *MS;
    UserMicroTileStruct *Tile;
    MvarMVStruct *DeformMV;
    TrivTVStruct *TVMap1;
    UserMicroParamStruct MSParam;
    UserMicroRegularParamStruct *MSRegularParam;

#ifdef DEBUG_TEST_DYN_MEM
    IritDynMemoryDbgCheckMark(&Start, &KeepStackStep, &TrackAllocID);
#endif /* DEBUG_TEST_DYN_MEM */

    IPSetFlattenObjects(FALSE);

    TVMap1 = TrivTVReadFromFile(defMap, &ErrStr, &ErrLine);
    if (TVMap1 == NULL) {
        fprintf(stderr, "Failed to load the deformation function\n");
	exit(1);
    }

    if (TrimFlag == 0) {
	/*	Surface list example.		*/
	SrfList = GenSrfList();
	Tile = UserMicroCreateTileFrmSrfs(SrfList);
    }
    else if (TrimFlag == 1){
	/*	Trimmed surface list example.    */
	TrimSrfList = GenTrimSrfList();
	Tile = UserMicroCreateTileFrmTrimSrfs(TrimSrfList);
    }
    else {
        /*      Load a tile.                     */
        IPObjectStruct
	    *PObjTile = IPGetDataFiles(&TileFileName, 1, TRUE, TRUE);

	if (PObjTile == NULL) {
	    fprintf(stderr, "Failed to load the tile\n");
	    exit(1);
	}
	Tile = UserMicroTileNew(PObjTile);
    }

    DeformMV = MvarCnvrtTVToMV(TVMap1);
    TrivTVFree(TVMap1);

    IRIT_ZAP_MEM(&MSParam, sizeof(UserMicroParamStruct));
    MSParam.DeformMV = DeformMV;
    MSParam.TilingType = USER_MICRO_TILE_REGULAR;
    MSParam.ShellThickness = 0;
    MSParam.ShellCapBits = 0;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = Tile;
    MSRegularParam -> TilingStepMode = TRUE;
    MSRegularParam -> C0DiscontWScale = 0;
    MSRegularParam -> MaxPolyEdgeLen = 0.1;
    MSRegularParam -> ApproxLowOrder = FALSE;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i] = 
	    (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i][0] = 1;
	MSRegularParam -> TilingSteps[i][1] = i + 2;
    }

    MS = UserMicroStructComposition(&MSParam);

    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i]);

    Handler = IPOpenDataFile(msFile, FALSE, 1);
    IPPutObjectToHandler(Handler, MS);
    IPFreeObject(MS);
    UserMicroTileFree(Tile);
    IPCloseStream(Handler, TRUE);

#ifdef DEBUG_TEST_DYN_MEM
    Start = 0;
    IritDynMemoryDbgCheckMark(&Start, &KeepStackStep, &TrackAllocID);
#endif /* DEBUG_TEST_DYN_MEM */
}

void main(int argc, char **argv)
{
  if ( argc != 3)
  printf(" ENTER DEFORMATION FIELD & tile FILE\n");
    snprintf(defMap,100,"%s",argv[1]);
    TileFileName = argv[2];
    snprintf(msFile,100, "XMSTEST0.itd");
    TestFunction(0);
    exit(0);
}
