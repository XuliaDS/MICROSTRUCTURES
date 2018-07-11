/*****************************************************************************
*   Constructs locally varying tiles in microstructure constructions using   *
* a call back function.							     *
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

typedef struct UserMicroLocalDataStruct { /* User specific data in CB funcs. */
    TrivTVStruct *DefMap;
    TrivTVStruct *DuDefMap, *DvDefMap, *DwDefMap;
} UserMicroLocalDataStruct;

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Maps a parameter between zero and one to a radius of the cone to use.    *
*                                                                            *
* PARAMETERS:                                                                *
*   Param:           Parameter between zero and one.                         *
*   MinRad, MaxRad:  Radii amplitudes to allow the sine wave in between.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdRType *:    Radius of cone.					     *
*****************************************************************************/
static CagdRType MapParamToAmplitude(CagdRType Param,
				     CagdRType MinRad,
				     CagdRType MaxRad)
{
    /* Map [0, 1] in Param to a sine period but start and end at MinRad: */
    return (sin(Param * 4 * M_PI - M_PI / 2) + 1.0) * 0.5
						* (MaxRad - MinRad) + MinRad;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Constructs a cone parallel to the Z axis from Z = 0 (radius Rad1) to     *
* Z = 1 (and radius Rad2).  Cone will have no bases.			     *
*                                                                            *
* PARAMETERS:                                                                *
*   Rad1, Rad2:  Radii of two cone bases.                                    *
*                                                                            *
* RETURN VALUE:                                                              *
*   CagdSrfStruct *:    The constructed cone.				     *
*****************************************************************************/
static CagdSrfStruct *MakeZConeTile(CagdRType Rad1, CagdRType Rad2)
{
    static const CagdPtStruct
	Center = { NULL, NULL, { 0.0, 0.0, 0.0 } };
    CagdCrvStruct *TCrv,
        *Circ1 = BspCrvCreatePCircle(&Center, Rad1),
        *Circ2 = BspCrvCreatePCircle(&Center, Rad2);
    CagdSrfStruct *Srf;
    IrtHmgnMatType Mat;

    MatGenMatTrans(0.0, 0.0, 1.0, Mat);
    TCrv = CagdCrvMatTransform(Circ2, Mat);
    CagdCrvFree(Circ2);
    Circ2 = TCrv;

    Srf = CagdRuledSrf(Circ1, Circ2, 2, 2);

    CagdCrvFree(Circ1);
    CagdCrvFree(Circ2);

    return Srf;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Constructs a tile as three open cones in [0, 1]^3 that locally fit the   *
* micro-structure in [UMin, UMax, VMin, VMax, WMIn, WMax], and assuming the  *
* entire trivariate mapping function domains' is [0, 1]^3.                   *
*   The tiles are sized so the entire structure will create a sine wave      *
* approximation.                                                             *
*                                                                            *
* PARAMETERS:                                                                *
*   UMin, UMax, VMin, VMax, WMin, WMax:  Position of this tile in the        *
*	deformation trivariate that is assumed to have domain [0, 1]^3.      *
*   MinRad, MaxRad:  Radii amplitudes to allow the sine wave in between.     *
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:    Tile with three ruled surfaces spanning [0, 1]^3 as *
*			 this local tile for this UVW location.		     *
*****************************************************************************/
static IPObjectStruct *MakeSineCrossTile(CagdRType UMin,
					 CagdRType UMax,
					 CagdRType VMin,
					 CagdRType VMax,
					 CagdRType WMin,
					 CagdRType WMax,
					 CagdRType MinRad,
					 CagdRType MaxRad)
{
    IrtHmgnMatType Mat1, Mat2;
    CagdSrfStruct *TSrf,
        *UCone = MakeZConeTile(MapParamToAmplitude(UMin, MinRad, MaxRad),
			       MapParamToAmplitude(UMax, MinRad, MaxRad)),
        *VCone = MakeZConeTile(MapParamToAmplitude(VMin, MinRad, MaxRad),
			       MapParamToAmplitude(VMax, MinRad, MaxRad)),
        *WCone = MakeZConeTile(MapParamToAmplitude(WMin, MinRad, MaxRad),
			       MapParamToAmplitude(WMax, MinRad, MaxRad));
    IPObjectStruct *PObj;

    /* Transform the cones to their correct position and orientation. */
    MatGenMatRotY1(M_PI / 2.0, Mat1);
    MatGenMatTrans(0.0, 0.5, 0.5, Mat2);
    MatMultTwo4by4(Mat1, Mat1, Mat2);
    TSrf = CagdSrfMatTransform(UCone, Mat1);
    CagdSrfFree(UCone);
    UCone = TSrf;

    MatGenMatRotX1(-M_PI / 2.0, Mat1);
    MatGenMatTrans(0.5, 0.0, 0.5, Mat2);
    MatMultTwo4by4(Mat1, Mat1, Mat2);
    TSrf = CagdSrfMatTransform(VCone, Mat1);
    CagdSrfFree(VCone);
    VCone = TSrf;

    MatGenMatTrans(0.5, 0.5, 0.0, Mat1);
    TSrf = CagdSrfMatTransform(WCone, Mat1);
    CagdSrfFree(WCone);
    WCone = TSrf;

    /* Build a tile from the three surfaces. */
    PObj = IPGenLISTObject(IPGenSRFObject(UCone));
    IPListObjectAppend(PObj, IPGenSRFObject(VCone));
    IPListObjectAppend(PObj, IPGenSRFObject(WCone));

    return PObj;
}

/*****************************************************************************
* DESCRIPTION:                                                               *
*   Called back data just before the tile is mapped - we will substitute in  *
* our local tile and free the old one.                                       *
*                                                                            *
* PARAMETERS:                                                                *
*   Tile:   Canonical tile as provided by the user.	                     *
*   CBData: Call back data as we provide in UserMicroStructComposition below.*
*                                                                            *
* RETURN VALUE:                                                              *
*   IPObjectStruct *:   New tile to use here.                                *
*****************************************************************************/
static IPObjectStruct *PreProcessSineTile(
				     IPObjectStruct *Tile,
				     UserMicroPreProcessTileCBStruct *CBData)
{
    CagdRType *R, UMin, UMax, VMin, VMax, WMin, WMax, Du, Dv, Dw;
    CagdPType DfDw;
    UserMicroLocalDataStruct
	*LclData = (UserMicroLocalDataStruct *) CBData -> CBFuncData;

    TrivTVDomain(LclData -> DefMap, &UMin, &UMax, &VMin, &VMax, &WMin, &WMax);
    R = TrivTVEvalMalloc(LclData -> DwDefMap,
		   CBData -> TileIdxsMin[0],
		   CBData -> TileIdxsMin[1],
		   CBData -> TileIdxsMin[2]);
    CagdCoerceToE3(DfDw, &R, -1, LclData -> DwDefMap -> PType);
    Du = UMax - UMin;
    Dv = UMax - UMin;
    Dw = UMax - UMin;

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

    IPFreeObject(Tile);                       /* Free the old (dummy) tile. */
    /* and create new tile as part of 3 orthogonal sine tubes, in [0, 1]^3: */
    Tile = MakeSineCrossTile((CBData -> TileIdxsMin[0] - UMin) / Du,
			     (CBData -> TileIdxsMax[0] - UMin) / Du,
			     (CBData -> TileIdxsMin[1] - VMin) / Dv,
			     (CBData -> TileIdxsMax[1] - VMin) / Dv,
			     (CBData -> TileIdxsMin[2] - WMin) / Dw,
			     (CBData -> TileIdxsMax[2] - WMin) / Dw,
			     0.1, 0.3);
    /* Map new tile to its proper place in deformation function UVW space. */
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
static void GenerateSineMicroStructure(void)
{
    char *ErrStr;
    int ErrLine, Handler, i;    
    CagdRType
        KeepStackStep = 1,
        TrackAllocID = 0;
    IPObjectStruct *MS;
    MvarMVStruct *DeformMV;
    TrivTVStruct *TVMap;
    UserMicroParamStruct MSParam;
    UserMicroRegularParamStruct *MSRegularParam;
    UserMicroLocalDataStruct LclData;

    TVMap = TrivTVReadFromFile("DefMap.itd", &ErrStr, &ErrLine);
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
    MSParam.ShellCapBits = USER_MICRO_BIT_CAP_UMIN |
                           USER_MICRO_BIT_CAP_UMAX |
                           USER_MICRO_BIT_CAP_VMIN |
                           USER_MICRO_BIT_CAP_VMAX |
                           USER_MICRO_BIT_SHELL_WMIN |
                           USER_MICRO_BIT_CAP_WMAX;
    MSParam.ShellThickness = 0.1;           /* In tile [0, 1]^3 space/size. */
    MSParam.DeformMV = DeformMV;

    MSRegularParam = &MSParam.U.RegularParam;
    MSRegularParam -> Tile = NULL;      /* Dummy tile - will not be used... */
    MSRegularParam -> TilingStepMode = TRUE;

    for (i = 0; i < 3; ++i) {
	MSRegularParam -> TilingSteps[i] = 
			      (CagdRType *) IritMalloc(sizeof(CagdRType) * 2);
	MSRegularParam -> TilingSteps[i][0] = 1;
	MSRegularParam -> TilingSteps[i][1] = 2;
    }
    MSRegularParam -> MaxPolyEdgeLen = 0.1;

    /* Call back function - will be called for each tile in the grid just   */
    /* before it is mapped through the deformation function, with the tile  */
    /* (that can be modified) and call back data.                           */
    MSRegularParam -> PreProcessCBFunc = PreProcessSineTile;
    MSRegularParam -> CBFuncData = &LclData;         /* The call back data. */
    
    MS = UserMicroStructComposition(&MSParam); /* Construct microstructure. */

    MvarMVFree(DeformMV);
    UserMicroTileFree(MSRegularParam -> Tile);
    for (i = 0; i < 3; ++i)
	IritFree(MSRegularParam -> TilingSteps[i]);

    Handler = IPOpenDataFile("MS2.itd", FALSE, 1);
    IPPutObjectToHandler(Handler, MS);
    IPFreeObject(MS);
    IPCloseStream(Handler, TRUE);

    for (i = 0; i < 3; i++) 
	IritFree(MSRegularParam -> TilingSteps[i]);

    /* Free the structure for the call back function. */
    TrivTVFree(TVMap);
    TrivTVFree(LclData.DuDefMap);
    TrivTVFree(LclData.DvDefMap);
    TrivTVFree(LclData.DwDefMap);
}

void main(int argc, char **argv)
{ 
    GenerateSineMicroStructure();
    
    exit(0);
}
