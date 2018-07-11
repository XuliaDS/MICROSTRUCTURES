/******************************************************************************
* BzrZrFct.c - Finding zeroes of a Bezier curve using factoring.		      *
*******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                 *
*******************************************************************************
* Written by Jinesh Machchhar, July 2015.				      *
******************************************************************************/

#include <ctype.h>
#include <stdio.h>
#include <string.h>

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/grap_lib.h"
#include "inc_irit/misc_lib.h"
#include "inc_irit/symb_lib.h"
#include "inc_irit/user_lib.h"
#include "inc_irit/mvar_lib.h"
#include "inc_irit/usertopo.h"



/* Create a vector of points. */
static UserTopoUnstrctGeomPtStruct *CreatePtVec(CagdPType Origin,
						CagdRType DistXYZ[3],
						int NumPtXYZ[3], 
						int FirstID)
{
    int x, y, z, 
	i = 0;
    CagdPType Pt;
    UserTopoUnstrctGeomPtStruct
	*PtVec = (UserTopoUnstrctGeomPtStruct *)
	IritMalloc(sizeof(UserTopoUnstrctGeomPtStruct) * (NumPtXYZ[0] * NumPtXYZ[1] * NumPtXYZ[2]));

    IRIT_PT_COPY(Pt, Origin);

    for (z = 0, Pt[2] = Origin[2]; z < NumPtXYZ[2]; z++, Pt[2] += DistXYZ[2]) {
	for (y = 0, Pt[1] = Origin[1]; y < NumPtXYZ[1]; y++, Pt[1] += DistXYZ[1]) {
	    for (x = 0, Pt[0] = Origin[0]; x < NumPtXYZ[0]; x++, Pt[0] += DistXYZ[0]) {
		IRIT_PT_COPY(PtVec[i].Pt.Pt, Pt);
		PtVec[i].ID = FirstID;
		PtVec[i].Attr = NULL;
		PtVec[i].Pt.Attr = NULL;
		FirstID++;
		i++;
	    }
	}
    }

    return PtVec;
}

/* Create a linear curve and add to the grid. */
static void AddCrvToGrid(UserTopoUnstrctGeomStruct *Ud,
		         int FirstPtId,
		         int NumPtXYZ[3],
		         int LineType,
		         int LineValue1,
		         int LineValue2)
{
    int Ind, NumPt, *PtIdVec,
	x = 0,
	y = 0,
	z = 0,
	CrvInd = 0;
    CagdCrvStruct *Crv;
    IPObjectStruct *IPObj;

    switch (LineType) {
	case 0: y = LineValue1;
	    z = LineValue2;	   
	    NumPt = NumPtXYZ[0];
	    Crv = BzrCrvNew(NumPt, CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (x = 0; x < NumPtXYZ[0]; x++) {
		Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) + y * NumPtXYZ[0] + x;
		PtIdVec[CrvInd++] = Ind;
	    }
	    break;
	case 1: x = LineValue1;
	    z = LineValue2;
	    NumPt = NumPtXYZ[1];
	    Crv = BzrCrvNew(NumPt, CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    Crv = BzrCrvNew(NumPt, CAGD_PT_E3_TYPE);
	    for (y = 0; y < NumPtXYZ[1]; y++) {
		Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) + y * NumPtXYZ[0] + x;
		PtIdVec[CrvInd++] = Ind;
	    }
	    break;
	case 2: x = LineValue1;
	    y = LineValue2;
	    NumPt = NumPtXYZ[2];
	    Crv = BzrCrvNew(NumPt, CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    Crv = BzrCrvNew(NumPt, CAGD_PT_E3_TYPE);
	    for (z = 0; z < NumPtXYZ[2]; z++) {
		Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) + y * NumPtXYZ[0] + x;
		PtIdVec[CrvInd++] = Ind;
	    }
	    break;
	default:
	    assert(0);
    }

    IPObj = IPGenCRVObject(Crv);    
    UserTopoAddCell(Ud, PtIdVec, NumPt, IPObj);
    IritFree(PtIdVec);
}


/* create a planar surface and add to the grid. */
static void AddSrfToGrid(UserTopoUnstrctGeomStruct *Ud,
		         int FirstPtId,
		         int NumPtXYZ[3],
		         int PlaneType,
		         int PlaneValue)
{
    int Ind, NumPt, *PtIdVec,
	x = 0, 
	y = 0, 
	z = 0,
	SrfInd = 0,
	xInc = 1,
	yInc = 1,
	zInc = 1;
    CagdSrfStruct *Srf;
    IPObjectStruct *IPObj;

    assert(PlaneValue >= 0 && PlaneValue < NumPtXYZ[PlaneType]);

    switch (PlaneType) {
	case 0:	x = PlaneValue;
	    xInc = 0;
	    NumPt = NumPtXYZ[1] * NumPtXYZ[2];
	    Srf = BzrSrfNew(NumPtXYZ[1], NumPtXYZ[2], CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (z = 0; z < NumPtXYZ[2]; z = z++) {
		for (y = 0; y < NumPtXYZ[1]; y = y++) {
		    Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) + y * NumPtXYZ[0] + x;
		    PtIdVec[SrfInd++] = Ind;
		}
	    }
	    break;
	case 1:	y = PlaneValue;
	    yInc = 0;
	    NumPt = NumPtXYZ[0] * NumPtXYZ[2];
	    Srf = BzrSrfNew(NumPtXYZ[0], NumPtXYZ[2], CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (z = 0; z < NumPtXYZ[2]; z = z++) {
		for (x = 0; x < NumPtXYZ[0]; x = x++) {
		    Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) + y * NumPtXYZ[0] + x;
		    PtIdVec[SrfInd++] = Ind;
		}
	    }
	    break;
	case 2: z = PlaneValue;
	    zInc = 0;
	    NumPt = NumPtXYZ[0] * NumPtXYZ[1];
	    Srf = BzrSrfNew(NumPtXYZ[0], NumPtXYZ[1], CAGD_PT_E3_TYPE);
	    PtIdVec = (int *) IritMalloc(sizeof(int) * NumPt);
	    for (y = 0; y < NumPtXYZ[1]; y = y++) {
		for (x = 0; x < NumPtXYZ[0]; x = x++) {
		    Ind = FirstPtId + z * (NumPtXYZ[0] * NumPtXYZ[1]) + y * NumPtXYZ[0] + x;
		    PtIdVec[SrfInd++] = Ind;
		}
	    }
	    break;
	default:
	    assert(0);
    }   
    
    IPObj = IPGenSRFObject(Srf);    
    UserTopoAddCell(Ud, PtIdVec, NumPt, IPObj);
    IritFree(PtIdVec);
}



/* Print all adjacency relations amongst entities of the grid. */
static void PrintAdjacencies(UserTopoUnstrctGeomStruct *Ud)
{
    int NumEntAdj, i, EntId, *AdjList;
    IPObjectStruct *Ent;

    for (Ent = Ud -> CellList; Ent != NULL; Ent = Ent -> Pnext) {
	EntId = UserTopoObjectToId(Ud, Ent);
	NumEntAdj = UserTopoEntitiesAdjacentToCell(Ud, EntId, &AdjList);
	printf("%d: ", EntId);
	for (i = 0; i < NumEntAdj; i++) {
	    printf("%d, ", AdjList[i]);
	}
/*	IritFree(AdjList);			dbgheap	*/
	printf("\n");
    }
    
}

/* Render all the entities of the grid onto the given handle. */
static void RenderGridEntities(const UserTopoUnstrctGeomStruct *Ud,
			       int ClearScreen,
			       int Handle)
{
    IPObjectStruct *Obj,
	*ClrObj = IPGenStrObject("command_", "clear", NULL);

    if (ClearScreen) {
	IPSocWriteOneObject(Handle, ClrObj);
    }

    assert(Ud);

    for (Obj = Ud -> CellList; Obj != NULL; Obj = Obj -> Pnext) {
	IPSocWriteOneObject(Handle, Obj);
    }

    IPFreeObject(ClrObj);
}

/* An example of a call-back function which keeps first 500 points and discards the rest. */
int *UserTopoFilterPtCBFunction(const UserTopoUnstrctGeomStruct *Ud)
{
    int i,
	*FilterPtIndx = (int *) IritMalloc(sizeof(int) * Ud -> NumPts);

    for (i = 0; i < Ud -> NumPts; i++) {
	if (i >= 0 && i < 500) {
	    FilterPtIndx[i] = 1;
	}
	else {
	    FilterPtIndx[i] = 0;
	}
    }
    return FilterPtIndx;
}


/* Test case for unstructured grid. */
static void TestCase1(int IOHandle) 
{
    CagdPType Origin;
    CagdRType DistXYZ[3];
    int NumPtXYZ[3], NumPt, NumEnt, *IdVec, *CloneMap,
	FirstPtId = 0;
    UserTopoUnstrctGeomPtStruct *PtVec, *PtVec2;
    UserTopoUnstrctGeomStruct *Ud1, *Ud2, *Ud12, *Ud12B2, *Ud12B1, *Ud12B12, *Ud12F;
    IPObjectStruct *Ent;

    Origin[0] = Origin[1] = Origin[2] = 0;
    DistXYZ[0] = DistXYZ[1] = DistXYZ[2] = 0.1;
    NumPtXYZ[0] = NumPtXYZ[1] = NumPtXYZ[2] = 10;
    NumPt = NumPtXYZ[0] * NumPtXYZ[1] * NumPtXYZ[2];

    PtVec = CreatePtVec(Origin, DistXYZ, NumPtXYZ, FirstPtId);
    Ud1 = UserTopoUnstrctGeomNew();
    UserTopoSetPoints(Ud1, PtVec, NumPt, FALSE, &CloneMap); 

    AddCrvToGrid(Ud1, FirstPtId, NumPtXYZ, 0, 1, 1);				    //Ent 0
    Ent = UserTopoIdToObject(Ud1, 0);
    UserTopoSetCellIntAttr(Ud1, Ent, "color", 1);

    AddCrvToGrid(Ud1, FirstPtId, NumPtXYZ, 1, NumPtXYZ[0] - 1, NumPtXYZ[2] - 1);    //Ent 1
    Ent = UserTopoIdToObject(Ud1, 1);
    UserTopoSetCellIntAttr(Ud1, Ent, "color", 1);

    AddCrvToGrid(Ud1, FirstPtId, NumPtXYZ, 2, NumPtXYZ[0] - 1, NumPtXYZ[1] - 1);    //Ent 2
    Ent = UserTopoIdToObject(Ud1, 2);
    UserTopoSetCellIntAttr(Ud1, Ent, "color", 3);

    AddSrfToGrid(Ud1, FirstPtId, NumPtXYZ, 2, 0);				    //Ent 3
    Ent = UserTopoIdToObject(Ud1, 3);
    UserTopoSetCellIntAttr(Ud1, Ent, "color", 10);

    UserTopoUnstrctGeomUpdate(&Ud1, IRIT_EPS, FALSE);
    printf("Ud1 adj\n");
    PrintAdjacencies(Ud1);  
    RenderGridEntities(Ud1, TRUE, IOHandle);

    Origin[2] = 0.9;
    NumPtXYZ[2] = 2;
    NumPt = NumPtXYZ[0] * NumPtXYZ[1] * NumPtXYZ[2];
    PtVec2 = CreatePtVec(Origin, DistXYZ, NumPtXYZ, FirstPtId);
    Ud2 = UserTopoUnstrctGeomNew();
    UserTopoSetPoints(Ud2, PtVec2, NumPt, FALSE, &CloneMap); 

    AddCrvToGrid(Ud2, FirstPtId, NumPtXYZ, 0, NumPtXYZ[1] - 1, 0);		//Ent 0, will become Ent 4    
    Ent = UserTopoIdToObject(Ud2, 0);
    UserTopoSetCellIntAttr(Ud2, Ent, "color", 4);
    
    UserTopoUnstrctGeomUpdate(&Ud2, IRIT_EPS, FALSE);
    printf("Ud2 adj\n");
    PrintAdjacencies(Ud2);  
    RenderGridEntities(Ud2, TRUE, IOHandle);
    
    Ud12 = UserTopoAppendUnstrctGeoms(Ud1, Ud2, IRIT_EPS, TRUE, &CloneMap);	//100 points common in Ud1 and Ud2
    printf("Ud12 adj\n");
    PrintAdjacencies(Ud12);
    RenderGridEntities(Ud12, TRUE, IOHandle);
    
    Ud12B1 = UserTopoCrvBndryFilter(Ud12);
    printf("Ud12B1 adj\n");
    PrintAdjacencies(Ud12B1);
    RenderGridEntities(Ud12B1, TRUE, IOHandle);

    Ud12B2 = UserTopoSrfBndryFilter(Ud12);
    printf("Ud12B2 adj\n");
    PrintAdjacencies(Ud12B2);
    RenderGridEntities(Ud12B2, TRUE, IOHandle);

    Ud12B12 = UserTopoAppendUnstrctGeoms(Ud12B2, Ud12B1, IRIT_EPS, TRUE, &CloneMap);
    printf("Ud12B12 adj\n");
    PrintAdjacencies(Ud12B12);
    RenderGridEntities(Ud12B12, TRUE, IOHandle);

    NumEnt = UserTopoGetCellAttrThreshold(Ud12, "color", 3, 4, &IdVec);
    IritFree(IdVec);

    UserTopoSetFilterGridCallBackFunc(UserTopoFilterPtCBFunction);
    Ud12F = UserTopoApplyFilterToGrid(Ud12, TRUE);
    printf("Ud12F adj\n");
    PrintAdjacencies(Ud12F);
    RenderGridEntities(Ud12F, TRUE, IOHandle);

    UserTopoUnstrctGeomFree(Ud1);
    UserTopoUnstrctGeomFree(Ud2);
    UserTopoUnstrctGeomFree(Ud12);
    UserTopoUnstrctGeomFree(Ud12B1);
    UserTopoUnstrctGeomFree(Ud12B2);
    UserTopoUnstrctGeomFree(Ud12B12);
    UserTopoUnstrctGeomFree(Ud12F);
}


void main(int argc, char **argv)
{
    int PrgmIO = 0;
    char *Program = getenv("IRIT_DISPLAY");
    
    IPSocSrvrInit();            /* Initialize the listen socket for clients. */

#ifdef __WINNT__
    if (Program == NULL)
	Program = "wntgdrvsD -s-";
#endif /* __WINNT__ */
#ifdef __UNIX__
    if (Program == NULL)
	Program = "x11drvs -s-";
#endif /* __UNIX__ */
     

    if ((PrgmIO = IPSocExecAndConnect(Program, getenv("IRIT_BIN_IPC") !=
	NULL)) >= 0)
    {
	CagdRType  
	    NumericTol = 1e-8,
	    SubdivTol = 1e-3;
	int i;

	TestCase1(PrgmIO);
	scanf("%d", &i);	
    }    
    
}



