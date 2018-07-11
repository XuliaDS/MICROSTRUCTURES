/*****************************************************************************
* Filter the data.							     *
******************************************************************************
* (C) Gershon Elber, Technion, Israel Institute of Technology                *
******************************************************************************
* Written by:  Gershon Elber				Ver 1.0, June 1995   *
*****************************************************************************/

#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"

static void AddAttributes(IPObjectStruct *PObj, IrtHmgnMatType Mat)
{
    switch (PObj -> ObjType) {
	case IP_OBJ_CURVE:
	    AttrSetIntAttrib(&PObj -> U.Crvs -> Attr, "IntAttr", 10000);
	    break;
	case IP_OBJ_SURFACE:
	    AttrSetIntAttrib(&PObj -> U.Srfs -> Attr, "IntAttr", 10001);
	    break;
	case IP_OBJ_TRISRF:
	    AttrSetIntAttrib(&PObj -> U.TriSrfs -> Attr, "IntAttr", 10002);
	    break;
	case IP_OBJ_TRIMSRF:
	    AttrSetIntAttrib(&PObj -> U.TrimSrfs -> Attr, "IntAttr", 10010);
	    AttrSetIntAttrib(&PObj -> U.TrimSrfs -> Srf -> Attr, "IntAttr",
			     10011);
	    AttrSetIntAttrib(&PObj -> U.TrimSrfs -> TrimCrvList -> Attr,
			     "IntAttr", 10012);
	    AttrSetIntAttrib(&PObj -> U.TrimSrfs -> TrimCrvList ->
							TrimCrvSegList -> Attr,
			     "IntAttr", 10013); 
	    AttrSetIntAttrib(&PObj -> U.TrimSrfs -> TrimCrvList ->
					       TrimCrvSegList -> UVCrv -> Attr,
			     "IntAttr", 10014);
	    break;
	case IP_OBJ_TRIVAR:
	    AttrSetIntAttrib(&PObj -> U.Srfs -> Attr, "IntAttr", 10020);
	    break;
	case IP_OBJ_MODEL:
	    AttrSetIntAttrib(&PObj -> U.Mdls -> Attr, "IntAttr", 10030);
	    AttrSetIntAttrib(&PObj -> U.Mdls -> TrimSegList -> Attr,
			     "IntAttr", 10031);
	    AttrSetIntAttrib(&PObj -> U.Mdls -> TrimSegList ->
							    UVCrvFirst -> Attr,
			     "IntAttr", 10032);
	    AttrSetIntAttrib(&PObj -> U.Mdls -> TrimSrfList -> Attr,
			     "IntAttr", 10033);
	    AttrSetIntAttrib(&PObj -> U.Mdls -> TrimSrfList -> Srf -> Attr,
			     "IntAttr", 10034);
	    AttrSetIntAttrib(&PObj -> U.Mdls -> TrimSrfList -> LoopList -> Attr,
			     "IntAttr", 10035);
	    break;
	case IP_OBJ_MULTIVAR:
	    AttrSetIntAttrib(&PObj -> U.MultiVars -> Attr, "IntAttr", 10040);
	    break;
    }
}

void main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IPSetFlattenObjects(FALSE);
	if ((Handler = IPOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    const char
		*TmpName = "/temp/AttrObjs.itd";
	    IPObjectStruct *PObj2, *PObj;
	    IrtHmgnMatType Mat;

	    PObj = IPGetObjects(Handler);

	    /* Done with file - close it. */
	    IPCloseStream(Handler, TRUE);

	    IPTraverseObjListHierarchy(PObj, Mat, AddAttributes);

	    /* Save data with attributes. */
	    IPPutObjectToFile3(TmpName, PObj, 0);

	    /* And reload. */
	    PObj2 = IPGetDataFiles(&TmpName, 1, TRUE, TRUE);
	    IPStdoutObject(PObj2, FALSE);

	    /* Now read the data. */
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: IOAttr geom.itd\n");
	exit(2);
    }

    exit(0);
}
