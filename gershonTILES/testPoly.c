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
#include "inc_irit/grap_lib.h"

void main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IPSetFlattenObjects(TRUE);
	if ((Handler = IPOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IPObjectStruct
		*PObj = IPGetObjects(Handler);

	    /* Done with file - close it. */
	    IPCloseStream(Handler, TRUE);

	    /* Process the geometry. */
	    if (IP_IS_POLY_OBJ(PObj)) {
	        IPPolygonStruct
		    *Pl = GMMergeClosedLoopHoles(PObj -> U.Pl,
						 PObj -> U.Pl -> Pnext);

		IPStderrObject(IPGenPOLYObject(Pl));
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: PolyArea geom.itd\n");
	exit(2);
    }

    exit(0);
}

void main2(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IPSetFlattenObjects(TRUE);
	if ((Handler = IPOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IPObjectStruct *PTmp,
		*PObj = IPGetObjects(Handler);

	    /* Done with file - close it. */
	    IPCloseStream(Handler, TRUE);

	    /* Filter the geometry. */
	    for (PTmp = PObj; PTmp != NULL; PTmp = PTmp -> Pnext) {
	        if (IP_IS_TRIMSRF_OBJ(PTmp)) {
		    IPObjectStruct
		        *E3Crv = AttrGetObjectObjAttrib(PTmp, "E3Crv");

		    //if (E3Crv != NULL)
			IPStdoutObject(PTmp, FALSE);
		}
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: PolyArea geom.itd\n");
	exit(2);
    }

    exit(0);
}

void main1(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
	IPSetFlattenObjects(TRUE);
	if ((Handler = IPOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IPObjectStruct
		*PObj = IPGetObjects(Handler);

	    /* Done with file - close it. */
	    IPCloseStream(Handler, TRUE);

	    if (IP_IS_CRV_OBJ(PObj)) {
		CagdCrvStruct
		    *MCrvs = TrimMergeTrimmingCurves2Loops2(PObj -> U.Crvs,
							    IRIT_EPS);
		CagdDbg(MCrvs);
	    }
	}
    }
}
