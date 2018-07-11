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
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

void main(int argc, char **argv)
{
    int Handler;

    if (argc == 2) {
        IPSetFlattenObjects(FALSE);
	if ((Handler = IPOpenDataFile(argv[1], TRUE, TRUE)) >= 0) {
	    IPObjectStruct *PTV, *PMdl, *PMdls,
		*PObj = IPGetObjects(Handler);
	    MdlModelStruct *AllMdls;

	    /* Done with file - close it. */
	    IPCloseStream(Handler, TRUE);

	    /* PObj should hold a list object with a TV and a Srf: */
	    PTV = IPListObjectGet(PObj, 0);
	    PMdl = IPListObjectGet(PObj, 1);

	    if (!IP_IS_TRIVAR_OBJ(PTV) || !IP_IS_MODEL_OBJ(PMdl)) {
	        fprintf(stderr, "Expects a trivar and a model\n");
		exit(1);
	    }

	    MdlBooleanSetTolerances(0.025, 1e-10, 0.01);

	    /* The key function - divide TSrfs at all knot planes of TV: */
	    AllMdls = UserDivideMdlsAtAllMVInteriorKnot(
		PMdl -> U.Mdls, MvarTVToMV(PTV -> U.Trivars)); /* Mem leak. */

	    PMdls = IPGenMODELObject(AllMdls);
	    IPStdoutObject(PMdls, FALSE);
	    IPFreeObject(PMdls);
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
