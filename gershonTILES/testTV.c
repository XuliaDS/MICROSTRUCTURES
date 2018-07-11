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
	    if (IP_IS_TRIVAR_OBJ(PObj)) {
		int IDir;
		CagdRType t;
		TrivTVDirType Dir;
		MvarMVStruct *MV;

		Dir = TRIV_NO_DIR;
		t = IRIT_INFNTY;
		if (TrivTVKnotHasC0Discont(PObj -> U.Trivars, &Dir, &t)) {
		    fprintf(stderr, "TV C0 Discont:  %d  %f: %smesh discont.\n",
			    Dir - TRIV_CONST_U_DIR, t,
			    TrivTVIsMeshC0DiscontAt(PObj -> U.Trivars, Dir, t)
								? "" : "no ");

		}
		else {
		    Dir = TRIV_NO_DIR;
		    t = IRIT_INFNTY;
		    if (TrivTVKnotHasC1Discont(PObj -> U.Trivars, &Dir, &t)) {
			fprintf(stderr, "TV C1 Discont:  %d  %f: %smesh discont.\n",
				Dir - TRIV_CONST_U_DIR, t,
				TrivTVIsMeshC1DiscontAt(PObj -> U.Trivars, Dir, t)
								? "" : "no ");
		    }
		}

		MV = MvarTVToMV(PObj -> U.Trivars);
		IDir = -1;
		t = IRIT_INFNTY;
		if (MvarMVKnotHasC0Discont(MV, &IDir, &t)) {
		    fprintf(stderr, "MV C0 Discont:  %d  %f: %smesh discont.\n",
			    IDir, t,
			    MvarMVIsMeshC0DiscontAt(MV, IDir, t)
								? "" : "no ");

		}
		else {
		    IDir = -1;
		    t = IRIT_INFNTY;
		    if (MvarMVKnotHasC1Discont(MV, &IDir, &t)) {
			fprintf(stderr, "MV C1 Discont:  %d  %f: %smesh discont.\n",
				IDir, t,
				MvarMVIsMeshC1DiscontAt(MV, IDir, t)
								? "" : "no ");
		    }
		}

	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: test tv.itd\n");
	exit(2);
    }

    exit(0);
}

void main4(int argc, char **argv)
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
	    if (IP_IS_CRV_OBJ(PObj)) {
	        int n;
		CagdRType R,
		    Tm = IritCPUTime(TRUE);
		CagdVType Pt;
		CagdCrvStruct
		    *Crv = PObj -> U.Crvs;
		void
		    *PrepPtCrv = SymbDistCrvPointPrep(Crv);

		Pt[2] = 0.0;
		for (n = 0; n < 100000; n++) {
		    CagdRType Pt[3];

		    Pt[0] = IritRandom(-2, 2);
		    Pt[1] = IritRandom(-2, 2);
		    R = SymbDistCrvPoint(Crv, PrepPtCrv, Pt,
					 TRUE, IRIT_EPS);
		}

		SymbDistCrvPointFree(PrepPtCrv);
		fprintf(stderr, "Time = %lf\n", IritCPUTime(FALSE));
	    }
	}
	else {
	    fprintf(stderr, "Failed to open file \"%s\"\n", argv[1]);
	    exit(1);
	}
    }
    else {
	fprintf(stderr, "Usage: test crv.itd\n");
	exit(2);
    }

    exit(0);
}

void main3(int argc, char **argv)
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
	    if (IP_IS_SRF_OBJ(PObj)) {
	        int n;
		CagdRType UMin, UMax, VMin, VMax, u, v, k1, k2;
		CagdVType D1, D2;
		CagdSrfStruct
		    *Srf = PObj -> U.Srfs;
		
		CagdSrfDomain(Srf, &UMin, &UMax, &VMin, &VMax);
		SymbEvalSrfCurvPrep(Srf, TRUE);
		for (n = 0; n < 1000000; n++) {
		    u = UMin + IritRandom(0, 1) * (UMax - UMin);
		    v = VMin + IritRandom(0, 1) * (VMax - VMin);

		    SymbEvalSrfCurvature(Srf, u, v, TRUE,
					 &k1, &k2, D1, D2);
		}
		SymbEvalSrfCurvPrep(Srf, FALSE);
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
