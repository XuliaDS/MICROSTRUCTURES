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

void main(int argc, char **argv)
{
    int Handler;

    /* START */
    int i;
    char *Name = "c:/temp/trim_srf.itd";
    IPObjectStruct 
	*trim_srf = IPGetDataFiles(&Name, 1, 1, 1);


    {
	IrtRType ri, rj, rk;
	ri = 1; rj = 30; rk = 0; IritDynMemoryDbgCheckMark(&ri, &rj, &rk);
    }

    for (i = 0; i < 100; i++) {
	CagdSrfStruct
	    *untrim_srfs_phys = TrimSrfCnvrt2BzrRglrSrf2(trim_srf -> U.TrimSrfs,
							 TRUE,TRUE,1.0e-15);
	CagdSrfFreeList(untrim_srfs_phys);
    }


    {
	IrtRType ri, rj, rk;
	ri = 0; rj = 30; rk = 0; IritDynMemoryDbgCheckMark(&ri, &rj, &rk);
    }
    IPFreeObject(trim_srf);

    exit(0);

    /* END */
}
