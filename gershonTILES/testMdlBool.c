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
#include "inc_irit/cagd_lib.h"
#include "inc_irit/mdl_lib.h"

int main(int argc, char **argv)
{
    CagdVType
        SphereCenter = { 0.0, 0.0, 0.0 },
        ConeCenter = { 0.2, 0.0, -2.0 };
    CagdSrfStruct
        *Srf1 = CagdPrimSphereSrf(SphereCenter, 1.0, FALSE),
        *Srf2 = CagdPrimConeSrf(ConeCenter, 1.0, 4, FALSE, TRUE);
    MdlModelStruct
        *Mdl1 = MdlCnvrtSrf2Mdl(Srf1),
        *Mdl2 = MdlCnvrtSrf2Mdl(Srf2);
    IPObjectStruct *Inter, *Union, *Subtr;

    CagdSrfFree(Srf1);
    CagdSrfFree(Srf2);

    Inter = MdlBooleanIntersection(Mdl1, Mdl2, NULL, NULL);
    Subtr = MdlBooleanSubtraction(Mdl1, Mdl2, NULL, NULL);
    Union = MdlBooleanUnion(Mdl1, Mdl2, NULL, NULL);

    MdlModelFree(Mdl1);
    MdlModelFree(Mdl2);

    IPPutObjectToFile3("InterMdlBool.itd", Inter, 0);
    IPPutObjectToFile3("SubtrMdlBool.itd", Subtr, 0);
    IPPutObjectToFile3("UnionMdlBool.itd", Union, 0);

    IPFreeObject(Inter);
    IPFreeObject(Subtr);
    IPFreeObject(Union);
}
