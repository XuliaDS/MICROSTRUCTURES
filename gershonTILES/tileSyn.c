#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/user_lib.h"

int main(int argc, char **argv)
{
    int i;
    UserMicroTileBndryPrmStruct
	UMinPrms, UMaxPrms, VMinPrms, VMaxPrms, WMinPrms, WMaxPrms;
    IPObjectStruct *Tile;
       
    IRIT_ZAP_MEM(&UMinPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&UMaxPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMinPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMaxPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMinPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMaxPrms, sizeof(UserMicroTileBndryPrmStruct));
    if ( argc > 1 ) { 
    UMinPrms.Circular = TRUE;
    VMinPrms.Circular = TRUE;
    WMaxPrms.Circular = TRUE;
    UMinPrms.OuterRadius = 0.12;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.12;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.12;
    Tile = UserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms);
    IPPutObjectToFile3("TileSolid1.itd", Tile, 0);
    IPFreeObject(Tile);
return 0;
   
     } 
    UMinPrms.Circular = TRUE;
    VMinPrms.Circular = TRUE;
    WMaxPrms.Circular = TRUE;
    WMaxPrms.Bndry = 0.1;
    WMinPrms.Bndry = 0.2;
    UMinPrms.OuterRadius = 0.12;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.12;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.12;
    Tile = UserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms);
    IPPutObjectToFile3("TileSolid1.itd", Tile, 0);
    IPFreeObject(Tile);

    WMaxPrms.Bndry = 0.0;
    WMinPrms.Bndry = 0.0;
    VMaxPrms.Bndry = 0.05;
    VMinPrms.Bndry = 0.14;
    UMinPrms.OuterRadius = 0.13;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.11;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.1;
    Tile = UserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms);
    IPPutObjectToFile3("TileSolid1diff.itd", Tile, 0);
    IPFreeObject(Tile);

    VMaxPrms.Bndry = 0.0;
    VMinPrms.Bndry = 0.0;
    UMaxPrms.Bndry = 0.15;
    UMinPrms.Bndry = 0.24;
    UMinPrms.OuterRadius = 0.12;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.12;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.12;
    UMinPrms.InnerRadius = 0.08;
    UMaxPrms.InnerRadius = 0.08;
    VMinPrms.InnerRadius = 0.08;
    VMaxPrms.InnerRadius = 0.08;
    WMinPrms.InnerRadius = 0.08;
    WMaxPrms.InnerRadius = 0.08;
    Tile = UserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms);
    IPPutObjectToFile3("TileHollowed1.itd", Tile, 0);
    IPFreeObject(Tile);

    UMaxPrms.Bndry = 0.0;
    UMinPrms.Bndry = 0.0;
    WMaxPrms.Bndry = 0.1;
    WMinPrms.Bndry = 0.12;
    UMinPrms.OuterRadius = 0.13;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.11;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.1;
    Tile = UserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				&VMinPrms, &VMaxPrms,
				&WMinPrms, &WMaxPrms);
    IPPutObjectToFile3("TileHollowed1diff.itd", Tile, 0);
    IPFreeObject(Tile);

    IRIT_ZAP_MEM(&UMinPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&UMaxPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMinPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&VMaxPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMinPrms, sizeof(UserMicroTileBndryPrmStruct));
    IRIT_ZAP_MEM(&WMaxPrms, sizeof(UserMicroTileBndryPrmStruct));
    UMinPrms.OuterRadius = 0.13;
    UMaxPrms.OuterRadius = 0.12;
    VMinPrms.OuterRadius = 0.11;
    VMaxPrms.OuterRadius = 0.12;
    WMinPrms.OuterRadius = 0.12;
    WMaxPrms.OuterRadius = 0.1;

#ifndef DEBUG
    fprintf(stderr, "starting 50000 solids...");
    for (i = 0; i < 50000; i++) {
	Tile = UserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				    &VMinPrms, &VMaxPrms,
				    &WMinPrms, &WMaxPrms);
        IPFreeObject(Tile);
    }
    fprintf(stderr, "done\n.");

    UMinPrms.InnerRadius = 0.1;
    UMaxPrms.InnerRadius = 0.11;
    VMinPrms.InnerRadius = 0.1;
    VMaxPrms.InnerRadius = 0.085;
    WMinPrms.InnerRadius = 0.09;
    WMaxPrms.InnerRadius = 0.08;
    UMinPrms.Bndry = 0.24;

    fprintf(stderr, "starting 10000 hollowed (and one bndry)...");
    for (i = 0; i < 10000; i++) {
	Tile = UserMicro3DCrossTile(&UMinPrms, &UMaxPrms,
				    &VMinPrms, &VMaxPrms,
				    &WMinPrms, &WMaxPrms);
        IPFreeObject(Tile);
    }
    fprintf(stderr, "done\n.");
#endif /* DEBUG */

    exit(0);
}
