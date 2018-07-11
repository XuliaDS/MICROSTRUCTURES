

#include <assert.h>

#include "egads.h"

#define  EPS06   1.0e-6
#define  HUGEQ   99999999.0

/* IRIT includes */
#include "inc_irit/irit_sm.h"
#include "inc_irit/iritprsr.h"
#include "inc_irit/allocate.h"
#include "inc_irit/attribut.h"
#include "inc_irit/cagd_lib.h"
#include "inc_irit/geom_lib.h"
#include "inc_irit/user_lib.h"

static int makeIRITsrf(ego eobj, CagdSrfStruct **surf);

int
EGADS_TO_IRIT ( ego  context, ego emodel )
{
  int     status = EGADS_SUCCESS;
  int     oclass, mtype, nchild, *senses, nface, iface, nloop, nedge, iedge;
  int     periodic, i;
  double  data[18], trange[2], tt, umin, umax, vmin, vmax;
  char    *ErrStr, fn[50], fnitd[50];
  ego      eref, esurf,  *ebodys,*efaces=NULL, *eloops, *eedges, epcurve;
  CagdSrfStruct **IRITsurf = NULL;
  TrivTVStruct   *TVMap;
  MvarMVStruct   *DeformMV;
  IRITsurf = malloc ( 6 * sizeof(CagdSrfStruct));
  if (IRITsurf == NULL ) return EGADS_MALLOC;

  /* check that Model was input that contains one Body */
  status = EG_getTopology(emodel, &eref, &oclass, &mtype,
			  data, &nchild, &ebodys, &senses);
  if ( status != EGADS_SUCCESS) {
      printf(" EG_getTopology = %d\n", status);
      return status;
  }
  if (oclass != MODEL) {
      printf(" udpExecute: expecting a Model\n");
      return  EGADS_NOTMODEL;
  } else if (nchild != 1) {
      printf(" udpExecute: expecting Model to contain one Bodys (not %d)\n", nchild);
      return EGADS_NOTBODY;
  }
  status = EG_getContext(emodel, &context);
  if ( status != EGADS_SUCCESS) {
      printf(" EG_getContext = %d\n", status);
      return status;
  }

  /* make sure that the first Body (the duct) contains 6 Faces */
  status = EG_getBodyTopos(ebodys[0], NULL, FACE, &nface, &efaces);
  if ( status != EGADS_SUCCESS) {
      printf(" EG_getBodyTopos = %d\n", status);
      return status;
  }

  if (nface != 6) {
      EG_free(efaces);
      printf(" DUCT  does not contain 6 Faces\n");
      status = EGADS_TOPOERR;
      return status;
  }

  /* make sure that the first 4 Faces all have exactly 4 Edges and
       are not trimmed (except along isoU or isoV lines) */
  printf("Checking input ...\n");
  // face 1 & 6 are top and bottom
  // face 2 & 4 are left right
  // face 3 & 5  are front back
  for (iface = 0; iface < 6; iface++) {
      status = EG_getTopology(efaces[iface], &esurf, &oclass, &mtype,
			      data, &nloop, &eloops, &senses);
      if ( status != EGADS_SUCCESS) {
	  printf(" EG_getTopology FACE %d = %d\n", iface, status);
	  return status;
      }
      if (nloop != 1) {
	  printf(" udpExecute: Face %d has more than one Loop\n", iface+1);
	  return EGADS_TOPOERR;
      }

      status = EG_getTopology(eloops[0], &eref, &oclass, &mtype,
			      data, &nedge, &eedges, &senses);
      if ( status != EGADS_SUCCESS) {
	  printf(" EG_getTopology LOOP = %d\n", status);
      }


      if (nedge != 4) {
	  printf(" udpExecute: Face %d is not bounded by 4 Edges\n", iface+1);
	  return EGADS_TOPOERR;
      }

      for (iedge = 0; iedge < 4; iedge++) {
	  umin = +HUGEQ;
	  umax = -HUGEQ;
	  vmin = +HUGEQ;
	  vmax = -HUGEQ;

	  status = EG_getRange(eedges[iedge], trange, &periodic);
	  if ( status != EGADS_SUCCESS) {
	      printf(" EG_getRange = %d\n", status);
	  }
	  status = EG_otherCurve(esurf, eedges[iedge], 0, &epcurve);
	  if ( status != EGADS_SUCCESS) {
	      printf(" EG_OtherCurve = %d\n", status);
	  }
	  for (i = 0; i < 51; i++) {
	      tt = trange[0] + (trange[1] - trange[0]) * i / 50.0;

	      status = EG_evaluate(epcurve, &tt, data);
	      if ( status != EGADS_SUCCESS) {
		  printf(" EG_evaluate = %d\n", status);
		  return status;
	      }
	      if (data[0] < umin) umin = data[0];
	      if (data[0] > umax) umax = data[0];
	      if (data[1] < vmin) vmin = data[1];
	      if (data[1] > vmax) vmax = data[1];
	  }
	  if (fabs(umax-umin) > EPS06 && fabs(vmax-vmin) > EPS06) {
	      printf(" Face %d has Edge %d that is not an isocline\n", iface+1, iedge+1);
	      printf( " umax %f umin %f    vmax %f  vmin %f\n ", umax, umin, vmax, vmin ) ;
	      //ereturn status;
	  }
	  else {
	      printf(" Face  %d EDGE %d  ISOCLINE\n", iface + 1, iedge + 1);
	  }
      }
  }

  /* set up the IRIT trivatiate for the duct (use Surfaces since in
     general not planar) */
  printf("Setting up IRIT duct...\n");


  for ( iface = 0; iface < 6; iface++ ) {
      printf(" Convert face %d to BSpline face \n", iface + 1);
      status = EG_convertToBSpline(efaces[iface], &esurf);
      if ( status != EGADS_SUCCESS) {
	  printf(" EG_covertToBSPline = %d\n", status);
	  return status;
      }

      status = makeIRITsrf(esurf, &IRITsurf[iface]);
      if ( status != EGADS_SUCCESS) {
	  printf(" EG_makeIRITsrf = %d\n", status);
	  return status;
      }
      sprintf ( fn   , "surf%d",iface+1);
      sprintf ( fnitd, "surf%d.itd",iface+1);
      BspSrfWriteToFile(IRITsurf[iface], fnitd, 0, fn, &ErrStr);
      if (ErrStr != NULL) printf(" udpExecute: surf1 ErrStr: %s\n", ErrStr);
  }

  TVMap    = MvarTrivarBoolSum3(IRITsurf[0],  IRITsurf[1],  IRITsurf[2],  IRITsurf[3],  IRITsurf[4], IRITsurf[5] ) ;
  DeformMV = MvarCnvrtTVToMV(TVMap);
  status = TrivBspTVWriteToFile(TVMap, "duct.itd", 0, "duct", &ErrStr);
  if (ErrStr != NULL) printf(" udpExecute: duct ErrStr: %s  (status=%d)\n", ErrStr, status);
  TrivTVFree(TVMap);
  for ( i = 0 ; i < 6; i++ )
    CagdSrfFree(IRITsurf[i]);
  free ( IRITsurf) ;
  EG_free(efaces);


  return status;

}


/*
 ************************************************************************
 *                                                                      *
 *   makeIRITsrf -- convert EGADS BSPINE into IRIT Surface              *
 *                                                                      *
 ************************************************************************
 */

static int
makeIRITsrf(ego           eobj,         /* (in)  Egads surface of Face */
	    CagdSrfStruct **Srf)        /* (in)  IRIT  surface (or NULL) */
{
  int     status = EGADS_SUCCESS;

  int     i, j, k, indx, oclass, mtype, *ivec=NULL, nloop, nedge, nnode, nchild, *senses, *senses2;
  double  *rvec=NULL, data[4], xyzsw[4], xyzse[4], xyznw[4], xyzne[4];
  ego     topRef, eprev, enext, eref, esurf, *eloops, *eedges, *enodes, *echilds;



  /* --------------------------------------------------------------- */

  /* default return */
  *Srf = NULL;

  /* get type of eobj */
  status = EG_getInfo(eobj, &oclass, &mtype, &topRef, &eprev, &enext);
  if ( status != EGADS_SUCCESS ) {
      printf( " EG_getInfo stat = %d\n", status);
      return status;
  }

  /* eobj is a Surface */
  printf(" OBJETC TYPE  %d (12 = surface 23 = face )  TYPE  %d  ( 8 = BSpline  1 = LINE  6 = TRIMMED ) \n", oclass, mtype ) ;
  if (oclass == SURFACE) {
      esurf = eobj;

      /* make sure esurf is a BSpline (and not a NURB and not periodic) */
      status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
      if ( status != EGADS_SUCCESS ) {
	  printf( " EG_getGeometry stat = %d\n", status);
	  return status;
      }

      if (mtype != BSPLINE) {
	  printf(" makeIRITsrf: esurf NOT a BSpline = %d\n", mtype);
	  status = EGADS_GEOMERR;
	  goto cleanup;
      } else if (ivec[0] != 0) {
	  printf(" makeIRITsrf: BSpline flags = %d\n", ivec[0]);
	  status = EGADS_GEOMERR;
	  goto cleanup;
      }
      /* start a new IRIT BSpline surface */
      printf(" Call IRIT Function  BspsSrfNew \n");
      *Srf = BspSrfNew(ivec[2], ivec[5], ivec[1]+1, ivec[4]+1, CAGD_PT_E3_TYPE);

      /* Populate the control points in Srf */
      k = ivec[3] + ivec[6];
      for (j = 0; j < ivec[5]; j++) {
	  for (i = 0; i < ivec[2]; i++) {
	      indx = CAGD_MESH_UV(*Srf, i, j);
	      (*Srf)->Points[1][indx] = rvec[k++];   /* X coefs */
	      (*Srf)->Points[2][indx] = rvec[k++];   /* Y coefs */
	      (*Srf)->Points[3][indx] = rvec[k++];   /* Z coefs */
	  }
      }

      /* Copy the knot vectors into Srf */
      for (i = 0; i < ivec[3]; i++) {
	  (*Srf)->UKnotVector[i] = rvec[i];
      }
      for (i = 0; i < ivec[6]; i++) {
	  (*Srf)->VKnotVector[i] = rvec[i+ivec[3]];
      }

      /* eobj is a Face */
  } else if (oclass == FACE) {
      printf( "\n\n ENTRA EN FACE\n");
      status = EG_getTopology(eobj, &eref, &oclass, &mtype, data,
			      &nloop, &eloops, &senses);
      if ( status != EGADS_SUCCESS ) {
	  printf( " EG_getTopology stat = %d\n", status);
	  return status;
      }



      status = EG_getTopology(eloops[0], &eref, &oclass, &mtype, data,
			      &nedge, &eedges, &senses);
      if ( status != EGADS_SUCCESS ) {
	  printf( " EG_getTopology stat = %d\n", status);
	  return status;
      }

      /* get the 4 corners as the beginning/end of each Edge */
      if (nedge == 4) {
	  status = EG_getTopology(eedges[0], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[0] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyzsw,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyzsw,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  status = EG_getTopology(eedges[1], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[1] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyzse,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyzse,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  status = EG_getTopology(eedges[2], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[2] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyzne,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyzne,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  status = EG_getTopology(eedges[3], &eref, &oclass, &mtype, data,
				  &nnode, &enodes, &senses2);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }
	  if (senses[3] == SFORWARD) {
	      status = EG_getTopology(enodes[0], &eref, &oclass, &mtype, xyznw,
				      &nchild, &echilds, &senses2);
	  } else {
	      status = EG_getTopology(enodes[1], &eref, &oclass, &mtype, xyznw,
				      &nchild, &echilds, &senses2);
	  }
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getTopology stat = %d\n", status);
	      return status;
	  }

	  /* make an untrimmed Bspline Surface */
	  *Srf = BspSrfNew(2, 2, 2, 2, CAGD_PT_E3_TYPE);

	  /* knot vectors */
	  (*Srf)->UKnotVector[0] = 0;
	  (*Srf)->UKnotVector[1] = 0;
	  (*Srf)->UKnotVector[2] = 1;
	  (*Srf)->UKnotVector[3] = 1;

	  (*Srf)->VKnotVector[0] = 0;
	  (*Srf)->VKnotVector[1] = 0;
	  (*Srf)->VKnotVector[2] = 1;
	  (*Srf)->VKnotVector[3] = 1;

	  /* control points */
	  indx = CAGD_MESH_UV(*Srf, 0, 0);
	  (*Srf)->Points[1][indx] = xyzsw[0];
	  (*Srf)->Points[2][indx] = xyzsw[1];
	  (*Srf)->Points[3][indx] = xyzsw[2];

	  indx = CAGD_MESH_UV(*Srf, 1, 0);
	  (*Srf)->Points[1][indx] = xyzse[0];
	  (*Srf)->Points[2][indx] = xyzse[1];
	  (*Srf)->Points[3][indx] = xyzse[2];

	  indx = CAGD_MESH_UV(*Srf, 0, 1);
	  (*Srf)->Points[1][indx] = xyznw[0];
	  (*Srf)->Points[2][indx] = xyznw[1];
	  (*Srf)->Points[3][indx] = xyznw[2];

	  indx = CAGD_MESH_UV(*Srf, 1, 1);
	  (*Srf)->Points[1][indx] = xyzne[0];
	  (*Srf)->Points[2][indx] = xyzne[1];
	  (*Srf)->Points[3][indx] = xyzne[2];

	  /* if Face does not have 4 Nodes but is a planar
           Face (hopefully with only a notch), so that we can
           convert Face to BSpline and then use it directly */
      } else {
	  status = EG_convertToBSpline(eobj, &esurf);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_convertToSPline stat = %d\n", status);
	      return status;
	  }

	  status = EG_getGeometry(esurf, &oclass, &mtype, &eref, &ivec, &rvec);
	  if ( status != EGADS_SUCCESS ) {
	      printf( " EG_getGeometry stat = %d\n", status);
	      return status;
	  }

	  /* start a new IRIT BSpline surface */
	  *Srf = BspSrfNew(ivec[2], ivec[5], ivec[1]+1, ivec[4]+1, CAGD_PT_E3_TYPE);

	  /* Populate the control points in Srf */
	  k = ivec[3] + ivec[6];
	  for (j = 0; j < ivec[5]; j++) {
	      for (i = 0; i < ivec[2]; i++) {
		  indx = CAGD_MESH_UV(*Srf, i, j);
		  (*Srf)->Points[1][indx] = rvec[k++];   /* X coefs */
		  (*Srf)->Points[2][indx] = rvec[k++];   /* Y coefs */
		  (*Srf)->Points[3][indx] = rvec[k++];   /* Z coefs */
	      }
	  }

	  /* Copy the knot vectors into Srf */
	  for (i = 0; i < ivec[3]; i++) {
	      (*Srf)->UKnotVector[i] = rvec[i];
	  }
	  for (i = 0; i < ivec[6]; i++) {
	      (*Srf)->VKnotVector[i] = rvec[i+ivec[3]];
	  }
      }

      /* eobj is an unknown type */
  } else {
      printf(" makeIRITsrf: eobj is neither Surface nor Face (oclass=%d)\n", oclass);
      status = EGADS_GEOMERR;
      goto cleanup;
  }

  cleanup:
  if (ivec != NULL) EG_free(ivec);
  if (rvec != NULL) EG_free(rvec);

  return status;
}


int main (int argc, char *argv[]){
  int  status;
  ego  context, model;
  if ( argc < 2 ) {
      printf(" Need EGADS file to convert.....\n\n");
      return EGADS_EMPTY;
  }
  /* dummy call to prevent compiler warnings */

  /* define a context */
  status = EG_open(&context);
  printf("EG_open -> status=%d\n", status);
  if (status < 0) exit(EXIT_FAILURE);
  if(argv[1])  status = EG_loadModel(context, 0, argv[1], &model);
  if ( status != EGADS_SUCCESS ) {
      printf(" EG_loadModel = %d !!\n", status ) ;
      EG_free(model);
      EG_free(context);
      return status;
  }
  status = EGADS_TO_IRIT(context, model);
  printf("udpExecute -> status=%d\n", status);
  EG_free(model);
  EG_free(context);
  return status;


}


















