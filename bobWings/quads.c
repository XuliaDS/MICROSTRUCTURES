

static void 
EG_setGrid ( int n[], double uv[][], double *grid ){
  int    i, j, k, i0, i1, i2, i3, ilast;
  int    ll, lr, ur, ul, j0, jm, ii, im, iv, sav;
  double corners[4 * 2], et, xi;
  if ( n[0] != n[2] || n[1] != n[3] ) {
      printf(" need same number of points in opposite directions\n");
      return;
  }
  corners[0] = uv[0][0];
  corners[1] = uv[0][1];
  corners[2] = uv[1][0];
  corners[3] = uv[1][1];
  corners[4] = uv[2][0];
  corners[5] = uv[2][1];
  corners[6] = uv[3][0];
  corners[7] = uv[3][1];
  printf(" %d -> %d -> %d  -> %d\n ", ll, lr, ur, ul);
  for (k = j = 0; j < n[1]; j++) {
      yj = ((double) (j+1)) / ((double) ny);
      grid[2 * k     ] = uv[3][2* ( n[3] - 1 - j )    ];
      grid[2 * k + 1 ] = uv[3][2* ( n[3] - 1 - j ) + 1];
      k++;
      for (i = 1; i < n[0] - 1; i++) {
	  if ( j == 0 ) {
	      grid[2 * k   ] = uv[0][2*i    ];
	      grid[2 * k + ] = uv[0][2*i + 1];
	      k++;
	  }
	  else if  ( j == n[1] - 1 ) {
	      grid[2 * k   ] = uv[2][2*i    ];
	      grid[2 * k + ] = uv[2][2*i + 1];
	      k++;
	  }
	  else {
	      xi          = ((double) (i+1)) / ((double) n[0]);
	      for (d = 0 ; d < 2; d++) {
		  grid[2 * k + d] = (1.0 - xi ) * uv[0][2 *         i + d  ] +
		      (      xi ) * uv[2][2* ( n[0] - i )   + d  ] +
		      (1.0 - yj ) * uv[1][2 *         j     + d  ] +
		      (      yj ) * uv[3][2 * (n[1] - j )   + d  ] -
		      (1.0 - xi ) * ( 1.0 - et ) * corners[0][d] -
		      (1.0 - xi ) * (       et ) * corners[3][d] -
		      (      xi ) * ( 1.0 - et ) * corners[1][d] -
		      (      xi ) * (        et) * corners[2][d];
	      }
	      k++;
	  }
      }
  }
}

