#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "buildNebs.h"

int nbLess( const void *nb1, const void *nb2 )
{
        //compare function for qsort of objects type NbData
        NbData *pnb1 = (NbData*) nb1;
        NbData *pnb2 = (NbData*) nb2;
        if ( pnb1->distance < pnb2->distance ) return -1;
        return 1;
}

int computeNeighbourList( int n_pts, float *rij, NbData **neighbours, float trySigma2 )
{
	double  weight, radius, sumR2, sumWeight; // sum of neighbour distances
	int     count;       // max number of neighbours
	int     i;           // a loop variable
	int     top_neb;
	double  min_contrib;
	NbData *nebs;
	
	//If there are more points then allow a smaller contribution before cutoff
    min_contrib = 1e-4 / n_pts;

    //make a list of points which can be sorted.
    nebs = (NbData *)malloc(n_pts * sizeof(NbData));
    for( i = 0 ; i < n_pts; i++ ){
       nebs[i].id       = i;
       nebs[i].distance = rij[i];
    }

	// Step 2:
	//   Sort in-place neighbours according to their distance in increasing order.
	qsort( nebs, n_pts, sizeof( NbData ), nbLess );

	// keep going out until new neighbours make negligible contribution to the weight
	sumWeight = 0.0;
	sumR2     = 0.0;
	top_neb   = n_pts - 2;
	for (i=1; i<n_pts; i++)	{ //loop from 1: not including self.
	    weight     = exp ( -nebs[i].distance * nebs[i].distance / trySigma2 );
	    if( weight / (sumWeight + weight) < min_contrib ){
	        top_neb = i - 1; //top_neb is always > 0 because min_contrib is < 1.
	        break;
	    }
	    sumWeight += weight;
	}
   *neighbours = (NbData *)malloc(top_neb * sizeof(NbData));
	for (i=1; i<=top_neb; i++)	{ 
	   (*neighbours)[i-1].id       = nebs[i].id;
	   (*neighbours)[i-1].distance = nebs[i].distance;
	}
	
	//don't need this workspace any more.
	free(nebs);
	
	//return number of neighbours assigned
	return top_neb;
}

void buildNeighbourGraph(KMMDCrd_t    *dh_coords, 
                         KMMDEneTrn_t *denes, 
                         int           n_snaps, 
                         int           n_dims,
                         KMMDCrd_t     test_sig2){
                                       
   //Function to estimate peak force in a given system of training data + sigma parameter
   //and verify if the simulation will be stable.
   

   int       i, j, ii, jj, d;
   float    *r_ij, dr2;
   NbData  **nebLists;
   int      *n_nebs;
   
   nebLists = (NbData  **)malloc(n_snaps * sizeof(NbData  *));
   n_nebs   = (int      *)malloc(n_snaps * sizeof(int));
   
   //build a distance list per-trainpt
   r_ij = (float *)malloc(n_snaps*sizeof(float));
   for( i = 0; i < n_snaps; i++ ){
      r_ij[i] = 0.0;
      ii      = i * n_dims;
      for( j = 0; j < n_snaps; j++ ){
         jj  = j * n_dims;
         dr2 = 0.0;
         for( d = 0; d < n_dims; d++ ){
            dr2 += (dh_coords[ii + d] - dh_coords[jj + d]) * (dh_coords[ii + d] - dh_coords[jj + d]);
         }
         r_ij[j] = sqrt(dr2);
      }
      //build a nearest neighbour list, saving pointId and L2 distance to point.
      n_nebs[i] = computeNeighbourList( n_snaps, r_ij, &(nebLists[i]), test_sig2 );
   }
   free(r_ij);
   
}


