#ifdef USE_MPI
#include <mpi.h>
#endif
#include "config.h"

#include <iostream>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
#include "emil.h"
#include "linearAssignment.h"

void hsc::MC(){


  for( int i = 0; i < numLiquids; i++ ){
    if( linearAssignment[i] != NULL && sim_nStep % assignmentEvery == 0 ){

        int    didJv;
        double dcost;
 
        double oldSumR2, checkSumR2;


        oldSumR2 = 0.0;
        for( int j = 0; j < chainsInLiquid[i]; j++){
              oldSumR2 += liquidLists[i][j]->rijsqW;
        }

        //MERR <<  "enteringAssignment: " <<  didJv 
        //         << " cost:   "    <<  oldSumR2    <<  endl ;

        didJv = linearAssignment[i]->refreshAssignment(&dcost);

        checkSumR2 = 0.0;
        for( int j = 0; j < chainsInLiquid[i]; j++){
            checkSumR2 += liquidLists[i][j]->rijsqW;
        }

        //MERR <<  "reassigned: " <<  didJv  
        //         << " newcost:  "  <<  checkSumR2 <<  endl ;

        //log the number of wells reassigned and the change in energy/k. */
        linearAssignment[i]->assignmentRate->accumulate( didJv / (double)chainsInLiquid[i] );
        linearAssignment[i]->assignmentDeRate->accumulate( checkSumR2-oldSumR2 );

    }
  }

  //do swap moves
  while( swapMoves < nSwap ){ 
    swapMove();  
    swapMoves++;
  }

  //only call relocs when the amber hamiltonian is fully off.
  if( mixFunc_mol(1.0 - amberLambda) == 0.0 && nReloc > 0){ 
      while( relocMoves < nReloc ){ //try turning off the reloc moves gradually as we turn on the molecular potential
              //Attempt to teleport a chain.  
              relocAccept += nonHastingsReloc(); //"long jump" - useful to fill voids in the system
              relocMoves++;
      }
  }

  if( nSwap > 0 ) {

    //log the acceptance rates
    swapAccRate->accumulate( swapAccept / (double)nSwap );
    swapAccept  = 0;
    swapMoves   = 0;
    relocAccept = 0;
    relocMoves  = 0;

   }
}

