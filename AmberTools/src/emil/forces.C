#include "config.h"
#ifdef USE_MPI 
#include "mpi.h"
#endif

#include <iostream>
#include <iomanip>
#include <sstream>

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
#include "emil.h"

#define _SSTI_FORCE_WARN  200.0
#define _FORCE_CEIL      145.0


void hsc::computeForces_mix(){
  
  double    dPhiDr, normFac, amberKt;
  double    scaleMolec, scaleAbs, scaleSoft;

  amberKt    = 1.0 / beta; // inside SSTI, forces are in units kT/Angstrom
                           // in amber, forces are kCal/mol / Angstrom
                           //... the amber kT is already in these units.
  
  
  //if( myTaskId == 0 )reportWellPart();
  //if( myTaskId == 0 )writeVtf();

  // *logFile<< "mixing forces, amberSoftCoring: " << amberSoftcoring << " emilSoftForce: " << emilSoftForce << endl; 
  // flushLog();


  //scale down the amber forces
  if( amberSoftcoring || !emilSoftForce ){
    scaleMolec = 1.0;
  }else{
    scaleMolec = mixFunc_mol( 1.0 - amberLambda );
  }
  if( emilSoftForce ){
    scaleSoft  = mixFunc_soft( restLambda );
  }else{
    scaleSoft  = 0.0;
  }
  scaleAbs   = mixFunc_abs( restLambda );
 *logFile << scientific;

  
  //scale down the amber forces
  if( scaleMolec != 1.0 ){
     
  //  *logFile<< "scaling molec: " << scaleMolec << endl; 
     for (int ii = 0; ii < myNatoms; ii++){

         int i = myAtomIds[ii];

         forces[3*i]     *= scaleMolec;
         forces[3*i + 1] *= scaleMolec;
         forces[3*i + 2] *= scaleMolec;
 
#ifdef SSTI_FORCE_WARN
 { 
   double F;
   F  = sqrt(forces[3*i]*forces[3*i]+forces[3*i+1]*forces[3*i+1]+forces[3*i+2]*forces[3*i+2]);
   if(F > SSTI_FORCE_WARN){
       cerr << "Step: " << sim_nStep << " mixed natural force: " << ii << " type: " << part[i].wellType 
            << " value: " << F << " gt warning threshold: " << SSTI_FORCE_WARN << endl;
   }
 }
#endif
 
      }
   }

   //add in the soft force
   softEnergy = 0.0;
   if( emilSoftForce ){
     
  //  *logFile<< "finding and scaling emilSoftForce: " << scaleSoft << endl; 
   
     for (int ii = 0; ii < myNatoms; ii++){
         int i = myAtomIds[ii];

       if( part[i].isRoot
        || part[i].wellType == HSC_POTENTIAL_EINSTEIN_CRYSTAL
        || part[i].wellType == HSC_POTENTIAL_WOBBLIUM ){
         softEnergy +=  softForceAll( &part[i], scaleSoft, &forces[3 * i] );
 
#ifdef SSTI_FORCE_WARN
 { 
   double F;
   F  = sqrt(forces[3*i]*forces[3*i]+forces[3*i+1]*forces[3*i+1]+forces[3*i+2]*forces[3*i+2]);
   if(F > SSTI_FORCE_WARN){
       cerr << "Step: " << sim_nStep << " mixed natural force + Soft: " << ii << " type: " << part[i].wellType 
            << " value: " << F << " gt warning threshold: " << SSTI_FORCE_WARN << endl;
   }
 }
#endif
 
 
        }
     }
#ifdef USE_MPI //collect the total soft energy
     {
         double softEnergyTmp;
  MPI_Allreduce ( &softEnergy, &softEnergyTmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, ourComm );
        softEnergy = softEnergyTmp;
     }
#endif
   }//end if emilSoftForce


  //do the EMIL Hamiltonian
  for (int ii = 0; ii < myNatoms; ii++){

      int i = myAtomIds[ii];
      
      //get force depending on well type (units of kT/Angstrom)
      dPhiDr = part[i].theWell->compute_dPhiOfDist( part[i].rijsqW );
      
      //add in the restraint-based force
      if( part[i].rijsqW > 0.0 ){//need to avoid an instability at r==0.
        
        //scale the abstract force, and
        //apply force as a vector to the well, accounting for PBCs
        normFac = amberKt / sqrt( part[i].rijsqW );
        dPhiDr *= normFac;
        
#ifdef SSTI_FORCE_WARN
      if(dPhiDr * scaleAbs > SSTI_FORCE_WARN){
        
        double phi; 
        phi =  part[i].theWell->computePhiOfDist( part[i].rijsqW );
        cerr << "Step: " << sim_nStep << " Well force: " << i << " type: " << part[i].wellType 
             << " scaled value: " << (dPhiDr * scaleAbs) << " gt warning threshold: " << SSTI_FORCE_WARN 
             << " Phi: " << phi << " r: " << sqrt(part[i].rijsqW) <<  endl;
             
       part[i].theWell->report(&cerr);   
       part[i].theWell->testMe(&cerr);      
      } 
#endif
      
        
        //add the abstract forces
        double scaledPhi;
        scaledPhi = dPhiDr * scaleAbs;
        
        forces[3 * i]      += part[i].rij[0] * scaledPhi;
        forces[3 * i + 1]  += part[i].rij[1] * scaledPhi;
        forces[3 * i + 2]  += part[i].rij[2] * scaledPhi;
        
        //conservative force due to particles with a relative-position restraint to this particle
        if( part[i].isRoot ){
          Particle *pp = part[i].strandAtom;
          while( pp ) {
            
            dPhiDr    = pp->theWell->compute_dPhiOfDist( pp->rijsqW );
            if( pp->rijsqW <= 0. ){
              pp = pp->strandAtom;
              continue;
            }
              
            scaledPhi = dPhiDr * scaleAbs * amberKt / sqrt( pp->rijsqW );
        
            forces[3 * i]      -= pp->rij[0] * scaledPhi;
            forces[3 * i + 1]  -= pp->rij[1] * scaledPhi;
            forces[3 * i + 2]  -= pp->rij[2] * scaledPhi;
            
            pp = pp->strandAtom;
            
          } 
        }
      }  
  }
}
 

void hsc::computeForces(){
  
  double dPhiDr, normFac, externalKt;
  int    atIndex;
  
  externalKt = 1.0 / beta;
  
  atIndex = 0;
  for (int i = 0; i < N; i++){
   
    //get force depending on well type (units of kT/Angstrom)
    if( part[i].isRoot ){
      dPhiDr = part[i].theWell->compute_dPhiOfDist( part[i].rijsqW );
    
      //apply force as a vector to the well, accounting for PBCs
      if( dPhiDr != 0.0 ){
        normFac = externalKt / sqrt( part[i].rijsqW );
        dPhiDr *= normFac;
        forces[atIndex++] = dPhiDr * part[i].rij[0];
        forces[atIndex++] = dPhiDr * part[i].rij[1];
        forces[atIndex++] = dPhiDr * part[i].rij[2]; 
      }
      else{
        forces[atIndex++] = 0.0;
        forces[atIndex++] = 0.0;
        forces[atIndex++] = 0.0; 
      }
    }
    else{//non-root particles have the force saved in rij already
        forces[atIndex++] = part[i].rij[0];
        forces[atIndex++] = part[i].rij[1];
        forces[atIndex++] = part[i].rij[2]; 
    }
  }
}
        
