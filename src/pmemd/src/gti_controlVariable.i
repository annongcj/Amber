#ifndef GTI_CONTROL
#define GTI_CONTROL

#ifdef C_COMPILER
#define STATIC_CONST_INT static const int
#else /* for FORTRAN */
#define STATIC_CONST_INT integer, parameter:: 
#endif

  STATIC_CONST_INT MaxNumberTIAtom = 750;
  STATIC_CONST_INT MaxNumberTIPair = 300;
  STATIC_CONST_INT MaxNumberNBPerAtom = 3072;

  STATIC_CONST_INT MaxNumberREAFAtom = 500;
  STATIC_CONST_INT MaxNumberRMSDAtom = 200;
  STATIC_CONST_INT MaxNumberRMSDRegion = 5;

  STATIC_CONST_INT GPUPotEnergyTerms = 54;
  STATIC_CONST_INT NumberTIEnergyTerms = 17;  
  STATIC_CONST_INT GPUKinEnergyTerms = 10;
  STATIC_CONST_INT TIEnergyBufferMultiplier = 9;
  STATIC_CONST_INT TIEnergyDLShift = 3;
  STATIC_CONST_INT TIEnergySCShift = 6;

  STATIC_CONST_INT TIExtraShift = 44;
  STATIC_CONST_INT TINumberExtraTerms = 6;

  STATIC_CONST_INT TISpecialShift = 52;
  STATIC_CONST_INT TINumberSpecialTerms = 2;

#ifdef C_COMPILER
  #define TCV_BEG typedef struct {
  #define TCV_END } TypeCtrlVar;
  #define TCV_DOUBLE double 
  #define TCV_INT int
  #define TCV_SP ;

#else /* for FORTRAN */

#define TCV_BEG type, bind(C) :: TypeCtrlVar
#define TCV_END end type TypeCtrlVar
#define TCV_DOUBLE double precision::
#define TCV_INT integer::

#endif

  TCV_BEG
    TCV_DOUBLE scalpha ;
    TCV_DOUBLE scbeta ;
    TCV_DOUBLE	 scgamma ;
    TCV_INT  klambda ;

    TCV_INT addSC ;
    TCV_INT autoAlpha ;
    TCV_INT scaleBeta ;
    TCV_INT batSC;

    TCV_INT eleExp ;
    TCV_INT eleGauss ;
    TCV_INT	eleSC ;

    TCV_INT vdwSC ;
    TCV_INT	vdwExp ;
    TCV_DOUBLE 	vdwCap ;

    TCV_INT 	tiCut ;
    TCV_INT 	cut_sc ;
    TCV_DOUBLE cut_sc0 ;
    TCV_DOUBLE cut_sc1 ;

  TCV_END

#endif
