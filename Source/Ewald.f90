MODULE Ewald
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Ewald
!     |_SUBROUTINE EwaldSetup
!     |_FUNCTION EwaldEnergy
!         |_FUNCTION EwaldRealSpace
!         |_FUNCTION EwaldRealSpaceSelfTerm
!         |_FUNCTION EwaldRecipSpace
!         |_FUNCTION EwaldRecipSpaceSpline
!         |_FUNCTION EwaldAveragePotentialTerm
!     |_FUNCTION EwaldForces
!         |_FUNCTION EwaldRealSpaceForces
!         |_FUNCTION EwaldRecipSpaceForces
!         |_FUNCTION EwaldRecipSpaceForcesSpline
!     |_FUNCTION EwaldStress
!         |_FUNCTION EwaldRealSpaceStress
!         |_FUNCTION EwaldRecipSpaceStress
!         |_FUNCTION EwaldRecipSpaceStressSpline
!         |_FUNCTION EwaldAveragePotentialTermStress
!
!
! DESCRIPTION:
!   This module contains the necessary functions and procedures to compute the
!   Ewald Summation (ion-ion interation) of the current cell, made
!   periodic in all 3 dimensions.  These procedures will work with any lattice
!   structure (such as asymmetric cells) though sub-optimally for very 
!   assymetric cells.
!
!   The Ewald sum splits up the slowly converging coulomb sum into two parts:
!   A real space part and a reciprocal space part.  The real-space sum takes
!   care of the short range interactions and the reciprocal space part takes
!   care of the long range interactions (see EwaldSum blurb for more details).
!
!   The typical Ewald routines scale as O(N^1.5) at best, but particle-mesh
!   ewald scales O(NlnN).  PME subroutines are indicated with "Spline".
!
! CONDITIONS AND ASSUMPTIONS:
!   1.  The overall charge on the system is 0
!   2.  There are no surfaces (surfaces will require additional correction
!       terms)
!   3.  Ion positions are in fractional coordinates
!   4.  The ions are not too close together.  The EWALD sum assumes point
!       charges, which is an assumption that breaks down at short range.
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   1.  Add surface term corrections [4]
!   2.  IMPORTANT:  After testing the error equations, they seem all right,
!       though I still don't completely trust them.  What's more troubling
!       is that when you change ETA (with converged real and recipcutoff)
!       we get different answers to something like the errorTolerance^(1/2).
!       So changing the eta actually changes the answer quite a bit...so if
!       errorTolerance is 1e-10, then our final answer can only be trusted
!       to about 1e-5 (i think).  More analysis on this is needed. It is 
!       still possible there is a bug.
!       ~Greg H, 11/24/2003
!   3.  The structure factor has many terms in common with the energy.
!       Consolidate some of them?
!    
! REFERENCES:
!   [1]  C. Kittel.  Introduction to Solid State Physics.  John Wiley & Sons,
!          Inc., New York, seventh edition, 1996.
!   [2]  A.Y. Toukmaji Jr. and J.A. Board.  Ewald Summation Techniques in 
!          Perspective: A survey.  Computer Physics Communications, 95:19,1996.
!   [3]  Naoki Karasawa and William A. Goddard III.  Acceleration of 
!          convergence for lattice sums.  Journal of Physical Chemistry, 
!          93(21):8, 1989.
!   [4]  In-Chun Yeh and Max Berkowitz.  Ewald summation for systems of slab
!          geometry, Journal of Chemical Physics, 111(3115):7, 1999.
!   [5]  Ulrich Essmann, Lalith Perera, and Max L. Berkowitz.  A smooth
!          particle mesh Ewald method, J. Chem. Phys., 103(19), 1995.
!   [6] M.T. Yin and M.L. Cohen, "Theory of lattice-dynamical properties of 
!       solids:  Application to Si and Ge," Phys. Rev. B 26, 3259 (1982)
!   [7] J. Ihm and M.L. Cohen, "Comment on 'Corretion to Fuch's Calculation of the 
!       Electrostatic Energy of a Wigner Solid'", Phys. Rev. B 21, 3754 (1980)
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  File Created (Greg Ho)
!   11/24/2003  First working version!! (Greg Ho)
!   02/04/2003  Forces added (Greg Ho), Setup tweaked
!   12/14/2007  Particle mesh Ewald subroutines added (Linda Hung)
!
!------------------------------------------------------------------------------

                                 !>> GLOBAL <<!

  USE Constants, ONLY: DP, PI, hartreeToeV, BOHR  
  USE CellInfo, ONLY: ion, element
  USE Fourier_NEW, ONLY: FFT_NEW
  USE Fourier_NEW, ONLY: FFT_STD_STATE
  USE MathFunctions, ONLY: Vecmul
  USE MathFunctions, ONLY: Norm
  USE MathFunctions, ONLY: Inverse
  USE MathFunctions, ONLY: Volume
  USE Timer
  USE OUTPUT, ONLY : WrtOut
  USE OutputFiles, ONLY : outputUnit
  USE CellInfo, ONLY: k1G, k2G, k3G, k3Goff
  USE CellInfo, ONLY: n1G, n2G, n3G, n3Goff
  USE CellInfo, ONLY: m1G, m2G, m3G, m123G
  USE MPI_Functions
#ifdef __USE_PARALLEL
  USE SYS, ONLY: numIonLoc ! Number of ions on the local processor
  USE SYS, ONLY: numIonInit ! Smallest number ion that is on local processor
#endif

  IMPLICIT NONE


  ! Here are the user-specified parameters and their default values. 
  REAL(KIND=DP) :: errorTolerance = 1.e-10_DP / hartreeToeV
  ! maximum error allowed in either the real or recip space energy (eV)
  !
  REAL(KIND=DP) :: maxRealCutoff = 12.0_DP / BOHR   
  ! Maximum interaction cutoff distance for the ions in real-space.  Each 
  ! ion sees nothing in real-space beyond a radius realCutoff of itself. 
  ! Note: this is not the realCutoff used in the program though, it's just the 
  ! maximum it can be. Generally, higher real cutoff = slower.
  !
  REAL(KIND=DP) :: etaIncrement = 0.01_DP * BOHR
  ! Increment used for optimizing eta. The larger this number, the harder to 
  ! converge on optimum number, but the faster the parameters are determined
  !
  REAL(KIND=DP) :: recipCutoffIncrement = 0.01_DP * BOHR
  ! Same as etaIncrement, except for recip cutoff
  !
  REAL(KIND=DP) :: realCutoffIncrement = 0.01_DP * BOHR
  ! Same as etaIncrement, except for recip cutoff


  REAL(KIND=DP), PRIVATE :: cellVol
  ! The volume of our calculational cell
  !
  REAL(KIND=DP), PRIVATE :: sumCharge
  ! The sum of all the charges of all ions in calc. cell
  !
  REAL(KIND=DP), PRIVATE :: sumChargeSquared
  ! The sum of squared charges of all ions in calc. cell
  !
  REAL(KIND=DP), PRIVATE :: sqrtPi
  !

  REAL(KIND=DP), PRIVATE, DIMENSION(3,3) :: cellRecip         
  ! 3x3 matrix containing reciprocal space lattice vectors.
                

  ! Here are parameters shared by the module that will be changed and edited.
  INTEGER, PRIVATE :: nArecip, nBrecip, nCrecip
  ! number of replications you need of the lattice 
  ! vector in order to form a shape that completely 
  ! encloses the sphere defined by recipCutoff.
  !
  INTEGER, PRIVATE :: nArecipSpl, nBrecipSpl, nCrecipSpl
  ! Same as above, but for use with particle mesh
  ! ewald approximation
  !
  INTEGER, PRIVATE :: numIon
  ! The number of ions in the calc. cell
  !
  INTEGER, PRIVATE :: numIonType
  ! The number of types of ions in the calc. cell
  !

  REAL(KIND=DP), PRIVATE :: eta
  ! Parameter which controls how much of the calculation
  ! is done in real space vs reciprocal space.  
  !
  REAL(KIND=DP), PRIVATE :: realCutoff
  ! The optimum real-space cutoff radius (calculated)
  !
  REAL(KIND=DP), PRIVATE :: recipCutoff
  ! The optimum reciprocal-space cutoff radius (calculated)
  !
  REAL(KIND=DP), PRIVATE :: realError
  ! The maximum error in real space
  !
  REAL(KIND=DP), PRIVATE :: recipError
  ! The maximum error in reciprocal space
  !
  REAL(KIND=DP), PRIVATE :: lattice1, lattice2, lattice3
  ! Cartesian length of 1st,2nd,3rd lattice vector
  !
  LOGICAL, PRIVATE :: ionScreen=.TRUE.  
  ! Whether efficient ion screening subroutines will be used
  !

  ! If memory is in free supply, subroutines can be changed to save values
  ! in arrays below instead of recalculating each time.

  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: tempRA
  ! An all-purpose work array for particle-mesh ewald
  !
  REAL(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: qIonTable         
  ! An array Q defined as equation 4.6 in [5]
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE, PRIVATE :: CCStructure       
  ! Approximate complex conjugate of the structure factor
  !

  REAL(KIND=DP) :: ionIonEnergy    
  ! The ion-ion energy, stored after calculation
  ! from every time the ions are moved

CONTAINS


SUBROUTINE EwaldSetup(cellReal, ionTable, elementTable, useSpline, sizeX, sizeY, sizeZ)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This subroutine calculates eta, the parameter that specifies how much of 
!   the Ewald sum is done in real space and reciprocal space, and also 
!   realCutoff and recipCutoff, the cutoff radius for the real and reciprocal
!   space calculations, respectively.
! 
!   This routine should be run every time the cell lattice parameters have 
!   changed.
!
!   For reasons of maximal efficiency, we want to do as much of this
!   calculation in reciprocal space O(N) as possible and as little of it in
!   real space O(N^2) as possible.  To this effect, we specify that we ONLY
!   ever want to sum in real space by interacting each ion in the 
!   calculational cell with only ions (including periodic images) within a
!   cutoff radius realCutoff of itself.  This cutoff, then, determines the
!   minimum value of eta for a given errorTolerance (the smaller eta is, the 
!   slower the real sum is to converge, and therefore the more error 
!   generated by a given cutoff).  
!
!   However, a bigger eta would make the real-space calculation go even faster,
!   so we want to maximize eta.  But, the greater eta is, the greater
!   recipCutoff needs to be.
!
!   This function first guesses at eta to make work in real and recip. space
!   equally loaded.  Then it sees if this eta is good enough for the realCutoff
!   specified (when sent to this routine, that variable contains the maximum
!   realCutoff allowed for the system.)  If so, it makes realCutoff smaller
!   until its optimized, if not it increases eta until its right at the border.
!   Then it calculates an appropriate recipCutoff for the error tolerance.
!
!   The relative speeds of the reciprocal and real space calculation should be
!   monitored closely.  Refer to Stuart's old OFDFT code for a more 
!   sophisticated way of calculating eta.
!
!   The functions to estimate the error were obtained from ref [3]. 
!
!   This function also initializes n(A,B,C)real and recip.  In order to do
!   that, (first we'll talk about the real case) we have to have some way 
!   (general for 
!   non-orthorhombic cells) to enclose for any given ion a sphere of radius 
!   realCutoff.  We do this by finding a shape centered on the sphere made 
!   by a set number of replications of all the lattice vectors that will 
!   enclose the sphere.  For example, for an orthorhombic cell we would get 
!   a cube that has dimensions equivalent to the diameter of the sphere.  
!   We then inspect each point in this cube to see whether it is within the 
!   cutoff radius.  This method is very wasteful for non-orthorhombic cells,
!   but I couldn't find any other way of being sure to get all the real-space
!   points within a cutoff radius.
!
!   The way we get this shape is by the following method.  For each lattice
!   vector, we consider the PLANE made by the other two lattice vectors.  
!   The lattice vector in consideration will have to be long enough to seperate
!   the two PLANES (there is one plane on either side of our lattice vector)
!   by an amount equal to the diameter of the sphere.  So, we take the NORMAL
!   with respect to the planes, and figure out (by a dot product) how many
!   of our lattice vectors we need to be able to fit a diameter realCutoff 
!   between these two points.  For instance, if Na is the number of lattice
!   vectors of a we need, b and c are the other lattice vectors, and ha is
!   the reciprocal lattice vector for a,
!
!   Na = 2*Rcutoff / (a dot ((b x c) / |b x c|))
!      = 2*Rcutoff * |b x c| / V   
!      = 2*Rcutoff * ha / 2 / pi = Rcutoff * ha / pi
!
!   We then divide this number by 2 and loop from -Na to Na, and then also we
!   add 1 to account for the fact that there is a distance between the two
!   interaction ions in the calculational cell that is not taken into account
!   here.  This process is simply repeated for all lattice vectors.
!
!   The only thing we do differently for the recip space sphere 
!    is that we take advantage of symmetry
!   in reciprocal space.  (we can't do this in real space because the ions are
!   displaced relative to on another).  To exploit inversion symmetry, we will
!   only sum over half the cell.  i_a are the number of reps of the first 
!   recip lattice vector, i_b are the number of reps of the second and i_c
!   are the number of reps of the third:   
!
!   i_a ranges over 0 TO na ONLY.                         
!   i_b ranges over 0 TO nb when i_a=0 and over -nb to nb OTHERWISE. 
!   i_c ranges over 1 TO nc when i_a = i_b = 0 and over -nc to nc OTHERWISE.
! 
!   ... and  we will multiply by two at the end.
!
! GLOBAL/MODULE VARIABLES CHANGED:
!   In This Module:
!     eta
!     realCutoff
!     recipCutoff
!     realError
!     recipError
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   I could use a reciprocal space vector to calculate the initial eta instead
!   of how I do a norm and a cross product now. 
!
!   I could looking at recipCutoff from a set number, like say start at 0.5.
!
!   I'M NOT SURE THAT THE FORMULAS FOR JUDGING THE MIN CUTOFFS ARE ACCURATE,
!   OR ACCURATELY IMPLEMENTED.  PLEASE CHECK.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------
  USE Mathfunctions, only : Metric3d

  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellReal
  ! The real-space lattice vectors
  !
  TYPE(ion), DIMENSION(:), INTENT(IN) :: ionTable
  ! Type + location of each ion, sorted by ion type
  ! We can get the number of ions (formerly numIon) by 
  ! the size of this array.
  !
  TYPE(element), DIMENSION(:), INTENT(IN):: elementTable     
  ! index of the FIRST ion of each type in iontable
  ! This array actually has one more element than
  ! numIonType, because the last entry here numIon+1
  ! to make it easy to use this to loop through stuff  
  !
  LOGICAL, INTENT(IN) :: useSpline
  ! Whether to use cardinal b-spline approximation
  !
  INTEGER, INTENT(IN), OPTIONAL :: sizeX, sizeY, sizeZ   
  ! Size of electron grid, used in particle mesh ewald
  !

                        !>> INTERNAL VARIABLES <<!

  INTEGER :: i_type
  ! a counter, counts through all the types of ions present
  !
  INTEGER :: numIonsInType     
  ! Used as a temporary variable for the number of ions in
  ! a specific type

  REAL(KIND=DP) :: averageChargeSquared
  REAL(KIND=DP) :: rmet(3,3)
  REAL(KIND=DP) :: gmet(3,3)

                          !>> INITIALIZATION <<!

  numIon = SIZE(ionTable)
  numIonType = SIZE(elementTable) - 1
  
  ! Compute metric tensor 
  CALL Metric3d(cellReal, rmet, gmet)

  sqrtPi = SQRT(pi)
  realCutoff = maxRealCutoff
 
  cellVol = Volume(cellReal) 
  cellRecip = TRANSPOSE(Inverse(cellReal)) * 2.0_DP * pi 
 
  ! Initialize sumCharge, sumChargeSquared
  sumCharge = 0.0_DP
  sumChargeSquared = 0.0_DP
  DO i_type = 1, numIonType  
    ! loop over # ion TYPES
    ! calculate the charge squared, multiply by the number of ions of this
    ! type, and add to sumChargeSquared
    numIonsInType = (elementTable(i_type+1)%firstIonID &
                     - elementTable(i_type)%firstIonID)
    sumChargeSquared = sumChargeSquared + &
      numIonsInType * elementTable(i_type)%charge**2 
    sumCharge = sumCharge + &
      numIonsInType * elementTable(i_type)%charge
  ENDDO

  averageChargeSquared = sumChargeSquared / REAL(numIon,KIND=DP)


                           !>> FUNCTION BODY <<!
    
  ! Chen's comment:
  ! I finally think we should just keep using one method for determining 
  ! eta. I removed the previous method by me and Greg, both method seems
  ! have some problem.
  !
  ! This method below is written by Linda, and is consistent with 
  ! the particle-mesh Ewald method, you should not change it unless you know 
  ! particle-mesh Ewald method well.
  !
  ! Linda's comment:
  ! This differs from above in that the Ewald reciprocal space grid is taken
  ! to be the same size and grid density as the reciprocal space grid for
  ! the electron density grid.  Perhaps in the future we will want to try a
  ! smaller grid?
  nArecipSpl = sizeX
  nBrecipSpl = sizeY
  nCrecipSpl = sizeZ

  ! Reciprocal cutoff is calculated from the grid dimensions
  recipCutoff = MIN((nArecipSpl-1)/Norm(cellReal(:,1)), &
                    (nBrecipSpl-1)/Norm(cellReal(:,2)), &
                    (nCrecipSpl-1)/Norm(cellReal(:,3))) * pi

  ! Pick a relatively large initial eta to allow as many computations in
  ! reciprocal space as possible.  From a few tests, 1.2 seems a decent
  ! upper bound, unless we are using extreme grid densities (>2000 eV).
  ! Given a KE cutoff, eta varies by < 0.5 depending on the number of ions.
  eta = 1.2_DP

  ! Now iterate eta downwards until desired error is reached, given that the
  ! reciprocal space cutoff is fixed.  From Eq. 38 of ref. 1.  
  DO 
    recipError = REAL(numIon,KIND=DP)**2 * averageChargeSquared / sqrtPi &
                 * eta * ERFC(recipCutoff/(2._DP*eta))
    IF (recipError < errorTolerance) EXIT
    eta = eta - etaIncrement
  END DO

  ! Pick a relatively small realCutoff, and adjust until desired accuracy is
  ! found.  A good lower bound is 6, which includes only the nearest
  ! neighbors for bulk equilibrium Al
  realCutoff = 6._DP

  DO
    realError = pi * REAL(numIon,KIND=DP)**2 * averageChargeSquared &
                / (cellVol * eta**2) * ERFC(realCutoff * eta)
    IF (realError < errorTolerance) EXIT   
    realCutoff = realCutoff + realCutoffIncrement
  END DO

  IF (.NOT. useSpline) then
    ! Now, with this real and recip cutoff, compute the sphere that we must
    ! perform all our calculations under.
    ! The goal here is to try to figure out how many lattice vectors we need
    ! for each of the 3 directions to make a shape that will always enclose our
    ! sphere, in recip. space.  See above blurb.
    nArecip = CEILING(recipCutoff * Norm(cellReal(:,1)) / pi / 2) 
    nBrecip = CEILING(recipCutoff * Norm(cellReal(:,2)) / pi / 2) 
    nCrecip = CEILING(recipCutoff * Norm(cellReal(:,3)) / pi / 2)
  ENDIF

 
  ! output some info about EwaldStep
  WRITE(outputUnit,'(/A)' ) " ------------------------------------------------------------------------------"
  WRITE(outputUnit,'(A)')   "                           SETUP EWALD "
  WRITE(outputUnit,'('' (Ewald) eta                             : '',G12.4)'), eta
  WRITE(outputUnit,'('' (Ewald) Real Space Cutoff               : '',G12.4,'' (Bohr)'')') realCutoff
  WRITE(outputUnit,'('' (Ewald) Reciprocal Cutoff               : '',G12.4,'' (1/Bohr)'')') recipCutoff

  RETURN

END SUBROUTINE EwaldSetup


FUNCTION EwaldEnergy(cellReal, ionTable, elementTable, useSpline) 
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function computes the ion-ion interaction energy (in eV) of current
!   system, assuming that the system is made periodic in three dimensions.  
!
!   This is the main program that runs all the right functions to calculate
!   the Ewald ENERGY, and then adds them together.
!
!   The Ewald sum splits up the slowly converging coulomb sum into two parts:
!   A real space part and a reciprocal space part.  The real-space sum takes
!   care of the short range interactions and the reciprocal space part takes
!   care of the long range interactions.  The real-space part is further
!   split up into 2 terms, a "self" term and a term for the system interacting
!   with its periodic images close to it. In order to perform the real-space
!   sum, we interact each ion with all the ions around (including periodic 
!   images) it within a spherical region whose radius is R_cut. (realCutoff).
!   The reciprocal space sum is also truncated.  In addition, there is an 
!   additional term to account for the uniform neutralizing charge background.
!
!   There are three main variables that control the convergence and accuracy
!   of these sums.  In the terminiology of [2], They are R_cut, the real-space
!   cutoff term, m_max, the maximum number of terms to sum over in reciprocal 
!   space, and eta, a term that controls how much of the sum is done in real 
!   space and how much is done in reciprocal space.
!
!   The subprograms in this main subprogram uses terminology found in [2], 
!   and basically seek to implement equations 3, 5, and 16 of [2] (we use eta 
!   instead of their alpha).  However, the units and constants we use are a 
!   bit different, and we use Atomic Units
!
!   The flow is as follows:  First we find an appropriate value of eta,
!   R_cut, and m_max, and then we run the real and reciprocal space routines
!   to compute the energy, and add everything up.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Maybe we don't have to run EwaldSetup also every time we do this.
!   
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------
  USE OUTPUT, ONLY: PrintEwald            
  ! Subroutine to do pretty printing of output
  !

  IMPLICIT NONE

                           !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP) :: EwaldEnergy               
  ! The answer!! (Ewald energy)
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellReal
  ! The real-space lattice vectors
  !
  TYPE(ion), DIMENSION(:), INTENT(IN) :: ionTable
  ! Type + location of each ion, sorted by ion type
  ! We can get the number of ions (formerly numIon) by 
  ! the size of this array.
  !
  TYPE(element), DIMENSION(:), INTENT(IN):: elementTable     
  ! index of the FIRST ion of each type in iontable
  ! This array actually has one more element than
  ! numIonType, because the last entry here numIon+1
  ! to make it easy to use this to loop through stuff
  !
  LOGICAL, INTENT(IN) :: useSpline
  ! Whether to use the particle-mesh ewald ONlnN method
  !

                        !>> INTERNAL VARIABLES <<!
  
  REAL(KIND=DP), DIMENSION(4) :: ewaldComponents   
  ! Array that holds all the values from the various components
  !
  REAL(KIND=DP), DIMENSION(4) :: ewaldTimes        
  ! Array holding the execution time of each procedure
  !
  TYPE(stopwatch):: watch ! Timer

                        !>> INITIALIZATION <<!

  CALL Title("Ewald::EwaldEnergy")
  CALL StartClock('EwaldEnergy')

                        !>> FUNCTION BODY <<!

  WRITE(outputUnit,*) " "
  WRITE(outputUnit,'(A)', ADVANCE="NO") " Calculating Ion-Ion Energy ..." 
  flush(outputUnit)
  watch=TimerStart()

  ! Here, we simply call the subroutines to calculate our energy and then exit,
  ! timng each procedure.

  ewaldComponents(1) = EwaldRealSpace(eta,realCutoff)
  ewaldTimes(1) = TimerStop(watch)

  watch = TimerStart()
  ewaldComponents(2) = EwaldRealSpaceSelfTerm(eta)
  ewaldTimes(2) = TimerStop(watch)

  watch = TimerStart()
  IF (useSpline .EQV. .TRUE.) THEN
    ewaldComponents(3) = EwaldRecipSpaceSpline(eta, recipCutoff)
  ELSE
    ewaldComponents(3) = EwaldRecipSpace(eta, recipCutoff)
  END IF

  ewaldTimes(3) = TimerStop(watch)

  watch = TimerStart()
  ewaldComponents(4) = EwaldAveragePotentialTerm(eta)
  ewaldTimes(4) = TimerStop(watch)
 
  EwaldEnergy = SUM(ewaldComponents)


  WRITE(message,'(A, F10.3, A)') "                   complete.", SUM(ewaldTimes), " s"
  CALL WrtOut(outputUnit,message)

  CALL PrintEwald(EwaldEnergy, eta, realCutoff, recipCutoff, &
                  realError, recipError, ewaldComponents, ewaldTimes)

  CALL StopClock('EwaldEnergy')

  RETURN

CONTAINS


FUNCTION EwaldRealSpace(eta, realCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function takes eta and a cutoff radius and outputs the real-space
!   contribution to the Ewald energy for this system via eqn. 3 in ref. 2.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   The way this is implemented is very wasteful for non-orthorhombic cells.
!   I can't figure out a better way to do it and still get all the cells,
!   however, but if anyone else has a way to get all the points within a
!   sphere of radius realCutoff around an ion, please implement it.
!
!   Also, another thing I might be able to do is to add the coordinates
!   of the lattice vector after each loop, so that we can eliminate the mult.
!   needed to compute the difference in icoord and jcoord.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  USE Mathfunctions, only: Metric3d

  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: eta
  ! The Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums 
  !
  REAL(KIND=DP), INTENT(IN) :: realCutoff   
  ! The optimum real-space cutoff radius
  !
  REAL(KIND=DP) :: EwaldRealSpace 
  ! the answer

                    !>> INTERNAL VARIABLES <<!
  CHARACTER(LEN=500) :: message

  INTEGER :: i_type, j_type
  ! counter, used to parse over all types of ions
  !
  INTEGER :: i, j      
  ! counter, used to parse over all ions of a given type
  !
  INTEGER :: i_a, i_b, i_c
  ! counter, used to parse over all reps of 1st, 2nd, 3rd lattice vector
  !
  REAL(KIND=DP) :: r1, r2, r3
  !
  REAL(KIND=DP) :: realCutoffSq
  ! realCutoff**2
  !
  REAL(KIND=DP) :: charge_ij
  ! The charge of ion i times the charge of ion j
  !
  REAL(KIND=DP) :: dist_ij
  ! The distance of ion i from ion j
  !
  INTEGER :: nr, newr
  ! 
  REAl(KIND=DP) :: fr(3)
  !
  REAL(KIND=DP) :: rmet(3,3), gmet(3,3)
  !
  TYPE(stopwatch):: watch
  ! Timer

                           !>> INITIALIZATION <<!

  watch = TimerStart()

  EwaldRealSpace = 0.0_DP
  realCutoffSq=realCutoff**2
  
  ! compute the metric tensor for realspace :rmet
  ! and for g-space : gmet
  CALL Metric3d(cellReal,rmet,gmet)

                           !>> FUNCTION BODY <<!

  DO i_type = 1, numIonType
#ifdef __USE_PARALLEL
    ! Skip this ion type if relevant ion calculations are not on this processor
    IF ((elementTable(i_type+1)%firstIonID-1<numIonInit) .OR. &
        (elementTable(i_type)%firstIonID>numIonInit+numIonLoc-1)) CYCLE
#endif
    DO j_type = 1, numIonType
    ! Calculates the product of the charges
    charge_ij = elementTable(i_type)%charge * elementTable(j_type)%charge
#ifdef __USE_PARALLEL
      DO i = MAX(elementTable(i_type)%firstIonID,numIonInit), &
             MIN(elementTable(i_type+1)%firstIonID-1, numIonInit+numIonLoc-1)
#else
      DO i = elementTable(i_type)%firstIonID, &
             elementTable(i_type+1)%firstIonID-1
#endif
        DO j = elementTable(j_type)%firstIonID, &
               elementTable(j_type+1)%firstIonID-1

          ! Determine the minimum fractional distance between ions,
          ! taking images into account, so that fractional distance in each
          ! dimension has absolute value less than 0.5
          fr(:) = ionTable(i)%coord(:)-ionTable(j)%coord(:)
          fr = fr - ANINT(fr)

          ! This loops over the geometrical shape that encloses the sphere   
          nr = 1
          DO ! loop over shells
             newr = 0
             
             DO i_a = -nr, nr 
               DO i_b = -nr, nr
                 DO i_c = -nr, nr
                   
                   ! we only loop over current shell, not inside shell
                   IF (ABS(i_a)==nr .OR. ABS(i_b)==nr .OR. ABS(i_c)==nr .OR. &
                       nr==1) THEN

                     r1 = REAL(i_a,KIND=DP) + fr(1)
                     r2 = REAL(i_b,KIND=DP) + fr(2)
                     r3 = REAL(i_c,KIND=DP) + fr(3)

                     ! Distance between ions, squared
                     dist_ij = rmet(1,1)*r1**2+rmet(2,2)*r2**2+rmet(3,3)*r3**2+&
                       2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r3*r1)

                     ! (1) Note: ERFC(8) is about 1.1e-29,  charge_ij is 
                     !     on the order of 1e1, dist_ij 's lower bound is 1e-1
                     !     bohr.  We assume two ions cannot get close to within
                     !     1e-1 bohr. Therefore we can neglect any contribution
                     !     from eta*dist_ij > 8
                     ! (2) dist_ij cannot = 0, to avoid zero in denominator
                     ! RECALL that dist_ij is actually storing the square of
                     ! the distance between ions at this point.
                     IF (dist_ij <= realCutoffSq .AND. dist_ij>1E-40) THEN
                       newr = 1
                       dist_ij=SQRT(dist_ij)
                       EwaldRealSpace = EwaldRealSpace + &
                                        charge_ij * ERFC(eta*dist_ij)/dist_ij
                     END IF

                  ENDIF ! on shell (not inside shell)
                  
                 END DO ! i_c
               END DO ! i_b
             END DO ! i_a

             IF (newr==0) EXIT

             nr = nr + 1

          ENDDO ! loop over nr
           
        END DO ! j
      END DO ! i
    END DO ! j_type
  END DO ! i_type

  ! Divide our final result by two, as per formula 3 of [2].
  EwaldRealSpace = EwaldRealSpace / 2._DP

#ifdef __USE_PARALLEL
  ! Use dist_ij as temporary storage for EwaldRealSpace, locally
  dist_ij = EwaldRealSpace
  CALL MPI_ALLREDUCE(dist_ij, EwaldRealSpace, 1, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN EwaldRealSpace***"
#endif

  RETURN

END FUNCTION EwaldRealSpace


FUNCTION EwaldRealSpaceSelfTerm(eta)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function will compute the real-space self-term contribution to the
!   Ewald sum (Eqn. 5 of [2]).  Uses units, of course, from [3].  
!   It does virtually nothing else.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: eta  
  ! or alpha, the Ewald convergence parameter that determines 
  ! relative rates of convergence between real and recip sums
  !
  REAL(KIND=DP) :: EwaldRealSpaceSelfTerm 
  ! the answer
  !

                    !>> INTERNAL VARIABLES <<!
                    !>> INITIALIZATION <<!
                    !>> FUNCTION BODY <<!
  EwaldRealSpaceSelfTerm = - eta/sqrtPi * sumChargeSquared

  RETURN

END FUNCTION EwaldRealSpaceSelfTerm


FUNCTION EwaldRecipSpace(eta, recipCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This program implements equation 14 of [2], using the units of [3].  
!   There's not much to say about it, except that we have the same issue of 
!   having to get a shape to enclose the sphere made by recipCutoff.  We
!   solve this problem the SAME WAY as above in EWALDREALSPACE, except we don't
!   add 1 at the end (because we don't have the displacement problem anymore).
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Again, this isn't all that efficient for non-orthonormal cells because it
!   is forced to sample more than it potentially has to.  Some smart guy can
!   come around and optimize it.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: eta
  ! eta, the Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums
  !
  REAL(KIND=DP), INTENT(IN) :: recipCutoff     
  ! the optimum reciprocal-space cutoff radius 
  !
  REAL(KIND=DP) :: EwaldRecipSpace 
  ! the answer
  !

                        !>> INTENRAL VARIABLES <<!
  INTEGER :: i_a, i_b, i_c
  ! Counter, goes over reps of first recip. latt. vector
  !
  INTEGER :: i_b_lowLimit
  ! the number of reps the i_b loop should start from.
  !  Its either 0 or -nb, depending on the value of i_a
  !
  INTEGER :: i_c_lowLimit
  ! the number of repsthe i_c loop should start from.  
  ! Its either 0 or -nc, depending on the value of i_a and
  !  i_b. (see above blurb)
  !
  INTEGER :: k
  ! counter to loop over all ions    
  !
  REAL(KIND=DP) :: structureFactorSin
  ! The SIN part of structure factor 
  ! (see eq. 11 of [3]).
  !
  REAL(KIND=DP) :: structureFactorCos
  ! The COS part of structure factor
  !
  REAL(KIND=DP) :: m_squared
  ! The dot product of m and m, where m is a 
  ! reciprocal lattice vector.
  !
  REAL(KIND=DP), DIMENSION(3) :: m 
  ! A 3-D reciprocal lattice vector
  !
  REAL(KIND=DP), DIMENSION(3) :: realCoord 
  ! The real-space coordinates of an ion, in Cartesian coordinates
  !
#ifdef __USE_PARALLEL
  INTEGER :: myRank, &     ! Rank for each processor
    mySize, &              ! Number of processors
    mpiErr, &              ! Error status from initializing, finalizing MPI
    countRecip             ! Count the reciprocal vector,

   REAL(KIND=DP) :: tmpE
#endif

                        !>> INITIALIZATION <<!
  EwaldRecipSpace = 0.0_DP

#ifdef __USE_PARALLEL 
  tmpE=0.0_DP
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, myRank, mpiErr)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mySize, mpiErr)
#endif

                         !>> FUNCTION BODY <<!
  ! This loops over the geometrical shape that encloses the sphere in
  ! reciprocal space.  We will try to exploit the symmetry of the reciprocal
  ! space lattice (see above blurb).
  ! This loops over the first recip. lattice vector
!  WRITE(*,*) "First recip vector ", nArecip
!  WRITE(*,*) "Second vector ", nBrecip
!  WRITE(*,*) "Third vecotr ", nCrecip
  
#ifdef __USE_PARALLEL
  countRecip=0
#endif
  
  DO i_a = 0, nArecip
    ! Sets the lower limit on i_b depending on value of i_a
    IF (i_a == 0) THEN
      i_b_lowLimit = 0
    ELSE
      i_b_lowLimit = -nBrecip
    ENDIF

    ! This loops over second recip. lattice vector
    DO i_b = i_b_lowLimit, nBrecip
     ! Sets the lower limit on i_c depending on value of i_a and i_b
      IF (i_a == 0 .AND. i_b == 0) THEN
        i_c_lowLimit = 1
      ELSE
        i_c_lowLimit = -nCrecip
      ENDIF

      ! This loops over the third recip. lattice vector
      DO 1000 i_c = i_c_lowLimit, nCrecip


#ifdef __USE_PARALLEL
        ! mohan add
        countRecip = countRecip+1
        IF( MOD(countRecip,mySize) .NE. myRank ) THEN 
          GO TO 1000
        END IF
#endif

        ! compute the reciprocal lattice vector m we're currently at
        m = i_a * cellRecip(:,1) + i_b * cellRecip(:,2) + i_c * cellRecip(:,3)

        ! If the length of m is within the cutoff, then add its contribution.
        IF (Norm(m) <= recipCutoff) THEN          

          ! This loops over all the ions, computes sturcture factor
          structureFactorSin = 0.0_DP
          structureFactorCos = 0.0_DP              
   
          DO k = 1, numIon
            realCoord = Vecmul(cellReal, ionTable(k)%coord)        
            structureFactorSin = structureFactorSin &
              + elementTable(ionTable(k)%elementID)%charge &
                * SIN(DOT_PRODUCT(m, realCoord))
            structureFactorCos = structureFactorCos &
              + elementTable(ionTable(k)%elementID)%charge &
                * COS(DOT_PRODUCT(m, realCoord))
          END DO
          ! Now put the formula together and add the contribution to the total
          m_squared = Dot_Product(m,m)
          EwaldRecipSpace = EwaldRecipSpace &
                            + exp(-(m_squared / eta**2 /4._DP)) / m_squared &
                              * (structureFactorSin**2 + structureFactorCos**2)
        ENDIF

1000    CONTINUE  ! i_c
    END DO ! i_b
  END DO ! i_a
  
  ! Multiply by 2*pi, which is the constant from eq. 7 of [3], and also by an
  ! extra factor of 2 to account for our symmetry exploitation.
  EwaldRecipSpace = EwaldRecipSpace * 4.0_DP * pi / cellVol

#ifdef __USE_PARALLEL
  CALL MPI_ALLREDUCE(EwaldRecipSpace, tmpE, 1, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  EwaldRecipSpace=tmpE
#endif

! test
!   IF(myRank == 0) THEN
!   WRITE(*,*) " EwaldRSpace=", EwaldRecipSpace
!   END IF
!  STOP

  RETURN

END FUNCTION EwaldRecipSpace


FUNCTION EwaldRecipSpaceSpline(eta, recipCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This program implements equation 4.7 of [5], using consistent units with
!   the remainder of the module.  It calculates ion-ion energy using particle-
!   mesh Ewald.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/09/2007  Function Created (Linda Hung)
!
!------------------------------------------------------------------------------
  USE CBSpline, ONLY : BSplineProduct
  USE CBSpline, ONLY : splineOrder
  USE IonElectronSpline, ONLY : FillQIonTable
  USE PlaneWave, ONLY: qVectors
  USE MathFunctions, ONLY: Norm

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: eta
  ! eta, the Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums
  !
  REAL(KIND=DP), INTENT(IN) :: recipCutoff     
  ! the optimum reciprocal-space cutoff radius 
  !
  REAL(KIND=DP) :: EwaldRecipSpaceSpline 
  ! the answer
  !

                        !>> INTENRAL VARIABLES <<!
  INTEGER :: i_a, i_b, i_c
  ! counter, goes over reps of 1st, 2nd, 3rd recip. latt. vector
  !
  INTEGER :: j, k
  ! counter to loop over ions
  !
  INTEGER :: allocateStatus
  !
  REAL(KIND=DP) :: m_squared
  ! The dot product of m and m, where m is a 
  ! reciprocal lattice vector.
  !
  INTEGER, DIMENSION(3) :: m
  ! A 3-D reciprocal lattice index
  !
  REAL(KIND=DP), DIMENSION(3) :: mCell
  ! A 3-D reciprocal lattice vector
  !
  REAL(KIND=DP) :: qNorm
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmpRecip

                        !>> INITIALIZATION <<!

  ALLOCATE(tempRA(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating c array for particle mesh Ewald"
  STOP
  END IF

  ALLOCATE(qIonTable(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating q convolution table for particle mesh Ewald"
    STOP
  END IF

  ALLOCATE(CCStructure(k1G,k2G,k3G),stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating CCStructure for particle mesh Ewald"
  STOP
  END IF

  ALLOCATE(tmpRecip(k1G,k2G,k3G),stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating tmpRecip for particle mesh Ewald"
  STOP
  END IF

  EwaldRecipSpaceSpline = 0.0_DP

                         !>> FUNCTION BODY <<!

  ! For each ion type, calculate qIonTable
  qIonTable = 0._DP
  DO j = 1, SIZE(elementTable) - 1
    CALL FillQIonTable(tempRA, ionTable(elementTable(j)%firstIonID: &
                                   (elementTable(j+1)%firstIonID-1)), &
                       splineOrder, m1G, m2G, m3G, n3Goff)

    qIonTable = qIonTable + tempRA*elementTable(j)%charge
  END DO

  ! FFT for the complex conjugate of the structure factor
  CALL FFT_NEW(FFT_STD_STATE, qIonTable, tmpRecip)
  CALL BSplineProduct(tmpRecip, CCStructure)

  DEALLOCATE(qIonTable)
  DEALLOCATE(tempRA)

  ALLOCATE(tempRA(k1G, k2G, k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating CCStructureNorm array for particle mesh Ewald"
  STOP
  END IF

  tempRA = CCStructure*CONJG(CCStructure)

  DEALLOCATE(CCStructure)

  DO i_c = 1, k3G 
    DO i_b = 1, k2G
      DO i_a = 1, k1G
        qNorm = Norm(qVectors(i_a,i_b,i_c,:))
        m_squared = qNorm*qNorm
        ! mohan add
        IF (m_squared <= recipCutoff**2 .AND. m_squared>0.d0) THEN
          IF (i_a==1) THEN ! This plane is its own reflection
            EwaldRecipSpaceSpline = EwaldRecipSpaceSpline + &
              EXP(-m_squared/eta**2/4._DP)/m_squared*tempRA(i_a,i_b,i_c)
          ELSE ! Add factor of 2 for symmetry for all other planes
            EwaldRecipSpaceSpline = EwaldRecipSpaceSpline + &
              2._DP*EXP(-m_squared/eta**2/4._DP)/m_squared*tempRA(i_a,i_b,i_c)
          END IF
        END IF
      END DO
    END DO
  END DO

  ! Finally multiply in conversion factors
  EwaldRecipSpaceSpline = EwaldRecipSpaceSpline * 2._DP * pi / cellVol &
                        * REAL(m123G,KIND=DP)**2 

#ifdef __USE_PARALLEL
  ! Sum contributions across processors using m_squared as temporary storage
  m_squared = EwaldRecipSpaceSpline
  CALL MPI_ALLREDUCE(m_squared, EwaldRecipSpaceSpline, 1, MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN EwaldRecipSpaceSpline***"
#endif

  DEALLOCATE(tempRA)
  DEALLOCATE(tmpRecip)

  RETURN

END FUNCTION EwaldRecipSpaceSpline


FUNCTION EwaldAveragePotentialTerm(eta)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   To be honest, I don't completely understand this term.  If anyone does, 
!   please feel free to insert your own explanation.
!
!   It is included in both the old OFDFT code and CASSTEP.  Basically, the
!   vague handwavy excuses they give is that its a term to "make the average
!   potential zero" and to deal w/ the "uniform charge-neutralizing background
!   and boundary conditions".  The best lead I can give to understand this term
!   is in the reference below.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/26/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: eta  
  ! or alpha, the Ewald convergence parameter that determines 
  ! relative rates of convergence between real and recip sums
  !
  REAL(KIND=DP) :: EwaldAveragePotentialTerm 
  ! the answer
  !

                    !>> INTERNAL VARIABLES <<!
                    !>> INITIALIZATION <<!
                    !>> FUNCTION BODY <<!

  EwaldAveragePotentialTerm = - 0.5_DP * pi / eta**2 / cellVol * sumCharge**2

  RETURN

END FUNCTION EwaldAveragePotentialTerm


END FUNCTION EwaldEnergy


FUNCTION EwaldForces(cellReal, ionTable, elementTable, useSpline)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function computes the ion-ion forces.  It outputs an array of forces
!   for each ion in each of the three lattice directions.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                           !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellReal
  ! The real-space lattice vectors
  !
  TYPE(ion), DIMENSION(:), INTENT(IN) :: ionTable
  ! Type + location of each ion, sorted by ion type
  ! We can get the number of ions (formerly numIon) by 
  ! the size of this array.
  !
  TYPE(element), DIMENSION(:), INTENT(IN):: elementTable
  ! index of the FIRST ion of each type in iontable
  ! This array actually has one more element than
  ! numIonType, because the last entry here numIon+1
  ! to make it easy to use this to loop through stuff
  !
  LOGICAL, INTENT(IN) :: useSpline
  ! Whether to use particle mesh ewald
  !
  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3) :: EwaldForces
  ! The answer!! (Ewald Forces)
  !
       
                        !>> INTERNAL VARIABLES <<!

#ifdef __USE_PARALLEL
  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3) :: tempLocForces 
  ! Temporary force storage used in MPI summations
  !
#endif

  INTEGER :: i
  ! A counter
  !
  REAL(KIND=DP), DIMENSION(2) :: ewaldTimes
  ! Array holding the execution time of each procedure
  !
  REAL(KIND=DP) :: realCutoffSq
  ! realCutoff**2
  !
  TYPE(stopwatch):: watch
  ! Timer
  !

                        !>> INITIALIZATION <<!

  CALL StartClock('EwaldForces')

  WRITE(outputUnit,*) " "
  WRITE(outputUnit,'(A)', ADVANCE="NO") " Calculating Ion-Ion Forces ..." 
  flush(outputUnit)
  watch=TimerStart()

  EwaldForces = 0._DP
  realCutoffSq = realCutoff**2

                         ! >> FUNCTION BODY <<!

#ifdef __USE_PARALLEL
  DO i = numIonInit, numIonInit+numIonLoc-1
#else
  DO i = 1, numIon
#endif
    EwaldForces(i,:)=EwaldRealSpaceForces( & 
        ionTable(i)%coord(:),elementTable(ionTable(i)%elementID)%charge,eta)
  END DO

#ifdef __USE_PARALLEL
  tempLocForces = EwaldForces
  CALL MPI_ALLREDUCE(tempLocForces, EwaldForces, SIZE(tempLocForces), &
                     MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN EwaldForces***"
#endif
  ewaldTimes(1) = TimerStop(watch)

  watch = TimerStart()
  IF (useSpline .EQV. .TRUE.) THEN
    !EwaldForces = EwaldForces + EwaldRecipSpaceForcesSpline(eta, recipCutoff)
    CALL EwaldRecipSpaceForcesSpline(eta, recipCutoff, EwaldForces)
  ELSE
    ! this sentence can not be compiled under PGI,
    ! so I change the function to subroutine.
    !EwaldForces = EwaldForces + EwaldRecipSpaceForces(eta, recipCutoff)
    CALL EwaldRecipSpaceForces(eta, recipCutoff, EwaldForces)
  END IF
  ewaldTimes(2) = TimerStop(watch)

!  WRITE(message,'(A, F7.1, F7.1)') "Ewald force times:", ewaldTimes
!  CALL WrtOut(6,message)

  WRITE(message,'(A, F10.3, A)') "                complete.", SUM(ewaldTimes), " s"
  CALL WrtOut(outputUnit,message)

  CALL StopClock('EwaldForces')

  RETURN

CONTAINS


FUNCTION EwaldRealSpaceForces(iCoord, iCharge, eta)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This is very simliar to EwaldEnergy version of EwaldRealSpaceEnergy
!   except we're implementing a different formula.  We only need to go through
!   half the ions because the forces for the other half are obtained by 
!   symmetry.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   I optimized for memory here insetad of speed.  Instead of calculating all
!   pairs of interactions I could STORE each pair of interaction and then
!   add them all up, thereby performing only 1/2 as many calculations.
!   However, this corresponds to a significant expense in regards to memory.
!   this method can be easily added, however.
!
!   The way this is implemented is very wasteful for non-orthorhombic cells.
!   I can't figure out a better way to do it and still get all the cells,
!   however, but if anyone else has a way to get all the points within a
!   sphere of radius realCutoff around an ion, please implement it.
!
!   Also, another thing I might be able to do is to add the coordinates
!   of the lattice vector after each loop, so that we can eliminate the mult.
!   needed to compute the difference in icoord and jcoord.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  USE Mathfunctions, only: Metric3d

  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(3), INTENT(IN) :: iCoord
  ! The coordinates of reference ion, in fractional coordinates
  !
  REAL(KIND=DP), INTENT(IN) :: eta
  ! The Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums 
  !
  REAL(KIND=DP), INTENT(IN) :: iCharge   
  ! The charge of the reference ion
  !
  REAL(KIND=DP), DIMENSION(3) :: EwaldRealSpaceForces 
  ! the answer
  !
             
                    !>> INTERNAL VARIABLES <<!
  INTEGER :: i_a, i_b, i_c
  ! counter, used to parse over all reps of 1st,2nd,3rd lattice vector
  !
  INTEGER :: j
  ! counter, used to parse over all ions
  !
  REAL(KIND=DP) :: r1, r2, r3
  !
  REAL(KIND=DP) :: dist_vec(3)
  ! Gives the distance between ion i and ion j+i_a+i_b+i_c
  !
  REAL(KIND=DP) :: dist_ijSq
  ! The distance of ion i from ion j
  !
  REAL(KIND=DP) :: etaDist
  ! eta * dist_ij
  !
  INTEGER ::  nr, newr
  ! 
  REAL(KIND=DP) :: oldfr1(3), oldfr2(3), fr1(3), fr2(3), fr(3)
  !
  REAL(KIND=DP) :: rmet(3,3), gmet(3,3)
  !
                           !>> INITIALIZATION <<!

  ! compute the metric tensor for realspace :rmet
  ! and for g-space : gmet
  call Metric3d(cellReal,rmet,gmet)

  EwaldRealSpaceForces = 0.0_DP

                           !>> FUNCTION BODY <<!
  ! Now we can loop through this shape for each ion in the cell
  ! and be sure we get all the cells in the sphere. 
  ! This loops over all the ions in the system
  DO j = 1, numIon

    ! Determine the minimum fractional distance between ions,
    ! taking images into account, so that fractional distance in each
    ! dimension has absolute value less than 0.5
    fr(:) = ionTable(i)%coord(:)-ionTable(j)%coord(:)
    fr = fr - ANINT(fr)

    ! This loops over the geometrical shape that encloses the sphere   
    nr = 1
    DO ! loop over shells
      newr = 0

      ! This loops over the geometrical shape that encloses the sphere   
      DO i_a = -nr,nr
        DO i_b = -nr,nr
          DO i_c = -nr,nr

             ! we only loop over current shell, not inside shell
             IF (ABS(i_a)==nr .OR. ABS(i_b)==nr .OR. ABS(i_c)==nr .OR. &
                 nr==1) THEN

               r1 = REAL(i_a,KIND=DP) + fr(1)
               r2 = REAL(i_b,KIND=DP) + fr(2)
               r3 = REAL(i_c,KIND=DP) + fr(3)
               dist_ijSq = rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3+&
                      2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r3*r1)

               ! (1) Note: ERFC(8) is about 1.1e-29,  charge_ij is 
               !     on the order of 1e1, dist_ij 's lower bound is 1e-1
               !     bohr.  We assume two ions cannot get close to within
               !     1e-1 bohr. Therefore we can neglect any contribution
               !     from eta*dist_ij > 8.
               !     Also exp(-64) is about 1.6e-28
               ! (2) dist_ij cannot = 0, to avoid zero in denominator
               !
               ! RECALL that dist_ij is actually storing the square of
               ! the distance between ions at this point.

               IF (dist_ijSq <= realCutoffSq .AND. dist_ijSq>1E-40) THEN
                 newr = 1
                 dist_vec(1) = cellReal(1,1)*r1+cellReal(1,2)*r2+&
                               cellReal(1,3)*r3
                 dist_vec(2) = cellReal(2,1)*r1+cellReal(2,2)*r2+&
                               cellReal(2,3)*r3
                 dist_vec(3) = cellReal(3,1)*r1+cellReal(3,2)*r2+&
                               cellReal(3,3)*r3
                 etaDist = eta * SQRT(dist_ijSq)
                 EwaldRealSpaceForces = EwaldRealSpaceForces &
                   + dist_vec * elementTable(ionTable(j)%elementID)%charge &
                              * (ERFC(etaDist) / etaDist**3 &
                   + 2._DP / sqrtPi * EXP(-(etaDist**2)) / etaDist**2)
               ENDIF
             ENDIF

          END DO ! i_c
        END DO ! i_b
      END DO ! i_a

      IF (newr == 0) EXIT

      nr = nr + 1

    ENDDO ! nr

  END DO ! j

  ! Multiply by the charge on the ion
  EwaldRealSpaceForces = EwaldRealSpaceForces * iCharge * eta**3

!  print *, 'realspace forces: ', EwaldRealSpaceForces

  RETURN

END FUNCTION EwaldRealSpaceForces


SUBROUTINE EwaldRecipSpaceForces(eta, recipCutoff, EwaldForces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This program implements equation 14 of [2], using the units of [3].  
!   There's not much to say about it, except that we have the same issue of 
!   having to get a shape to enclose the sphere made by recipCutoff.  We
!   solve this problem the SAME WAY as above in EWALDREALSPACE, except we don't
!   add 1 at the end (because we don't have the displacement problem anymore).
!
!   The only thing we do differently here is that we take advantage of symmetry
!   in reciprocal space.  (we can't do this in real space because the ions are
!   displaced relative to on another).  To exploit inversion symmetry, we will
!   only sum over half the cell.  i_a are the number of reps of the first 
!   recip lattice vector, i_b are the number of reps of the second and i_c
!   are the number of reps of the third:   
!
!   i_a ranges over 0 TO na ONLY.                         
!   i_b ranges over 0 TO nb when i_a=0 and over -nb to nb OTHERWISE. 
!   i_c ranges over 1 TO nc when i_a = i_b = 0 and over -nc to nc OTHERWISE.
! 
!   ... and  we will multiply by two at the end.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Again, this isn't all that efficient for non-orthonormal cells because it
!   is forced to sample more than it potentially has to.  Some smart guy can
!   come around and optimize it.
!   We might also find a way to store the structure factor.
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: eta
  ! eta, the Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums
  !
  REAL(KIND=DP), INTENT(IN) :: recipCutoff    
  ! the optimum reciprocal-space cutoff radius 
  !
  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3), INTENT(INOUT) :: EwaldForces 
  ! the answer
  !

                        !>> INTENRAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3) :: tmpForces 
  ! used for tmp storage 
  !
  INTEGER :: i_a, i_b, i_c
  ! counter, goes over reps of first recip. latt. vector
  !
  INTEGER :: i_b_lowLimit
  ! the number of reps the i_b loop should start from.
  !  Its either 0 or -nb, depending on the value of i_a
  !
  INTEGER :: i_c_lowLimit
  ! the number of repsthe i_c loop should start from.  
  !  Its either 0 or -nc, depending on the value of i_a and
  !  i_b. (see above blurb)
  !
  INTEGER ::  k
  ! counter to loop over all ions    
  !
  REAL(KIND=DP) :: structureFactorSin
  ! The SIN part of structure factor 
  ! (see eq. 11 of [3]).
  !
  REAL(KIND=DP) ::  structureFactorCos
  ! The COS part of structure factor
  !
  REAL(KIND=DP) ::  m_squared
  ! The dot product of m and m, where m is a 
  ! reciprocal lattice vector.
  !
  REAL(KIND=DP) ::  dotProd
  ! Temporary variable to store dot products
  !
  REAL(KIND=DP) ::  inverseNeg4eta2
  !
  REAL(KIND=DP), DIMENSION(3) :: m
  ! A 3-D reciprocal lattice vector
  !
  REAL(KIND=DP), DIMENSION(3) :: realCoord 
  ! The real-space coordinates of an ion, in Cartesian coordinates
  !

                        !>> INITIALIZATION <<!
  ! not used anymore
  ! EwaldRecipSpaceForces = 0.0_DP
  tmpForces = 0.0_DP
  inverseNeg4eta2 = -.25_DP/eta**2

                         !>> FUNCTION BODY <<!
  ! This loops over the geometrical shape that encloses the sphere in
  ! reciprocal space.  We will try to exploit the symmetry of the reciprocal
  ! space lattice (see above blurb).
  ! This loops over the first recip. lattice vector
  DO i_a = 0, nArecip
    ! Sets the lower limit on i_b depending on value of i_a
    IF (i_a == 0) THEN
      i_b_lowLimit = 0
    ELSE
      i_b_lowLimit = -nBrecip
    ENDIF

    ! This loops over second recip. lattice vector
    DO i_b = i_b_lowLimit, nBrecip
      ! Sets the lower limit on i_c depending on value of i_a and i_b
      IF (i_a == 0 .AND. i_b == 0) THEN
        i_c_lowLimit = 1
      ELSE
        i_c_lowLimit = -nCrecip
      ENDIF

      ! This loops over the third recip. lattice vector
      DO i_c = i_c_lowLimit, nCrecip
      
        ! compute the reciprocal lattice vector m we're currently at
        m = i_a * cellRecip(:,1) + i_b * cellRecip(:,2) + i_c * cellRecip(:,3)

        ! If the length of m is within the cutoff, then add its contribution.
        IF (Norm(m) <= recipCutoff) THEN          

          ! This loops over all the ions, computes sturcture factor
          structureFactorSin = 0.0_DP
          structureFactorCos = 0.0_DP              

          DO k = 1, numIon
            realCoord = Vecmul(cellReal, ionTable(k)%coord)
            dotProd = DOT_PRODUCT(m, realCoord)
            structureFactorSin = structureFactorSin &
              + elementTable(ionTable(k)%elementID)%charge * SIN(dotProd)
            structureFactorCos = structureFactorCos &
              + elementTable(ionTable(k)%elementID)%charge * COS(dotProd)
          END DO

          DO k = 1, numIon
            realCoord = Vecmul(cellReal, ionTable(k)%coord)
            dotProd = DOT_PRODUCT(m, realCoord)
            m_squared = DOT_PRODUCT(m,m)
            tmpForces(k,:) = tmpForces(k,:) &
                 + elementTable(ionTable(k)%elementID)%charge * &
                  EXP(m_squared * inverseNeg4eta2) / m_squared * &
                  m * ( SIN(dotProd) * structureFactorCos &
                       -COS(dotProd) * structureFactorSin)
          END DO

        ENDIF

      END DO ! i_c
    END DO ! i_b
  END DO ! i_a

  ! Multiply by 4*pi, which is the constant from eq. 13 of [3], and also by an
  ! extra factor of 2 to account for our symmetry exploitation. 
  tmpForces = tmpForces * 8.0_DP * pi / cellVol

  ! mohan add
  EwaldForces = EwaldForces + tmpForces

  RETURN

END SUBROUTINE EwaldRecipSpaceForces


SUBROUTINE EwaldRecipSpaceForcesSpline(eta, recipCutoff, EwaldForces)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This program implements equation 4.9 of [5], calculating forces using
!   particle-mesh Ewald.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/09/2007  Function Created (Linda Hung)
!
!------------------------------------------------------------------------------
  USE CBSpline, ONLY : splineOrder
  USE CBSpline, ONLY : BSplineProduct
  USE IonElectronSpline, ONLY: FillQIonTable
  USE IonElectronSpline, ONLY: CalculateSplineForces

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP), INTENT(IN) :: eta
  ! eta, the Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums
  !
  REAL(KIND=DP), INTENT(IN) :: recipCutoff    
  ! the optimum reciprocal-space cutoff radius 
  !
  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3), INTENT(INOUT) :: EwaldForces 
  ! the answer
  !
                        !>> INTENRAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(SIZE(ionTable),3) :: tmpForces 
  ! the answer
  !
  INTEGER :: allocateStatus
  ! flag in case of allocation error
  !
  INTEGER :: i_a, i_b, i_c
  ! counter, goes over reps of 1st, 2nd, 3rd recip. latt. vector
  !
  INTEGER :: k
  ! counter to loop over all ions    
  !
  REAL(KIND=DP) :: m_squared
  ! The dot product of m and m, where m is a 
  ! reciprocal lattice vector.
  !
  INTEGER, DIMENSION(3) :: m
  ! A 3-D reciprocal lattice index
  !
  REAL(KIND=DP), DIMENSION(3) :: mCell
  ! A 3-D reciprocal lattice vector
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmpRecip
  !

                        !>> INITIALIZATION <<!

  ! First use temp array to store qIonTable information for each ion type
  ALLOCATE(tempRA(n1G, n2G ,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating temporary array for particle mesh Ewald"
    STOP
  END IF

  ALLOCATE(qIonTable(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating q convolution table for particle mesh Ewald"
    STOP
  END IF

  ALLOCATE(CCStructure(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating CCStructure for particle mesh Ewald"
  STOP
  END IF

  ALLOCATE(tmpRecip(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating tmpRecip for particle mesh Ewald"
  STOP
  END IF

  !EwaldRecipSpaceForcesSpline = 0.0_DP
  tmpForces = 0.0_DP

                         !>> FUNCTION BODY <<!

  ! For each ion type, calculate qIonTable
  qIonTable = 0._DP
  DO k = 1, SIZE(elementTable) - 1
    CALL FillQIonTable(tempRA, ionTable(elementTable(k)%firstIonID: &
                                   (elementTable(k+1)%firstIonID-1)), &
                       splineOrder, nArecipSpl, nBrecipSpl, nCrecipSpl, &
                       n3GOff)

    qIonTable = qIonTable + tempRA*elementTable(k)%charge
  END DO

  ! FFT for complex conjugate of structure factor
  ! Fourier transform adds a normalization factor that is removed at end of
  ! this subroutine
  CALL FFT_NEW(FFT_STD_STATE, qIonTable, tmpRecip)
  CALL BSplineProduct(CONJG(tmpRecip),CCStructure, .TRUE.)

  DEALLOCATE(qIonTable)
  DEALLOCATE(tempRA)

  ! Now multiply contributions from the m-grid
#ifdef __USE_PARALLEL
  DO i_c= 1, k3G
    m(3) = i_c + k3GOff - 1
    IF (m(3) > m3G/2) m(3) = m(3) - m3G
    DO i_b = 1, k2G
      m(2) = i_b - 1 
      IF (m(2) > m2G/2) m(2) = m(2) - m2G
#else
  DO i_c = 1, nCrecipSpl
    m(3) = i_c - 1
    IF (m(3) > nCrecipSpl/2) m(3) = m(3) - nCrecipSpl
    DO i_b = 1, nBrecipSpl
      m(2) = i_b - 1
      IF (m(2) > nBrecipSpl/2) m(2) = m(2) - nBrecipSpl
#endif
      DO i_a = 1, nArecipSpl/2+1
        m(1) = i_a - 1
        IF (ALL(m==0)) THEN
          CCStructure(i_a,i_b,i_c) = (0._DP,0._DP)
        ELSE
          ! compute the reciprocal lattice vector m we're currently at
          mCell = Vecmul(cellRecip, REAL(m,KIND=DP))
          m_squared = DOT_PRODUCT(mCell,mCell)

          ! If the length of m is within the cutoff, then add its contribution.
          IF (m_squared <= recipCutoff**2) THEN
            CCStructure(i_a,i_b,i_c) = CCStructure(i_a,i_b,i_c) &
                                     * EXP(-m_squared/eta**2/4._DP)/m_squared
          ELSE
            CCStructure(i_a,i_b,i_c) = (0._DP,0._DP)
          END IF
        END IF
      END DO !i_a
    END DO !i_b (or i_c for parallel)
  END DO !i_c (or i_b for parallel)

  ! FFT to return to real space
#ifdef __USE_PARALLEL
  ALLOCATE(tempRA(nArecipSpl,nBrecipSpl,n3G), stat=allocateStatus)
#else
  ALLOCATE(tempRA(nArecipSpl,nBrecipSpl,nCrecipSpl), stat=allocateStatus)
#endif
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating temporary array for particle mesh Ewald"
    STOP
  END IF

  CALL FFT_NEW(FFT_STD_STATE, CCStructure, tempRA)

  ! Finally, add contributions to the derivative Q array
  DO k=1,SIZE(elementTable)-1
    CALL CalculateSplineForces( &
                 tmpForces(elementTable(k)%firstIonID:&
                                           elementTable(k+1)%firstIonID-1,:), &
                 ionTable(elementTable(k)%firstIonID: &
                          elementTable(k+1)%firstIonID-1), &
                 tempRA*elementTable(k)%charge, &
                 splineOrder, cellRecip/2._DP/pi, &
                 n3GOff, nCrecipSpl)
  END DO

  DEALLOCATE(tempRA)
  DEALLOCATE(CCStructure)
  DEALLOCATE(tmpRecip)


  ! Multiply by 4*pi, which is the constant from eq. 13 of [3], and also by an
  ! extra factor of 2 to account for our symmetry exploitation. 
  tmpForces = tmpForces &
                          * 4._DP * pi / cellVol &
                          * REAL(nArecipSpl,KIND=DP) &
                          * REAL(nBrecipSpl,KIND=DP) &
                          * REAL(nCrecipSpl,KIND=DP)

  ! For parallel code, ion-ion interactions need to be summed across processors
#ifdef __USE_PARALLEL
  tempLocForces = tmpForces
  CALL MPI_ALLREDUCE(tempLocForces, tmpForces, &
                     SIZE(tempLocForces), MPI_REAL8, MPI_SUM, &
                     MPI_COMM_WORLD, mpiErr)
  IF (mpiErr /= MPI_SUCCESS) STOP &
    "**PROBLEMS WITH MPI_ALLREDUCE IN CALCULATEFORCES***"
#endif

  EwaldForces = EwaldForces + tmpForces

  RETURN

END SUBROUTINE EwaldRecipSpaceForcesSpline


END FUNCTION EwaldForces


FUNCTION EwaldStress(cellReal, ionTable, elementTable, useSpline)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   This function calculates the Ewald stress
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   02/05/2004  Function created.  (Vincent Ligneres)
!
!------------------------------------------------------------------------------

  USE Constants, ONLY: auToGPa

  IMPLICIT NONE

                      !>> EXTERNAL VARIABLES <<!

  REAL(KIND=DP),DIMENSION(3,3) :: EwaldStress
  ! The answer
  !
  REAL(KIND=DP), DIMENSION(3,3), INTENT(IN) :: cellReal 
  ! The real-space lattice vectors
  !
  TYPE(ion), DIMENSION(:), INTENT(IN) :: ionTable
  ! Type + location of each ion, sorted by ion type
  ! We can get the number of ions (formerly numIon) by 
  ! the size of this array.
  !
  TYPE(element), DIMENSION(:), INTENT(IN):: elementTable
  ! index of the FIRST ion of each type in iontable
  ! This array actually has one more element than
  ! numIonType, because the last entry here numIon+1
  ! to make it easy to use this to loop through stuff
  !
  LOGICAL, INTENT(IN) :: useSpline
  ! Whether to use particle mesh ewald
  !

                        !>> INTERNAL VARIABLES <<!

  REAL(KIND=DP), DIMENSION(2) :: ewaldTimes
  ! Array holding the execution time of each procedure
  !
  TYPE(stopwatch):: watch 
  ! Timer
  !
  
                           !>> INITIALIZATION <<!

  CALL StartClock("EwaldStress")

  WRITE(outputUnit,*) " "
  WRITE(outputUnit,'(A)', ADVANCE="NO") " Calculating Ion-Ion Stress ..." 
  flush(outputUnit)

                           !>> FUNCTION BODY <<!

  watch = TimerStart()
  EwaldStress = EwaldRealSpaceStress(eta) + &
                EwaldAveragePotentialTermStress(eta)
  ewaldTimes(1) = TimerStop(watch)

  watch = TimerStart()
  IF (useSpline .EQV. .TRUE.) THEN
    EwaldStress = EwaldStress + &
                  EwaldRecipSpaceStressSpline(eta, recipCutoff)
  ELSE
    EwaldStress = EwaldStress + &
                  EwaldRecipSpaceStress(eta, recipCutoff)
  END IF
  ewaldTimes(2) = TimerStop(watch)

  WRITE(message,'(A, F10.3, A)') "                complete.", SUM(ewaldTimes), " s"
  CALL WrtOut(outputUnit,message)

  !WRITE(*,*) " stress for Ewald (Gpa): "
  !WRITE(*,'(3F12.6)') EwaldStress(1,:)*auToGPa
  !WRITE(*,'(3F12.6)') EwaldStress(2,:)*auToGPa
  !WRITE(*,'(3F12.6)') EwaldStress(3,:)*auToGPa

  CALL StopClock("EwaldStress")

  RETURN

CONTAINS


FUNCTION EwaldRealSpaceStress(eta)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   The real-space contribution to the Ewald stress
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  USE Mathfunctions, only: Metric3d

  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: eta
  ! The Ewald convergence parameter that determines relative 
  !
  REAL(KIND=DP), DIMENSION(3,3) :: EwaldRealSpaceStress 
  ! the answer
  !

                    !>> INTERNAL VARIABLES <<!

#ifdef __USE_PARALLEL
  REAL(KIND=dp), DIMENSION(3,3) :: tempLocStress 
  ! Local buffer for stress
  !
#endif

  INTEGER :: a, b
  ! The two indices that define the stress component
  !
  INTEGER :: i_type, j_type
  ! counter, used to parse over all types of ions
  !
  INTEGER :: i, j
  ! counter, used to parse over all ions in the system
  !
  INTEGER :: i_a, i_b, i_c
  ! counter, used to parse over all reps of 1st,2nd,3rd lattice vector
  !
  REAL(KIND=DP) :: r1, r2, r3
  !
  REAL(KIND=DP) :: charge_ij
  ! The charge of ion i times the charge of ion j
  !
  REAL(KIND=DP) :: dist_vec(3)
  !
  !
  REAL(KIND=DP) :: dist_ij
  ! The distance of ion i from ion j
  !
  REAL(KIND=DP) :: temp
  ! Temporary variable for calculational efficiency
  !
  REAL(KIND=DP) :: realCutoffSq
  !
  INTEGER ::  nr, newr
  ! 
  REAL(KIND=DP) :: fr(3)
  !
  REAL(KIND=DP) :: rmet(3,3), gmet(3,3)
  !

                           !>> INITIALIZATION <<!

  ! compute the metric tensor for realspace :rmet
  ! and for g-space : gmet
  call Metric3d(cellReal,rmet,gmet)
  EwaldRealSpaceStress = 0.0_DP
  realCutoffSq=realCutoff**2

                           !>> FUNCTION BODY <<!

  ! Now we can loop through this shape for each ion in the cell
  ! and be sure we get all the cells in the sphere. 
  DO i_type = 1, numIonType
#ifdef __USE_PARALLEL
    ! Skip this ion type if relevant ion calculations are not on this processor
    IF ((elementTable(i_type+1)%firstIonID-1<numIonInit) .OR. &
        (elementTable(i_type)%firstIonID>numIonInit+numIonLoc-1)) CYCLE
#endif
    DO j_type = 1, numIonType
    ! Calculates the product of the charges
    charge_ij = elementTable(i_type)%charge * elementTable(j_type)%charge
#ifdef __USE_PARALLEL
      DO i = MAX(elementTable(i_type)%firstIonID,numIonInit), &
             MIN(elementTable(i_type+1)%firstIonID-1, numIonInit+numIonLoc-1)
#else
      DO i = elementTable(i_type)%firstIonID, &
             elementTable(i_type+1)%firstIonID-1
#endif
        DO j = elementTable(j_type)%firstIonID, &
               elementTable(j_type+1)%firstIonID-1

          ! Determine the minimum fractional distance between ions,
          ! taking images into account, so that fractional distance in each
          ! dimension has absolute value less than 0.5
          fr(:) = ionTable(i)%coord(:)-ionTable(j)%coord(:)
          fr = fr - ANINT(fr)

          ! This loops over the geometrical shape that encloses the sphere   
          nr = 1
          DO ! loop over shells
             newr = 0

             ! This loops over the geometrical shape that encloses the sphere   
             DO i_a = -nr,nr
               DO i_b = -nr,nr
                 DO i_c = -nr,nr

                   ! we only loop over current shell, not inside shell
                   IF (ABS(i_a)==nr .OR. ABS(i_b)==nr .OR. ABS(i_c)==nr .OR. &
                       nr==1) THEN

                     r1 = REAL(i_a,KIND=DP) + fr(1)
                     r2 = REAL(i_b,KIND=DP) + fr(2)
                     r3 = REAL(i_c,KIND=DP) + fr(3)

                     ! Distance between ions, squared
                     dist_ij = rmet(1,1)*r1**2+rmet(2,2)*r2**2+rmet(3,3)*r3**2+&
                       2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r3*r1)

                     ! (1) Note: ERFC(8) is about 1.1e-29,  charge_ij is 
                     !     on the order of 1e1, dist_ij 's lower bound is 1e-1
                     !     bohr.  We assume two ions cannot get close to within
                     !     1e-1 bohr. Therefore we can neglect any contribution
                     !     from eta*dist_ij > 8
                     ! (2) dist_ij cannot = 0, to avoid zero in denominator
                     ! RECALL that dist_ij is actually storing the square of
                     ! the distance between ions at this point.
                     IF (dist_ij <= realCutoffSq .AND. dist_ij>1E-40) THEN
                       newr = 1
                       dist_ij=SQRT(dist_ij)
                       temp = charge_ij &
                           * ( ERFC(eta * dist_ij) / dist_ij &
                               + 2._DP * eta * EXP( -(eta * dist_ij)**2) &
                               / sqrtPi &
                              ) / (dist_ij**2)
                       dist_vec(1) = cellReal(1,1)*r1+cellReal(1,2)*r2+&
                                     cellReal(1,3)*r3
                       dist_vec(2) = cellReal(2,1)*r1+cellReal(2,2)*r2+&
                                     cellReal(2,3)*r3
                       dist_vec(3) = cellReal(3,1)*r1+cellReal(3,2)*r2+&
                                     cellReal(3,3)*r3

                       DO a=1, 3
                         DO b = a, 3
                           EwaldRealSpaceStress(a,b) = &
                             EwaldRealSpaceStress(a,b) &
                             + temp*dist_vec(a)*dist_vec(b)
                         END DO ! a
                       END DO  ! b
                     ENDIF ! Within cutoff, not identical

                   ENDIF ! On shell, not inside of sphere
                 END DO ! i_c
               END DO ! i_b
             END DO ! i_a

           IF (newr==0) EXIT

           nr = nr + 1

           ENDDO  ! nr

        END DO ! j
      END DO ! i
    END DO ! j_type
  END DO ! i_type
  
  EwaldRealSpaceStress = - EwaldRealSpaceStress / 2._DP / cellVol

#ifdef __USE_PARALLEL
  tempLocStress = EwaldRealSpaceStress
  CALL MPI_ALLREDUCE(tempLocStress, EwaldRealSpaceStress, 9, MPI_REAL8, &
                     MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN EwaldRealSpaceStress***"
#endif

!  print *, 'Ewald Real stress = ', EwaldRealSpaceStress

END FUNCTION EwaldRealSpaceStress


FUNCTION EwaldRecipSpaceStress(eta, recipCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   The reciprocal-space contribution to the Ewald stress
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   10/31/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: eta
  ! eta, the Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums
  !
  REAL(KIND=DP), INTENT(IN) :: recipCutoff
  ! the optimum reciprocal-space cutoff radius 
  !
  REAL(KIND=DP), DIMENSION(3,3) :: EwaldRecipSpaceStress 
  ! the answer
  !

                        !>> INTENRAL VARIABLES <<!
  INTEGER :: a, b
  ! The two indices that define the stress component
  !
  INTEGER :: i_a, i_b, i_c
  ! counter, goes over reps of 1st,2nd,3rd recip. latt. vector
  !
  INTEGER :: i_b_lowLimit
  ! the number of reps the i_b loop should start from.
  !  Its either 0 or -nb, depending on the value of i_a
  !
  INTEGER :: i_c_lowLimit
  ! the number of repsthe i_c loop should start from.  
  ! Its either 0 or -nc, depending on the value of i_a and
  ! i_b. (see above blurb)
  INTEGER :: k
  ! counter to loop over all ions    
  !
  REAL(KIND=DP) :: structureFactorSin
  ! The SIN part of structure factor 
  ! (see eq. 11 of [3]).
  !
  REAL(KIND=DP) :: structureFactorCos
  ! The COS part of structure factor
  !
  REAL(KIND=DP) :: m_squared
  ! The dot product of m and m, where m is a 
  ! reciprocal lattice vector.
  !
  REAL(KIND=DP) :: temp1, temp2
  ! Temporary variables
  !
  REAL(KIND=DP), DIMENSION(3) :: m
  ! A 3-D reciprocal lattice vector
  !
  REAL(KIND=DP), DIMENSION(3) :: realCoord 
  ! The real-space coordinates of an ion, in Cartesian coordinates
  !

                        !>> INITIALIZATION <<!
  EwaldRecipSpaceStress = 0.0_DP

                         !>> FUNCTION BODY <<!
  ! This loops over the geometrical shape that encloses the sphere in
  ! reciprocal space.  We will try to exploit the symmetry of the reciprocal
  ! space lattice (see above blurb).
  ! This loops over the first recip. lattice vector
  DO i_a = 0, nArecip
    ! Sets the lower limit on i_b depending on value of i_a
    IF (i_a == 0) THEN
      i_b_lowLimit = 0
    ELSE
      i_b_lowLimit = -nBrecip
    ENDIF

    ! This loops over second recip. lattice vector
    DO i_b = i_b_lowLimit, nBrecip
      ! Sets the lower limit on i_c depending on value of i_a and i_b
      IF (i_a == 0 .AND. i_b == 0) THEN
        i_c_lowLimit = 1
      ELSE
        i_c_lowLimit = -nCrecip
      ENDIF

      ! This loops over the third recip. lattice vector
      DO i_c = i_c_lowLimit, nCrecip
      
        ! compute the reciprocal lattice vector m we're currently at
        m = i_a * cellRecip(:,1) + i_b * cellRecip(:,2) + i_c * cellRecip(:,3)

        ! If the length of m is within the cutoff, then add its contribution.
        IF (Norm(m) <= recipCutoff) THEN          

          ! This loops over all the ions, computes sturcture factor
          structureFactorSin = 0.0_DP
          structureFactorCos = 0.0_DP              
   
          DO k = 1, numIon
            realCoord = Vecmul(cellReal, ionTable(k)%coord)        
            structureFactorSin = structureFactorSin &
              + elementTable(ionTable(k)%elementID)%charge &
                * SIN(DOT_PRODUCT(m, realCoord))
            structureFactorCos = structureFactorCos &
              + elementTable(ionTable(k)%elementID)%charge &
                * COS(DOT_PRODUCT(m, realCoord)) 
          END DO

          ! Now put the formula together and add the contribution to the total
          m_squared = Dot_Product(m,m)

          temp1 = exp(-(m_squared / eta**2 / 4._DP)) / m_squared &
                  * (structureFactorSin**2 + structureFactorCos**2)

          temp2 = temp1 * 2._DP / m_squared &
                  * (1._DP + m_squared / eta**2 / 4._DP)

          DO a = 1, 3
            DO b = a, 3
              IF (a == b) THEN
                EwaldRecipSpaceStress(a,b) = EwaldRecipSpaceStress(a,b) + &
                  temp2 * m(a) * m(b) - temp1 
              ELSE
                EwaldRecipSpaceStress(a,b) = EwaldRecipSpaceStress(a,b) + &
                  temp2 * m(a) * m(b) 
              END IF              
            END DO
          END DO

        ENDIF

      END DO ! i_c
    END DO ! i_b
  END DO ! i_a
  
  ! Multiply by an extra factor of 2 to account for our symmetry exploitation.
  EwaldRecipSpaceStress = EwaldRecipSpaceStress * 4.0_DP * pi / cellVol**2

  RETURN
 
END FUNCTION EwaldRecipSpaceStress


FUNCTION EwaldRecipSpaceStressSpline(eta, recipCutoff)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   The reciprocal-space contribution to the Ewald stress with particle mesh
!   Ewald approximation.  (Analogous to PME reciprocal energy calculation.)
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/16/2007 Function Created (Linda Hung)
!
!------------------------------------------------------------------------------
  USE CBSpline, ONLY: splineOrder
  USE CBSpline, ONLY: BSplineProduct
  USE IonElectronSpline, ONLY : FillQIonTable

  IMPLICIT NONE

                        !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: eta
  ! eta, the Ewald convergence parameter that determines relative 
  ! rates of convergence between real and recip sums
  REAL(KIND=DP), INTENT(IN) :: recipCutoff     
  ! the optimum reciprocal-space cutoff radius 
  !
  REAL(KIND=DP), DIMENSION(3,3) :: EwaldRecipSpaceStressSpline 
  ! the answer
  !

                        !>> INTERNAL VARIABLES <<!
  INTEGER :: allocateStatus
  ! flag in case of allocation error
  !
  INTEGER :: a, b
  ! The two indices that define the stress component
  !
  INTEGER :: i_a, i_b, i_c
  ! counter, goes over reps of 1st,2nd,3rd recip. latt. vector
  !
  INTEGER :: k 
  ! counter to loop over all ion types
  !
  REAL(KIND=DP) :: m_squared
  ! The dot product of m and m, where m is a 
  ! reciprocal lattice vector.
  !
  REAL(KIND=DP) :: temp1, temp2
  ! Temporary variables
  !
  INTEGER, DIMENSION(3) :: m 
  ! A 3-D reciprocal lattice index
  !
  REAL(KIND=DP), DIMENSION(3) :: mCell
  ! A 3-D reciprocal lattice vector
  !
  COMPLEX(KIND=DP), DIMENSION(:,:,:), ALLOCATABLE :: tmpRecip
  !

#ifdef __USE_PARALLEL
  REAL(KIND=DP), DIMENSION(3,3) :: tempLocStress 
  ! Total stress from local processor - used with MPI
  !
#endif

                        !>> INITIALIZATION <<!

  ALLOCATE(tempRA(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating c array for particle mesh Ewald"
  STOP
  END IF

  ALLOCATE(qIonTable(n1G,n2G,n3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating q convolution table for particle mesh Ewald"
    STOP
  END IF

  ALLOCATE(CCStructure(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating CCStructure for particle mesh Ewald"
  STOP
  END IF

  ALLOCATE(tmpRecip(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating tmpRecip for particle mesh Ewald"
  STOP
  END IF


  EwaldRecipSpaceStressSpline = 0.0_DP

                         !>> FUNCTION BODY <<!

  ! For each ion type, calculate qIonTable
  qIonTable = 0._DP
  DO k = 1, SIZE(elementTable) - 1
    CALL FillQIonTable(tempRA, ionTable(elementTable(k)%firstIonID: &
                                        (elementTable(k+1)%firstIonID-1)), &
                       splineOrder, nArecipSpl, nBrecipSpl, nCrecipSpl, &
                       n3Goff)

    qIonTable = qIonTable + tempRA*elementTable(k)%charge
  END DO

  ! FFT for complex conjugate of structure factor
  ! Fourier transform adds a normalization factor that is removed at end of
  ! this subroutine
  CALL FFT_NEW(FFT_STD_STATE, qIonTable, tmpRecip)
  CALL BSplineProduct(tmpRecip, CCStructure)

  DEALLOCATE(qIonTable)
  DEALLOCATE(tempRA)
  ALLOCATE(tempRA(k1G,k2G,k3G), stat=allocateStatus)
  IF (allocateStatus/=0) THEN
    WRITE(*,*) "Error allocating c array for particle mesh Ewald"
    STOP
  END IF

  ! And now determine |S(m)|**2
  tempRA = CCStructure*CONJG(CCStructure)

  DEALLOCATE(CCStructure)

  ! Now sum up contributions from the m-grid
  DO i_c = 1, k3G 
    m(3) = i_c + k3GOff - 1
    IF (m(3) > m3G/2) m(3) = m(3) - m3G
    DO i_b = 1, k2G 
      m(2) = i_b - 1 
      IF (m(2) > m2G/2) m(2) = m(2) - m2G 
      DO i_a = 1, k1G 
        m(1) = i_a - 1
        IF (ANY(m/=0)) THEN
        ! compute the reciprocal lattice vector m we're currently at
          mCell = Vecmul(cellRecip, REAL(m,KIND=DP))
          m_squared = DOT_PRODUCT(mCell,mCell)

          ! If the length of m is within the cutoff, then add its contribution.
          IF (m_squared <= recipCutoff**2) THEN

            temp1 = EXP(-(m_squared / eta**2 / 4._DP)) / m_squared &
                    * tempRA(i_a,i_b,i_c)

            ! Now add in factor for symmetry in non-m(1)=0 plane
            IF (i_a/=1) temp1 = 2*temp1
                  

            temp2 = temp1 * 2._DP &
                    * (1._DP / m_squared + 1._DP / eta**2 / 4._DP)

            DO a = 1, 3
              DO b = a, 3
                IF (a == b) THEN
                  EwaldRecipSpaceStressSpline(a,b) = &
                    EwaldRecipSpaceStressSpline(a,b) + &
                    temp2 * mCell(a) * mCell(b) - temp1
                ELSE
                  EwaldRecipSpaceStressSpline(a,b) = &
                    EwaldRecipSpaceStressSpline(a,b) + &
                    temp2 * mCell(a) * mCell(b) 
                END IF 
              END DO
            END DO

          ENDIF

        END IF
      END DO
    END DO
  END DO
  
  ! Multiply for grid normalization factor (Lost during FFT)
  EwaldRecipSpaceStressSpline  = EwaldRecipSpaceStressSpline &
                               * 2._DP * pi / cellVol**2 & ! Is this right?
                               * REAL(m123G,KIND=DP)**2 

#ifdef __USE_PARALLEL
  ! Sum contributions across processors using m_squared as temporary storage
  tempLocStress = EwaldRecipSpaceStressSpline
  CALL MPI_ALLREDUCE(tempLocStress, EwaldRecipSpaceStressSpline, 9, MPI_REAL8, &
                     MPI_SUM, MPI_COMM_WORLD, mpiErr)
  IF (mpiErr/=MPI_SUCCESS) STOP &
    "***PROBLEMS WITH MPI_ALLREDUCE IN EwaldRecipSpaceStressSpline***"
#endif

  DEALLOCATE(tempRA)
  DEALLOCATE(tmpRecip)

  RETURN

END FUNCTION EwaldRecipSpaceStressSpline



FUNCTION EwaldAveragePotentialTermStress(eta)
!------------------------------------------------------------------------------
! DESCRIPTION:
!   See explanation for the Average Potential Energy Term above.  This is the
!   stress contribution from that.
!
! CONDITIONS AND ASSUMPTIONS:
! 
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/26/2003  Function Created (Greg Ho)
!
!------------------------------------------------------------------------------
  IMPLICIT NONE

                    !>> EXTERNAL VARIABLES <<!
  REAL(KIND=DP), INTENT(IN) :: eta
  ! or alpha, the Ewald convergence parameter that determines 
  ! relative rates of convergence between real and recip sums
  !
  REAL(KIND=DP), DIMENSION(3,3) :: EwaldAveragePotentialTermStress 
  ! the answer
  !

                    !>> INTERNAL VARIABLES <<!
  INTEGER :: a 
  ! The two indices that define the stress component we're looking at.
  !

                   !>> INITIALIZATION <<!

  EwaldAveragePotentialTermStress = 0._DP

                    !>> FUNCTION BODY <<!

  DO a=1,3
    EwaldAveragePotentialTermStress(a,a) = &
          0.5_DP * pi / eta**2 / cellVol**2 * sumCharge**2
  END DO

  RETURN

END FUNCTION EwaldAveragePotentialTermStress


END FUNCTION EwaldStress


END MODULE Ewald
