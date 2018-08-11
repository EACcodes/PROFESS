MODULE Constants
!------------------------------------------------------------------------------
! STRUCTURE OF MODULE:
!   MODULE Constants
!     |_FUNCTION AtomicMass
!
! DESCRIPTION:
!   Just a bunch of universal constants (pi, speed of light, conversion
!   factors, etc.)  They will never ever change. Supposedly.  Also, this
!   program contains a function to look up the atomic mass of an element.
!   
!   The constants here were chosen in general to try to conform with CASTEP.
!   Of course, CASTEP has the nasty habit of changing their own constants from
!   version to version, so this may have to be checked periodically.
!
! CONDITIONS AND ASSUMPTIONS:
!
! FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!   Add more elements to the AtomicMass table!
!    
! REFERENCES:
! [1] Edwin R. Williams, Richard L. Steiner, David B. Newell, and Paul T. Olsen
!     Physical Review Letters, v.81, p. 2404, 21 September 1998. (value of h.) 
! [2] The Fundamental Physical Constants, by E. Richard Cohen and 
!     Barry N. Taylor, Physics Today, August 1995, page BG9i
! [3] Moher and Taylor, Rev. Mod. Phys. 72 April 2000, p. 447
!------------------------------------------------------------------------------
! REVISION LOG:
!   11/20/2003  File Created (Greg Ho)
!   11/24/2003  Added more accuracy to mu0
!   12/04/2003  Changed Planck's constant to the ref. 2 value. (Vince Ligneres)
!   12/12/2003  Fixed everything to be "_DP" at the end.  This fixes a problem
!               in the vW term that makes the value fluxuate depending on the
!               value in bohr!
!   09/30/2004  Changed the value of mu0 from fixed constant to derived
!               Added a few more digits of pi (which is basically worthless)
!               changed qElec to conform with ref. [3]
!               Changed the value of bohr to conform with CASTEP II
!               Pretty much, these changes are just to get it to conformity 
!               with CASSTEP II
!------------------------------------------------------------------------------
                                 !>> GLOBAL <<!

  INTEGER, PARAMETER :: &
  DP = 8, &                ! Double Precision, ~14 significant digits
  systemNameLen = 80       ! The maximum number of characters we can have
                           ! in the system name.  This is determined by how
                           ! much we can fit into our output file! :)  

  REAL(kind=DP), PARAMETER :: & 
  ! Defined Constants
  tiny = 1.E-12_DP, &                   ! tiny
  pi = 3.141592653589793238462643_DP, & ! you know, pi.
  qElec = 1.602176462E-19_DP, & ! C.    ! Electron charge [3]
  light = 299792458._DP, &      ! m/s.  ! Speed of light in vacuum [2]
  bohr = 0.529177208607388_DP, &! Ang   ! Bohr length.
  mToA = 1.E10_DP, &                    ! Conversion for m to angstroms
  golden = 0.38196601125010515179541316563436_DP, & ! The Golden Ratio
  rydbergToeV = 13.60569193, &  ! eV    ! Rydberg (energy unit)
  emass = 9.1093826E-28, &      ! g     ! Electron mass
  mol = 6.02214179E23, &                ! Avagdro's number


  ! Derived Constants
  boltzmann = 3.16682E-6_dp, &                     ! in unit of [Ha/K] ! Boltzmann constant
  fundamentalTime = 2.418884326505E-17_DP, &       ! [second / a.u.]
  auPressure = 1._dp/2.9421912e13_dp  , &          ! [a.u./Pa], convert from Pa to a.u. of pressure
  mu0 = 4._DP*pi*1.E-7_DP, &                       ! N/A^2. ! Permeability of vacuum [2]
  eps0 = 1._DP/(mu0*light**2), &                   ! Permittivity of vacuum
                                                   ! Conv, factor from Hartree to eV
  hartreeToeV = qElec * mToA / (4._DP * pi * eps0 * bohr), &  

  ! Conversion from AU of pressure to GPa
  auToGPa = hartreeToeV / bohr**3 * qElec * mToA**3 / 1.E9_DP

  COMPLEX(KIND=DP), PARAMETER :: &
  imag = CMPLX(0.0_DP, 1.0_DP)               ! The square root of -1 

  REAL(KIND=DP) :: machPrec  ! the machine precision

CONTAINS


FUNCTION AtomicMass(symbol)
!------------------------------------------------------------------------------
!  DESCRIPTION:  
!    Gets the atomic mass of the atomic element given by the two character
!    symbol in atomic units.
!
!  GLOBAL/MODULE VARIABLES CHANGED:
!
!  CONDITIONS AND ASSUMPTIONS:
!
!  FUTURE OPTIMIZATIONS AND IMPROVEMENTS:
!    Eventually we want all the important elements used in OFDFT reprsented
!    in this table
!------------------------------------------------------------------------------
!  REVISION LOG:
!  09-27-12 Mohan add new function, the input name of element don't need
!  to have 3 characters, it can also has only one or three characters.
!------------------------------------------------------------------------------
  IMPLICIT NONE

                         !>> EXTERNAL VARIABLES <<!

  CHARACTER(LEN=3), INTENT(IN) :: symbol

  REAL(KIND=DP) :: AtomicMass

                          !>> INTERNAL VARIABLES <<!
                            !>> INITIALIZATION <<!
                             !>> FUNCTION BODY <<!

  SELECT CASE(TRIM(symbol)) 

    ! 1
    CASE("H")
      AtomicMass = 1.00794_DP/(mol*emass)

    ! 1 Deuterium
    CASE("D")
      AtomicMass = 2.014102_DP/(mol*emass)

    ! 1 Deuteron (Nucleus of deuterium)
!    CASE("D")
!      AtomicMass = 2.01355_DP/(mol*emass)

    ! 1 Tritium
    CASE("T")
      AtomicMass = 3.0160492_DP/(mol*emass)

    ! 2
    CASE("He")
      AtomicMass = 4.002602_DP/(mol*emass)

    ! 3
    CASE("Li")
      AtomicMass = 6.9412_DP/(mol*emass)

    ! 4
    CASE("Be")
      AtomicMass = 9.012182_DP/(mol*emass)

    ! 5
    CASE("B")
      AtomicMass = 10.811_DP/(mol*emass)

    ! 6
    CASE("C")
      AtomicMass = 12.0107_DP/(mol*emass)

    ! 7
    CASE("N")
      AtomicMass = 14.00674_DP/(mol*emass)

    ! 8
    CASE("O")
      AtomicMass = 15.9994_DP/(mol*emass)

    ! 9
    CASE("F")
      AtomicMass = 18.9984032_DP/(mol*emass)

    ! 10
    CASE("Ne")
      AtomicMass = 20.1797_DP/(mol*emass)

    ! 11
    CASE("Na")
      AtomicMass = 22.98976928_DP/(mol*emass)
      
    ! 12
    CASE("Mg")
      AtomicMass = 24.3050_DP/(mol*emass)
    
    ! 13
    CASE("Al")
      AtomicMass = 26.981539_DP/(mol*emass)
    
    ! 14
    CASE("Si")
      AtomicMass = 28.0855_dp/(mol*emass)
      
    ! 15
    CASE("P")
      AtomicMass = 30.973762_DP/(mol*emass)

    ! 16
    CASE("S")
      AtomicMass = 32.066_DP/(mol*emass)

    ! 17
    CASE("Cl")
      AtomicMass = 35.4527_DP/(mol*emass)

    ! 18
    CASE("Ar")
      AtomicMass = 39.948_DP/(mol*emass)
    
    ! 19
    CASE("K")
      AtomicMass = 39.0983_DP/(mol*emass)

    ! 20
    CASE("Ca")
      AtomicMass = 40.078_DP/(mol*emass)
 
    ! 21
    CASE("Sc")
      AtomicMass = 44.955910_DP/(mol*emass)
    
    ! 22
    CASE("Ti")
      AtomicMass = 47.867_DP/(mol*emass)

    ! 23
    CASE("V")
      AtomicMass = 50.9415_DP/(mol*emass)

    ! 24
    CASE("Cr")
      AtomicMass = 51.9961_DP/(mol*emass)

    ! 25
    CASE("Mn")
      AtomicMass = 54.938049_DP/(mol*emass)

    ! 26
    CASE("Fe")
      AtomicMass = 55.845_DP/(mol*emass)

    ! 27
    CASE("Co")
      AtomicMass = 58.933200_DP/(mol*emass)

    ! 28
    CASE("Ni")
      AtomicMass = 58.6934_DP/(mol*emass)

    ! 29
    CASE("Cu")
      AtomicMass = 63.546_DP/(mol*emass)

    ! 30
    CASE("Zn")
      AtomicMass = 65.39_DP/(mol*emass)

    ! 31
    CASE("Ga")
      AtomicMass = 69.723_DP/(mol*emass)

    ! 32
    CASE("Ge")
      AtomicMass = 72.61_DP/(mol*emass)

    ! 33
    CASE("As")
      AtomicMass = 74.92160_DP/(mol*emass)

    ! 34
    CASE("Se")
      AtomicMass = 78.96_DP/(mol*emass)

    ! 35
    CASE("Br")
      AtomicMass = 79.904_DP/(mol*emass)

    ! 36
    CASE("Kr")
      AtomicMass = 83.80_DP/(mol*emass)

    ! 37
    CASE("Rb")
      AtomicMass = 85.4678_DP/(mol*emass)

    ! 38
    CASE("Sr")
      AtomicMass = 87.62_DP/(mol*emass)

    ! 39
    CASE("Y")
      AtomicMass = 88.90585_DP/(mol*emass)

    ! 40
    CASE("Zr")
      AtomicMass = 91.224_DP/(mol*emass)

    ! 41
    CASE("Nb")
      AtomicMass = 92.90638_DP/(mol*emass)

    ! 42
    CASE("Mo")
      AtomicMass = 95.94_DP/(mol*emass)

    ! 43
    CASE("Tc")
      AtomicMass = 98_DP/(mol*emass)

    ! 44
    CASE("Ru")
      AtomicMass = 101.07_DP/(mol*emass)

    ! 45
    CASE("Rh")
      AtomicMass = 102.90550_DP/(mol*emass)

    ! 46
    CASE("Pd")
      AtomicMass = 106.42_DP/(mol*emass)

    ! 47
    CASE("Ag")
      AtomicMass = 107.8682_DP/(mol*emass)

    ! 48
    CASE("Cd")
      AtomicMass = 112.411_DP/(mol*emass)

    ! 49
    CASE("In")
      AtomicMass = 114.818_DP/(mol*emass)

    ! 50
    CASE("Sn")
      AtomicMass = 118.710_DP/(mol*emass)

    ! 51
    CASE("Sb")
      AtomicMass = 121.760_DP/(mol*emass)

    ! 52
    CASE("Te")
      AtomicMass = 127.60_DP/(mol*emass)

    ! 53
    CASE("I")
      AtomicMass = 126.90447_DP/(mol*emass)

    ! 54
    CASE("Xe")
      AtomicMass = 131.29_DP/(mol*emass)

    ! 55
    CASE("Cs")
      AtomicMass = 132.90545_DP/(mol*emass)

    ! 56
    CASE("Ba")
      AtomicMass = 137.327_DP/(mol*emass)

    ! 57
    CASE("La")
      AtomicMass = 138.9055_DP/(mol*emass)

    ! 58
    CASE("Ce")
      AtomicMass = 140.116_DP/(mol*emass)

    ! 59
    CASE("Pr")
      AtomicMass = 140.90765_DP/(mol*emass)

    ! 60
    CASE("Nd")
      AtomicMass = 144.24_DP/(mol*emass)

    ! 61
    CASE("Pm")
      AtomicMass = 145_DP/(mol*emass)

    ! 62
    CASE("Sm")
      AtomicMass = 150.36_DP/(mol*emass)

    ! 63
    CASE("Eu")
      AtomicMass = 151.964_DP/(mol*emass)

    ! 64
    CASE("Gd")
      AtomicMass = 157.25_DP/(mol*emass)

    ! 65
    CASE("Tb")
      AtomicMass = 158.92534_DP/(mol*emass)

    ! 66
    CASE("Dy")
      AtomicMass = 162.50_DP/(mol*emass)

    ! 67
    CASE("Ho")
      AtomicMass = 164.93032_DP/(mol*emass)

    ! 68
    CASE("Er")
      AtomicMass = 167.26_DP/(mol*emass)

    ! 69
    CASE("Tm")
      AtomicMass = 168.93421_DP/(mol*emass)

    ! 70
    CASE("Yb")
      AtomicMass = 173.04_DP/(mol*emass)

    ! 71
    CASE("Lu")
      AtomicMass = 174.967_DP/(mol*emass)

    ! 72
    CASE("Hf")
      AtomicMass = 178.49_DP/(mol*emass)

    ! 73
    CASE("Ta")
      AtomicMass = 180.9479_DP/(mol*emass)

    ! 74
    CASE("W")
      AtomicMass = 183.84_DP/(mol*emass)

    ! 75
    CASE("Re")
      AtomicMass = 186.207_DP/(mol*emass)

    ! 76
    CASE("Os")
      AtomicMass = 190.23_DP/(mol*emass)

    ! 77
    CASE("Ir")
      AtomicMass = 192.217_DP/(mol*emass)
      
    ! 78
    CASE("Pt")
      AtomicMass = 195.078_DP/(mol*emass)

    ! 79
    CASE("Au")
      AtomicMass = 196.96655_DP/(mol*emass)

    ! 80
    CASE("Hg")
      AtomicMass = 200.59_DP/(mol*emass)
    
    ! 81
    CASE("Tl")
      AtomicMass = 204.3833_DP/(mol*emass)

    ! 82
    CASE("Pb")
      AtomicMass = 207.2_DP/(mol*emass)

    ! 83
    CASE("Bi")
      AtomicMass = 208.98038_DP/(mol*emass)
      
    ! 84
    CASE("Po")
      AtomicMass = 209_DP/(mol*emass)

    ! 85
    CASE("At")
      AtomicMass = 210_DP/(mol*emass)

    ! 86
    CASE("Rn")
      AtomicMass = 222_DP/(mol*emass)

    ! 87
    CASE("Fr")
      AtomicMass = 223_DP/(mol*emass)

    ! 88
    CASE("Ra")
      AtomicMass = 226_DP/(mol*emass)
      
    ! 89
    CASE("Ac")
      AtomicMass = 227_DP/(mol*emass)

    ! 90
    CASE("Th")
      AtomicMass = 232.0381_DP/(mol*emass)

    ! 91
    CASE("Pa")
      AtomicMass = 231.03588_DP/(mol*emass)

    ! 92
    CASE("U")
      AtomicMass = 238.0289_DP/(mol*emass)

    ! 93
    CASE("Np")
      AtomicMass = 237_DP/(mol*emass)
      
    ! 94
    CASE("Pu")
      AtomicMass = 244_DP/(mol*emass)

    ! 95
    CASE("Am")
      AtomicMass = 243_DP/(mol*emass)
    
    ! 96
    CASE("Cm")
      AtomicMass = 247_DP/(mol*emass)

    ! 97
    CASE("Bk")
      AtomicMass = 247_DP/(mol*emass)

    ! 98
    CASE("Cf")
      AtomicMass = 251_DP/(mol*emass)
      
    ! 99
    CASE("Es")
      AtomicMass = 252_DP/(mol*emass)

    ! 100
    CASE("Fm")
      AtomicMass = 257_DP/(mol*emass)

    ! 101
    CASE("Md")
      AtomicMass = 258_DP/(mol*emass)

    ! 102
    CASE("No")
      AtomicMass = 259_DP/(mol*emass)

    ! 103
    CASE("Lr")
      AtomicMass = 262_DP/(mol*emass)
      
    ! 104
    CASE("Rf")
      AtomicMass = 261_DP/(mol*emass)

    ! 105
    CASE("Db")
      AtomicMass = 262_DP/(mol*emass)

    ! 106
    CASE("Sg")
      AtomicMass = 263_DP/(mol*emass)

    ! 107
    CASE("Bh")
      AtomicMass = 262_DP/(mol*emass)

    ! 108
    CASE("Hs")
      AtomicMass = 265_DP/(mol*emass)
      
    ! 109
    CASE("Mt")
      AtomicMass = 266_DP/(mol*emass)

    ! 110
    CASE("Uun")
      AtomicMass = 269_DP/(mol*emass)
    
    ! 111
    CASE("Uuu")
      AtomicMass = 272_DP/(mol*emass)

    ! 112
    CASE("Uub")
      AtomicMass = 277_DP/(mol*emass)

    CASE DEFAULT
      PRINT *, "Atom name input=>",symbol
      STOP "ERROR: Mass not defined for this atom.  Add to Constants.f90 before proceeding!"
!      AtomicMass = 1._DP*mol/emass

  END SELECT

END FUNCTION AtomicMass


END MODULE Constants
