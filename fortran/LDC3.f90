MODULE LDC3

! ==============================================================================
! =============================     LDC3 v1.0     ==============================
! ==============================================================================
!
! AUTHOR: David Kipping
!         Columbia University, Dept. of Astronomy
!         Please report any problems to: d.kipping@columbia.edu
!
! CITATION: If using this code, please cite:
!           Kipping, D. M., 2015, 'Efficient, uninformative sampling of limb 
!           darkening coefficients for a three-parameter law', MNRAS, accepted
!
! DESCRIPTION: LDC3 is a module containing three subroutines: "forward", 
!              "inverse" and "criteriatest". These subroutines perform tasks
!              related to identifying physically plausible limb darkening 
!              coefficients (LDCs) in the case of the 3-parameter limb
!              darkening law proposed by Sing et al. (2009), A&A, 505, 891:
!              I(mu)/I(1) = 1 - c_2*[1-mu] - c_3*[1-mu^(3/2)] - c_4*[1-mu^2].
!              The expressions used in LDC3 were derived in Kipping (2015),
!              and we recommend the reader review this paper before continuing.
!              See each subroutine for a description of its purpose. Note that
!              alpha_t denotes alpha_theta throughout this code.
!
! HISTORY: v1.0 Initial version released

 implicit none

CONTAINS
      
!===============================================================================
SUBROUTINE forward(alpha_h,alpha_r,alpha_t,c_2,c_3,c_4)

! DESCRIPTION: forward takes a set of LDCs in the alpha-parameterization and
!              converts them to the c-parameterization. A likely common use
!              of this subroutine would be when fitting LDCs. One would treat
!              the alpha-parameters as the free parameters, but then convert 
!              them to the more standard c-parameters before calling whatever 
!              astronomical model describes one's observations.
!
! INPUTS: alpha_h, alpha_r, alpha_t, all of which are bound between 0 and 1 and 
!         the latter is a wrap-around parameter.
!
! OUTPUTS: c_2, c_3, c_4

 implicit none

 REAL(8), INTENT(IN) :: alpha_h, alpha_r, alpha_t
 REAL(8), INTENT(OUT) :: c_2, c_3, c_4
 REAL(8), PARAMETER :: third = 0.3333333333333333D0
 REAL(8), PARAMETER :: twopi = 6.283185307179586D0
 REAL(8), PARAMETER :: P1 = 4.5008417723138905D0
 REAL(8), PARAMETER :: P2 = 17.14213562373095D0
 REAL(8), PARAMETER :: Q1 = 7.99682547780603D0
 REAL(8), PARAMETER :: Q2 = 8.566161603278331D0

 c_2 = (alpha_h**third)*( P1 + 0.25D0*DSQRT(alpha_r)*( &
       -6.0D0*DCOS(twopi*alpha_t) + P2*DSIN(twopi*alpha_t) ) )

 c_3 = (alpha_h**third)*( -Q1 - Q2*DSQRT(alpha_r)*DSIN(twopi*alpha_t) )

 c_4 = (alpha_h**third)*( P1 + 0.25D0*DSQRT(alpha_r)*( &
       6.0D0*DCOS(twopi*alpha_t) + P2*DSIN(twopi*alpha_t) ) )

END SUBROUTINE forward

!===============================================================================
SUBROUTINE inverse(c_2,c_3,c_4,alpha_h,alpha_r,alpha_t)

! DESCRIPTION: inverse performs the reverse transformation of subroutine 
!              forward. It therefore converts c-parameters into 
!              alpha-parameters.
!
! INPUTS: c_2, c_3, c_4 (= LDCs)
!
! OUTPUTS: alpha_h, alpha_r, alpha_t

 implicit none

 REAL(8), INTENT(IN) :: c_2, c_3, c_4
 REAL(8), INTENT(OUT) :: alpha_h, alpha_r, alpha_t
 REAL(8), PARAMETER :: third = 0.3333333333333333D0
 REAL(8), PARAMETER :: twopi = 6.283185307179586D0
 REAL(8), PARAMETER :: F1 = 0.9997221357486548D0
 REAL(8), PARAMETER :: F2 = 1.0002947195445628D0
 REAL(8), PARAMETER :: G1 = 51.396969619669996D0
 REAL(8), PARAMETER :: G2 = 51.426406871192850D0
 REAL(8), PARAMETER :: G3 = 0.5098605589396862D0
 REAL(8), PARAMETER :: G4 = 67.195959492893320D0
 REAL(8), PARAMETER :: G5 = 75.639610306789280D0
 REAL(8), PARAMETER :: H1 = 0.4666386075895370D0
 REAL(8), PARAMETER :: H2 = 0.5252750715749255D0

 alpha_h = ( F1*c_2 + F2*c_3 + F1*c_4 )**3

 alpha_r = G3/( G1*c_2 + G2*c_3 + G1*c_4 )**2
 alpha_r = alpha_r*( 576.0D0*(c_2 - c_4)**2 + ( G4*c_2 + G5*c_3 + G4*c_4 )**2 )

 alpha_t = DATAN2( -H1*c_2 - H2*c_3 - H1*c_4 , third*0.5D0*(c_4-c_2) )
 alpha_t = ( alpha_t - twopi*FLOOR( alpha_t/twopi ) )/twopi

END SUBROUTINE inverse
!===============================================================================

!===============================================================================
SUBROUTINE criteriatest(usemod,c_2,c_3,c_4,passed)

! DESCRIPTION: criteriatest tests whether a set of LDCs (with the 
!              c-parameterization) satisfy the seven analytic criteria defined
!              in Kipping (2015). Passing this test indicates that LDCs
!              correspond to a physically allowed limb darkened intensity
!              profile, as defined in Kipping (2015). One may choose to use
!              modified or unmodified versions of the criteria with the usemod
!              logical control. The modified versions are slightly more
!              conservative, cropping ~5% of the allowed parameter volume, but
!              yielding a more symmetric volume.
!
! INPUTS: usemod (logical; true = use modified criteria; false = use unmodified 
!                          criteria)
!         c_2, c_3, c_4 (= LDCs)
!
! OUTPUTS: passed (logical; true = test passed, false = test failed)

 implicit none

 LOGICAL, INTENT(IN) :: usemod
 REAL(8), INTENT(IN) :: c_2, c_3, c_4
 LOGICAL, INTENT(OUT) :: passed

 passed = .TRUE.

 ! Criteria A
 IF( c_2 + c_3 + c_4 .GT. 1.0D0 ) THEN
   passed = .FALSE.
 END IF

 ! Criteria B
 IF( 2.0D0*c_2 + 3.0D0*c_3 + 4.0D0*c_4 .LT. 0.0D0 ) THEN
   passed = .FALSE.
 END IF

 ! Criteria C
 IF( c_2 .LT. 0.0D0 ) THEN
   passed = .FALSE.
 END IF

 ! Criteria D
 IF( c_2 + c_3 + c_4 .LT. 0.0D0 ) THEN
   passed = .FALSE.
 END IF

 ! Criteria E
 IF( c_3 .GT. 0.0D0 ) THEN
   passed = .FALSE.
 END IF

 ! Criteria F
 IF( usemod ) THEN
   ! Modified Criterion F
   IF( c_4 .LT. 0.0D0 ) THEN
     passed = .FALSE.
   END IF
 ELSE
   ! Unmodified Criterion F
   IF( c_4 .LT. -1.0D0 ) THEN
     passed = .FALSE.
   END IF
 END IF

 ! Criteria G
 IF( usemod ) THEN
   ! Modified Criterion F
   IF( 32.0D0*c_2*c_4 .LT. 9.0D0*c_3*c_3 ) THEN
     passed = .FALSE.
   END IF
 ELSE
   ! Unmodified Criterion F
   IF( (0.375D0*c_3/c_4) .LE. 0.0D0 .AND. (0.375D0*c_3/c_4) .GE. -1.0D0 ) THEN
     IF( 32.0D0*c_2*c_4 .LT. 9.0D0*c_3*c_3 ) THEN
       passed = .FALSE.
     END IF
   END IF
 END IF

END SUBROUTINE criteriatest
!===============================================================================

END MODULE LDC3
