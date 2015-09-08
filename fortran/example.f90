PROGRAM example

use LDC3

 implicit none

 REAL(8), DIMENSION(3) :: alphas, c, alphas_inv
 LOGICAL :: passed

 ! Generate random alpha values
 call random_flat(3,alphas)
 write(*,*) ' '
 write(*,*) 'Generating random alpha vector:'
 write(*,*) 'alpha_h = ',alphas(1)
 write(*,*) 'alpha_r = ',alphas(2)
 write(*,*) 'alpha_t = ',alphas(3)
 write(*,*) ' '

 ! Calculate the corresponding LDCs
 write(*,*) 'Converting alpha vector -> LDCs vector:'
 call forward(alphas(1),alphas(2),alphas(3),c(1),c(2),c(3))
 write(*,*) 'c_2 = ',c(1)
 write(*,*) 'c_3 = ',c(2)
 write(*,*) 'c_4 = ',c(3)

 ! Test if the LDCs satisfy the seven analytic criteria
 call criteriatest(.TRUE.,c(1),c(2),c(3),passed)
 IF ( passed ) THEN
   write(*,*) 'LDCs satisfy the 7 analytic criteria => physically valid'
   write(*,*) ' '
 ELSE
   write(*,*) 'LDCs violate one or more of the 7 analytic criteria'
   write(*,*) ' '
 END IF

 ! Invert LDCs back to alphas
 write(*,*) 'Converting LDCs vector -> alpha vector:'
 call inverse(c(1),c(2),c(3),alphas_inv(1),alphas_inv(2),alphas_inv(3))
 write(*,*) 'alpha_h_inv = ',alphas_inv(1)
 write(*,*) 'alpha_r_inv = ',alphas_inv(2)
 write(*,*) 'alpha_t_inv = ',alphas_inv(3)
 write(*,*) ' '

CONTAINS

! ==============================================================================
SUBROUTINE random_flat(n,r)

 INTEGER i, n
 DOUBLE PRECISION seeda, r(n)

 call random_seed()
 DO i=1,n
   call random_number(seeda)
   r(i)=seeda
 END DO

END SUBROUTINE random_flat
! ==============================================================================

END PROGRAM example
