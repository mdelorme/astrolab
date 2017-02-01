SUBROUTINE Fitter(n,t0,acc,data,Rin,Pin,incin,Rg,Pg,ing)
  ! Subroutine to Find the orbital distance, period and observer inclination angle.  
  ! This subroutine makes use of 'Nonlinear Least Squares Fitting'  
  ! Information on this method can be found at: http://mathworld.wolfram.com/NonlinearLeastSquaresFitting.html  
  ! The user needs to supply the observed data in nx5 array (where n is the number of  
  !  observational data points), the value of t0, the value of n, the required accuracy  
  !   of the fit and initial guesses of R, P and inclination angle  
  !
  ! The Columns of the data array should be (in this order):
  !   X-position, Z-position, t and R'
  !
  ! List of Inputs:
  !  data = User observed data
  !  n = number of data points
  !  t0 = time of first observation
  !  acc = desired accuracy of fit    
  !  Rin = starting R guess value    
  !  Pin = starting R guess value    
  !  incin = starting R guess value    
  !  res = array containing the residual for each guess of R, P and inc angle    
  !  Rg = Current R guess    
  !  Pg = Current P guess    
  !  ing = Current inclination angel guess    
  !    
  ! The program calling this Subroutine will need to define four functions:    
  !  func = equation defining R' as a function of t,t0,R,P and inclination angle    
  !  df\_dR = equation defining dR'/dR as a function of t,t0,R,P and inclination angle    
  !  df\_dP = equation defining dR'/dP as a function of t,t0,R,P and inclination angle    
  !  df\_dinc = equation defining dR'/d(inclination angle) as a function of t,t0,R,P and inclination angle    
  !    
  ! \#\#\#\#\#\#\# NB: Ensure that initial guesses are as accurate as possible \#\#\#\#\#\#\#    
  ! \#\#\#\#\#\#\# NB: ENSURE THAT P AND t ARE IN THE SAME UNITS \#\#\#\#\#\#    
  ! \#\#\#\#\#\#\# NB: This Subroutine will need to be compiled with ifort or gfortran \#\#\#\#\#\#    
  !---------------------------------------------

  IMPLICIT NONE    

  INTEGER,INTENT(IN)::n    
  REAL(KIND=8),INTENT(IN)::t0,acc,Rin,Pin,incin    
  REAL(KIND=8),DIMENSION(:,:),INTENT(IN),ALLOCATABLE::data    
  REAL(KIND=8),INTENT(OUT)::Rg,Pg,ing

  INTEGER::i,j,cnt
  INTEGER,DIMENSION(:),ALLOCATABLE::indx    
  REAL(KIND=8)::S,d    
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE::res,Atres,dlam,temp    
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE::A,At,AtA,Id

  write(*,*) '------------------------------------------------------------------------'

  Rg = Rin                                      !----- Sets Rin as the first guess of R     
  Pg = Pin                                      !----- Sets Pin as the first guess of P     
  ing = incin                                   !----- Sets incin as the first guess of inclination angle

  !$<$$<$$<$$<$$<$----- Allocates Arrays needed for Calculation -----$>$$>$$>$$>$!    
  write(*,*) 'Now Allocating Required Arrays'    
  allocate(res(n))                              !----- Allocates Array for the residual    
  allocate(A(n,3))                              !----- Allocates Array for     
  allocate(At(3,n))                             !----- Allocates Array for Transpose of A    
  allocate(AtA(3,3))                            !----- Allocates Array for Transpose of A*A    
  allocate(Atres(3))                            !----- Allocates Array for Transpose of A*res    
  allocate(indx(3))                             !----- Allocates Array for index (See Numerical Recipes)    
  allocate(Id(3,3))                             !----- Allocates Array for Indenity/Inverse Matrix    
  allocate(temp(3))                             !----- Allocates Array for a temperarry array    
  allocate(dlam(3))                             !----- Allocates Array for Storing change in d/di,d/dR and d/dP

  write(*,*) 'All Arrays Allocated'    
  write(*,*) ''

  ! $<$$<$$<$$<$$<$----- Nonlinear Least Squares Fitting -----$>$$>$$>$$>$!    
  S = 100.0d0                                   !----- Sets S to an arbitarily large value to begin loop    
  cnt = 0                                       !----- Set Iteration count to zero

  write(*,*) 'Beginning Nonlinear Least Squares Fitting'    
  write(*,*) ' This may take some time please be patient'

  do while (abs(S).gt.acc)                      !----- Exit when required accuracy achieved    
     S = 0.0d0                                  !----- Sets S to zero to begin iterations       
     do i=1, n
        res(i) = data(i,4) - func(data(i,3),t0,Rg,Pg,ing) !----- Calculates the Residual for each data point          
        S = S + (res(i)**2.0d0)                 !----- Calculates the square of the residuals and adds to S          
        A(i,1) = df\_dR(data(i,3),t0,Rg,Pg,ing)    !----- Computes ith dR'(t)/dR          
        A(i,2) = df\_dP(data(i,3),t0,Rg,Pg,ing)    !----- Computes ith dR'(t)/dP          
        A(i,3) = df\_dinc(data(i,3),t0,Rg,Pg,ing)   !----- Computes ith dR'(t)/dinclination          
     end do

     At = transpose(A)                          !----- Calculates Transpose of A       
     AtA = matmul(At,A)                         !----- Calculates Transpose of A*A       
     Atres = matmul(At,res)                     !----- Calculates Transpose of A*residual       
     Id = 0.0d0                                 !----- Zeros all elements in the Id array

     !$<$$<$$<$$<$$<$----- Calculate Inverse Array -----$>$$>$$>$$>$!       
     !\#\#\#\#\# DO NOT MODIFY THIS SECTION OR INVERSE ARRAY MAY NOT BE CALCULATED CORRECTLY \#\#\#\#\#!       
     do i=1,3       
        Id(i,i)=1.0d0                           !----- Assigns values to the diagonal element of Id (Creates Identity Matrix)       
     end do

     call ludcmp(AtA,3,indx,d)                  !----- Converts(Destroys) AtA       

     do i=1, 3       
        temp(:) = Id(:,i)                       !----- Moves Column i of Idensity Matrix to temp array       
        call lubksb(AtA,3,indx,temp)            !----- Calculates the column i of the Inverse of AtA         
        Id(:,i) = temp(:)                       !----- Moves temp to column i in Id array       
     end do

     !\#\#\#\#\# END OF NO-MODIFICATION ZONE \#\#\#\#\#!

     !$<$$<$$<$$<$$<$---- Calculates Changes New Guesses for next iteration -----$>$$>$$>$!       
     dlam = matmul(Id,Atres)                    !----- Calculates the needed changes to Rg, Pg and ing       
     Rg = Rg + dlam(1)                          !----- Calculates the new value of Rg       
     Pg = Pg + dlam(2)                          !----- Calculates the new value of Pg       
     ing = ing + dlam(3)                        !----- Calculates the new value of ing       
     cnt = cnt + 1                              !----- Increases the Iteration count by 1    
  end do

  deallocate(res,A,At,AtA,Atres,indx,Id,temp,dlam)

  write(*,'(A43)') 'The Values of R, P and Inc have been found'    
  write(*,'(A8,I10,A11)') 'It took ',cnt,' iterations'    
  write(*,*) '------------------------------------------------------------------------'    
  write(*,*) ''

END SUBROUTINE Fitter
!------------------------------------------------------- \\


SUBROUTINE lubksb(a,n,indx,b)  
  ! Subroutine to Carry out backward and forward subsitution.    
  ! This subroutine, when combined with ludcmp is used to find    
  !  the inverse of of a    
  !    
  ! Full details of this Subroutine can be found on page 39    
  !  of Numerical Recipes in Fortran 77    
  !    
  ! The user is advised to read through Numerical Recpipes in Fortran 77    
  !  pages 33 - 40
  !
  ! List of Inputs:    
  !  a = 2D (NxM) array to which the inverse is needed    
  !  n = number of rows in the array    
  !  indx = Index of perumation created by ludcmp    
  !  b = 1D output array containing a single column of the inverse of a    
  !    
  ! Solves the set of n linear equations A $\sum$ X = B. Here a is input, not as the    
  !  matrix A but rather as its LU decomposition, determined by the routine ludcmp.    
  !  indx is input as the permutation vector returned by ludcmp. b(n) is input as    
  !  the right-hand side vector B, and returns with the solution vector X. a, n    
  !  and indx are not modified by this routine and can be left in place for successive    
  !  calls with different right-hand sides b. This routine takes into account the    
  !  possibility that b will begin with many zero elements, so it is efficient for use    
  !  in matrix inversion.    
  !    
  ! \#\#\#\#\#\#\#\#\# NB: This subroutine should NOT be Modified \#\#\#\#\#\#\#\#\#    
  !---------------------------------------------

  IMPLICIT NONE

  INTEGER,INTENT(IN)::n    
  INTEGER,INTENT(IN),DIMENSION(:),ALLOCATABLE::indx    
  REAL(KIND=8),DIMENSION(:,:),INTENT(IN),ALLOCATABLE::a    
  REAL(KIND=8),DIMENSION(:),INTENT(INOUT),ALLOCATABLE::b    
  INTEGER:: i,ii,j,ll    
  REAL(KIND=8)::sum

  ii=0

  do i=1,n    
     ll=indx(i)       
     sum=b(ll)       
     b(ll)=b(i)

     if (ii.ne.0)then       
        do j=ii,i-1          
           sum=sum-a(i,j)*b(j)             
        end do
     else if (sum.ne.0.) then       
        ii=i          
     end if

     b(i)=sum       
  end do

  do i=n,1,-1    
     sum=b(i)

     do j=i+1,n       
        sum=sum-a(i,j)*b(j)          
     end do

     b(i)=sum/a(i,i)       
  end do

  return

END SUBROUTINE lubksb
!------------------------------------------------------- \\


SUBROUTINE ludcmp(a,n,indx,d)  
  ! Subroutine to carry out LU Decomposiontion.        
  ! This subroutine, when combined with ludksb is used to find    
  !  the inverse of of a   
  !   
  ! Full details of this Subroutine can be found on page 38-39    
  !  of Numerical Recpipes in Fortran 77   
  !   
  ! The user is advised to read through Numerical Recpipes in Fortran 77    
  !  pages 33 - 40 
  ! 
  ! Litst of Inputs:    
  !  a = 2D (NxM) array to which the inverse is needed    
  !  n = number of rows in the array    
  !  indx = Index of perumation    
  !  d = Record odd or even row interchanges   
  !   
  ! Given a matrix a(1:n,1:n), with physical dimension np by np, this routine    
  !  replaces it by the LU decomposition of a rowwise permutation of itself.    
  !  a and n are input. a is output, arranged as in equation (2.3.14) above;    
  !  indx(1:n) is an output vector that records the row permutation effected    
  !  by the partial pivoting; d is output as Â±1 depending on whether the number    
  !  of row interchanges was even or odd, respectively. This routine is used in    
  !  combination with lubksb to solve linear equations or invert a matrix.    
  !  
  ! \#\#\#\#\#\#\#\#\# NB: This subroutine should NOT be Modified \#\#\#\#\#\#\#\#\#  
  !---------------------------------------------

  IMPLICIT NONE

  INTEGER,INTENT(IN)::n    
  INTEGER,INTENT(INOUT),DIMENSION(:),ALLOCATABLE::indx    
  REAL(KIND=8),DIMENSION(:,:),INTENT(INOUT),ALLOCATABLE::a    
  REAL(KIND=8),INTENT(OUT)::d    
  INTEGER::NMAX,i,imax,j,k    
  REAL(KIND=8)::TINY,aamax,dum,sum    
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE::vv    
  PARAMETER (NMAX=500,TINY=0.0d0) !Largest expected n, and a small number.

  ALLOCATE(VV(NMAX))    
  d=1.

  !$<$$<$$<$$<$$<$----- Finds the Implicite scaling for Each Row -----$>$$>$$>$$>$$>$!    
  do i=1,n    
     aamax=0.       
     do  j=1,n       
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))          
     end do

     if (aamax.eq.0.) pause 'singular matrix in ludcmp'        
     vv(i)=1./aamax       
  end do

  !$<$$<$$<$$<$$<$----- Crout's Method -----$>$$>$$>$$>$$>$!    
  do j=1,n    
     do i=1,j-1       
        sum=a(i,j)          
        do k=1,i-1          
           sum=sum-a(i,k)*a(k,j)             
        end do
        a(i,j)=sum          
     end do

     aamax=0.0     

     do i=j,n       
        sum=a(i,j)          
        do k=1,j-1          
           sum=sum-a(i,k)*a(k,j)             
        end do
        a(i,j)=sum          
        dum=vv(i)*abs(sum)          
        if (dum.ge.aamax) then          
           imax=i             
           aamax=dum             
        end if
     end do

     if (j.ne.imax)then       
        do k=1,n          
           dum=a(imax,k)             
           a(imax,k)=a(j,k)             
           a(j,k)=dum             
        end do

        d=-d          
        vv(imax)=vv(j)          
     end if

     indx(j)=imax

     if (a(j,j).eq.0.) a(j,j)=TINY

     if (j.ne.n)then       
        dum=1./a(j,j)          
        do i=j+1,n
           a(i,j)=a(i,j)* dum             
        end do
     end if
  end do

  return

END SUBROUTINE ludcmp
!-------------------------------------------------------
