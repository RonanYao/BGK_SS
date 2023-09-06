!----------------------------------------------------------------------!
!                         BGK - SS - Common                            !
!----------------------------------------------------------------------!
!                    Fundamental Function Module                       !
!----------------------------------------------------------------------!
! contents:                                                            !
!----------                                                            !
!  ProjsInit                                                           !
!  ProjsInit                                                           !
!  Get_proc_task                                                       !
!  create_rand_matrix  :                                               !                                      
!  create_hutch_samplesv                                               !
!  quad_ell_trap                                                       !
!  Readmat                                                             !
!----------------------------------------------------------------------! 
#include "_BGK_SS.h"
!------------------------------------------------------------------------!
! ProjsInit && Proj_resInit                                              !
!------------------------------------------------------------------------! 
!                                                                        ！
! Initialise and allocate the memories for the projection matrix         !           
! Note: It's better not to use when the matrix is too large.             ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! ptr  =  The global Indicator                                           ！
! nrow  = Row dimension of the projection matrix                         ！
! rank  = Column dimension of the projection matrix                      ！
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! projs/Proj_resInit =  The projection matrix                            ！
!                                                                        ！
!------------------------------------------------------------------------！    
subroutine ProjsInit(ptr, nrow, ncol, projs)
   implicit none
   integer, intent(in)               :: nrow, ncol
   type(Indicator), intent(in)       :: ptr
   MATRIX_TYPE, intent(out), pointer :: projs(:,:)
   
   allocate(ptr%BGK_SS_Fortran(projS)(nrow, ncol))
   ptr%BGK_SS_Fortran(projS) = ZERO
   projs => ptr%BGK_SS_Fortran(projS)
   
end subroutine ProjsInit   

subroutine Proj_resInit(ptr, nrow, ncol, proj_res)
   implicit none
   integer, intent(in)               :: nrow, ncol
   type(Indicator), intent(in)       :: ptr
   MATRIX_TYPE, intent(out), pointer :: proj_res(:,:)
   
   allocate(ptr%BGK_SS_Fortran(proj_res)(nrow, ncol))
   ptr%BGK_SS_Fortran(proj_res) = ZERO
   proj_res => ptr%BGK_SS_Fortran(proj_res)
   
   
end subroutine Proj_resInit   
!------------------------------------------------------------------------!
! Get_proc_task                                                          !
!------------------------------------------------------------------------! 
!                                                                        ！
! Assign the task and get taskindex on every processes                   !           
! Note: default set is 1 and N*L respectively.                           ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! ptr  =  The global Indicator                                           ！
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! taskstart =  The initial index of Current process                      !
!                                                                        ！
! tasknum =  The total num of tasks                                      ！
!                                                                        ！
!------------------------------------------------------------------------！
subroutine Get_proc_task(ptr, taskstart, tasknum)
    implicit none
    type(Indicator), intent(in)  :: ptr
    integer, intent(out)         :: taskstart, tasknum
    !------------local vairables-------------
    integer :: Mpi_Common, rank, size, ierr
    integer :: Num
    
    Mpi_Common = ptr%Mpi_Common
    Num = ptr%N * ptr%L
    call MPI_COMM_RANK(Mpi_Common, rank, ierr)
    call MPI_COMM_SIZE(Mpi_Common, size, ierr)
    
    taskstart = rank*(Num/size) + 1 + min(rank, mod(Num, size))
    tasknum   = Num/size
    if ( mod(Num, size) > rank ) tasknum = tasknum + 1


end subroutine Get_proc_task

!------------------------------------------------------------------------!
! create_rand_matrix                                                     !
!------------------------------------------------------------------------! 
!                                                                        ！
! To initialise the inital rand block matrix                             !           
! Note: 0 or 1 is even distributed in Every entry of this matrix         !                         ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! nrow = row-dimension of V                                              i
! ncol = column dimension of V                                           i
! rank  = maximum order of V allowed.                                    ！
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! V =  The initial matrix                                                !
!                                                                        ！
! ierr  = integer error indicator:                                       !
!         ierr .eq. 0 means normal reture                                ！
!                                                                        ！
!------------------------------------------------------------------------！-     
#ifdef REALMAT
  
  subroutine create_rand_matrix(V, nrow, ncol, rank, ierr)
    implicit none
    integer, intent(in)    :: nrow, ncol, rank
    REAL_TYPE, intent(out) :: V(nrow,ncol)
    integer, intent(out)   :: ierr
    !-----------local scalar-------------
    integer :: iseed(4)
    
    iseed(1) = modulo(rank-2*4096, 4096) ! must be between 0 and 4095
    iseed(2) = modulo(rank-4096, 4096)   ! must be between 0 and 4095
    iseed(3) = modulo(rank, 4096)        ! must be between 0 and 4095
    iseed(4) = 1                         ! must be between 0 and 4095 and odd
    ierr     = -1 
    if(iseed(1) < 0) return
    if(iseed(2) < 0) return
    if(iseed(3) < 0) return
    
    call SLARNV(2, iseed, nrow*ncol, tmpV)
    ierr = 0
  end subroutine

#else

  subroutine create_rand_matrix(V, nrow, ncol, rank,ierr)
    implicit none
    integer, intent(in)       :: nrow, ncol, rank
    COMPLEX_TYPE, intent(out) :: V(nrow,ncol)
    integer, intent(out)      :: ierr
    !-----------local scalar-------------
    integer      :: iseed(4)
    COMPLEX_TYPE :: tmpV(nrow, ncol)
    
    iseed(1) = modulo(rank-2*4096, 4096) ! must be between 0 and 4095
    iseed(2) = modulo(rank-4096, 4096)   ! must be between 0 and 4095
    iseed(3) = modulo(rank, 4096)        ! must be between 0 and 4095
    iseed(4) = 1                         ! must be between 0 and 4095 and odd
    if(iseed(1) < 0) return
    if(iseed(2) < 0) return
    if(iseed(3) < 0) return
    
    ierr     = -1 
    call BGK_SS_Fortran(LARNV) &
         (2, iseed, nrow*ncol, tmpV)
    V(:,1:ncol) = tmpV(:,1:ncol)
    ierr = 0
    
    end subroutine
    
#endif
!------------------------------------------------------------------------!
! Creat hutch samples                                                    !
!------------------------------------------------------------------------! 
!                                                                        ！
! To initialise the inital rand block matrix                             !           
! Note: -1 or 1 is even distributed in Every entry of this matrix        !                         ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! nrow = row-dimension of V                                              i
! ncol = column dimension of V                                           i
! rank  = maximum order of V allowed.                                    ！
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! V =  The initial matrix                                                !
!                                                                        ！
! ierr  = integer error indicator:                                       !
!         ierr .eq. 0 means normal reture                                ！
!                                                                        ！
!------------------------------------------------------------------------！ 
subroutine create_hutch_samples(V, nrow, ncol, ierr)
    implicit none
    integer, intent(in)      :: nrow, ncol
    MATRIX_TYPE, intent(out) :: V(nrow,ncol)
    integer, intent(out)     :: ierr
    !----------------local variables-----------------------------------
    integer                  :: i, j
    double precision         :: ONE_R, ZERO_R
    parameter(ONE_R = 1D0, ZERO_R = 0D0)
    
    call create_rand_matrix(V, nrow, ncol, 0, ierr)
    if(ierr .ne. 0) return
    do j = 1, ncol
       do i = 1, nrow
          V(i,j) = cmplx(sign(ONE_R, real(V(i,j), kind(ZERO_R))), ZERO_R, kind(ZERO_R))
       end do
    end do
    ierr = 0
  end subroutine
!------------------------------------------------------------------------!
! create_hutch_vectors                                                   !
!------------------------------------------------------------------------! 
!                                                                        ！
! To initialise the inital rand block matrix                             !           
! Note: -1 or 1 is even distributed in Every entry of this matrix        !                         ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! nrow = row-dimension of input matrix                                   i
! rhs  = col-dimension of initial block matrix V                         i
! work  = The working array, dimention(Iwork)                            !
! Iwork = The demension of the working array                             !
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! work  = work(nrow+1:rhs*nrow) storing the hutch matrix                 !
!                                                                        ！
! ierr  = integer error indicator:                                       !
!         ierr .eq. 0 means normal reture                                ！
!                                                                        ！
!------------------------------------------------------------------------！ 
subroutine create_hutch_vectors(nrow, rhs, work, Iwork, ierr)
    implicit none
    integer, intent(in)        :: nrow, rhs, Iwork
    MATRIX_TYPE, intent(inout) :: work(Iwork)
    integer, intent(out)       :: ierr
    !----------------local variables-----------------------------------
    integer                  :: i, j
    MATRIX_TYPE              :: V(nrow, rhs)
    double precision         :: ONE_R, ZERO_R
    parameter(ONE_R = 1D0, ZERO_R = 0D0)
    
    call create_rand_matrix(V, nrow, rhs, 0, ierr)
    if(ierr .ne. 0) return
    do i = 1, rhs
       do j = 1, nrow
          work(i*nrow+j) = cmplx(sign(ONE_R, real(V(j,i), kind(ZERO_R))), ZERO_R, kind(ZERO_R))
       end do
    end do
    ierr = 0
  end subroutine   
!------------------------------------------------------------------------!
! quad_ell_trap                                                          !
!------------------------------------------------------------------------! 
!                                                                        ！
! To calculate the Sum up nodes on every processor during the numerical  !
! integral                             !                                 !
! Note: Left and right is adapted for better extraction of interval      !                         ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! ptr  =  The global Indicator                                           ！
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! In ptr:                                                                !
!                                                                        ！
! zeta   =  The exponent of natural constant about coefficient z         !                                      !
! weight =  the weight of numerical integral                             ！
! z      =  The coefficient of matrix B                                  ！
!                                                                        ！
!------------------------------------------------------------------------！     
  subroutine quad_ell_trap(ptr)
    implicit none
    type(Indicator), intent(inout) :: ptr
    !---------------local scalar------------------------
    integer      :: i, N
    COMPLEX_TYPE :: centre
    REAL_TYPE    :: t, pi, radius, ONE_R, ZERO_R, asp_ratio
    PARAMETER(pi = 3.14159265358979323846D0)
    parameter(ONE_R = 1.0D0, ZERO_R = 0.0D0)

    N = ptr%N
    centre = ptr%centre
    radius = ptr%radius
    asp_ratio = ptr%asp_ratio
    allocate(ptr%zeta(N))
    allocate(ptr%weight(N))
    allocate(ptr%z(N))
    
    do i = 1, N
        t = (ONE_R+ONE_R) * pi / N * ((i - 1) + 1/(ONE_R+ONE_R))
        ptr%zeta(i) = cmplx(cos(t), asp_ratio * sin(t), 8)
        ptr%weight(i) = radius * (cmplx(asp_ratio * cos(t), sin(t), 8)) / N
        ptr%z(i) = centre + radius * ptr%zeta(i)
    enddo
    
    end subroutine
    
!------------------------------------------------------------------------!
! Readmat                                                                !
!------------------------------------------------------------------------! 
!                                                                        ！
! To read dense matrix from .mat                                         !                                 !
! Note: Left and right is adapted for better extraction of interval      !                         ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! path       = Filepath of the .mat file                                 ！
! filename   = The filename of .mat file                                 ！
! nrow      = The row of matrix.                                         ！
! ncol      = The col of matrix.                                         ！
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! Amat/Bmat = The dense matrix                                           ！
!                                                                        ！
!------------------------------------------------------------------------！
#include "fintrf.h"
  subroutine Readmat(path, filename, nrow, ncol, Amat, Bmat)
    implicit none
    integer, intent(in)     :: nrow, ncol
    character, intent(in)   :: path(:), filename(:)
    MATRIX_TYPE,intent(out) :: Amat(nrow,ncol), Bmat(nrow,ncol)
    !-----------------mw type---------------
    mwPointer matOpen, matGetVariable, matGetNextVariable
    mwPointer matGetVariableInfo, matGetNextVariableInfo
    mwPointer mp, dir, adir(100), pa, mxGetM, mxGetN  
    mwPointer mxCopyPtrToComplex16,mxCopyPtrToComplex8
    mwPointer mxcopyptrtoreal4, mxcopyptrtoreal8
    mwpointer mxGetComplexDoubles, mxGetComplexSingles
    mwpointer mxGetSingles, mxGetDoubles
    mwPointer i
    mwPointer mrows, ncols
    mwSize size
    !-----------------local scalar---------------
    integer   matClose
    integer   stat
    integer*4 temp
    integer*4 status
    character*32 names(100), name


    mp = matOpen(path//filename, 'r')
    if (mp .eq. 0) then
        write(6,*) 'Can''t open ''matdemo.mat''.'
        stop
    end if
    
      write(6,*) 'Getting Header info from first array.'
      pa = matGetVariableInfo(mp, names(1))
      write(6,*) 'Retrieved ', names(1)
      mrows = mxGetM(pa)
      ncols = mxGetN(pa)
      size = mrows*ncols
      write(6,*) '  With size ', mrows, '-by-', ncols 
      pa = matGetVariable(mp, names(1))
      
#ifdef SINGLE
#ifdef REALMAT
      status = mxcopyptrtoreal4(mxGetSingles(pa),Amat,size)
#else
      status = mxCopyPtrToComplex8(mxGetComplexSingles(pa),Amat,size)
#endif
#else
#ifdef REALMAT
      status = mxcopyptrtoreal8(mxGetDoubles(pa),Amat,size)
#else
      status = mxCopyPtrToComplex16(mxGetComplexDoubles(pa),Amat,size)
#endif
#endif
      if (status .ne. 1) then
          write(6,*) 'Could not copy pa to Complex*16 matrix A .'
      endif
      call mxDestroyArray(pa)
      
      write(6,*) 'Getting Header info from next array.'
      pa = matGetNextVariableInfo(mp, name)
      write(6,*) 'Retrieved ', name
      mrows = mxGetM(pa)
      ncols = mxGetN(pa)
      size = mrows*ncols
      write(6,*) '  With size ', mrows, '-by-', ncols 
      pa = matGetVariable(mp, name)
      
#ifdef SINGLE
#ifdef REALMAT
      status = mxcopyptrtoreal4(mxGetSingles(pa),Bmat,size)
#else
      status = mxCopyPtrToComplex8(mxGetComplexSingles(pa),Bmat,size)
#endif
#else
#ifdef REALMAT
      status = mxcopyptrtoreal8(mxGetDoubles(pa),Bmat,size)
#else
      status = mxCopyPtrToComplex16(mxGetComplexDoubles(pa),Bmat,size)
#endif
#endif
      if (status .ne. 1) then
           write(6,*) 'Could not copy pa to Complex*16 matrix B .'
      endif
      call mxDestroyArray(pa)
      
      write(6,*) 'Getting rest of array contents:'
      pa = matGetNextVariable(mp, name)
      do while (pa .ne. 0)       
          i=mxGetM(pa)
          write(*, *) i 
          write(6,*) 'Retrieved ', name
          write(6,*) '  With size ', mxGetM(pa), '-by-', mxGetN(pa)
          call mxDestroyArray(pa)
          pa = matGetNextVariable(mp, name)
      end do

      stat = matClose(mp)
      if (stat .ne. 0) then
         write(6,*) 'Error closing ''matdemo.mat''.'
         stop
      end if 
    end subroutine
    