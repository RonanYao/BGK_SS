#include "_BGK_SS.h"
subroutine BGK_SS_Fortran(_VecSum) &
    (Vec, n, ierr, mpi_comm)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    MATRIX_TYPE, intent(inout) :: Vec(n)
    integer, intent(in) :: n
    integer, intent(out) :: ierr
    integer, intent(in), optional :: mpi_comm
!---------------local variables--------------
    integer :: comm


#ifdef MPI
    if ( present(mpi_comm) ) then
       comm = mpi_comm       
    else
       comm = MPI_COMM_SELF
    end if
    call MPI_ALLREDUCE(MPI_IN_PLACE, vec, n, MPI_TYPE, MPI_SUM, comm, ierr)
#else
    ierr = 0
#endif
    end subroutine
    
subroutine BGK_SS_Fortran(_ModifyGS_QR) &
           (ptr, nrow, ncol, work, IWORK, R)
    implicit none
    type(Indicator), intent(in) :: ptr
    integer, intent(in)         :: IWORK, nrow, ncol
    MATRIX_TYPE, intent(inout)  :: work(IWORK)
    MATRIX_TYPE, intent(out)    :: R(ncol,ncol)
!---------------local variables--------------
    MATRIX_TYPE :: ddot, ZDOTC
    REAL_TYPE   :: DNRM2, DZNRM2
    integer     :: info, i, j, L, indx
    
    L = ptr%L
    do i = 1, ncol 
         do j = 1, i-1
#ifdef REALMAT
              R(j,i) =  ddot &
#else
              R(j,i) =  ZDOTC &
#endif
             (nrow, work((j+L+1)*nrow+1:(j+L+2)*nrow), 1, work((i+L+1)*nrow+1:(i+L+2)*nrow),1)
              work((i+L+1)*nrow+1:(i+L+2)*nrow) = work((i+L+1)*nrow+1:(i+L+2)*nrow) - R(j,i)*&
                  work((j+L+1)*nrow+1:(j+L+2)*nrow)
         enddo
         do j = 1, ncol
              if ( i > j ) then
                 R(i,j) = ZERO
              end if
         enddo
#ifdef REALMAT
        R(i,i) =  DNRM2 &
#else
        R(i,i) =  DZNRM2 &
#endif
        (nrow, work((i+L+1)*nrow+1:(i+L+2)*nrow), 1)
        work((i+L+1)*nrow+1:(i+L+2)*nrow) = work((i+L+1)*nrow+1:(i+L+2)*nrow)/R(i,i)
    enddo
  end subroutine
    
  subroutine BGK_SS_Fortran(_SVD) &
  (ptr, NLRB, nrow, work, IWORK, sigma, U, VT, num_rank, info)
    implicit none
    type(Indicator), intent(in) :: ptr
    character, intent(in)       :: NLRB
    integer, intent(in)         :: nrow, IWORK
    MATRIX_TYPE, intent(in)     :: work(IWORK)
    integer, intent(out)        :: num_rank, info
    MATRIX_TYPE, intent(out)    :: U(:,:), VT(:,:)
    REAL_TYPE, intent(out)      :: sigma(:)
!---------------local variables--------------
    character   :: jobU, jobVT
    integer     :: i, L, ncol, sig_size, lwork, infola, LDU, LDVT
    REAL_TYPE   :: sigma_max, num_thres, optlwork
    MATRIX_TYPE,allocatable :: mat(:,:), cwork(:)
#ifndef REALMAT
    REAL_TYPE, allocatable :: rwork(:)
    allocate(rwork(5*max(nrow,ncol)))
#endif

    if ( NLRB == 'N' ) then
       jobU  = 'N'
       jobVT = 'N'
       LDU   =  1
       LDVT  =  1
    else if ( NLRB == 'L' ) then
       jobU  = 'O'
       jobVT = 'N'
       LDU   = nrow
       LDVT  =  1
    else if ( NLRB == 'R' ) then
       jobU  = 'N'
       jobVT = 'O'
       LDU   = 1
       LDVT  = nrow
    else if ( NLRB == 'B' ) then
       jobU  = 'S'
       jobVT = 'S'
       LDU   = nrow
       LDVT  = nrow
    end if
    
    L = ptr%L
    ncol = ptr%L * ptr%M
    num_thres = ptr%num_thres
    sig_size = min(nrow, ncol)
    allocate(mat(nrow, ncol))
    
    do i = 1, ncol
        mat(1:nrow, i) = work((i+L)*nrow+1:(i+L+1)*nrow)
    enddo
    
    call BGK_SS_Fortran(GESVD) &
    (jobU, jobVT, nrow, ncol, mat, nrow, sigma, U, LDU, VT, LDVT, optlwork, -1 &
#ifdef REALMAT
    , info)
#else
    , rwork, info)
#endif
    lwork = int(optlwork)
    allocate(cwork(lwork))

    call BGK_SS_Fortran(GESVD) &
    (jobU, jobVT, nrow, ncol, mat, nrow, sigma, U, LDU, VT, LDVT, cwork, lwork &
#ifdef REALMAT
    , info)
#else
    , rwork, info)
#endif 
        
    sigma_max = sigma(1)
    do num_rank = 1, sig_size
       if ( sigma(num_rank) < num_thres*sigma_max ) then
          exit
       end if
    end do
    
    num_rank = num_rank - 1
    deallocate(mat, cwork)
#ifndef REALMAT
    deallocate(rwork)
#endif
    info = 0 
    end subroutine
    
subroutine BGK_SS_Fortran(_spec_eig)&
           (bat, issymtric, isHermit, order, A, LDA, B, LDB, work, Iwork, eigval, info)
    implicit none
    character, intent(in)       :: bat
    integer, intent(in)         :: order, LDA, LDB, Iwork
    logical, intent(in)         :: issymtric, isHermit
    integer, intent(out)        :: info
    MATRIX_TYPE, intent(inout)  :: A(LDA, *), B(LDB, *), work(Iwork)
    MATRIX_TYPE, intent(out)    :: eigval(*)
!---------------local variables--------------
    integer      :: lwork, infola
    REAL_TYPE    :: optlwork
    MATRIX_TYPE,allocatable  :: cwork(:)
    REAL_TYPE,allocatable    :: rwork(:)
    COMPLEX_TYPE,allocatable :: VL(:,:), VR(:,:)

    
    infola = 0
    if(bat == "G")then
        if(issymtric) then
            

            
        elseif(isHermit) then
            
            

        elseif((.not. issymtric) .or. (.not. isHermit))then
            
            
        endif
        
    elseif(bat == "S") then
        
        if(issymtric) then
            
            
        elseif(isHermit) then
            
            
        elseif((.not. issymtric) .or. (.not. isHermit)) then
            
            allocate(VR(order, order), rwork(2*order), VL(1,1))
            call BGK_SS_Fortran(GEEV) &
                ('N', 'V', order, A, LDA, eigval, VL, 1, VR, order, optlwork, -1, rwork, infola)
            lwork = int(optlwork)
            call BGK_SS_Fortran(GEEV) &
                ('N', 'V', order, A, LDA, eigval, VL, 1, VR, order, work(1:lwork), lwork, rwork, infola)
            A(1:order,1:order) = VR(1:order,1:order)
            deallocate(VR, VL, rwork)
            if(infola .ne. 0) info = -1
            
        endif
    endif
            
    end subroutine
    
subroutine BGK_SS_Fortran(_eigenvec_rotation) &
  (trans, nrow, ncol, k, mat, work, Iwork, ierr)
    implicit none
    character, intent(in)      :: trans
    integer, intent(in)        :: nrow, ncol, k, Iwork
    MATRIX_TYPE, intent(in)    :: mat(k, ncol)
    MATRIX_TYPE, intent(inout) :: work(Iwork)
    integer, intent(out)       :: ierr
 !---------------local variables--------------   
    integer     :: i, jdx, j, l
    MATRIX_TYPE :: sum
    

    if(k > nrow) then
        write(*,*) "Please choose a apt k or nrow."
        ierr = -1
    elseif(Iwork .le. max(nrow, ncol, k))then
        write(*,*) "The size of work is too small"
        ierr = -1
    endif
    
    if(trans == "N") then
        
        
        do i = 1, nrow
            do l = 1, k
                sum = 0
                do j = 1, ncol
                    jdx = (j+L)*nrow
                    sum = work(jdx + i) * mat(k, j) + sum
                enddo
                work(l) = sum
            enddo
            do j = 1, k
                jdx = (j+L)*nrow
                work(jdx + i) = work(j)
            enddo
        enddo
        ierr = 0
        
    elseif(trans == "T") then
        
        
    endif
end subroutine