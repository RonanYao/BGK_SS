#include "_BGK_SS.h"
!------------------------------------------------------------------------!
! BGK_SS_DG                                                              !
!------------------------------------------------------------------------! 
!                                                                        ！
! To solve the general eigenvalue problem of dense matrix                !           
! Note: This reverse communication interface is only used for debug      !                         ！
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! nrow = row-dimension of A and B                                        !
! ncol = column dimension of A and B                                     !
! A  = The dense matrix A                                                ！
! B  = The dense matrix B                                                !                       ！
!                                                                        ！
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! eigval = the general eigenvalue of dense matriX pencil (A,B)           !
!                                                                        !
!                                                                        ！
!------------------------------------------------------------------------！-
subroutine BGK_SS_DG(ptr, nrow, ncol, A, B, eigval)
    implicit none
    type(Indicator), intent(inout)  :: ptr
    integer, intent(in)          :: nrow, ncol
    MATRIX_TYPE, intent(in)      :: A(nrow, ncol), B(nrow, ncol)
    MATRIX_TYPE, intent(out)     :: eigval(*)
    !------------local vairables-------------
    integer :: i, L, M, Iwork, &
               Mpi_Common, rank, size, ierr
    MATRIX_TYPE, allocatable :: work(:)
    
    if(nrow .ne. ncol) then
        write(*,*) "This method only provides the eigensolver for square matrix."
        return
    endif
    
    L = ptr%L
    M = ptr%M
    Mpi_Common = ptr%Mpi_Common
#ifdef MPI
    call MPI_COMM_RANK(Mpi_Common, rank, ierr)
    call MPI_COMM_SIZE(Mpi_Common, size, ierr)
#else
    rank = 0
    size = 1
#endif

!        %-------------------------------------------------------------%
!        | Pointer into WORKL for address of temporal vector, V, S,    |
!        | etc... and the remaining workspace.                         |
!        | Also update pointer to be used on output.                   |
!        | Memory is laid out as follows:                              |
!        | workl(1:nrow) := temporal vector                            |
!        | workl(nrow+1:(L+1)*nrow) := the initial Block matrix V      |
!        | workl((L+1)*nrow+1:(L*(M+1)+1)*nrow):= projection matrix S  |
!        | In SS-CAA method,                                           |
!        | workl((L+1)*nrow+1:(L*(M+2)+1)*nrow):= projection matrix S  |
!        %-------------------------------------------------------------%
    select case(ptr%Method)
    case(SS_RaylaignRitzs)
        
    case(SS_Arnoldi)
        
    case(SS_ComAvoArnoldi)
        
        Iwork = (L*(2+m)+1)*nrow
        
    end select
    
    allocate(work(Iwork))
    work = 0
    call create_hutch_vectors(nrow, L, work, Iwork, ierr)
    call DenseLinearSolver(ptr, nrow, ncol, A, B, work, IWORK)
    call DenseSVD(ptr, nrow, work, IWORK)
    call DenseEigenSolver(ptr, nrow, work, IWORK, eigval, ierr)
    
    if(ierr .ne. 0) then
    write(*,*) "The eigensolver of BGK_SS_DG had a mistake"
    return
    endif
    
    deallocate(work)
end subroutine BGK_SS_DG  
    