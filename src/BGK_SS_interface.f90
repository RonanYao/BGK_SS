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
    MATRIX_TYPE, allocatable :: work(:), V(:, :) 
    
    if(nrow .ne. ncol) then
        write(*,*) "This method only provides the eigensolver for square matrix."
        return
    endif
    
    L = ptr%L
    M = ptr%M
    Mpi_Common = ptr%Mpi_Common
    Iwork = (L*(1+m)+1)*nrow
    
    allocate(work(Iwork), V(nrow, L))
    call get_rank_and_size(Mpi_Common, rank, size)
    call create_hutch_samples(V, nrow, L, rank, ierr)
    do i = 1, L
        work((1+i-1)*nrow+1:(i+1)*nrow) = V(1:nrow, i)
    enddo

    call DenseLinearSolver(ptr, nrow, ncol, A, B, work, IWORK, Mpi_Common)
    call DenseSVD(ptr, nrow, work, IWORK)
    call DenseEigenSolver(ptr, nrow, work, IWORK, eigval, ierr)
    
    if(ierr .ne. 0) then
    write(*,*) "The eigensolver of BGK_SS_DG had a mistake"
    return
    endif
    
    deallocate(work, V)
end subroutine BGK_SS_DG  
    