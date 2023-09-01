module BGK_SS
  implicit none

  
  integer, parameter ::  SS_RaylaignRitzs = 1
  integer, parameter ::  SS_Arnoldi       = 2
  integer, parameter ::  SS_ComAvoArnoldi = 3
  
  integer, parameter ::  Standard_ell     = 1
  integer, parameter ::  Users_rule       = 0
  
  integer, parameter ::  BGK_INIT         = 0
  integer, parameter ::  BGK_LinearSolver = 1
  integer, parameter ::  BGK_SVD          = 2
  
  
  
  type timer
     double precision :: start_time
     double precision :: t = 0.0D0
  end type timer

  !> This derived type is used to store optional parameters and context for BGK_SS routines
  type Indicator

     ! Eigenvalue Solver inputs
     integer :: Method   =  SS_ComAvoArnoldi  !< Type of Eigen method
     integer :: Integral =  Standard_ell      !< Type of Integral/quadrature rule
     integer :: N = 32                        !< Number of integral points
     integer :: M = 10                        !< Degree of moment
     integer :: Lstep = 8
     integer :: Lmax = 32                     !< Maximum number of columns of the source matrix
     integer :: n_orth = 3                    !< Number of itaration for orthonormalization
     
     logical :: calc_res  = .true.            !< Flag for calculating residuals
     logical :: calc_vec  = .true.            !< Flag for calculating right eigenvector
     logical :: issymtric = .false.           !< Whether the matrix is symtric
     logical :: isHermit  = .false.           !< Whether the matrix is Hermition
     logical :: trim_out = .true.             !< Discard eigenvalues located outside of the ellipse
     logical :: trim_res = .false. 
     logical :: trim_spu = .true.             !< Discard spurious eigenvalues
     logical :: user_path = .false.
     logical :: sym_contour = .false.
     logical :: force_inner = .false.
     
     double precision :: num_thres = 1d-12 !< Threshold value for numerical rank
     double precision :: asp_ratio = 1d0   !< Aspect ratin of the ellipse
     double precision :: tol = 1d-14       !< Tolerance for residual
     double precision :: spu_thres = 1d-7  !< Thresfold for trimming away the spirious eigenpairs
     double precision :: radius = 1.0D0
     
     complex*16 :: centre = cmplx(0.0D0, 0.0D0, 8)
     complex*16, pointer :: zeta(:)   => NULL()
     complex*16, pointer :: weight(:) => NULL()
     complex*16, pointer :: z(:)      => NULL()
     !Linear Solver Inputs
     
     
     !Mpi inputs
     integer :: Mpi_Common
     
     
     
     !Other inputs
     integer :: Write_Unit = 6 
     integer :: Print_Level = 0 

     ! output
     integer :: itask = BGK_INIT          !< Reverse communication task
     integer :: nev                       !< Estimation of eigenvalue count
     integer :: ws                        !< starting columns of work arrays
     integer :: xs                        !< starting columns of X array
     integer :: nc                        !< number of columns for solve or matvec
     integer :: quad_idx                  !< Index of quadrature points
     integer :: iter = 0                  !< Number of iteration
     integer :: num_basis                 !< Number of basis vectors
     double precision, pointer :: sig_val(:)   => NULL()
     double precision, pointer :: indi_spu(:)  => NULL()

     integer :: imisc(64)
     double precision :: dmisc(64)

     type(timer) :: timer_Init


     integer :: ierr_init
     ! input/output
     integer :: L = 16 !< Number of columns of the source matrix

     ! private

     integer :: x_offset !< User should not touch this variable
     integer :: x_ncol !< User should not touch this variable


     real, pointer             :: sprojS(:,:) => NULL(), sproj_res(:,:) => NULL()
     complex, pointer          :: cprojS(:,:) => NULL(), cproj_res(:,:) => NULL()
     double precision, pointer :: dprojS(:,:) => NULL(), dproj_res(:,:) => NULL()
     complex*16, pointer       :: zprojS(:,:) => NULL(), zproj_res(:,:) => NULL()
     
  end type Indicator
  
  private :: check_inputs, start_timer, stop_timer
    
    interface VecSum
        module procedure s_VecSum
        module procedure d_VecSum
        module procedure c_VecSum
        module procedure z_VecSum
    end interface
    
    interface ModifyGS_QR
        module procedure c_ModifyGS_QR
        module procedure d_ModifyGS_QR
        module procedure s_ModifyGS_QR
        module procedure z_ModifyGS_QR
    end interface
    
    interface SVD
        module procedure s_SVD
        module procedure d_SVD
        module procedure c_SVD
        module procedure z_SVD
    end interface
    
    interface spec_eig
        module procedure s_spec_eig
        module procedure d_spec_eig
        module procedure c_spec_eig
        module procedure z_spec_eig
    end interface
    
    
    contains
  

#define REALMAT
#include "BGK_SS_DBLA.f90"       !real double routine 
#define SINGLE
#include "BGK_SS_DBLA.f90"       !real single routine
#undef REALMAT
#include "BGK_SS_DBLA.f90"       !complex single routine
#undef SINGLE
#include "BGK_SS_DBLA.f90"       !complex double routine
#undef

#include "BGK_SS_DM.f90"
    

#include "BGK_SS_Comm.f90"
#include "BGK_SS_interface.f90"

    
  subroutine Indicator_init(ptr)
    implicit none
#ifdef MPI
    include 'mpif.h'
#endif
    type(Indicator) :: ptr, new_ptr
    
    ptr = new_ptr
#ifdef MPI
    Indicator%Mpi_Common = MPI_COMM_SELF    
#endif


  end subroutine Indicator_init

  integer function Indicator_get_ncv(ptr)
    implicit none
    type(Indicator) :: ptr
    
    if( ptr%Method == SS_ComAvoArnoldi)then
        Indicator_get_ncv = ptr%L*(ptr%M+1) 
    else
        Indicator_get_ncv = ptr%L*ptr%M
    endif
  end function Indicator_get_ncv

  subroutine Indicator_finalize(ptr)
    implicit none
    type(Indicator) :: ptr

    if ( associated(ptr%zeta) ) then
       deallocate(ptr%zeta)
    end if
    if ( associated(ptr%weight) ) then
       deallocate(ptr%weight)
    end if
    if ( associated(ptr%z) ) then
       deallocate(ptr%z)
    end if
    if ( associated(ptr%sprojS) ) then
       deallocate(ptr%sprojS)
    end if
    if ( associated(ptr%cprojS) ) then
       deallocate(ptr%cprojS)
    end if
    if ( associated(ptr%dprojS) ) then
       deallocate(ptr%dprojS)
    end if
    if ( associated(ptr%zprojS) ) then
       deallocate(ptr%zprojS)
    end if 
    
     if ( associated(ptr%sproj_res) ) then
       deallocate(ptr%sproj_res)
    end if
    if ( associated(ptr%cproj_res) ) then
       deallocate(ptr%cproj_res)
    end if
    if ( associated(ptr%dproj_res) ) then
       deallocate(ptr%dproj_res)
    end if
    if ( associated(ptr%zproj_res) ) then
       deallocate(ptr%zproj_res)
    end if 

  end subroutine Indicator_finalize 
  
logical function check_inputs(ptr)
    type(Indicator), intent(inout) :: ptr
    integer :: unit
    logical :: v

    unit = ptr%write_unit
    check_inputs = .false.
    v = ptr%Print_Level 

 end function

  subroutine start_timer(t)
    implicit none
    type(timer), intent(inout) :: t
#ifndef MPI
    integer                    :: itc
    
    call system_clock(itc)
    t%start_time = dble(itc)
#else
    include 'mpif.h'
    t%start_time = mpi_wtime()
#endif
  end subroutine start_timer

  subroutine stop_timer(t)
    implicit none
    type(timer), intent(inout) :: t 
#ifdef MPI  
    include 'mpif.h'
    
    t%t = t%t + mpi_wtime() - t%start_time
#else
    integer          :: itc, t_rate, t_max
    double precision :: dtc, diff
    
    call system_clock(itc, t_rate, t_max)
    dtc = dble(itc)
    if ( dtc < t%start_time ) then
       diff = dble(t_max) - t%start_time + dtc
    else
       diff = dtc - t%start_time
    end if
    diff = diff / dble(t_rate)
    t%t = t%t + diff
#endif
  end subroutine stop_timer

  subroutine clear_timer(t)
    implicit none
    type(timer), intent(inout) :: t
    t%t = 0d0
    
    end subroutine clear_timer
  
end module BGK_SS
