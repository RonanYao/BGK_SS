!------------------------------------------------------------------------!
!                               BGK__SS__DM                              !
!------------------------------------------------------------------------!
!                     Dense Matrix Arithmetic Module                     ！
!------------------------------------------------------------------------!
! contents:                                                              !
!----------                                                              !
!  DenseLinearSolver                                                     !
!  DenseStochastic                                                       !        !
!------------------------------------------------------------------------! 
#include "_BGK_SS.h"  
!------------------------------------------------------------------------!
! DenseLinearSolver                                                      !
!------------------------------------------------------------------------! 
!                                                                        ！
! To solve the general eigenvalue problem of dense matrix                !   
! Note: This reverse communication interface is only used for debug      !   
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! ptr = The global Indicator                                             ！
! nrow = row-dimension of A and B                                        !
! ncol = column dimension of A and B                                     !
! A  = The dense matrix A                                                ！
! B  = The dense matrix B                                                !
! Iwork = The demension of the working array                             !
!                                                                        ！  
! on return:                                                             ！
!----------                                                              ！
!                                                                        ！
! work = The working array, the dimention is equal (L+1+LM)*nrow         !
!        The Top (1:nrow) workspace is for calculation                   !
!        The nrow+1: (L+1)*nrow th is used for storing the initial block !
!        matrix                                                          !
!        The (L+1)*nrow+1: (L+1+LM)*nrow th is used for storing the      ！
!        projection matrix S                                             ！
!                                                                        ！
!        More details about work please see BGK_SS_interface             ！
!                                                                        ！
!------------------------------------------------------------------------！-
subroutine DenseLinearSolver(ptr, nrow, ncol, A, B, work, IWORK)
  implicit none
  type(Indicator), intent(inout), target :: ptr
  integer, intent(in)                    :: nrow, ncol, IWORK
  MATRIX_TYPE, intent(in)                :: A(nrow, ncol), B(nrow, ncol)
  MATRIX_TYPE, intent(inout)             :: work(IWORK)
  !---------------local variables--------------
  integer      :: i, j, jx, L, M, N, IPIV(nrow), mpi_comm, &
                  taskstart, tasknum, Index,  ierr, subs
  COMPLEX_TYPE :: z, weight, zeta, const, temp(ncol)
  REAL_TYPE    :: ratio
  MATRIX_TYPE  :: C(nrow, ncol)
  
  if(nrow .ne. ncol)then
      write(*, *) 'The Matrix is not square !'
      return
  end if
  
  L = ptr%L
  Index = 0
  M = ptr%M
  N = ptr%N
  mpi_comm = ptr%Mpi_Common
  select case(ptr%Method)
      
    case(SS_RaylaignRitzs)
      subs = M
    case(SS_ComAvoArnoldi)
      subs = M + 1
  end select
    
  
  select case(ptr%Integral)
      
    case(Standard_ell)   
      call quad_ell_trap(ptr) 
    case(Users_rule)
      !call users_rule
    end select

  call Get_proc_task(ptr, taskstart, tasknum)
  
  do i = taskstart, (taskstart + tasknum)
      ratio = (i-1)/L
      if(Index .ne. floor(ratio))then
          Index = floor(ratio) + 1
          z  = ptr%z(Index)
          weight = ptr%weight(Index)
          zeta = ptr%zeta(Index)
          C = z*B - A
      endif
      ratio = (i-1)/N
      temp = work((floor(ratio)+1)*ncol+1 : (floor(ratio)+2)*ncol)
      call BGK_SS_Fortran(GETRF) &
           (nrow, ncol, C, nrow, ipiv, ierr) 
      call BGK_SS_Fortran(GETRS) &
           ('N', nrow, 1, C, nrow, ipiv, temp, ncol, ierr) 
      do j = 1, subs
          jx = L*J+1 + floor(ratio)
          work(jx*ncol+1 : (jx+1)*ncol) =  weight*zeta**(j-1) * &
          work(jx*ncol+1 : (jx+1)*ncol)
      enddo
  enddo
  call BGK_SS_Fortran(_VecSum)&
       (work((L+1)*ncol+1:(L*m+L+1)*ncol), L*M*ncol, ierr, mpi_comm)
end subroutine
    
!------------------------------------------------------------------------!
! DenseStochastic                                                        !
!------------------------------------------------------------------------! 
!                                                                        ！
! To approximate the numerical eigenvalue number within the region       !   
! Note: This subroutine is not necessarily called for users              !   
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! ptr   = The global Indicator                                           ！
! nrow  = row-dimension of A and B                                       ！
! work  = The working array, the dimention is equal Iwork                !
!         More details about work please see BGK_SS_interface            ！
! Iwork = The demension of the working array                             !
!                                                                        ！  
! on return:                                                             ！
!----------                                                              ！
! m_num = The approximated number of selected region                     ！
!                                                                        ！
!------------------------------------------------------------------------！- 
subroutine DenseStochastic(ptr, nrow, work, IWORK, m_num)
  implicit none
  type(Indicator), intent(inout), target :: ptr
  integer, intent(in)                    :: nrow, IWORK
  MATRIX_TYPE, intent(in)                :: work(IWORK) 
  MATRIX_TYPE, intent(out)               :: m_num
  !---------------local variables--------------
  integer      :: i, j, L, nev, Mpi_Common, ierr
  MATRIX_TYPE  :: dot(ptr%L)
  MATRIX_TYPE, allocatable :: tmp_m1d(:)

  L = ptr%L

  allocate(tmp_m1d(L))
  
  do i = 1, L
      do j = 1, nrow
#ifdef REALMAT
          dot(i) = dot(i) + work(i*nrow+j)*work((i+L)*nrow+j)
#else
          dot(i) = dot(i) + conjg(work(i*nrow+j))*work((i+L)*nrow+j)
#endif
      enddo
  enddo
  m_num = ZERO
  do i = 1, L
      m_num = m_num + tmp_m1d(i)
  end do
  deallocate(tmp_m1d)

  m_num = m_num / L
  nev = nint(abs(m_num))
  ptr%nev = nev

end subroutine
!------------------------------------------------------------------------!
! DenseSVD                                                               !
!------------------------------------------------------------------------! 
!                                                                        ！
! Singular value decomposition of projection matrix S                    !   
! Note: SVD is adapted for lowrank approximation of projection matrix S  !
!       This subroutine is differ with _SVD for various uses             ！   
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! ptr   = The global Indicator                                           ！
! nrow  = row-dimension of A and B                                       ！
! work  = The working array, dimention(Iwork)                            !
! Iwork = The demension of the working array                             !
!                                                                        ！  
! on return:                                                             ！
!----------                                                              ！
! work  = The working array, dimention(Iwork)                            !
!         More details about work please see BGK_SS_interface            ！
!                                                                        ！
!------------------------------------------------------------------------！
subroutine DenseSVD(ptr, nrow, work, IWORK)
  implicit none
  type(Indicator), intent(inout), target :: ptr
  integer, intent(in)                    :: nrow, IWORK
  MATRIX_TYPE, intent(inout)             :: work(IWORK) 
  !---------------local variables--------------
  integer                  :: i, L, LM, num_basis, infola
  REAL_TYPE, allocatable   :: sigma(:)
  MATRIX_TYPE, allocatable :: tmp_mat(:,:), R(:,:), &
                              H0(:,:), H1(:,:)
  MATRIX_TYPE, pointer     :: projs(:,:), proj_res(:,:)
  
  L = ptr%L
  LM = L * ptr%M
  num_basis = ptr%num_basis
  
  select case(ptr%Method)
  case(SS_RaylaignRitzs)
      
     
  case(SS_ComAvoArnoldi)
      
     
     print *, LM+L, LM

     allocate(R(LM+L,LM+L),tmp_mat(LM+L,LM+L))
     pause
     tmp_mat = 0D0
     allocate(tmp_mat(LM+L,LM+L),sigma(LM))
     call ModifyGS_QR(ptr, nrow, LM + L, work, IWORK, R)

     call ProjsInit(ptr, LM, 2*LM, projs)
     !R_m = R_m+1(1:LM,1:LM)
     R = tmp_mat(1:LM,1:LM)
     
     call BGK_SS_Fortran(_SVD)(ptr, 'B', LM, work, IWORK, sigma, &
          projs(:,1:LM), projs(:,LM+1:2*LM), num_basis, infola)

     ptr%sig_val(1:LM) = sigma(1:LM)

     do i = 1, LM
        projs(1:num_basis,LM+i) = projs(1:num_basis,LM+i) /sigma(1:num_basis)
     end do
     deallocate(R, sigma)
     if(ptr%calc_res) then
         
         allocate(H0(LM,LM), H1(num_basis,LM))        
         H1 = projs(1:num_basis,LM+1: 2*LM) 
         call BGK_SS_Fortran(GEMM)('C', 'C', LM, LM , num_basis, &
              ONE, H1, num_basis, projs(1,1), LM, ZERO, H0, LM)
         
         call Proj_resInit(ptr, LM+L, 2*LM, proj_res)
         proj_res(:,:) = 0
  
         call BGK_SS_Fortran(GEMM)('N', 'N', LM+L, LM, LM &
          , ONE, tmp_mat(:,L+1:LM+L), LM+L, H0,  LM, ZERO, proj_res(1,1), LM+L)
         deallocate(H0,H1)
         
     endif
     allocate(R(LM,num_basis))

     call BGK_SS_Fortran(GEMM)('N', 'C', LM, num_basis, LM &
          , ONE, tmp_mat(1:LM,L+1:LM+L), LM, projs(1,LM+1), LM, ZERO, R, LM)

     tmp_mat = ZERO
     call BGK_SS_Fortran(GEMM)('C', 'N', num_basis, num_basis, LM &
          , ONE, projs(1,1), LM, R, LM, ZERO, tmp_mat(1:LM,L+1:LM+L), LM)

     projs(1:num_basis,LM+1:LM+num_basis)= tmp_mat(1:num_basis,L+1:LM+num_basis)
     deallocate(tmp_mat,R)
     
     end select 

    end subroutine
    
!------------------------------------------------------------------------!
! DenseEigenSolver                                                       !
!------------------------------------------------------------------------! 
!                                                                        ！
! Calculate the eigenvalue problem of small matrix                       !   
! Note: Different method will generate a total different matrix, thus    !
!       It's better not to modify this subroutine.                       ！   
!------------------------------------------------------------------------!
! on entry:                                                              ！
!---------                                                               ！
!                                                                        ！
! ptr   = The global Indicator                                           ！
! nrow  = row-dimension of A and B                                       ！
! work  = The working array, dimention(Iwork)                            !
! Iwork = The demension of the working array                             !
!                                                                        ！  
! on return:                                                             ！
!----------                                                              ！
! work  = The working array, dimention(Iwork)                            !
!         More details about work please see BGK_SS_interface            ！
!                                                                        ！
!------------------------------------------------------------------------！ 
subroutine DenseEigenSolver(ptr, nrow, work, IWORK, eigval, info)
  implicit none
  type(Indicator), intent(inout), target :: ptr
  integer, intent(in)                    :: nrow, IWORK
  MATRIX_TYPE, intent(inout)             :: work(IWORK) 
  integer, intent(out)                   :: info
  MATRIX_TYPE, intent(out)               :: eigval(*)
  !---------------local variables--------------
  integer      :: i, num_basis, ierr, &
                  L, M, LM
  logical      :: issymtric, isHermit
  MATRIX_TYPE  :: hmat(nrow, ptr%L), dot(ptr%L)
  MATRIX_TYPE, allocatable :: Q(:,:), R(:,:), Ut(:,:), &
                              tmpvec(:), tmpmat(:,:)
  MATRIX_TYPE, pointer     :: projS(:,:) => NULL()
  
  select case(ptr%Method)
  case(SS_RaylaignRitzs)
      
     
  case(SS_ComAvoArnoldi)
      
     num_basis = ptr%num_basis
     L = ptr%L
     M = ptr%M
     LM = ptr%L*ptr%M
     issymtric = ptr%issymtric
     isHermit = ptr%isHermit
     projs => ptr%BGK_SS_Fortran(projS)
     allocate(R(LM, LM), Ut(num_basis,num_basis), tmpvec(LM), tmpmat(1,1))
     
     !call start_timer(prm%timer_reduced_eig)
     R = projs(1:LM,LM+1:LM+LM)
     
     call spec_eig('S', issymtric, isHermit, num_basis, R, LM, tmpmat, 1, work, Iwork, eigval, info) 
     
     call BGK_SS_Fortran(GEMM)('N','N', num_basis, num_basis, num_basis &
          , ONE, projs(1:num_basis,1:num_basis), num_basis ,R(1:num_basis,1:num_basis), num_basis, ZERO, Ut, num_basis) 

     projs(1:num_basis,LM+1:LM+num_basis) = Ut
     
     do i = 1, num_basis
           ptr%indi_spu(i) = sum(abs(R(:,i))**2) / sum(abs(R(:,i))**2/ptr%sig_val(1:num_basis))
     end do
     ptr%indi_spu(1:num_basis) = ptr%indi_spu(1:num_basis) / maxval(ptr%indi_spu(1:num_basis))
     
     
     !call start_timer(prm%timer_sub_rot)
     call BGK_SS_Fortran(_eigenvec_rotation)&
          ("N", nrow, num_basis, num_basis, Ut, work, Iwork, ierr)
     !call stop_timer(prm%timer_sub_rot)
     
     
     select case(ptr%Integral)
     case(Standard_ell)
         
        eigval(1:num_basis) = ptr%centre + ptr%radius*eigval(1:num_basis) 
        
     case(Users_rule)

         
     endselect
     !print *,  eigval(1:num_basis)
     !print *, X(:,1),eigval(1)

     deallocate(R, tmpvec, tmpmat, Ut)
     
  end select
  
end subroutine

!------------------------------------------------------------------------£¡
! Dgre
!----------------------------------------------------------------------- 
!
! To initialise the begin block matrix
! Note: -1 or 1 is even distributed in Every entry of this matrix 
!-----------------------------------------------------------------------
! on entry:
!---------
!
! nrow	= row-dimension of V
! ncol	= column dimension of V
! rank  = maximum order of V allowed.
!
! on return:
!---------- 
! 
! V =  The initial matrix
!
! ierr	= integer error indicator: 
!         ierr .eq. 0 means normal retur
!         ierr .eq.-1 means that the the code stopped 
!----------------------------------------------------------------------- 
    
!subroutine DenseCalRes(ptr, nrow, work, IWORK, eigval, info)
!  implicit none
!  type(Indicator), intent(inout), target :: ptr
!  integer, intent(in)                    :: nrow, IWORK
!  MATRIX_TYPE, intent(in)                :: work(IWORK) 
!  integer, intent(out)                   :: info
!  MATRIX_TYPE, intent(out)               :: eigval(*)
!  !---------------local variables--------------
!  integer      :: i, , num_basis, &
!                  L, M, LM
!  MATRIX_TYPE  :: hmat(nrow, L), dot(L)
!  MATRIX_TYPE, allocatable :: Q(:,:), R(:,:), Ut(:,:), &
!                              tmpvec(:), tmpmat(:,:)
!  MATRIX_TYPE, pointer     :: projS(:,:) => NULL()
!  
!  
!  select case(ptr%Method)
!      
!  case(SS_RaylaignRitzs)
!      
!      
!  case(SS_Arnoldi)
!      
!      
!  case(SS_ComAvoArnoldi)
!      
!      allocate(tmp_mat(LM+L,LM),dwork(num_basis))
!       select case(prm%quad_type)
!         case(ZPARES_QUAD_ELL_TRAP)
!             call calc_center_radius(left, right, center, radius)
!             dwork(:) = (eigval(1:num_basis) - center)/radius
!         case(ZPARES_QUAD_USER)
!             if ( present(set_rule) ) then
!                 do i = 1, num_basis
!                     call set_rule(ZPARES_QU_MODE_BACK, dummy1, dummy2, left, right, dummy3, dummy4, dummy5, eigval(i))
!                 end do
!             end if
!         end select
!         !dwork(:) = eigval(1:num_basis)
!         call zpares_rci_sub_get_Mu(prm, Mu)
!         tmp_mat = ZERO_M
!         do i = 1, LM
!             tmp_mat(i,i) =  ONE_M
!         enddo
!         tmp_mat = dwork(ws) * tmp_mat
!         !Mu´æ´¢(H - theta*I)
!         Mu(:,LM+1:2*LM) =  Mu(:,1:LM) -  tmp_mat
!         itask = ZPARES_TASK_MULT_C
!         call start_timer(prm%timer_mult_C)
!         deallocate(tmp_mat,dwork)
!         return
!  end select 
!  end subroutine 
!------------------------------------------------------------------------£¡
! Dgre
!----------------------------------------------------------------------- 
!
! To initialise the begin block matrix
! Note: -1 or 1 is even distributed in Every entry of this matrix 
!-----------------------------------------------------------------------
! on entry:
!---------
!
! nrow	= row-dimension of V
! ncol	= column dimension of V
! rank  = maximum order of V allowed.
!
! on return:
!---------- 
! 
! V =  The initial matrix
!
! ierr	= integer error indicator: 
!         ierr .eq. 0 means normal retur
!         ierr .eq.-1 means that the the code stopped 
!----------------------------------------------------------------------- 
!subroutine EigenCheck(ptr, nrow, work, IWORK, eigval, info)
!  implicit none
!  type(Indicator), intent(inout), target :: ptr
!  integer, intent(in)                    :: nrow, IWORK
!  MATRIX_TYPE, intent(in)                :: work(IWORK) 
!  integer, intent(out)                   :: info
!  MATRIX_TYPE, intent(out)               :: eigval(*)
!  
!   do i = 1, m
!       z = (eigval(i) - center) / radius
!       if ( real(z,kind(ZERO_R))**2 + aimag(z)**2/(asp_ratio**2) <= ONE_R ) then
!          flags(i) = .true.
!          num_true = num_true + 1
!       else
!          flags(i) = .false.
!       end if
!   end do
!   
!     res_max = ZERO_R
!        do i = 1, num_basis
!           if (flags(i)) then
!              res_max = max(res_max, res(i))
!           end if
!        end do
!        flags(1:num_basis) = flags(1:num_basis) .and. (prm%indi_spu(1:num_basis) > prm%spu_thres)
!        counter = count(flags)
!        num_ev = num_basis  
!     if ( prm%trim_out .and. prm%quad_type == ZPARES_QUAD_ELL_TRAP ) then
!        allocate(flags(num_ev))
!        call inside_ellipse(left, right, asp_ratio, num_basis, eigval, flags, j)
!        call packing(num_ev, flags, eigval, X, nrow_local, res, prm%indi_spu)
!        num_ev = j
!        deallocate(flags)
!     end if
!     if ( prm%trim_spu ) then
!        allocate(flags(num_ev))
!        flags(1:num_ev) = prm%indi_spu(1:num_ev) > prm%spu_thres
!        call packing(num_ev, flags, eigval, X, nrow_local, res, prm%indi_spu)
!        num_ev = count(flags)
!        deallocate(flags)
!     end if
!     if ( prm%trim_res ) then
!        allocate(flags(num_ev))
!        flags(1:num_ev) = res(1:num_ev) > prm%tol
!        call packing(num_ev, flags, eigval, X, nrow_local, res, prm%indi_spu)
!        num_ev = count(flags)
!        deallocate(flags)
!     end if
!     
!     end subroutine