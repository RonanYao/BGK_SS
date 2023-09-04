!  ex1.f90 
!
!  FUNCTIONS:
!  ex1 - general eigenvalue problem with a general matrix pencil
!
!****************************************************************************


    program ex1
  

  use BGK_SS
  include 'mpif.h'
  type(Indicator) :: ptr
  integer, parameter :: nrow = 500
  integer            :: i, j, L, M, N, ncv, ierr, &
                        rank, infola, iopt, error
  double precision   :: radius
  complex(kind(0d0)) :: A(nrow, nrow), B(nrow, nrow), &
                        centre
  complex(kind(0d0)), allocatable :: eigval(:)

  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)


  do i = 1, nrow
    do j = 1, nrow
        if ( i == j ) then
            A(i,j) = (2d0,0d0)
        elseif ( abs(i-j) == 1 ) then
          A(i,j) = (1d0,0d0)
        endif
    enddo
  enddo
  do i = 1, nrow
   B(i,i) = (1d0,0d0)
  enddo

  radius = 0.01D0
  centre = (3.99d0,-0.01d0)

  L = 3
  N = 32
  M = 10
  
  !------------------------------------------------------
  call Indicator_init(ptr)
  ptr%L = L
  ptr%N = N
  ptr%M = M
  ptr%centre = centre
  ptr%radius = radius
  ptr%Mpi_Common = MPI_COMM_SELF

  ncv = Indicator_get_ncv(ptr)
  allocate(eigval(ncv))
  CALL cpu_time(system_begin)

  call BGK_SS_DG(ptr, nrow, nrow, A, B, eigval)
  
  CALL cpu_time(system_end)

  if ( rank == 0 ) then
     write(*,'(A23,f15.8,A9)') 'Time of operation was ', system_end -  system_begin,' seconds'
  end if

  deallocate(eigval)
  call Indicator_finalize(ptr)
    end program ex1

