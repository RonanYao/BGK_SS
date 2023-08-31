#define package_name "BGKSS"
#define package_version "0.1"
#define package_string "BGKSS-0.1"

#undef REAL_TYPE
#undef COMPLEX_TYPE
#undef ZERO
#undef ONE
#undef MATRIX_TYPE
#undef MPI_TYPE
#undef head
#undef Sparse

#ifdef  SINGLE
#define REAL_TYPE    real
#define COMPLEX_TYPE complex
#ifdef  REALMAT
#define MATRIX_TYPE  real
#define ZERO         0.0
#define ONE          1.0
#define MPI_TYPE     MPI_REAL
#define head         c
#else
#define MATRIX_TYPE  complex
#define ZERO         (0.0,0.0)
#define ONE          (1.0,0.0)
#define MPI_TYPE      MPI_COMPLEX
#define head         s
#endif
#else
#define REAL_TYPE    double precision
#define COMPLEX_TYPE complex(kind(0d0))
#define ZERO         (0d0,0d0)
#define ONE          (1d0,0d0)
#ifdef  REALMAT
#define MATRIX_TYPE  double precision
#define MPI_TYPE     MPI_DOUBLE_PRECISION
#define head         d
#else
#define MATRIX_TYPE  complex(kind(0d0))
#define MPI_TYPE     MPI_DOUBLE_COMPLEX
#define head         z
#endif
#endif

#define BGK_SS_Object(A) A
#define BGK_SS_Fortran(A) BGK_SS_Object(head)A

#define SS_RaylaignRitzs 1
#define SS_Arnoldi       2
#define SS_ComAvoArnoldi 3

#define Standard_ell     1
#define Users_rule       0

#define BGK_None         0
#define BGK_INIT         1
#define BGK_LinearSolver 2