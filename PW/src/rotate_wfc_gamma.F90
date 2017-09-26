!
! Copyright (C) 2003-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_gpu
#else
#define MY_ROUTINE(x)  x##_cpu
#endif
! for calling cublasdger()
#ifdef USE_GPU
  SUBROUTINE cu_dger(m, n, alpha, x, incx, y, incy, A, lda)
    USE kinds, ONLY : DP
    USE cudafor
    USE cublas
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: m, n, incx, incy, lda
    REAL(DP), INTENT(IN) :: alpha
    REAL(DP), DEVICE, INTENT(IN) :: x(incx,m), y(incx,m)
    REAL(DP), DEVICE, INTENT(OUT) :: A(m,n)
    call cublasDger( m, n, alpha, x, incx, y, incy, A, lda )
    !
    RETURN
  END SUBROUTINE cu_dger
#endif
!----------------------------------------------------------------------------
SUBROUTINE MY_ROUTINE( rotate_wfc_gamma )( npwx, npw, nstart, gstart, nbnd, &
                             psi, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Serial version of rotate_wfc for Gamma-only calculations
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : gamma_only 
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_sum 
  USE cpu_gpu_interface
#ifdef USE_GPU
  USE cudafor
  USE gpu_routines
#endif
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd, gstart, ibnd, i, j, istat
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
#ifdef USE_GPU
  ATTRIBUTES( DEVICE ) :: psi, evc, e
#endif
  !
  ! ... auxiliary variables:
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), en(:)
  REAL(DP),    ALLOCATABLE :: hr_h(:,:), sr_h(:,:), vr_h(:,:), en_h(:)
  !
#ifdef USE_GPU
  ATTRIBUTES( DEVICE ) :: aux, hr, sr, vr, en
#endif
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hr( nstart, nstart ), hr_h( nstart, nstart ))    
  ALLOCATE( sr( nstart, nstart ), sr_h( nstart, nstart ) )    
  ALLOCATE( vr( nstart, nstart ), vr_h( nstart, nstart  ) )  
  ALLOCATE( en( nstart ), en_h( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) &
     psi(1,1:nstart) = CMPLX( DBLE( psi(1,1:nstart) ), 0.D0 ,kind=DP)
  !
  CALL h_psi( npwx, npw, nstart, psi, aux )
  !
#ifdef USE_GPU
  istat = cublasDgemm( cublasH, CUBLAS_OP_C, CUBLAS_OP_N, nstart, nstart, 2 * npw, 2.D0 , psi,  2 * npwx, aux, 2 * npwx, 0.D0, hr, nstart )
#else
  CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi,  2 * npwx, aux, 2 * npwx, 0.D0, hr, nstart )
#endif  
  !  
  IF ( gstart == 2 ) THEN
#ifdef USE_GPU
     call cu_dger( nstart, nstart, -1.D0, psi, 2 * npwx, aux, &
                2 * npwx, hr, nstart )
#else
     call DGER( nstart, nstart, -1.D0, psi, 2 * npwx, aux, &
                2 * npwx, hr, nstart )
#endif
  END IF
  !     
  CALL mp_sum(  hr , intra_bgrp_comm )
  !     
  IF ( overlap ) THEN 
     ! 
     CALL s_psi( npwx, npw, nstart, psi, aux )
     !
#ifdef USE_GPU
     istat = cublasDgemm( cublasH, CUBLAS_OP_C, CUBLAS_OP_N, nstart, nstart, 2 * npw, 2.D0 , psi,  2 * npwx, aux, 2 * npwx, 0.D0, sr, nstart )
#else
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0 , psi,  2 * npwx, aux, 2 * npwx, 0.D0, sr, nstart )
#endif
     !            
     IF ( gstart == 2 ) THEN
#ifdef USE_GPU
        CALL cu_dger( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   aux, 2 * npwx, sr, nstart )
#else
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   aux, 2 * npwx, sr, nstart )
#endif
     END IF
     !              
  ELSE
     !
#ifdef USE_GPU
     istat = cublasDgemm( cublasH, CUBLAS_OP_C, CUBLAS_OP_N, nstart, nstart, 2 * npw, 2.D0, psi,  2 * npwx, psi, 2 * npwx, 0.D0, sr, nstart )
#else
     CALL DGEMM( 'T', 'N', nstart, nstart, 2 * npw, 2.D0, psi,  2 * npwx, psi, 2 * npwx, 0.D0, sr, nstart )
#endif
     !
     IF ( gstart == 2 ) THEN
#ifdef USE_GPU
        CALL cu_dger( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   psi, 2 * npwx, sr, nstart )
#else
        CALL DGER( nstart, nstart, -1.D0, psi, 2 * npwx, &
                   psi, 2 * npwx, sr, nstart )
#endif
     END IF
     !
  END IF
  !
  CALL mp_sum(  sr , intra_bgrp_comm )
  !
  ! ... Diagonalize
  !
  ! FIXME ===========================================
  hr_h = hr
  sr_h = sr 
  en_h = en 
  vr_h = vr
  CALL rdiaghg( nstart, nbnd, hr, sr, nstart, en, vr )
  hr = hr_h
  sr = sr_h 
  en = en_h 
  vr = vr_h
  ! ===========================================
  !
#ifdef USE_GPU
!$cuf kernel do(1) <<<*,*>>>
#endif
  DO i=1,nbnd
    e(i) = en(i)
  END DO
  ! e(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
#ifdef USE_GPU
  istat = cublasDgemm( cublasH, CUBLAS_OP_N, CUBLAS_OP_N, 2 * npw, nbnd, nstart, 1.D0, psi, 2 * npwx,  vr, nstart, 0.D0, aux, 2 * npwx )
#else
  CALL DGEMM( 'N', 'N', 2 * npw, nbnd, nstart, 1.D0, psi, 2 * npwx,  vr, nstart, 0.D0, aux, 2 * npwx )
#endif
  !   
  ! evc(:,:) = aux(:,1:nbnd)
#ifdef USE_GPU
!$cuf kernel do(1) <<<*,*>>>
#endif
  DO i=1, nbnd
     DO j=1,npwx
       evc(j,i) = aux(j,i)
     END DO
  END DO
  !
  DEALLOCATE( en, en_h )
  DEALLOCATE( vr, vr_h )
  DEALLOCATE( sr, sr_h )
  DEALLOCATE( hr, hr_h )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE MY_ROUTINE( rotate_wfc_gamma )
!
!
#ifndef USE_GPU
!----------------------------------------------------------------------------
SUBROUTINE protate_wfc_gamma( npwx, npw, nstart, gstart, nbnd, psi, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Parallel version of rotate_wfc for Gamma-only calculations
  ! ... Subroutine with distributed matrices, written by Carlo Cavazzoni
  ! ... This version assumes real wavefunctions (k=0) with only
  ! ... half plane waves stored: psi(-G)=psi*(G), except G=0
  !
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only 
  USE mp_bands,         ONLY : intra_bgrp_comm, nbgrp
  USE mp_diag,          ONLY : ortho_comm, np_ortho, me_ortho, ortho_comm_id,&
                               leg_ortho, ortho_parent_comm, ortho_cntx
  USE descriptors,      ONLY : la_descriptor, descla_init
  USE parallel_toolkit, ONLY : dsqmsym
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
  USE cpu_gpu_interface
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER :: npw, npwx, nstart, nbnd, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
  LOGICAL :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP) :: psi(npwx,nstart), evc(npwx,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP) :: e(nbnd)
    ! eigenvalues
  !
  ! ... auxiliary variables:
  !
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  REAL(DP),    ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), en(:)
  !
  TYPE(la_descriptor) :: desc
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  TYPE(la_descriptor), ALLOCATABLE :: desc_ip( :, : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  !
  Integer :: ibnd
  !

  ALLOCATE( desc_ip( np_ortho(1), np_ortho(2) ) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ) )
  !
  CALL desc_init( nstart, desc, desc_ip )
  !
  ALLOCATE( aux(  npwx, nstart ) )    
  ALLOCATE( hr( nx, nx ) )    
  ALLOCATE( sr( nx, nx ) )    
  ALLOCATE( vr( nx, nx ) )    
  ALLOCATE( en( nstart ) )
  !
  ! ... Set up the Hamiltonian and Overlap matrix on the subspace :
  !
  ! ...      H_ij = <psi_i| H |psi_j>     S_ij = <psi_i| S |psi_j>
  !
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  !
  IF ( gstart == 2 ) &
     psi(1,1:nstart) = CMPLX( DBLE( psi(1,1:nstart) ), 0.D0 ,kind=DP)
  !
  CALL h_psi( npwx, npw, nstart, psi, aux )
  !
  CALL compute_distmat( hr, psi, aux )
  !
  IF ( overlap ) THEN
     !
     CALL s_psi( npwx, npw, nstart, psi, aux )
     CALL compute_distmat( sr, psi, aux )
     !              
  ELSE
     !
     CALL compute_distmat( sr, psi, psi )
     !
  END IF
  !
  ! ... Diagonalize
  !
  CALL prdiaghg( nstart, hr, sr, nx, en, vr, desc )
  !
  e(:) = en(1:nbnd)
  !
  ! ... update the basis set
  !
  CALL refresh_evc( )
  !   
  evc(:,:) = aux(:,1:nbnd)
  !
  DEALLOCATE( desc_ip )
  DEALLOCATE( rank_ip )
  DEALLOCATE( en )
  DEALLOCATE( vr )
  DEALLOCATE( sr )
  DEALLOCATE( hr )
  DEALLOCATE( aux )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE desc_init( nsiz, desc, desc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     TYPE(la_descriptor), INTENT(OUT) :: desc
     TYPE(la_descriptor), INTENT(OUT) :: desc_ip(:,:)
     INTEGER :: i, j, rank
     INTEGER :: coor_ip( 2 )
     ! 
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
     ! 
     nx = desc%nrcx
     !
     DO j = 0, desc%npc - 1
        DO i = 0, desc%npr - 1
           coor_ip( 1 ) = i
           coor_ip( 2 ) = j
           CALL descla_init( desc_ip(i+1,j+1), desc%n, desc%nx, np_ortho, coor_ip, ortho_comm, ortho_cntx, 1 )
           CALL GRID2D_RANK( 'R', desc%npr, desc%npc, i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc%active_node > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init
  !
  !
  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm 
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = 0.0d0
     !
     DO ipc = 1, desc%npc !  loop on column procs 
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = desc_ip( ipr, ipc )%nr
           ir = desc_ip( ipr, ipc )%ir
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, 2*npw, 2.D0 ,  v(1,ir), 2*npwx, w(1,ic), 2*npwx, 0.D0, work, nx )

           IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), 2*npwx, w(1,ic), 2*npwx, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO

     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL dsqmsym( nstart, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = desc_ip( 1, ipc )%nc
        ic = desc_ip( 1, ipc )%ic
        !
        IF( ic <= nbnd ) THEN
           !
           nc = min( nc, nbnd - ic + 1 )
           !
           beta = 0.0d0

           DO ipr = 1, desc%npr
              !
              nr = desc_ip( ipr, ipc )%nr
              ir = desc_ip( ipr, ipc )%ir
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vr(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', 2*npw, nc, nr, 1.D0,  psi(1,ir), 2*npwx, vr, nx, beta, aux(1,ic), 2*npwx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', 2*npw, nc, nr, 1.D0,  psi(1,ir), 2*npwx, vtmp, nx, beta, aux(1,ic), 2*npwx )
              END IF
              ! 

              beta = 1.0d0

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc
  !
END SUBROUTINE protate_wfc_gamma
#endif