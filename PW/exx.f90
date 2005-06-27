!
! Copyright (C) 2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

#include "f_defs.h"

module exx

  USE klist,    ONLY :  nks 
  USE kinds,    ONLY : DP
  implicit none

  logical:: lexx=.true. ! if .true. exx is used

  integer :: currentk
  real (kind=DP):: exxalfa=0.d0 ! 1 if exx, 0 elsewhere
  logical:: exxstart=.false. !1 if initialited
  integer:: iunexx
  integer :: exx_nwordwfc
  !
  ! variables defining the auxiliary k-point grid used in X BZ integration
  !
  integer :: nq1=1, nq2=1, nq3=1 ! integers defining the X integration mesh
  integer :: nqs                 ! number of points in the q-gridd
  integer :: nkqs                ! total number of different k+q
  real (kind=DP), allocatable :: &
             xkq(:,:)               ! xkq(3,nkqs) the auxiliary k+q set
  !
  ! let xk(:,ik) + xq(:,iq) = xkq(:,ikq) = S(isym)*xk(ik') + G
  ! 
  !     index_xkq(ik,iq) = ikq
  !     index_xk(ikq)    = ik'
  !     index_sym(ikq)   = isym
  !
  integer, allocatable :: index_xkq(:,:) ! index_xkq(nks,nqs) 
  integer, allocatable :: index_xk(:)    ! index_xk(nkqs)  
  integer, allocatable :: index_sym(:)   ! index_sym(nkqs)

  real (kind=DP):: exxdiv ! 1 if exx, 0 elsewhere

contains
  !------------------------------------------------------------------------
  subroutine exx_grid_init()
  !------------------------------------------------------------------------
  USE symme,     ONLY : nsym, s
  USE cell_base, ONLY : bg, at
  USE lsda_mod,  ONLY : nspin
  USE klist

  integer :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
  integer, allocatable :: temp_index_xk(:), temp_index_sym(:)
  integer, allocatable :: temp_index_ikq(:), new_ikq(:)
  real (kind=DP), allocatable :: temp_xkq(:,:)
  logical:: xk_not_found
  real (kind=DP) :: sxk(3), dxk(3), xk_cryst(3)
  real (kind=DP) :: eps, dq1, dq2, dq3

  call start_clock ('exx_grid')
  eps = 1.d-6

  !
  ! set a safe limit as the maximum number of auxiliary points we may need
  ! and allocate auxiliary arrays
  !
  max_nk = nks * min(48, 2 * nsym)
  allocate ( temp_index_xk(max_nk), temp_index_sym(max_nk) )
  allocate ( temp_index_ikq(max_nk), new_ikq(max_nk) )
  allocate ( temp_xkq(3,max_nk) )
  !
  ! find all k-points equivalent by symmetry to the points in the k-list
  !
  temp_nkqs = 0
  do isym=1,nsym
     do ik =1, nks
        xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
        sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                 s(:,2,isym)*xk_cryst(2) + &
                 s(:,3,isym)*xk_cryst(3)
        ! add sxk to the auxiliary list if it is not already present
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) = isym 
        end if

        sxk(:) = - sxk(:)
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:) - temp_xkq(:,ikq) - nint(sxk(:) - temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) =-isym 
        end if

     end do
  end do

  !
  ! define the q-mesh step-sizes
  !
  dq1= 1.d0/dble(nq1)
  dq2= 1.d0/dble(nq2)
  dq3= 1.d0/dble(nq3)
  !
  ! allocate and fill the array index_xkq(nks,nqs)
  !
  nqs = nq1 * nq2 * nq3
  if ( nspin == 2 ) then
     allocate ( index_xkq(2*nks,nqs) )
  else
     allocate ( index_xkq(nks,nqs) )
  end if
  nkqs = 0
  new_ikq(:) = 0
  do ik=1,nks
     xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1

            xk_not_found = .true.
            do ikq=1, temp_nkqs
               if ( xk_not_found ) then
                  dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                  if ( abs(dxk(1)).le.eps .and. &
                       abs(dxk(2)).le.eps .and. &
                       abs(dxk(3)).le.eps ) then

                       xk_not_found = .false.

                       if( new_ikq(ikq) == 0) then
                           nkqs = nkqs + 1
                           temp_index_ikq(nkqs) = ikq
                           new_ikq(ikq) = nkqs
                       end if
                       index_xkq(ik,iq) = new_ikq(ikq)

                  end if
               end if
            end do
            if (xk_not_found) call errore('exx_grid_init', &
                              ' k + q is not an S*k ', (ik-1) * nqs + iq )
         end do
       end do
     end do

  end do
  !
  ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
  !
  if ( nspin == 2 ) then
     allocate ( xkq(3,2*nkqs), index_xk(2*nkqs), index_sym(2*nkqs) )
  else
     allocate ( xkq(3,nkqs), index_xk(nkqs), index_sym(nkqs) )
  end if

  do ik =1, nkqs
     ikq = temp_index_ikq(ik)
     xkq(:,ik) = bg(:,1)*temp_xkq(1,ikq) + &
                 bg(:,2)*temp_xkq(2,ikq) + &
                 bg(:,3)*temp_xkq(3,ikq)
     index_xk(ik)  = temp_index_xk(ikq)
     index_sym(ik) = temp_index_sym(ikq)
  end do
  !
  ! clean up
  !
  deallocate (temp_index_xk, temp_index_sym, temp_index_ikq, new_ikq, temp_xkq)
  !
  ! check that everything is what it should be
  !
  call exx_grid_check
  call stop_clock ('exx_grid')

  return
  end subroutine exx_grid_init

  !------------------------------------------------------------------------
  subroutine exx_grid_check
  !------------------------------------------------------------------------
  USE symme,     ONLY : nsym, s
  USE cell_base, ONLY : bg, at
  USE lsda_mod,  ONLY : nspin
  USE klist
  real (kind=DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
  integer :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
  real (kind=DP) :: eps, dq1, dq2, dq3
  eps = 1.d-6
  dq1= 1.d0/dble(nq1)
  dq2= 1.d0/dble(nq2)
  dq3= 1.d0/dble(nq3)

  do ik =1, nks
     xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1
            
            ikq  = index_xkq(ik,iq) 
            ikk  = index_xk(ikq)
            isym = index_sym(ikq)

            xkk_cryst(:)=at(1,:)*xk(1,ikk)+at(2,:)*xk(2,ikk)+at(3,:)*xk(3,ikk)
            if (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
            isym = abs (isym)
            dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                     s(:,2,isym)*xkk_cryst(2) + &
                     s(:,3,isym)*xkk_cryst(3) - sxk(:)
            dxk(:) = dxk(:) - nint(dxk(:))
            if ( .not. ( abs(dxk(1)).le.eps .and. &
                         abs(dxk(2)).le.eps .and. &
                         abs(dxk(3)).le.eps )   ) then
                 write(*,*) ik,iq
                 write(*,*) ikq,ikk,isym
                 write(*,*) dxk(:)
                 call errore('exx_grid_check', &
                             'something wrong', 1 )
            end if

         end do
       end do
     end do
  end do

  write (*,*) ' EXX GRID CHECK SUCCESSFUL '

  return

  end subroutine exx_grid_check

  !------------------------------------------------------------------------
  subroutine exxinit()
  !------------------------------------------------------------------------

    !This subroutine is run before the first H_psi() of each iteration.
    !It saves the wavefunctions for the right density matrix. in real space
    !It saves all the wavefunctions in a single file called prefix.exx
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc
    USE io_files,             ONLY : prefix
    USE io_files,             ONLY : tmp_dir, iunwfc, iunigk
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
    USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                     nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
    USE parser,               ONLY : find_free_unit
    USE symme,                ONLY : nsym, s, ftau

    integer :: ios, ik,ibnd, i, j, k, ir, ri, rj, rk, isym, ikq
    COMPLEX(KIND=DP),allocatable :: temppsic(:), psic(:), tempevc(:,:)
    logical, allocatable :: present(:)
    logical :: exst
    integer, allocatable :: rir(:,:)

    call start_clock ('exxinit')

    allocate(present(nsym),rir(nrxxs,nsym))
    allocate(temppsic(nrxxs), psic(nrxxs), tempevc( npwx, nbnd ))

    exx_nwordwfc=2*nrxxs
    if (.not.exxstart) then 
       iunexx = find_free_unit()
       call diropn(iunexx,'exx', exx_nwordwfc, exst) 

       exxdiv = exx_divergence(nq1) 
       exxalfa = 1.d0

       exxstart=.true.
    endif


    IF ( nks > 1 ) REWIND( iunigk )

    present(1:nsym) = .false.
    do ikq =1,nkqs
       isym = abs(index_sym(ikq))
       if (.not. present(isym) ) then
          present(isym) = .true.
          if ( mod (s (2, 1, isym) * nr1s, nr2s) .ne.0 .or. &
               mod (s (3, 1, isym) * nr1s, nr3s) .ne.0 .or. &
               mod (s (1, 2, isym) * nr2s, nr1s) .ne.0 .or. &
               mod (s (3, 2, isym) * nr2s, nr3s) .ne.0 .or. &
               mod (s (1, 3, isym) * nr3s, nr1s) .ne.0 .or. &
               mod (s (2, 3, isym) * nr3s, nr2s) .ne.0 ) then
             call errore ('exxinit',' EXX + smooth grid is not working',isym)
          end if
          do ir=1,nrxxs
             rir(ir,isym) = ir
          end do
          do k = 1, nr3s
             do j = 1, nr2s
                do i = 1, nr1s
                   call ruotaijk (s(1,1,isym), ftau(1,isym), i, j, k, &
                                  nr1s, nr2s, nr3s, ri, rj , rk )
                   ir =   i + ( j-1)*nrx1s + ( k-1)*nrx1s*nrx2s
                   rir(ir,isym) = ri + (rj-1)*nrx1s + (rk-1)*nrx1s*nrx2s
                end do
             end do
          end do

       end if
    end do

    DO ik = 1, nks
       call davcio (tempevc, nwordwfc, iunwfc, ik, -1)
       IF ( nks > 1 ) READ( iunigk ) npw, igk

       do ibnd =1, nbnd     
          temppsic(:) = ( 0.D0, 0.D0 )
          temppsic(nls(igk(1:npw))) = tempevc(1:npw,ibnd)
          CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )

          do ikq=1,nkqs
             if (index_xk(ikq) .ne. ik) cycle

             isym = abs(index_sym(ikq) )
             psic(rir(1:nrxxs,isym)) = temppsic(1:nrxxs)
             if (index_sym(ikq) < 0 ) psic(1:nrxxs) = conjg(psic(1:nrxxs))

             CALL davcio( psic, exx_nwordwfc, iunexx, (ikq-1)*nbnd+ibnd, 1 )
          end do
       end do
    end do
    deallocate(temppsic, psic, tempevc)
    deallocate(present,rir)

    call stop_clock ('exxinit')  

  end subroutine exxinit

  !-----------------------------------------------------------------------
  subroutine vexx(lda, n, m, psi, hpsi)
  !-----------------------------------------------------------------------
    !This routine calculates V_xx \Psi

    ! ... This routine computes the product of the Hamiltonian
    ! ... matrix with m wavefunctions contained in psi
    !
    ! ... input:
    ! ...    lda   leading dimension of arrays psi, spsi, hpsi
    ! ...    n     true dimension of psi, spsi, hpsi
    ! ...    m     number of states psi
    ! ...    psi
    !
    ! ... output:
    ! ...    hpsi  Vexx*psi
    !
    USE cell_base, ONLY : alat, omega, bg, at
    USE symme,     ONLY : nsym, s
    USE gvect,     ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, ngm, gstart
    USE gsmooth,   ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                           nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg, et
    USE klist,     ONLY : xk,wk
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl

    INTEGER          :: lda, n, m, kpsi
    COMPLEX(KIND=DP) :: psi(lda,m) 
    COMPLEX(KIND=DP) :: hpsi(lda,m)

    ! local variables
    COMPLEX(KIND=DP), allocatable :: tempphic(:), temppsic(:), result(:)
    COMPLEX(KIND=DP), allocatable :: rhoc(:), vc(:)
    real (kind=DP),   allocatable :: fac(:)
    integer          :: ibnd, ik, im , ig, ir,  ikq, iq, isym
    real(kind=DP), parameter  :: fpi = 4.d0 * 3.14159265358979d0, e2  = 2.d0
    real(kind=DP) :: tpiba2, qq, xk_cryst(3), sxk(3), xkq(3)

    call start_clock ('vexx')

    allocate (tempphic(nrxxs), temppsic(nrxxs), result(nrxxs), &
              rhoc(nrxxs), vc(nrxxs), fac(ngm) )

    tpiba2 = (fpi / 2.d0 / alat) **2
  
    ! write (*,*) exx_nwordwfc,lda,n,m, lda*n

    do im=1,m !for each band of psi (the k cycle is outside band)
       temppsic(:) = ( 0.D0, 0.D0 )
       temppsic(nls(igk(1:npw))) = psi(1:npw,im)
       CALL cft3s( temppsic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2 )
       
       result(:) = (0.d0,0.d0)

       do iq = 1, nqs
          ikq  = index_xkq(currentk,iq)
          ik   = index_xk(ikq)
          isym = abs(index_sym(ikq))

          xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
          if (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                   s(:,2,isym)*xk_cryst(2) + &
                   s(:,3,isym)*xk_cryst(3) 
          xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

          do ig=1,ngm
             qq = (xkq(1)-xk(1,currentk)+g(1,ig))**2 + &
                  (xkq(2)-xk(2,currentk)+g(2,ig))**2 + &
                  (xkq(3)-xk(3,currentk)+g(3,ig))**2
             if (qq.gt.1.d-8) then
                fac(ig)=e2*fpi/tpiba2/qq
             else
                fac(ig)= - exxdiv ! or rather something else (see F.Gygi)
!                fac(ig)= 0.d0 ! or rather something else (see F.Gygi)
             end if
          end do
          do ibnd=1,nbnd !for each band of psi
             if ( abs(wg(ibnd,ik)) < 1.d-6) cycle
             !
             !loads the phi from file
             !
             CALL davcio(tempphic,exx_nwordwfc,iunexx,(ikq-1)*nbnd+ibnd,-1)
             !calculate rho in real space
             rhoc(:)=conjg(tempphic(:))*temppsic(:) / omega
             !brings it to G-space
             CALL cft3s( rhoc,nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -1 )
   
             vc(:) = ( 0.D0, 0.D0 )
             do ig=1,ngm
                vc(nls(ig)) = fac(ig) * rhoc(nls(ig))
             end do
             vc = vc * wg (ibnd, ik) / wk(ik) / nqs

             !brings back v in real space
             CALL cft3s( vc, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 1 ) 

             do ir=1,nrxxs
                !accumulates over bands and k points
                result(ir)=result(ir)+vc(ir)*tempphic(ir)
             end do
          end do
       end do

       !brings back result in G-space
       CALL cft3s( result, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2 )
       !adds it to hpsi
       hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
    end do

    deallocate (tempphic, temppsic, result, rhoc, vc, fac )

    call stop_clock ('vexx')

     end subroutine vexx

  function exxenergy ()
    ! This function is called to correct the deband value and have 
    ! the correct energy 
    USE io_files,   ONLY : iunigk,iunwfc, nwordwfc
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg
    USE wavefunctions_module, ONLY : evc
    USE lsda_mod,   ONLY : lsda, current_spin, isk

    implicit none
    REAL (KIND=DP)   :: exxenergy,  energy
    INTEGER          :: ibnd, ik
    COMPLEX(KIND=DP) :: vxpsi ( npwx, nbnd ), psi(npwx,nbnd)
    COMPLEX(KIND=DP) :: ZDOTC

    call start_clock ('exxenergy')

    energy=0.d0
    IF ( nks > 1 ) REWIND( iunigk )
    do ik=1,nks
       currentk = ik
       IF ( nks > 1 ) READ( iunigk ) npw, igk
       IF ( lsda ) current_spin = isk(ik)
       call davcio (psi, nwordwfc, iunwfc, ik, -1)
       vxpsi(:,:) = (0.d0, 0.d0)
!!    subroutine vexx(lda, n, m, psi, hpsi)
       call vexx(npwx,npw,nbnd,psi,vxpsi)
       do ibnd=1,nbnd
          energy = energy + &
                   wg(ibnd,ik) * ZDOTC(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1)
       end do
    end do
    exxenergy = energy

    call stop_clock ('exxenergy')
  end function exxenergy

  function exx_divergence (nq)

     USE cell_base, ONLY : bg, alat, omega
     USE gvect,     ONLY : ngm, g, ecutwfc

     real(kind=DP) :: exx_divergence
     integer       :: nq, nqq

     ! local variables
     integer :: iq1,iq2,iq3, ig
     real(kind=DP) :: div, dq1, dq2, dq3, xq(3), qq, tpiba2, alpha
     real(kind=DP), parameter  :: fpi = 4.d0 * 3.14159265358979d0, e2  = 2.d0

     call start_clock ('exx_div')

     tpiba2 = (fpi / 2.d0 / alat) **2

     alpha  = 10.d0 * tpiba2 / ecutwfc

     dq1= 1.d0/dble(nq) ! dq1= 1.d0/dble(nq1)
     dq2= 1.d0/dble(nq) ! dq1= 1.d0/dble(nq2)
     dq3= 1.d0/dble(nq) ! dq1= 1.d0/dble(nq3)

     nqq = nq**3
 
     div = 0.d0
     do iq1=1,nq
        do iq2=1,nq
           do iq3=1,nq
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3 
              do ig=1,ngm
                 qq = ( xq(1) + g(1,ig) )**2 + &
                      ( xq(2) + g(2,ig) )**2 + &
                      ( xq(3) + g(3,ig) )**2
                 if ( qq.gt.1.d-8 ) then
                    div = div + exp( -alpha * qq) / qq
                 else
                    div = div - alpha ! or maybe something else
                 end if
              end do
           end do
        end do
     end do
!     div = div * e2 * fpi / tpiba2 / nqs
     div = div * e2 * fpi / tpiba2 / nqq

     alpha = alpha / tpiba2

     div = div - e2*omega/sqrt(alpha*0.25d0*fpi)

     exx_divergence = div * nqs

     write (*,'(a,i4,a,3f12.4)') 'EXX divergence (',nq,')= ', &
                                  div, alpha

     call stop_clock ('exx_div')

     call print_clock ('exx_div')

     return
  end function exx_divergence 

end module exx
