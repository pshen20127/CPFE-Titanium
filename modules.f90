!***************
 module NumKind
!***************

    ! This module defines the kind of integer and real numbers.
    ! Every module, subroutine or func must use this module.
    ! To change the precision from double to single,
    ! only this module needs to be changed.
    implicit none

    integer(kind(1)),parameter :: Ikind=kind(1),Rkind=kind(0.0d0)

    real(Rkind), parameter		:: pzero=0.d0,pone=1.d0,ptwo=2.d0,pthree=3.d0
    real(Rkind), parameter		:: pfour=4.d0,pfive=5.d0,psix=6.d0,pseven=7.d0
    real(Rkind), parameter		:: peight=8.d0,pnine=9.d0,pten=10.d0
    real(Rkind), parameter		:: pthird= 0.333333333333333d0, phalf = 0.5d0, ptwothrd= 0.666666666666667d0    
    real(Rkind), parameter      :: sqr23 = 0.816496580927726d0
    real(Rkind), parameter      :: sqr32 = 1.224744871391589d0
    real(Rkind), parameter      :: sqr2  = 1.414213562373095d0
    real(Rkind), parameter      :: sqr3  = 1.732050807568877d0 
    real(Rkind), parameter      :: Ident2nd(9) = (/1.d0, 0.d0, 0.d0,                &
                                                   0.d0, 1.d0, 0.d0,                &
                                                   0.d0, 0.d0, 1.d0/)
!                                              
    real(Rkind), parameter      ::   Ident(6) = (/1.d0, 1.d0, 1.d0, 0.d0,  0.d0,  0.d0/)
    real(Rkind), parameter      ::   factoTwo(6) = (/1.d0, 1.d0, 1.d0, 2.d0,  2.d0,  2.d0/)
    real(Rkind), parameter      ::   factoHalf(6) = (/1.d0, 1.d0, 1.d0, 0.5d0, 0.5d0, 0.5d0/)                                
!
end module NumKind
!----------------------------------------------------------------------
module AbaData
    ! print the results at element NO. numqpt_aba integration point at NO. numel_aba
	integer, parameter  :: numel_aba=-1, numqpt_aba=-1
! 
end
!----------------------------------------------------------------------
!
module timing
  integer c1, c2, cm, cr
  real(kind=8) :: longest_step,printed_step
  integer :: most_psis,printed_psis
  integer :: most_updates,printed_updates
  logical :: printing

end module
!-----------------------------------------------------------------------
!
module WorkDir
!
!-----------------------------------------------------------------------
!----Define the file path and root name
    character(len=255) :: jobname, outdir, outdir1, outdir2, outdir3
    integer :: lenjobname, lenoutdir
    character(len=255), parameter :: fileroot="test"
!   NumMat indicated number of materials usedd in the simulation, current code support maximum 3 different material
!   and need to expand IterPar, ElaPar, PlaPar, ... accordingly if more materal type presented 
    integer, parameter :: NumMat=1
!   
end module WorkDir
!-----------------------------------------------------------------------
!
module DataType
    use numkind
!-----------------------------------------------------------------------
    type EProp
        real*8  fCeDev(5,5)
        real*8  fCeiDev(5,5)
        real*8  fCeDevVol(5)
        real*8  fCeVol
    end type EProp        
!-----------------------------------------------------------------------
    type XtalPar
        !real*8, allocatable ::  matProp(20,48)
        !real*8, allocatable ::  hardmtx(48,48)
        integer numgrn
        integer numslip
        integer numvtx
        real*8  kappa0(48)
        real*8  rhofor0(48)
        real*8  rhodeb0(48)
        real*8  rhofwd0(48)
        real*8  rhorevp0(48)
        real*8  rhorevn0(48)
        real*8  rss0(48)
        real*8  rho0s(48)
        real*8  matProp(40,48)
        real*8  hardmtx(48,48)  
        real*8  sigfs(5, 241)   
        real(kind=rkind)    overstress(48)            
    end type XtalPar  
!-----------------------------------------------------------------------
    type SlipSys
        !real(kind=rkind), allocatable ::  zBar0(3, :, :)
        !real(kind=rkind), allocatable ::  pBar0(:, :, :), qBar0(:, :, :)
        !real(kind=rkind), allocatable ::  pBar0Vec(:, :), qBar0Vec(:, :)
        !real(kind=rkind), allocatable ::  ppTBar0(:, :, :)
        !real(kind=rkind), allocatable ::  tauSlip(:)
        real(kind=rkind)    zBar0(3, 3, 48)
        real(kind=rkind)    pBar0(3, 3, 48), qBar0(3, 3, 48)
        real(kind=rkind)    pBar0Vec(5, 48), qBar0Vec(3, 48)
        real(kind=rkind)    ppTBar0(5, 5, 48)
        real(kind=rkind)    tauSlip(48)  
    end type SlipSys 
!-----------------------------------------------------------------------
    type OriData
        integer kODF, kODFout
        real*8  angles(3, 10000)
        real*8  euler(3, 1000)
        real*8  gcrot0(3, 3, 1000)
        integer numor, seed           
    end type OriData 
!-----------------------------------------------------------------------
    type xtalVars
        !real(kind=rkind), allocatable ::  gstress    (:, :)
        !real(kind=rkind), allocatable ::  gestran    (:, :)
        !real(kind=rkind), allocatable ::  gkappa     (:, :)
        !real(kind=rkind), allocatable ::  gstatev    (:, :)
        !real(kind=rkind), allocatable ::  geqvalues  (:, :)
        !real(kind=rkind), allocatable ::  ggamdot    (:, :)
        !real(kind=rkind), allocatable ::  gcrot      (:, :)
        !real(kind=rkind), allocatable ::  grrot      (:, :)   
        real*8  gstress    (5, 1000)
        real*8  gestran    (5, 1000)
        real*8  gkappa     (48, 1000)
        real*8  grhofor     (48, 1000)
        real*8  grhodeb     (48, 1000)
        real*8  grhofwd     (48, 1000)
        real*8  grhorevp     (48, 1000)
        real*8  grhorevn     (48, 1000)
        real*8  grss     (48, 1000)
        real*8  grho0     (48, 1000)

        real*8  gstatev    (5, 1000)
        real*8  geqvalues  (8, 1000)
        real*8  ggamdot    (48, 1000)
        real*8  gcrot      (3, 3, 1000)
        real*8  grrot      (3, 3, 1000)              
    end type xtalVars
!-----------------------------------------------------------------------
    type xtalVars_n
        !real(kind=rkind), allocatable ::  gstress_n    (:, :)
        !real(kind=rkind), allocatable ::  gestran_n    (:, :)
        !real(kind=rkind), allocatable ::  gkappa_n     (:, :)
        !real(kind=rkind), allocatable ::  gstatev_n    (:, :)
        !real(kind=rkind), allocatable ::  gcrot_n      (:, :)
        !real(kind=rkind), allocatable ::  grrot_n      (:, :)  
        real*8  gstress_n    (5, 1000)
        real*8  gestran_n    (5, 1000)
        real*8  gkappa_n     (48, 1000)
        real*8  grhofor_n     (48, 1000)
        real*8  grhodeb_n     (48, 1000)
        real*8  grhofwd_n     (48, 1000)
        real*8  grhorevp_n     (48, 1000)
        real*8  grhorevn_n     (48, 1000)
        real*8  grss_n     (48, 1000)
        real*8  grho0_n     (48, 1000)

        real*8  gstatev_n    (5, 1000)
        real*8  gcrot_n      (3, 3, 1000)
        real*8  grrot_n      (3, 3, 1000)                
    end  type xtalVars_n
!-----------------------------------------------------------------------    
    type  IterData
        integer maxIterstate, MaxIterNewt
        real*8  tolerState, tolerNewt 
    end type  IterData
end module DataType 
!-----------------------------------------------------------------------
!
module IterPar
    use datatype
    type(IterData) iterP1, iterP2, iterP3
end module IterPar
!-----------------------------------------------------------------------
module ElaPar
    use datatype
    type(EProp) Elap1, Elap2, Elap3
end module ElaPar
!-----------------------------------------------------------------------
module PlaPar
    use datatype
    type(XtalPar) Plap1, Plap2, Plap3
end module PlaPar
!-----------------------------------------------------------------------
module SlipGeo
    use datatype
    type(SlipSys) SlipG1, SlipG2, SlipG3
end module SlipGeo
!-----------------------------------------------------------------------
module OriPar
    use datatype
    type(OriData) OriP1, Orip2, Orip3
end module OriPar
!-----------------------------------------------------------------------
module CPVars
!
    use datatype
    type(xtalVars) CPVars1, CPVars2, CPVars3
!
end module CPVars
!-----------------------------------------------------------------------
module CPVars_n
!
    use datatype
    type(xtalVars_n) CPVars_n1, CPVars_n2, CPVars_n3
!
end module CPVars_n
!-----------------------------------------------------------------------
module TransData
!
    real*8  fDevMat5x6(5, 6), fDevMat6x5(6, 5)
    real*8  fMatTId5x6(5, 6)
!
end module TransData
!-----------------------------------------------------------------------
