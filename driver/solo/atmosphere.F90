!***********************************************************************
!*                   GNU Lesser General Public License
!*
!* This file is part of the FV3 dynamical core.
!*
!* The FV3 dynamical core is free software: you can redistribute it
!* and/or modify it under the terms of the
!* GNU Lesser General Public License as published by the
!* Free Software Foundation, either version 3 of the License, or
!* (at your option) any later version.
!*
!* The FV3 dynamical core is distributed in the hope that it will be
!* useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!* of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!* See the GNU General Public License for more details.
!*
!* You should have received a copy of the GNU Lesser General Public
!* License along with the FV3 dynamical core.
!* If not, see <http://www.gnu.org/licenses/>.
!***********************************************************************

module atmosphere_mod

!-----------------------------------------------------------------------
!
!    interface for FV dynamical core with Held-Suarez forcing
!
!-----------------------------------------------------------------------
!$ser verbatim use mpi
!$ser verbatim use m_serialize, ONLY: fs_is_serialization_on
use constants_mod, only: grav, kappa, cp_air, pi, rdgas, rvgas, SECONDS_PER_DAY
use fms_mod,       only: error_mesg, FATAL,                &
                         check_nml_error, stdlog, stdout,  &
                         write_version_number,             &
                         mpp_pe, mpp_root_pe, &
                         mpp_error, FATAL, NOTE
use fms2_io_mod,   only: file_exists
use mpp_mod,       only: input_nml_file
use time_manager_mod, only: time_type, get_time, set_time, operator(+)
use mpp_domains_mod,  only: domain2d
use mpp_mod,          only: input_nml_file, mpp_sync_self, mpp_sync, &
                                  mpp_set_current_pelist, mpp_npes, &
                                  mpp_get_current_pelist
!------------------
! FV specific codes:
!------------------
use fv_arrays_mod,      only: fv_atmos_type
use fv_control_mod,     only: fv_control_init, fv_end, ngrids
use fv_phys_mod,        only: fv_phys, fv_nudge, fv_phys_init
use fv_diagnostics_mod, only: fv_diag_init, fv_diag, fv_time, eqv_pot
use fv_timing_mod,      only: timing_on, timing_off, timing_init, timing_prt
use fv_restart_mod,     only: fv_restart
use fv_dynamics_mod,    only: fv_dynamics
use fv_nesting_mod,     only: twoway_nesting
use gfdl_mp_mod,        only: gfdl_mp_init, gfdl_mp_end
use fv_nwp_nudge_mod,   only: fv_nwp_nudge_init, fv_nwp_nudge_end, do_adiabatic_init
use field_manager_mod,  only: MODEL_ATMOS
use tracer_manager_mod, only: get_tracer_index
!-----------------------------------------------------------------------

implicit none
private

public   atmosphere_init, atmosphere,  atmosphere_end, atmosphere_domain
public   mygrid

!-----------------------------------------------------------------------
!---- private data ----

type        (time_type) :: Time_step_atmos
real                    :: dt_atmos
integer :: sec
integer days, seconds

logical :: cold_start      = .false.       ! read in initial condition
integer :: mytile = 1
integer :: p_split = 1
real, allocatable:: lprec(:,:), fprec(:,:), f_land(:,:)

type(fv_atmos_type), allocatable, target :: Atm(:)

logical, allocatable :: grids_on_this_pe(:)
integer :: mygrid = 1 !not used yet
integer, allocatable :: pelist(:)
integer :: axes(4)
integer:: isd, ied, jsd, jed, ngc
!-----------------------------------------------------------------------

! version number of this module
! Include variable "version" to be written to log file.
#include<file_version.h>

contains

!#######################################################################

  subroutine atmosphere_init ( Time_init, Time, Time_step )

    type (time_type), intent(in) :: Time_step
    type (time_type), intent(in) :: Time_init
    type (time_type), intent(in) :: Time

    ! local:
    integer isc, iec, jsc, jec
    real:: zvir
    integer :: n, theta_d

    integer :: nlunit = 9999
    character (len = 64) :: fn_nml = 'input.nml'

   call timing_init
   call timing_on('ATMOS_TOTAL')
   call timing_on('ATMOS_INIT')

    allocate(pelist(mpp_npes()))
    call mpp_get_current_pelist(pelist)
  !----- write version and namelist to log file -----

    call write_version_number ( 'SOLO/ATMOSPHERE_MOD', version )

  !---- compute physics/atmos time step in seconds ----

    Time_step_atmos = Time_step
    call get_time (Time_step_atmos, sec)
    dt_atmos = real(sec)

  !----- initialize FV dynamical core -----
    cold_start = (.not.file_exists('INPUT/fv_core.res.nc') .and. .not.file_exists('INPUT/fv_core.res.tile1.nc'))

    call fv_control_init(Atm, dt_atmos, mygrid, grids_on_this_pe, p_split)  ! allocates Atm components

    call mpp_set_current_pelist(Atm(mygrid)%pelist, no_sync=.TRUE.)

   call timing_on('FV_RESTART')
    call fv_restart(Atm(mygrid)%domain, Atm, dt_atmos, seconds, days, cold_start, &
         Atm(mygrid)%flagstruct%grid_type, mygrid)
   call timing_off('FV_RESTART')

     fv_time = time

     isc = Atm(mygrid)%bd%isc
     iec = Atm(mygrid)%bd%iec
     jsc = Atm(mygrid)%bd%jsc
     jec = Atm(mygrid)%bd%jec
     isd = Atm(mygrid)%bd%isd
     ied = Atm(mygrid)%bd%ied
     jsd = Atm(mygrid)%bd%jsd
     jed = Atm(mygrid)%bd%jed

     Atm(mygrid)%flagstruct%moist_phys = .false. ! need this for fv_diag calendar
     call fv_diag_init(Atm(mygrid:mygrid), axes, Time, Atm(mygrid)%npx, Atm(mygrid)%npy, Atm(mygrid)%npz, Atm(mygrid)%flagstruct%p_ref)

     !   if ( Atm(n)%flagstruct%adiabatic .or. Atm(n)%flagstruct%do_Held_Suarez ) then
     zvir = 0.
     if ( Atm(mygrid)%flagstruct%adiabatic ) then
        Atm(mygrid)%flagstruct%moist_phys = .false.
     else
        zvir = rvgas/rdgas - 1.
        Atm(mygrid)%flagstruct%moist_phys = .true.
     endif
     if (.not. Atm(mygrid)%flagstruct%adiabatic .or. Atm(mygrid)%flagstruct%fv_sg_adj > 0.) then
        call fv_phys_init(isc,iec,jsc,jec,Atm(mygrid)%npz,Atm(mygrid)%flagstruct%nwat, Atm(mygrid)%ts, Atm(mygrid)%pt(isc:iec,jsc:jec,:),   &
                          Time, axes, Atm(mygrid)%gridstruct%agrid(isc:iec,jsc:jec,2))
     endif

     if (.not. Atm(mygrid)%flagstruct%adiabatic) call gfdl_mp_init (input_nml_file, stdlog(), Atm(mygrid)%flagstruct%hydrostatic)


     if ( Atm(mygrid)%flagstruct%nudge )    &
          call fv_nwp_nudge_init( Time, axes, Atm(mygrid)%npz, zvir, Atm(mygrid)%ak, Atm(mygrid)%bk, Atm(mygrid)%ts, &
          Atm(mygrid)%phis, Atm(mygrid)%gridstruct, Atm(mygrid)%ks, Atm(mygrid)%npx, Atm(mygrid)%neststruct, Atm(mygrid)%bd)

     if ( Atm(mygrid)%flagstruct%make_nh ) then
        Atm(mygrid)%w(:,:,:) = 0.
     endif

     if ( Atm(mygrid)%flagstruct%na_init>0 ) then
        call adiabatic_init(zvir,mygrid)
     endif

     theta_d = get_tracer_index (MODEL_ATMOS, 'theta_d')
     if ( theta_d > 0 ) then
        call eqv_pot(Atm(mygrid)%q(isc:iec,jsc:jec,:,theta_d), Atm(mygrid)%pt, Atm(mygrid)%delp,    &
             Atm(mygrid)%delz, Atm(mygrid)%peln, Atm(mygrid)%pkz, Atm(mygrid)%q(isd,jsd,1,1), isc, iec, jsc, jec, Atm(mygrid)%ng,   &
             Atm(mygrid)%npz,  Atm(mygrid)%flagstruct%hydrostatic, Atm(mygrid)%flagstruct%moist_phys)
     endif

     Atm(mygrid)%time_init = Time_init

   call timing_off('ATMOS_INIT')

  end subroutine atmosphere_init

 subroutine adiabatic_init(zvir, n)
   real, allocatable, dimension(:,:,:):: u0, v0, t0, dp0
   real, intent(in):: zvir
   integer, intent(in) :: n
   real, parameter:: wt = 1.5  !  2.
   real:: xt, esl
   integer:: isc, iec, jsc, jec, npz
   integer:: m, i,j,k

   character(len=80) :: errstr

   xt = 1./(1.+wt)
   if ( Atm(n)%flagstruct%moist_phys ) then
        esl = zvir
   else
        esl = 0.
   endif

   write(errstr,'(A, I4, A)') 'Performing adiabatic init',  Atm(n)%flagstruct%na_init, ' times'
   call mpp_error(NOTE, errstr)

    npz = Atm(n)%npz

    isc = Atm(n)%bd%isc
    iec = Atm(n)%bd%iec
    jsc = Atm(n)%bd%jsc
    jec = Atm(n)%bd%jec

    ngc = Atm(n)%ng
    isd = isc - ngc
    ied = iec + ngc
    jsd = jsc - ngc
    jed = jec + ngc

     do_adiabatic_init = .true.

     allocate ( u0(isc:iec,  jsc:jec+1, npz) )
     allocate ( v0(isc:iec+1,jsc:jec,   npz) )
     allocate ( t0(isc:iec,jsc:jec, npz) )
     allocate (dp0(isc:iec,jsc:jec, npz) )
     call p_adi(npz, Atm(n)%ng, isc, iec, jsc, jec, Atm(n)%ptop,  &
                Atm(n)%delp, Atm(n)%ps, Atm(n)%pe,     &
                Atm(n)%peln, Atm(n)%pk, Atm(n)%pkz, Atm(n)%flagstruct%hydrostatic)

!$omp parallel do default(shared)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                u0(i,j,k) = Atm(n)%u(i,j,k)
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                v0(i,j,k) = Atm(n)%v(i,j,k)
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec
!               t0(i,j,k) = Atm(n)%pt(i,j,k)*(1.+esl*Atm(n)%q(i,j,k,1))*(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))
                t0(i,j,k) = Atm(n)%pt(i,j,k)
               dp0(i,j,k) = Atm(n)%delp(i,j,k)
             enddo
          enddo
       enddo

     do m=1,Atm(n)%flagstruct%na_init
! Forward call
    call fv_dynamics(Atm(n)%npx, Atm(n)%npy, npz,  Atm(n)%ncnst, Atm(n)%ng, dt_atmos, 0.,      &
                     Atm(n)%flagstruct%fill, Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(n)%ptop, Atm(n)%ks, Atm(n)%ncnst, Atm(n)%flagstruct%n_split,        &
                     Atm(n)%flagstruct%q_split, Atm(n)%u0, Atm(n)%v0, Atm(n)%u,       &
                     Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%flagstruct%hydrostatic,  &
                     Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,                     &
                     Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz, Atm(n)%phis,      &
                     Atm(n)%q_con, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc, &
                     Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy,                    &
                     Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0, Atm(n)%flagstruct%hybrid_z,    &
                     Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct,         &
                     Atm(n)%thermostruct, Atm(n)%idiag, Atm(n)%bd, Atm(n)%parent_grid,&
                     Atm(n)%domain, Atm(n)%inline_mp, Atm(n)%heat_source, Atm(n)%diss_est)
! Backward
    call fv_dynamics(Atm(n)%npx, Atm(n)%npy, npz,  Atm(n)%ncnst, Atm(n)%ng, -dt_atmos, 0.,      &
                     Atm(n)%flagstruct%fill, Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(n)%ptop, Atm(n)%ks, Atm(n)%ncnst, Atm(n)%flagstruct%n_split,        &
                     Atm(n)%flagstruct%q_split, Atm(n)%u0, Atm(n)%v0, Atm(n)%u,       &
                     Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%flagstruct%hydrostatic,  &
                     Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,                     &
                     Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz, Atm(n)%phis,      &
                     Atm(n)%q_con, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc, &
                     Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy,                    &
                     Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0, Atm(n)%flagstruct%hybrid_z,    &
                     Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct,         &
                     Atm(n)%thermostruct, Atm(n)%idiag, Atm(n)%bd, Atm(n)%parent_grid,  &
                     Atm(n)%domain, Atm(n)%inline_mp, Atm(n)%heat_source, Atm(n)%diss_est)
! Nudging back to IC
!$omp parallel do default(shared)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                Atm(n)%u(i,j,k) = xt*(Atm(n)%u(i,j,k) + wt*u0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                Atm(n)%v(i,j,k) = xt*(Atm(n)%v(i,j,k) + wt*v0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec
                Atm(n)%delp(i,j,k) = xt*(Atm(n)%delp(i,j,k) + wt*dp0(i,j,k))
             enddo
          enddo
       enddo

     call p_adi(npz, Atm(n)%ng, isc, iec, jsc, jec, Atm(n)%ptop,  &
                Atm(n)%delp, Atm(n)%ps, Atm(n)%pe,     &
                Atm(n)%peln, Atm(n)%pk, Atm(n)%pkz, Atm(n)%flagstruct%hydrostatic)
!$omp parallel do default(shared)
       do k=1,npz
          do j=jsc,jec
             do i=isc,iec
!               Atm(n)%pt(i,j,k) = xt*(Atm(n)%pt(i,j,k)+wt*t0(i,j,k)/((1.+esl*Atm(n)%q(i,j,k,1))*(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))))
                Atm(n)%pt(i,j,k) = xt*(Atm(n)%pt(i,j,k)+wt*t0(i,j,k))
             enddo
          enddo
       enddo

! Backward
    call fv_dynamics(Atm(n)%npx, Atm(n)%npy, npz,  Atm(n)%ncnst, Atm(n)%ng, -dt_atmos, 0.,      &
                     Atm(n)%flagstruct%fill, Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(n)%ptop, Atm(n)%ks, Atm(n)%ncnst, Atm(n)%flagstruct%n_split,        &
                     Atm(n)%flagstruct%q_split, Atm(n)%u0, Atm(n)%v0, Atm(n)%u,       &
                     Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%flagstruct%hydrostatic,  &
                     Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,                     &
                     Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz, Atm(n)%phis,      &
                     Atm(n)%q_con, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc, &
                     Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy,                    &
                     Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0, Atm(n)%flagstruct%hybrid_z,    &
                     Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct,         &
                     Atm(n)%thermostruct, Atm(n)%idiag, Atm(n)%bd, Atm(n)%parent_grid,  &
                     Atm(n)%domain, Atm(n)%inline_mp, Atm(n)%heat_source, Atm(n)%diss_est)
! Forwardward call
    call fv_dynamics(Atm(n)%npx, Atm(n)%npy, npz,  Atm(n)%ncnst, Atm(n)%ng, dt_atmos, 0.,      &
                     Atm(n)%flagstruct%fill, Atm(n)%flagstruct%reproduce_sum, kappa, cp_air, zvir,  &
                     Atm(n)%ptop, Atm(n)%ks, Atm(n)%ncnst, Atm(n)%flagstruct%n_split,        &
                     Atm(n)%flagstruct%q_split, Atm(n)%u0, Atm(n)%v0, Atm(n)%u,       &
                     Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%flagstruct%hydrostatic,  &
                     Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps,                     &
                     Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz, Atm(n)%phis,      &
                     Atm(n)%q_con, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc, &
                     Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy,                    &
                     Atm(n)%cx, Atm(n)%cy, Atm(n)%ze0, Atm(n)%flagstruct%hybrid_z,    &
                     Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct,         &
                     Atm(n)%thermostruct, Atm(n)%idiag, Atm(n)%bd, Atm(n)%parent_grid,  &
                     Atm(n)%domain, Atm(n)%inline_mp, Atm(n)%heat_source, Atm(n)%diss_est)
! Nudging back to IC
!$omp parallel do default(shared)
       do k=1,npz
          do j=jsc,jec+1
             do i=isc,iec
                Atm(n)%u(i,j,k) = xt*(Atm(n)%u(i,j,k) + wt*u0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec+1
                Atm(n)%v(i,j,k) = xt*(Atm(n)%v(i,j,k) + wt*v0(i,j,k))
             enddo
          enddo
          do j=jsc,jec
             do i=isc,iec
                Atm(n)%delp(i,j,k) = xt*(Atm(n)%delp(i,j,k) + wt*dp0(i,j,k))
             enddo
          enddo
       enddo

     call p_adi(npz, Atm(n)%ng, isc, iec, jsc, jec, Atm(n)%ptop,  &
                Atm(n)%delp, Atm(n)%ps, Atm(n)%pe,     &
                Atm(n)%peln, Atm(n)%pk, Atm(n)%pkz, Atm(n)%flagstruct%hydrostatic)

!$omp parallel do default(shared)
       do k=1,npz
          do j=jsc,jec
             do i=isc,iec
!               Atm(n)%pt(i,j,k) = xt*(Atm(n)%pt(i,j,k)+wt*t0(i,j,k)/((1.+esl*Atm(n)%q(i,j,k,1))*(Atm(n)%peln(i,k+1,j)-Atm(n)%peln(i,k,j))))
                Atm(n)%pt(i,j,k) = xt*(Atm(n)%pt(i,j,k)+wt*t0(i,j,k))
             enddo
          enddo
       enddo
     enddo

     deallocate ( u0 )
     deallocate ( v0 )
     deallocate ( t0 )
     deallocate (dp0 )

     do_adiabatic_init = .false.

 end subroutine adiabatic_init

!#######################################################################

  subroutine atmosphere (Time)
    type(time_type), intent(in) :: Time

    real:: zvir
    real:: time_total
    integer :: n, sphum, p, nc
    integer :: psc ! p_split counter

    !$ser verbatim integer :: mpi_rank,ier
    !$ser verbatim real :: bdt
    !$ser verbatim logical :: ser_on
    !$ser verbatim  call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank,ier)

    call timing_on('ATMOS_DYNAMICS')

    fv_time = Time + Time_step_atmos
    call get_time (fv_time, seconds,  days)

    time_total = days*SECONDS_PER_DAY + seconds

    !$ser verbatim save_step=fs_is_serialization_on()

    do psc=1,abs(p_split)

      !$ser verbatim if ((psc == abs(p_split)) .and. save_step) then
       !$ser verbatim save_step = .true.
       !$ser on
       !$ser savepoint Grid-Info
       !$ser data is_=Atm(n)%bd%is ie=Atm(n)%bd%ie isd=Atm(n)%bd%isd ied=Atm(n)%bd%ied js=Atm(n)%bd%js je=Atm(n)%bd%je jsd=Atm(n)%bd%jsd jed=Atm(n)%bd%jed
       !$ser data npx=npx npy=npy npz=npz nested=Atm(n)%gridstruct%nested grid_type=Atm(n)%gridstruct%grid_type dya=Atm(n)%gridstruct%dya
       !$ser data dxa=Atm(n)%gridstruct%dxa dxc=Atm(n)%gridstruct%dxc dyc=Atm(n)%gridstruct%dyc rdxc=Atm(n)%gridstruct%rdxc rdx=Atm(n)%gridstruct%rdx
       !$ser data rdyc=Atm(n)%gridstruct%rdyc  rdxa=Atm(n)%gridstruct%rdxa  rdya=Atm(n)%gridstruct%rdya rdy=Atm(n)%gridstruct%rdy
       !$ser data cosa_u=Atm(n)%gridstruct%cosa_u cosa_v=Atm(n)%gridstruct%cosa_v  cosa_s=Atm(n)%gridstruct%cosa_s sina_v=Atm(n)%gridstruct%sina_v
       !$ser data sina_u=Atm(n)%gridstruct%sina_u rsin_u=Atm(n)%gridstruct%rsin_u rsin_v=Atm(n)%gridstruct%rsin_v rsin2=Atm(n)%gridstruct%rsin2
       !$ser data sin_sg=Atm(n)%gridstruct%sin_sg cos_sg=Atm(n)%gridstruct%cos_sg area=Atm(n)%gridstruct%area dy=Atm(n)%gridstruct%dy
       !$ser data rarea=Atm(n)%gridstruct%rarea  rarea_c=Atm(n)%gridstruct%rarea_c area_64=Atm(n)%gridstruct%area_64 rsina=Atm(n)%gridstruct%rsina
       !$ser data cosa=Atm(n)%gridstruct%cosa  dx=Atm(n)%gridstruct%dx fC=Atm(n)%gridstruct%fC da_min=Atm(n)%gridstruct%da_min
       !$ser data da_min_c=Atm(n)%gridstruct%da_min_c del6_u=Atm(n)%gridstruct%del6_u del6_v=Atm(n)%gridstruct%del6_v f0=Atm(n)%gridstruct%f0
       !$ser data divg_u=Atm(n)%gridstruct%divg_u divg_v=Atm(n)%gridstruct%divg_v stretched_grid=Atm(n)%gridstruct%stretched_grid
       !$ser data agrid1=Atm(n)%gridstruct%agrid(:,:,1) agrid2=Atm(n)%gridstruct%agrid(:,:,2) bgrid1=Atm(n)%gridstruct%grid(:,:,1)
       !$ser data bgrid2=Atm(n)%gridstruct%grid(:,:,2) edge_w=Atm(n)%gridstruct%edge_w(jsc:jec+1) edge_e=Atm(n)%gridstruct%edge_e(jsc:jec+1)
       !$ser data edge_s=Atm(n)%gridstruct%edge_s(isc:iec+1) edge_n=Atm(n)%gridstruct%edge_n(isc:iec+1) a11=Atm(n)%gridstruct%a11
       !$ser data a12=Atm(n)%gridstruct%a12  a21=Atm(n)%gridstruct%a21  a22=Atm(n)%gridstruct%a22 ak=Atm(n)%ak bk=Atm(n)%bk ptop=Atm(n)%ptop ks=Atm(n)%ks
       !$ser data ee1=Atm(n)%gridstruct%ee1 ee2=Atm(n)%gridstruct%ee2 ew1=Atm(n)%gridstruct%ew(:,:,:,1) ew2=Atm(n)%gridstruct%ew(:,:,:,2)
       !$ser data  es1=Atm(n)%gridstruct%es(:,:,:,1) es2=Atm(n)%gridstruct%es(:,:,:,2)
       !$ser data vlon=Atm(n)%gridstruct%vlon vlat=Atm(n)%gridstruct%vlat
       !$ser data edge_vect_w=Atm(n)%gridstruct%edge_vect_w edge_vect_e=Atm(n)%gridstruct%edge_vect_e edge_vect_s=Atm(n)%gridstruct%edge_vect_s edge_vect_n=Atm(n)%gridstruct%edge_vect_n
      !$ser verbatim else
       !$ser verbatim save_step = .false.
       !$ser off
      !$ser verbatim endif

    n=mygrid

    call mpp_set_current_pelist(Atm(n)%pelist, no_sync=.TRUE.)

       if ( Atm(n)%flagstruct%nudge_ic )     &
            call  fv_nudge(Atm(n)%npz, Atm(n)%bd%isc, Atm(n)%bd%iec, Atm(n)%bd%jsc, Atm(n)%bd%jec, Atm(n)%ng, &
            Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz, Atm(n)%delp, Atm(n)%pt, dt_atmos/real(abs(p_split)), Atm(n)%flagstruct%hydrostatic )

    !---- call fv dynamics -----
    !if ( Atm(n)%flagstruct%adiabatic .or. Atm(n)%flagstruct%do_Held_Suarez ) then
    if ( Atm(n)%flagstruct%adiabatic ) then
       zvir = 0.         ! no virtual effect
    else
       zvir = rvgas/rdgas - 1.
    endif

       call timing_on('FV_DYNAMICS')
       call fv_dynamics(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%ncnst, Atm(n)%ng,   &
            dt_atmos/real(abs(p_split)), Atm(n)%flagstruct%consv_te, Atm(n)%flagstruct%fill, &
            Atm(n)%flagstruct%reproduce_sum, kappa,   &
            cp_air, zvir, Atm(n)%ptop, Atm(n)%ks, Atm(n)%ncnst, &
            Atm(n)%flagstruct%n_split, Atm(n)%flagstruct%q_split, &
            Atm(n)%u0, Atm(n)%v0, Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%delz, &
            Atm(n)%flagstruct%hydrostatic, Atm(n)%pt, Atm(n)%delp, Atm(n)%q, Atm(n)%ps, &
            Atm(n)%pe, Atm(n)%pk, Atm(n)%peln, Atm(n)%pkz,             &
            Atm(n)%phis, Atm(n)%q_con, Atm(n)%omga, Atm(n)%ua, Atm(n)%va, Atm(n)%uc, Atm(n)%vc,  &
            Atm(n)%ak, Atm(n)%bk, Atm(n)%mfx, Atm(n)%mfy, Atm(n)%cx, Atm(n)%cy,    &
            Atm(n)%ze0, Atm(n)%flagstruct%hybrid_z, Atm(n)%gridstruct, Atm(n)%flagstruct, &
            Atm(n)%neststruct, Atm(n)%thermostruct, Atm(n)%idiag, Atm(n)%bd, Atm(n)%parent_grid, Atm(n)%domain, &
            Atm(n)%inline_mp, Atm(n)%heat_source, Atm(n)%diss_est, time_total=time_total)
       call timing_off('FV_DYNAMICS')

    if (ngrids > 1 .and. (psc < p_split .or. p_split < 0)) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir, fv_time, mygrid)
       call timing_off('TWOWAY_UPDATE')
    endif

    end do !p_split

    if(Atm(n)%npz /=1 .and. (.not. Atm(n)%flagstruct%adiabatic .or. Atm(mygrid)%flagstruct%fv_sg_adj > 0.)) then

           call timing_on('FV_PHYS')
           call fv_phys(Atm(n)%npx, Atm(n)%npy, Atm(n)%npz, Atm(n)%bd%isc, Atm(n)%bd%iec, &
                        Atm(n)%bd%jsc, Atm(n)%bd%jec, Atm(n)%ng, Atm(n)%ncnst,                &
                        Atm(n)%u, Atm(n)%v, Atm(n)%w, Atm(n)%pt, Atm(n)%q, Atm(n)%pe,   &
                        Atm(n)%delp, Atm(n)%peln, Atm(n)%pkz, dt_atmos,                 &
                        Atm(n)%ua, Atm(n)%va, Atm(n)%phis, Atm(n)%gridstruct%agrid,     &
                        Atm(n)%ptop, Atm(n)%ak, Atm(n)%bk, Atm(n)%ks, Atm(n)%ps, Atm(n)%pk, &
                        Atm(n)%u_srf, Atm(n)%v_srf,  Atm(n)%ts, Atm(n)%delz,            &
                        Atm(n)%flagstruct%hydrostatic, Atm(n)%oro, .false.,             &
                        Atm(n)%flagstruct%p_ref, Atm(n)%flagstruct%do_Held_Suarez,      &
                        Atm(n)%gridstruct, Atm(n)%flagstruct, Atm(n)%neststruct,        &
                        Atm(n)%thermostruct, Atm(n)%flagstruct%nwat, Atm(n)%bd,         &
                        Atm(n)%domain, fv_time, Atm(n)%phys_diag, Atm(n)%nudge_diag,    &
                        Atm(n)%sg_diag, time_total)
           call timing_off('FV_PHYS')
    endif

    if (ngrids > 1 .and. p_split > 0) then
       call timing_on('TWOWAY_UPDATE')
       call twoway_nesting(Atm, ngrids, grids_on_this_pe, zvir, fv_time, mygrid)
       call timing_off('TWOWAY_UPDATE')
    endif


  !---- diagnostics for FV dynamics -----


   !For correct diagnostics (may need to be changed for moist Held-Suarez)
   if ( Atm(n)%flagstruct%adiabatic .or. Atm(n)%flagstruct%do_Held_Suarez ) then
       zvir = 0.         ! no virtual effect
    else
       zvir = rvgas/rdgas - 1.
    endif

    call timing_on('FV_DIAG')
    call fv_diag(Atm(n:n), zvir, fv_time, Atm(n)%flagstruct%print_freq)
    call timing_off('FV_DIAG')

    call timing_off('ATMOS_DYNAMICS')
    call mpp_set_current_pelist()

 end subroutine atmosphere


 subroutine atmosphere_end

   integer n

    call timing_on('ATMOS_END')

    call get_time (fv_time, seconds,  days)

    if ( Atm(mygrid)%flagstruct%moist_phys .and. Atm(mygrid)%flagstruct%nwat==6 ) call gfdl_mp_end

    call fv_end(Atm, mygrid)
    deallocate(Atm)

    call timing_off('ATMOS_END')
    call timing_off('ATMOS_TOTAL')
    call timing_prt( mpp_pe() )
    call mpp_set_current_pelist()

  end subroutine atmosphere_end

 subroutine atmosphere_domain ( fv_domain )
 type(domain2d), intent(out) :: fv_domain

!  returns the domain2d variable associated with the coupling grid
!  note: coupling is done using the mass/temperature grid with no halos

   fv_domain = Atm(mygrid)%domain

 end subroutine atmosphere_domain

 subroutine p_adi(km, ng, ifirst, ilast, jfirst, jlast, ptop,   &
                  delp, ps, pe, peln, pk, pkz, hydrostatic)

! Given (ptop, delp) computes (ps, pk, pe, peln, pkz)
! Input:
   integer,  intent(in):: km, ng
   integer,  intent(in):: ifirst, ilast            ! Longitude strip
   integer,  intent(in):: jfirst, jlast            ! Latitude strip
   logical, intent(in)::  hydrostatic
   real, intent(in):: ptop
   real, intent(in):: delp(ifirst-ng:ilast+ng,jfirst-ng:jlast+ng, km)
! Output:
   real, intent(out) ::   ps(ifirst-ng:ilast+ng, jfirst-ng:jlast+ng)
   real, intent(out) ::   pk(ifirst:ilast, jfirst:jlast, km+1)
   real, intent(out) ::   pe(ifirst-1:ilast+1,km+1,jfirst-1:jlast+1) ! Ghosted Edge pressure
   real, intent(out) :: peln(ifirst:ilast, km+1, jfirst:jlast)    ! Edge pressure
   real, intent(out) ::  pkz(ifirst:ilast, jfirst:jlast, km)
! Local
   real pek
   integer i, j, k

   pek = ptop ** kappa

!$OMP parallel do default(none) shared(ifirst,ilast,jfirst,jlast,km,ptop,pek,pe,pk, &
!$OMP                                  ps,delp,peln,hydrostatic,pkz)
   do j=jfirst,jlast
      do i=ifirst,ilast
         pe(i,1,j) = ptop
         pk(i,j,1) = pek
      enddo

      do k=2,km+1
         do i=ifirst,ilast
            pe(i,k,j) = pe(i,k-1,j) + delp(i,j,k-1)
            peln(i,k,j) = log(pe(i,k,j))
            pk(i,j,k) = exp( kappa*peln(i,k,j) )
         enddo
      enddo

      do i=ifirst,ilast
         ps(i,j) = pe(i,km+1,j)
      enddo

      if ( hydrostatic ) then
         do k=1,km
            do i=ifirst,ilast
               pkz(i,j,k) = (pk(i,j,k+1)-pk(i,j,k))/(kappa*(peln(i,k+1,j)-peln(i,k,j)))
            enddo
         enddo
      endif
   enddo

 end subroutine p_adi
end module atmosphere_mod
