!=======================================================================

! Aerosol tracer within sea ice
!
! authors Marika Holland, NCAR
!         David Bailey, NCAR

      module icepack_aerosol

      use icepack_kinds
      use icepack_parameters, only: c0, c1, c2, p5, puny, rhoi, rhos, hs_min
      use icepack_parameters, only: hi_ssl, hs_ssl, hs_ssl_min
      use icepack_tracers, only: max_aero, nilyr, nslyr, nblyr, ntrcr, nbtrcr, n_aero
      use icepack_warnings, only: warnstr, icepack_warnings_add
      use icepack_warnings, only: icepack_warnings_setabort, icepack_warnings_aborted

      use icepack_zbgc_shared, only: kscavz

      implicit none

      private
      public :: update_aerosol, update_snow_bgc

!=======================================================================

      contains

!=======================================================================

!  Deposition and vertical cycling of aerosol in ice or snow
!  Called from icepack_step_therm1 when tr_aero=T (not used for zbgc tracers)

      subroutine update_aerosol(dt,                   &
                                meltt,    melts,      &
                                meltb,    congel,     &
                                snoice,               &
                                fsnow,                &
                                aerosno,  aeroice,    &
                                aice_old,             &
                                vice_old, vsno_old,   &
                                vicen, vsnon, aicen,  &
                                faero_atm, faero_ocn)

      real (kind=dbl_kind), intent(in) :: &
         dt,       & ! time step
         meltt,    & ! thermodynamic melt/growth rates
         melts,    &
         meltb,    &
         congel,   &
         snoice,   &
         fsnow,    &
         vicen,    & ! ice volume (m)
         vsnon,    & ! snow volume (m)
         aicen,    & ! ice area fraction
         aice_old, & ! values prior to thermodynamic changes
         vice_old, &
         vsno_old

      real (kind=dbl_kind), dimension(:), intent(in) :: &
         faero_atm   ! aerosol deposition rate (W/m^2 s)

      real (kind=dbl_kind), dimension(:), intent(inout) :: &
         faero_ocn   ! aerosol flux to ocean (W/m^2 s)

      real (kind=dbl_kind), dimension(:,:), intent(inout) :: &
         aerosno,  aeroice    ! kg/m^2

      !  local variables
      integer (kind=int_kind) :: k

      real (kind=dbl_kind) :: &
         dzssl,  dzssl_new,      & ! snow ssl thickness
         dzint,                  & ! snow interior thickness
         dzssli, dzssli_new,     & ! ice ssl thickness
         dzinti,                 & ! ice interior thickness
         dznew,                  & ! tracks thickness changes
         hs, hi,                 & ! snow/ice thickness (m)
         dhs_evap, dhi_evap,     & ! snow/ice thickness change due to evap
         dhs_melts, dhi_meltt,   & ! ... due to surface melt
         dhs_snoice, dhi_snoice, & ! ... due to snow-ice formation
         dhi_congel, dhi_meltb,  & ! ... due to bottom growth, melt
         hslyr, hilyr,           & ! snow, ice layer thickness (m)
         hslyr_old, hilyr_old,   & ! old snow, ice layer thickness (m)
         hs_old, hi_old,         & ! old snow, ice thickness (m)
         sloss1, sloss2,         & ! aerosol mass loss (kg/m^2)
         ar                        ! 1/aicen(i,j)

      real (kind=dbl_kind), dimension(max_aero) :: &
         kscav, kscavsi       ! scavenging by melt water

      real (kind=dbl_kind), dimension(n_aero) :: &
         aerotot, aerotot0, & ! for conservation check
         focn_old             ! for conservation check

      ! echmod:  this assumes max_aero=6
      data kscav   / .03_dbl_kind, .20_dbl_kind, .02_dbl_kind, &
                     .02_dbl_kind, .01_dbl_kind, .01_dbl_kind /
      data kscavsi / .03_dbl_kind, .20_dbl_kind, .02_dbl_kind, &
                     .02_dbl_kind, .01_dbl_kind, .01_dbl_kind /

     character(len=*),parameter :: subname='(update_aerosol)'

    !-------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------
      focn_old(:) = faero_ocn(:)

      hs_old    = vsno_old/aice_old
      hi_old    = vice_old/aice_old
      hslyr_old = hs_old/real(nslyr,kind=dbl_kind)
      hilyr_old = hi_old/real(nilyr,kind=dbl_kind)

      dzssl  = min(hslyr_old/c2, hs_ssl)
      dzssli = min(hilyr_old/c2, hi_ssl)
      dzint  = hs_old - dzssl
      dzinti = hi_old - dzssli

      if (aicen > c0) then
         ar = c1/aicen
      else ! ice disappeared during time step
         ar = c1/aice_old
      endif

      hs = vsnon*ar
      hi = vicen*ar
      ! fluxes were divided by aice for thermo, not yet multiplied by aice
      dhs_melts  = -melts
      dhi_snoice = snoice
      dhs_snoice = dhi_snoice*rhoi/rhos
      dhi_meltt  = -meltt
      dhi_meltb  = -meltb
      dhi_congel = congel

      dhs_evap = hs - (hs_old + dhs_melts - dhs_snoice &
                              + fsnow/rhos*dt)
      dhi_evap = hi - (hi_old + dhi_meltt + dhi_meltb &
                              + dhi_congel + dhi_snoice)

      do k = 1, n_aero
         aerotot0(k) = aerosno(k,2) + aerosno(k,1) &
                     + aeroice(k,2) + aeroice(k,1)
      enddo

    !-------------------------------------------------------------------
    ! evaporation
    !-------------------------------------------------------------------
      dzint  = dzint  + min(dzssl  + dhs_evap, c0)
      dzinti = dzinti + min(dzssli + dhi_evap, c0)
      dzssl  = max(dzssl  + dhs_evap, c0)
      dzssli = max(dzssli + dhi_evap, c0)

    !-------------------------------------------------------------------
    ! basal ice growth
    !-------------------------------------------------------------------
      dzinti = dzinti + dhi_congel

    !-------------------------------------------------------------------
    ! surface snow melt
    !-------------------------------------------------------------------
      if (-dhs_melts > puny) then
         do k = 1, n_aero
            sloss1 = c0
            sloss2 = c0
            if (dzssl > puny)  &
                 sloss1 = kscav(k)*aerosno(k,1)  &
                                  *min(-dhs_melts,dzssl)/dzssl
            aerosno(k,1) = aerosno(k,1) - sloss1
            if (dzint > puny)  &
                 sloss2 = kscav(k)*aerosno(k,2) &
                                  *max(-dhs_melts-dzssl,c0)/dzint
            aerosno(k,2) = aerosno(k,2) - sloss2
            faero_ocn(k) = faero_ocn(k) + (sloss1+sloss2)/dt
         enddo  ! n_aero

         ! update snow thickness
         dzint=dzint+min(dzssl+dhs_melts, c0)
         dzssl=max(dzssl+dhs_melts, c0)

         if ( dzssl <= puny ) then ! ssl melts away
            aerosno(:,2) = aerosno(:,1) + aerosno(:,2)
            aerosno(:,1) = c0
            dzssl = max(dzssl, c0)
         endif
         if (dzint <= puny ) then  ! all snow melts away
            aeroice(:,1) = aeroice(:,1) &
                         + aerosno(:,1) + aerosno(:,2)
            aerosno(:,:) = c0
            dzint = max(dzint, c0)
         endif
      endif

    !-------------------------------------------------------------------
    ! surface ice melt
    !-------------------------------------------------------------------
      if (-dhi_meltt > puny) then
         do k = 1, n_aero
            sloss1 = c0
            sloss2 = c0
            if (dzssli > puny)  &
                 sloss1 = kscav(k)*aeroice(k,1)  &
                                  *min(-dhi_meltt,dzssli)/dzssli
            aeroice(k,1) = aeroice(k,1) - sloss1
            if (dzinti > puny)  &
                 sloss2 = kscav(k)*aeroice(k,2)  &
                                  *max(-dhi_meltt-dzssli,c0)/dzinti
            aeroice(k,2) = aeroice(k,2) - sloss2
            faero_ocn(k) = faero_ocn(k) + (sloss1+sloss2)/dt
         enddo

         dzinti = dzinti + min(dzssli+dhi_meltt, c0)
         dzssli = max(dzssli+dhi_meltt, c0)
         if (dzssli <= puny) then   ! ssl ice melts away
            do k = 1, n_aero
               aeroice(k,2) = aeroice(k,1) + aeroice(k,2)
               aeroice(k,1) = c0
            enddo
            dzssli = max(dzssli, c0)
         endif
         if (dzinti <= puny) then   ! all ice melts away
            do k = 1, n_aero
               faero_ocn(k) = faero_ocn(k)  &
                            + (aeroice(k,1)+aeroice(k,2))/dt
               aeroice(k,:)=c0
            enddo
            dzinti = max(dzinti, c0)
         endif
      endif

    !-------------------------------------------------------------------
    ! basal ice melt.  Assume all aero lost in basal melt
    !-------------------------------------------------------------------
      if (-dhi_meltb > puny) then
         do k=1,n_aero
            sloss1=c0
            sloss2=c0
            if (dzssli > puny)  &
                 sloss1 = max(-dhi_meltb-dzinti, c0)  &
                        *aeroice(k,1)/dzssli
            aeroice(k,1) = aeroice(k,1) - sloss1
            if (dzinti > puny)  &
                 sloss2 = min(-dhi_meltb, dzinti)  &
                        *aeroice(k,2)/dzinti
            aeroice(k,2) = aeroice(k,2) - sloss2
            faero_ocn(k) = faero_ocn(k) + (sloss1+sloss2)/dt
         enddo

         dzssli = dzssli + min(dzinti+dhi_meltb, c0)
         dzinti = max(dzinti+dhi_meltb, c0)
      endif

    !-------------------------------------------------------------------
    ! snowfall
    !-------------------------------------------------------------------
      if (fsnow > c0) dzssl = dzssl + fsnow/rhos*dt

    !-------------------------------------------------------------------
    ! snow-ice formation
    !-------------------------------------------------------------------
      if (dhs_snoice > puny) then
         do k = 1, n_aero
            sloss1 = c0
            sloss2 = c0
            if (dzint > puny)  &
                 sloss2 = min(dhs_snoice, dzint)  &
                        *aerosno(k,2)/dzint
            aerosno(k,2) = aerosno(k,2) - sloss2
            if (dzssl > puny)  &
                 sloss1 = max(dhs_snoice-dzint, c0)  &
                        *aerosno(k,1)/dzssl
            aerosno(k,1) = aerosno(k,1) - sloss1
            aeroice(k,1) = aeroice(k,1) &
                         + (c1-kscavsi(k))*(sloss2+sloss1)
            faero_ocn(k) = faero_ocn(k) &
                         + kscavsi(k)*(sloss2+sloss1)/dt
         enddo
         dzssl  = dzssl - max(dhs_snoice-dzint, c0)
         dzint  = max(dzint-dhs_snoice, c0)
         dzssli = dzssli + dhi_snoice
      endif

    !-------------------------------------------------------------------
    ! aerosol deposition
    !-------------------------------------------------------------------
      if (aicen > c0) then
         hs = vsnon * ar
      else
         hs = c0
      endif
      if (hs > hs_min) then    ! should this really be hs_min or 0?
         ! should use same hs_min value as in radiation
         do k=1,n_aero
            aerosno(k,1) = aerosno(k,1) &
                         + faero_atm(k)*dt*aicen
         enddo
      else
         do k=1,n_aero
            aeroice(k,1) = aeroice(k,1) &
                         + faero_atm(k)*dt*aicen
         enddo
      endif

    !-------------------------------------------------------------------
    ! redistribute aerosol within vertical layers
    !-------------------------------------------------------------------
      if (aicen > c0) then
         hs = vsnon  * ar     ! new snow thickness
         hi = vicen  * ar     ! new ice thickness
      else
         hs = c0
         hi = c0
      endif
      if (dzssl <= puny) then   ! nothing in SSL
         do k=1,n_aero
            aerosno(k,2) = aerosno(k,2) + aerosno(k,1)
            aerosno(k,1) = c0
         enddo
      endif
      if (dzint <= puny) then   ! nothing in Snow Int
         do k = 1, n_aero
            aeroice(k,1) = aeroice(k,1) + aerosno(k,2)
            aerosno(k,2) = c0
         enddo
      endif
      if (dzssli <= puny) then  ! nothing in Ice SSL
         do k = 1, n_aero
            aeroice(k,2) = aeroice(k,2) + aeroice(k,1)
            aeroice(k,1) = c0
         enddo
      endif

      if (dzinti <= puny) then  ! nothing in Ice INT
         do k = 1, n_aero
            faero_ocn(k) = faero_ocn(k) &
                         + (aeroice(k,1)+aeroice(k,2))/dt
            aeroice(k,:)=c0
         enddo
      endif

      hslyr      = hs/real(nslyr,kind=dbl_kind)
      hilyr      = hi/real(nilyr,kind=dbl_kind)
      dzssl_new  = min(hslyr/c2, hs_ssl)
      dzssli_new = min(hilyr/c2, hi_ssl)

      if (hs > hs_min) then
         do k = 1, n_aero
            dznew = min(dzssl_new-dzssl, c0)
            sloss1 = c0
            if (dzssl > puny) &
                 sloss1 = dznew*aerosno(k,1)/dzssl ! not neccesarily a loss
            dznew = max(dzssl_new-dzssl, c0)
            if (dzint > puny) &
                 sloss1 = sloss1 + aerosno(k,2)*dznew/dzint
            aerosno(k,1) = aerosno(k,1) + sloss1
            aerosno(k,2) = aerosno(k,2) - sloss1
         enddo
      else
         aeroice(:,1) = aeroice(:,1)  &
                      + aerosno(:,1) + aerosno(:,2)
         aerosno(:,:) = c0
      endif

      if (vicen > puny) then ! may want a limit on hi instead?
         do k = 1, n_aero
            sloss2 = c0
            dznew = min(dzssli_new-dzssli, c0)
            if (dzssli > puny) &
                 sloss2 = dznew*aeroice(k,1)/dzssli
            dznew = max(dzssli_new-dzssli, c0)
            if (dzinti > puny) &
                 sloss2 = sloss2 + aeroice(k,2)*dznew/dzinti
            aeroice(k,1) = aeroice(k,1) + sloss2
            aeroice(k,2) = aeroice(k,2) - sloss2
         enddo
      else
         faero_ocn(:) = faero_ocn(:) + (aeroice(:,1)+aeroice(:,2))/dt
         aeroice(:,:) = c0
      endif

    !-------------------------------------------------------------------
    ! check conservation
    !-------------------------------------------------------------------
      do k = 1, n_aero
         aerotot(k) = aerosno(k,2) + aerosno(k,1) &
                    + aeroice(k,2) + aeroice(k,1)
         if ((aerotot(k)-aerotot0(k)) &
              - (   faero_atm(k)*aicen &
              - (faero_ocn(k)-focn_old(k)) )*dt  > puny) then

            write(warnstr,*) subname, 'aerosol tracer:  ',k
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'aerotot-aerotot0 ',aerotot(k)-aerotot0(k)
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'faero_atm-faero_ocn      ', &
                 (faero_atm(k)*aicen-(faero_ocn(k)-focn_old(k)))*dt
            call icepack_warnings_add(warnstr)
         endif
      enddo

    !-------------------------------------------------------------------
    ! check for negative values
    !-------------------------------------------------------------------

!echmod:  note that this does not test or fix all aero tracers
      if (aeroice(1,1) < -puny .or. &
          aeroice(1,2) < -puny .or. &
          aerosno(1,1) < -puny .or. &
          aerosno(1,2) < -puny) then

         write(warnstr,*) subname, 'aerosol negative in aerosol code'
         call icepack_warnings_add(warnstr)

         aeroice(1,1) = max(aeroice(1,1), c0)
         aeroice(1,2) = max(aeroice(1,2), c0)
         aerosno(1,1) = max(aerosno(1,1), c0)
         aerosno(1,2) = max(aerosno(1,2), c0)

      endif

      end subroutine update_aerosol

!=======================================================================

!  Aerosol in snow for vertical biogeochemistry with mushy thermodynamics
!  Called from icepack_algae.F90 when z_tracers=T (replaces update_aerosol)

      subroutine update_snow_bgc(dt,                   &
                                meltt,    melts,       &
                                meltb,    congel,      &
                                snoice,   fsnow,       &
                                trcrn,    bio_index,   &
                                aice_old, zbgc_snow,   &
                                vice_old, vsno_old,    &
                                vicen,    vsnon,       &
                                aicen,    flux_bio_atm,&
                                zbgc_atm, flux_bio,    &
                                bio_index_o)

      integer (kind=int_kind), dimension (nbtrcr), intent(in) :: &
         bio_index,  &
         bio_index_o   ! provides index of scavenging (kscavz) data array

      real (kind=dbl_kind), intent(in) :: &
         dt                    ! time step

      real (kind=dbl_kind), intent(in) :: &
         meltt,    & ! thermodynamic melt/growth rates
         melts,    &
         meltb,    &
         congel,   &
         snoice,   &
         fsnow,    &
         vicen,    & ! ice volume (m)
         vsnon,    & ! snow volume (m)
         aicen,    & ! ice area fraction
         aice_old, & ! values prior to thermodynamic changes
         vice_old, &
         vsno_old

      real (kind=dbl_kind), dimension(nbtrcr), intent(out) :: &
         zbgc_snow, & ! aerosol contribution from snow to ice
         zbgc_atm     ! and atm to ice concentration * volume (kg or mmol/m^3*m)

      real (kind=dbl_kind),dimension(nbtrcr), intent(inout) :: &
         flux_bio     ! total ocean tracer flux (mmol/m^2/s)

      real (kind=dbl_kind), dimension(nbtrcr), intent(in) :: &
         flux_bio_atm   ! aerosol deposition rate (kg or mmol/m^2 s)

      real (kind=dbl_kind), dimension(ntrcr), intent(inout) :: &
         trcrn       ! ice/snow tracer array

      !  local variables

      integer (kind=int_kind) ::  k, n

      real (kind=dbl_kind) :: &
         dzssl,  dzssl_new,      & ! snow ssl thickness
         dzint,  dzint_new,      & ! snow interior thickness
         dz,                     & !
         hi,                     & ! ice  thickness (m)
         hilyr,                  & ! ice layer thickness (m)
         hs,                     & ! snow thickness (m)
         dhs_evap,               & ! snow thickness change due to evap
         dhs_melts,              & ! ... due to surface melt
         dhs_snoice,             & ! ... due to snow-ice formation
         hslyr,                  & ! snow layer thickness (m)
         hslyr_old,              & ! old snow layer thickness (m)
         hs_old,                 & ! old snow thickness (m)
         dznew,                  & ! change in the snow sl (m)
         sloss1, sloss2,         & ! aerosol mass loss (kg/m^2)
         ar                        ! 1/aicen(i,j)

      real (kind=dbl_kind), dimension(nbtrcr) :: &
         aerotot, aerotot0, & ! for conservation check (mmol/m^3)
         aero_cons        , & ! for conservation check (mmol/m^2)
         flux_bio_o           ! initial ocean tracer flux (mmol/m^2/s)

      real (kind=dbl_kind), dimension(nbtrcr,2) :: &
         aerosno,  & ! kg/m^2
         aerosno0    ! for diagnostic prints

      character(len=*),parameter :: subname='(update_snow_bgc)'

    !-------------------------------------------------------------------
    ! initialize
    !-------------------------------------------------------------------
         aerosno (:,:) = c0
         aerosno0(:,:) = c0
         aero_cons(:) = c0
         zbgc_snow(:) = c0
         zbgc_atm(:) = c0

         hs_old    = vsno_old/aice_old
         if (aice_old .gt. puny) then
            hs_old    = vsno_old/aice_old
         else
            hs_old = c0
         end if
         hslyr_old = hs_old/real(nslyr,kind=dbl_kind)

         dzssl  = hslyr_old/c2
         dzint  = hs_old - dzssl

         if (aicen > c0) then
            ar = c1/aicen
            hs = vsnon*ar
            hi = vicen*ar
         else ! ice disappeared during time step
            ar = c1
            hs = c0
            hi = c0
            if (aice_old > c0) hs = vsnon/aice_old
         endif
         hilyr = hi/real(nblyr,kind=dbl_kind)
         hslyr = hs/real(nslyr,kind=dbl_kind)
         dzssl_new  = hslyr/c2
         dhs_melts  = -melts
         dhs_snoice = snoice*rhoi/rhos
         dhs_evap = hs - (hs_old + dhs_melts - dhs_snoice &
                                 + fsnow/rhos*dt)

         ! trcrn() has units kg/m^3

      if (dzssl_new .lt. hs_ssl_min) then ! Put atm BC/dust flux directly into the sea ice
         do k=1,nbtrcr
            flux_bio_o(k) = flux_bio(k)
            if (hilyr .lt. hs_ssl_min) then
               flux_bio(k) = flux_bio(k) +  &
                         (trcrn(bio_index(k)+ nblyr+1)*dzssl+ &
                          trcrn(bio_index(k)+ nblyr+2)*dzint)/dt
               flux_bio(k) = flux_bio(k) + flux_bio_atm(k)
            else
               zbgc_snow(k) = zbgc_snow(k) +  &
                         (trcrn(bio_index(k)+ nblyr+1)*dzssl+ &
                          trcrn(bio_index(k)+ nblyr+2)*dzint)
               zbgc_atm(k) = zbgc_atm(k) &
                                + flux_bio_atm(k)*dt
            end if
            trcrn(bio_index(k) + nblyr+1) = c0
            trcrn(bio_index(k) + nblyr+2) = c0
         enddo

      else

         do k=1,nbtrcr
            flux_bio_o(k) = flux_bio(k)
            aerosno (k,1) = trcrn(bio_index(k)+ nblyr+1) * dzssl
            aerosno (k,2) = trcrn(bio_index(k)+ nblyr+2) * dzint
            aerosno0(k,:) = aerosno(k,:)
            aerotot0(k)   = aerosno(k,2) + aerosno(k,1)
         enddo

    !-------------------------------------------------------------------
    ! evaporation
    !-------------------------------------------------------------------
         dzint  = dzint  + min(dzssl  + dhs_evap, c0)
         dzssl  = max(dzssl  + dhs_evap, c0)
         if (dzssl <= puny) then
           do k = 1,nbtrcr
              aerosno(k,2)  = aerosno(k,2) + aerosno(k,1)
              aerosno(k,1) = c0
           end do
         end if
         if (dzint <= puny) then
           do k = 1,nbtrcr
              zbgc_snow(k) = zbgc_snow(k) + (aerosno(k,2) + aerosno(k,1))
              aerosno(k,2) = c0
              aerosno(k,1) = c0
           end do
         end if

    !------------------------------------------------------------------
    ! snowfall
    !-------------------------------------------------------------------
         if (fsnow > c0) then
           sloss1 = c0
           dz = min(fsnow/rhos*dt,dzssl)
           do k = 1, nbtrcr
              if (dzssl > puny) &
                 sloss1 = aerosno(k,1)*dz/dzssl
               aerosno(k,1) = max(c0,aerosno(k,1) - sloss1)
               aerosno(k,2) = aerosno(k,2) + sloss1
           end do
           dzssl = dzssl - dz + fsnow/rhos*dt
           dzint = dzint + dz
         end if
         if (dzssl <= puny) then
           do k = 1,nbtrcr
              aerosno(k,2)  = aerosno(k,2) + aerosno(k,1)
              aerosno(k,1) = c0
           end do
         end if
         if (dzint <= puny) then
           do k = 1,nbtrcr
              zbgc_snow(k) = zbgc_snow(k) + (aerosno(k,2) + aerosno(k,1))
              aerosno(k,2) = c0
              aerosno(k,1) = c0
           end do
         end if

    !-------------------------------------------------------------------
    ! surface snow melt
    !-------------------------------------------------------------------
         if (-dhs_melts > puny) then
            do k = 1, nbtrcr
               sloss1 = c0
               sloss2 = c0
               if (dzssl > puny) &
                  sloss1 = kscavz(bio_index_o(k))*aerosno(k,1)  &
                                    *min(-dhs_melts,dzssl)/dzssl
               aerosno(k,1) = max(c0,aerosno(k,1) - sloss1)
               if (dzint > puny) &
                   sloss2 = kscavz(bio_index_o(k))*aerosno(k,2) &
                                     *max(-dhs_melts-dzssl,c0)/dzint
               aerosno(k,2) = max(c0,aerosno(k,2) - sloss2)
               zbgc_snow(k) = zbgc_snow(k) + (sloss1+sloss2)  ! all not scavenged ends in ice
            enddo

            ! update snow thickness
            dzint=dzint+min(dzssl+dhs_melts, c0)
            dzssl=max(dzssl+dhs_melts, c0)

            if ( dzssl .le. puny ) then ! ssl melts away
               do k = 1,nbtrcr
                  aerosno(k,2) = aerosno(k,1) + aerosno(k,2)
                  aerosno(k,1) = c0
               end do
               dzssl = max(dzssl, c0)
            endif
            if (dzint .le. puny ) then  ! all snow melts away
               do k = 1,nbtrcr
                  zbgc_snow(k) = zbgc_snow(k) &
                                + aerosno(k,1) + aerosno(k,2)
                  aerosno(k,:) = c0
               enddo
               dzint = max(dzint, c0)
            endif
         endif ! -dhs_melts > puny

    !-------------------------------------------------------------------
    ! snow-ice formation
    !-------------------------------------------------------------------
         if (dhs_snoice > puny) then
            do k = 1, nbtrcr
               sloss1 = c0
               sloss2 = c0
               if (dzint > puny .and. aerosno(k,2) > c0) &
                   sloss2 = min(dhs_snoice, dzint)  &
                            *aerosno(k,2)/dzint
               aerosno(k,2) = max(c0,aerosno(k,2) - sloss2)
               if (dzssl > puny .and. aerosno(k,1) > c0)  &
                  sloss1 = max(dhs_snoice-dzint, c0)  &
                           *aerosno(k,1)/dzssl

               aerosno(k,1) = max(c0,aerosno(k,1) - sloss1)
               flux_bio(k)  = flux_bio(k) &
                               + kscavz(bio_index_o(k)) * (sloss2+sloss1)/dt
               zbgc_snow(k) = zbgc_snow(k) &
                               + (c1-kscavz(bio_index_o(k)))*(sloss2+sloss1)
            enddo
            dzssl  = max(c0,dzssl - max(dhs_snoice-dzint, c0))
            dzint  = max(dzint-dhs_snoice, c0)
         endif  ! dhs_snowice > puny

    !-------------------------------------------------------------------
    ! aerosol deposition
    !-------------------------------------------------------------------
         ! Spread out the atm dust flux in the snow interior for small snow surface layers
         if (dzssl .ge. hs_ssl*p5)  then

            do k=1,nbtrcr
               aerosno(k,1) = aerosno(k,1) &
                                + flux_bio_atm(k)*dt
            enddo
         else
            dz = (hs_ssl*p5 - dzssl)/(hs_ssl*p5)
            do k=1,nbtrcr
               aerosno(k,1) = aerosno(k,1) &
                                + flux_bio_atm(k)*dt*(c1-dz)
               aerosno(k,2) = aerosno(k,2) &
                                + flux_bio_atm(k)*dt*dz
            enddo
         endif

    !-------------------------------------------------------------------
    ! redistribute aerosol within vertical layers
    !-------------------------------------------------------------------
         if (aicen > c0) then
            hs = vsnon  * ar     ! new snow thickness
         else
            hs = c0
         endif
         if (dzssl <= puny) then   ! nothing in SSL
            do k=1,nbtrcr
               aerosno(k,2) = aerosno(k,2) + aerosno(k,1)
               aerosno(k,1) = c0
            enddo
         endif
         if (dzint <= puny) then   ! nothing in Snow Int
            do k = 1, nbtrcr
               zbgc_snow(k) = zbgc_snow(k) + max(c0,aerosno(k,2)+aerosno(k,1))
               aerosno(k,1) = c0
               aerosno(k,2) = c0
            enddo
         endif

         hslyr      = hs/real(nslyr,kind=dbl_kind)
         dzssl_new  = hslyr/c2
         dzint_new  = max(c0,hs - dzssl_new)

         if (hs > hs_min) then
            do k = 1, nbtrcr
               dznew = min(dzssl_new-dzssl, c0)
               sloss1 = c0
               if (dzssl > puny .and. aerosno(k,1) > c0) &
                  sloss1 = dznew*aerosno(k,1)/dzssl ! not neccesarily a loss
               dznew = max(dzssl_new-dzssl, c0)
               if (dzint > puny  .and. aerosno(k,2) > c0) &
                  sloss1 = aerosno(k,2)*dznew/dzint
               aerosno(k,1) = max(c0,aerosno(k,1) + sloss1)
               aerosno(k,2) = max(c0,aerosno(k,2) - sloss1)
            enddo
         else
            zbgc_snow(:) = zbgc_snow(:)  &
                             + aerosno(:,1) + aerosno(:,2)
            aerosno(:,:) = c0
         endif

    !-------------------------------------------------------------------
    ! check conservation
    !-------------------------------------------------------------------
         do k = 1, nbtrcr
            aerotot(k) = aerosno(k,2) + aerosno(k,1) &
                       + zbgc_snow(k) + zbgc_atm(k)
            aero_cons(k) = aerotot(k)-aerotot0(k) &
                          - (    flux_bio_atm(k) &
                          - (flux_bio(k)-flux_bio_o(k))) * dt
            if (aerotot0(k) > aerotot(k) .and. aerotot0(k) > c0) then
                 aero_cons(k) = aero_cons(k)/aerotot0(k)
            else if (aerotot(k) > c0) then
                 aero_cons(k) = aero_cons(k)/aerotot(k)
            end if
            if (aero_cons(k)  > puny .or. zbgc_snow(k) + zbgc_atm(k) < c0) then
               write(warnstr,*) subname, 'Conservation failure: aerosols in snow'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'test aerosol 1'
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'aerosol tracer:  ',k
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'aero_cons(k),puny:', aero_cons(k),puny
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'aerotot,aerotot0 ',aerotot(k),aerotot0(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, ' aerosno(k,2),aerosno(k,1) ', aerosno(k,2),aerosno(k,1)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'flux_bio_atm(k)*aicen*dt', &
                    flux_bio_atm(k)*aicen*dt
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zbgc_snow(k)', &
                    zbgc_snow(k)
               call icepack_warnings_add(warnstr)
               write(warnstr,*) subname, 'zbgc_atm(k)', &
                    zbgc_atm(k)
               call icepack_warnings_add(warnstr)
            endif
         enddo

    !-------------------------------------------------------------------
    ! reload tracers
    !-------------------------------------------------------------------
         if (dzssl_new > puny .and. dzint_new > puny .and. vsnon > puny) then
         do k = 1,nbtrcr
             trcrn(bio_index(k)+nblyr+1)=aerosno(k,1)/dzssl_new
             trcrn(bio_index(k)+nblyr+2)=aerosno(k,2)/dzint_new
         enddo
         else
         do k = 1,nbtrcr
            trcrn(bio_index(k)+nblyr+1)= c0
            trcrn(bio_index(k)+nblyr+2)= c0
         enddo
         endif

    !-------------------------------------------------------------------
    ! check for negative values
    !-------------------------------------------------------------------
         if (minval(aerosno(:,1)) < -puny  .or. &
            minval(aerosno(:,2)) < -puny) then

            write(warnstr,*) subname, 'Snow aerosol negative in update_snow_bgc'
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'aicen= '      ,aicen
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'vicen= '      ,vicen
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'vsnon= '      ,vsnon
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'viceold= '    ,vice_old
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'vsnoold= '    ,vsno_old
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'melts= '      ,melts
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'meltt= '      ,meltt
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'meltb= '      ,meltb
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'congel= '     ,congel
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'snoice= '     ,snoice
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'aero evap snow= '  ,dhs_evap
            call icepack_warnings_add(warnstr)
            write(warnstr,*) subname, 'fsnow= '      ,fsnow
            call icepack_warnings_add(warnstr)
            do k = 1, nbtrcr
              write(warnstr,*) subname, 'NBTRCR value k = ', k
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'aero snowssl (k)= '    ,aerosno0(k,1)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'aero new snowssl (k)= ',aerosno(k,1)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'aero snowint (k)= '    ,aerosno0(k,2)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'aero new snowint(k)= ',aerosno(k,2)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'flux_bio_atm(k)= ' , flux_bio_atm(k)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'zbgc_snow(k)= '  ,zbgc_snow(k)
              call icepack_warnings_add(warnstr)
              write(warnstr,*) subname, 'zbgc_atm(k)= '  ,zbgc_atm(k)
              call icepack_warnings_add(warnstr)

              do n = 1,2
                trcrn(bio_index(k)+nblyr+n)=max(trcrn(bio_index(k)+nblyr+n), c0)
              enddo
            enddo
         endif
        endif

      end subroutine update_snow_bgc

!=======================================================================

      end module icepack_aerosol

!=======================================================================
