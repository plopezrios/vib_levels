PROGRAM vib_levels
  !-------------------------------------------------------------!
  ! VIB_LEVELS                                                  !
  ! ==========                                                  !
  ! Given a one-dimensional potential V(x) at various values of !
  ! x, interpolate the potential and obtain the first few       !
  ! eigenfunctions and corresponding energy levels.             !
  ! PLR 05.2019                                                 !
  !-------------------------------------------------------------!
  IMPLICIT NONE

  call main()


CONTAINS


  SUBROUTINE main()
    !--------------!
    ! Main driver. !
    !--------------!
    IMPLICIT NONE
    ! Data in file.
    INTEGER nxy
    DOUBLE PRECISION, ALLOCATABLE :: xy_x(:), xy_y(:), xy_dy(:)
    ! System definitions and variables for fit to potential.
    INTEGER norder_v
    DOUBLE PRECISION ucentre, omega, rmass, fit_x0, fit_y0, xleft, xright, &
       &xmin, ymin
    DOUBLE PRECISION, ALLOCATABLE :: vcoeff_nat(:), vcoeff(:)
    ! Variables defining ground state of Hamiltonian.
    INTEGER, PARAMETER :: MAX_NORDER = 200 ! NB, dsyev* hang at norder=202
    DOUBLE PRECISION, PARAMETER :: VIRIAL_TOL = 1.d-9
    INTEGER norder
    DOUBLE PRECISION e0, vratio
    DOUBLE PRECISION, ALLOCATABLE :: orbcoeff(:)
    ! Random sampling and objects to evaluate.
    INTEGER irandom, iexpval
    INTEGER, PARAMETER :: nsample = 100, nexpval = 8
    CHARACTER(16), PARAMETER :: expval_name(nexpval) = (/&
       &'Delta0         ',&
       &'Delta1         ',&
       &'Delta2         ',&
       &'Delta3         ',&
       &'Virial ratio   ',&
       &'Orb. exp. order',&
       &'r of minimum   ',&
       &'Minimum energy '/)
    DOUBLE PRECISION mean, var
    DOUBLE PRECISION, ALLOCATABLE :: sample_y(:), fmc(:,:), fatx(:,:), &
       &fplot(:,:), fmeanx(:), w_one(:)
    ! Variables for plotting potential and eigenfunctions.
    INTEGER, PARAMETER :: PLOT_NPOINT = 200
    DOUBLE PRECISION, PARAMETER :: CDF_TOL = epsilon(1.d0)
    INTEGER ipoint, ixy, ipoly, iplot
    DOUBLE PRECISION xp, xl, xr, up, ul, ur, fp, fp_CI(4)
    DOUBLE PRECISION psix, plot_orb_scale_factor
    DOUBLE PRECISION, ALLOCATABLE :: all_eigval(:), all_eigvec(:,:), hbasis(:)
    ! Misc local variables.
    CHARACTER(2048) line, fname, pname
    INTEGER ncolumn, ierr
    DOUBLE PRECISION t1
    INTEGER, PARAMETER :: io=10

    ! Get name of data file.
    write(6,'(a)') 'Enter name of data file containing {r,E,dE} data:'
    read(5,'(a)',iostat=ierr) fname
    if (ierr/=0) call quit()
    write(6,'()')

    ! Load data.
    call check_file (fname, nxy, ncolumn)
    select case(nxy)
    case(-1)
      call quit('Could not open file "'//trim(fname)//'".')
    case(-2)
      call quit('Problem reading file "'//trim(fname)//'".')
    case(-3)
      call quit('Column count problem in file "'//trim(fname)//'".')
    case(0)
      call quit('File "'//trim(fname)//'" contains no useful data.')
    case default
      write(6,'(a)')'File "'//trim(fname)//'" contains '//trim(i2s(nxy))//&
         &' data lines.'
      write(6,'()')
    end select
    if (ncolumn<3) call quit( 'File "'//trim(fname)//'" must contain &
       &at least 3 data columns.')
    allocate (xy_x(nxy), xy_y(nxy), xy_dy(nxy))
    call read_file (fname, nxy, xy_x, xy_y, xy_dy, ierr)
    fit_x0 = xy_x(sum(minloc(xy_y)))
    fit_y0 = xy_y(sum(minloc(xy_y)))
    xleft = minval(xy_x)
    xright = maxval(xy_x)

    ! Get effective mass.
    write(6,'(a)') 'Effective mass in Da [e.g., 0.5 for H2; empty for 1 a.u.]:'
    read(5,'(a)',iostat=ierr) line
    if (ierr/=0) call quit()
    read(line,*,iostat=ierr) rmass
    if (ierr/=0) then
      rmass=1.d0
    else
      if (rmass<=0.d0) call quit ('Mass must be positive.')
      ! Convert to a.u.
      rmass = rmass*1822.888486209d0
    endif
    write(6,'()')

    ! Get expansion order.
    write(6,'(a)') 'Polynomial expansion order for potential [empty for 4 &
       &(= quartic)]:'
    read(5,'(a)',iostat=ierr) line
    if (ierr/=0) call quit()
    write(6,'()')
    read(line,*,iostat=ierr) norder_v
    if (ierr/=0) norder_v=4
    if (norder_v+1>nxy) call quit('Cannot have more parameters than data &
       &points.')
    if (norder_v+1==nxy) then
      write(6,'(a)') 'WARNING: using the same number of parameters as points &
         &is discouraged.'
      write(6,'()')
    endif
    allocate (vcoeff_nat(0:norder_v), vcoeff(0:norder_v))

    ! Get name of output plot file.
    write(6,'(a)') 'Enter name *root* for plot files [empty to skip plot]:'
    read(5,'(a)',iostat=ierr) pname
    if (ierr/=0) pname=''
    write(6,'()')

    ! Perform Monte Carlo resample.
    allocate ( fmc(nsample,nexpval), fatx(nsample,nxy), &
       &fplot(nsample,PLOT_NPOINT) )
    do irandom = 1, nsample
      ! Generate random instance of data.
      sample_y = xy_y + gaussian_random_number(xy_dy)
      ! Perform fit.
      call perform_fit (nxy, norder_v+1, xy_x-fit_x0, sample_y-fit_y0, &
         &vcoeff_nat, ierr)
      if (ierr/=0) call quit ('Could not perform fit.')
      ! Locate minimum of fit.
      call find_fit_minimum (norder_v+1, vcoeff_nat, xleft-fit_x0, 0.d0, &
         &xright-fit_x0, xmin, ymin)
      ! Regularize problem and transform the potential into Hermite polynomial.
      call obtain_ucentre_omega (rmass, norder_v, vcoeff_nat, 40, ucentre, &
         &omega)
      call transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
      ! Converge trial ground-state wave function with expansion order.
      do norder = max(norder_v,2), max(norder_v,2,MAX_NORDER)
        if (allocated(orbcoeff)) deallocate(orbcoeff, all_eigval, all_eigvec)
        allocate (orbcoeff(0:norder), all_eigval(0:norder), &
           &all_eigvec(0:norder,0:norder))
        call get_ground_state (rmass, norder, norder_v, vcoeff, e0, orbcoeff, &
           &vratio, all_eigval, all_eigvec)
        if (abs(vratio-1.d0)<VIRIAL_TOL) exit
      enddo ! norder
      if (norder>MAX_NORDER) call quit &
         &('failed to converge virial ratio to target accuracy.')
      norder = min(norder,MAX_NORDER)
      ! Store data for expectation value evaluation.
      fmc(irandom,1) = all_eigval(0)*omega - ymin
      fmc(irandom,2) = all_eigval(1)*omega - ymin
      fmc(irandom,3) = all_eigval(2)*omega - ymin
      fmc(irandom,4) = all_eigval(3)*omega - ymin
      fmc(irandom,5) = vratio
      fmc(irandom,6) = dble(norder)
      fmc(irandom,7) = xmin + fit_x0
      fmc(irandom,8) = ymin + fit_y0
      ! Store data for mean-curve fit.
      do ixy = 1, nxy
        fatx(irandom,ixy) = fit_y0 + &
           &eval_poly (norder_v+1, vcoeff_nat, xy_x(ixy)-fit_x0)
      enddo ! ixy
      ! Store data for plot.
      do ipoint = 1, PLOT_NPOINT
        xp = xleft-fit_x0 + dble(ipoint-1)/dble(PLOT_NPOINT-1) * (xright-xleft)
        fp = eval_poly (norder_v+1, vcoeff_nat, xp)
        fplot(irandom,ipoint) = fp+fit_y0
      enddo ! ipoint
    enddo ! irandom

    ! Dummy weight vector.
    allocate(w_one(nsample))
    w_one=1.d0

    ! Report expvals.
    write(6,'(a)')'Expectation values:'
    do iexpval = 1, nexpval
      call characterize_dist(nsample,fmc(1,iexpval),w_one,mean,var=var)
      write(6,'(2x,a16,":  ",es20.12,1x,es20.12)') expval_name(iexpval), &
         &mean, sqrt(var)
    enddo ! iexpval
    write(6,'()')

    ! Perform fit to mean V(r).
    allocate(fmeanx(nxy))
    do ixy = 1, nxy
      call characterize_dist(nsample,fatx(1,ixy),w_one,fmeanx(ixy))
    enddo ! ixy
    call perform_fit (nxy, norder_v+1, xy_x-fit_x0, fmeanx-fit_y0, &
       &vcoeff_nat, ierr)
    if (ierr/=0) call quit ('Could not perform fit.')
    ! Report.
    write(6,'(a)')'Parameters reproducing mean (offset by x,y of min in data):'
    do ipoly = 1, norder_v+1
      write(6,'(2x,"c_",a,":  ",es20.12,1x,es20.12)') trim(i2s(ipoly-1)),&
         &vcoeff_nat(ipoly-1)
    enddo ! ipoly
    write(6,'()')

    ! Solve in mean V(r).
    call obtain_ucentre_omega (rmass, norder_v, vcoeff_nat, 40, ucentre, omega)
    call transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
    ! Converge trial ground-state wave function with expansion order.
    do norder = max(norder_v,2), max(norder_v,2,MAX_NORDER)
      if (allocated(orbcoeff)) deallocate(orbcoeff, all_eigval, all_eigvec)
      allocate (orbcoeff(0:norder), all_eigval(0:norder), &
         &all_eigvec(0:norder,0:norder))
      call get_ground_state (rmass, norder, norder_v, vcoeff, e0, orbcoeff, &
         &vratio, all_eigval, all_eigvec)
      if (abs(vratio-1.d0)<VIRIAL_TOL) exit
    enddo ! norder
    if (norder>MAX_NORDER) call quit &
       &('failed to converge virial ratio to target accuracy.')
    norder = min(norder,MAX_NORDER)

    ! Find range of x where wave function is non-negligible.
    xl = locate_quantile(norder,orbcoeff,CDF_TOL)
    xr = locate_quantile(norder,orbcoeff,1.d0-CDF_TOL)
    ul = xl/sqrt(omega) - ucentre
    ur = xr/sqrt(omega) - ucentre

    ! Make plot data file.
    open(unit=io,file=trim(pname)//'.dat',status='replace',iostat=ierr)
    if (ierr/=0) call quit('Could not create '//trim(pname)//'.dat.')
    ! Plot original data.
    do ixy = 1, nxy
      write(io,'(3(1x,es20.12))') xy_x(ixy), xy_y(ixy), xy_dy(ixy)
    enddo ! ixy
    ! gnuplot "index" separator.
    write(io,'()')
    write(io,'()')
    ! Plot fit with confidence intervals.
    do ipoint = 1, PLOT_NPOINT
      xp = xleft + dble(ipoint-1)/dble(PLOT_NPOINT-1) * (xright-xleft)
      call get_conf_intervals (nsample, fplot(1,ipoint), fp, fp_CI)
      write(io,'(6(1x,es20.12))') xp, fp, fp_CI
    enddo ! ipoint
    ! gnuplot "index" separator.
    write(io,'()')
    write(io,'()')
    ! Plot mean fit.
    do ipoint = 1, PLOT_NPOINT
      xp = xleft + dble(ipoint-1)/dble(PLOT_NPOINT-1) * (xright-xleft)
      write(io,'(6(1x,es20.12))') xp, &
         &fit_y0 + eval_poly (norder_v+1, vcoeff_nat, xp-fit_x0)
    enddo ! ipoint
    ! gnuplot "index" separator.
    write(io,'()')
    write(io,'()')
    ! Plot first four states.
    allocate (hbasis(0:norder))
    plot_orb_scale_factor = 0.5d0*(all_eigval(1)-all_eigval(0))
    do iplot = 0, 3
      ! Plot eigenvalue.
      write(io,'(2(1x,es20.12))') fit_x0 + xl/sqrt(omega) - ucentre, &
         &fit_y0 + all_eigval(iplot)*omega
      write(io,'(2(1x,es20.12))') fit_x0 + xr/sqrt(omega) - ucentre, &
         &fit_y0 + all_eigval(iplot)*omega
      ! gnuplot "index" separator.
      write(io,'()')
      write(io,'()')
      ! Plot eigenfunction.
      do ipoint = 1, PLOT_NPOINT
        xp = xl + dble(ipoint-1)/dble(PLOT_NPOINT-1) * (xr-xl)
        up = xp/sqrt(omega) - ucentre
        call eval_hermite_poly_norm (norder, xp, hbasis)
        t1 = exp(-0.5d0*xp*xp)
        psix = t1*sum(all_eigvec(0:norder,iplot)*hbasis(0:norder))
        write(io,'(2(1x,es20.12))') fit_x0 + up, &
           &fit_y0 + all_eigval(iplot)*omega + plot_orb_scale_factor*omega*psix
      enddo ! ipoint
      ! gnuplot "index" separator.
      write(io,'()')
      write(io,'()')
    enddo ! iplot
    deallocate (hbasis)
    close(io)

    ! Create gnuplot file.
    open(unit=io,file=trim(pname)//'.gpi',status='replace',iostat=ierr)
    if (ierr/=0) call quit('Could not create '//trim(pname)//'.gpi.')
    write(io,'(a)')"set lmargin 9"
    write(io,'(a)')"set border lw 3"
    write(io,'(a)')"set xlabel '$r$ (a.u.)'"
    write(io,'(a)')"set ylabel '$E(r)$ (a.u.)'"
    write(io,'(a)')"set xtics autofreq 1"
    write(io,'(a)')"set mxtics 2"
    write(io,'(a)')"set ytics autofreq 0.1"
    write(io,'(a)')"set mytics 2"
    write(io,'(a)')"plot &
       &'plot.dat' index 3 using 1:2 w l lw 2 lc rgb '#888888' t '', &
       &'plot.dat' index 4 using 1:2 w l lw 2 lc rgb '#008800' t '', &
       &'plot.dat' index 5 using 1:2 w l lw 2 lc rgb '#888888' t '', &
       &'plot.dat' index 6 using 1:2 w l lw 2 lc rgb '#008800' t '', &
       &'plot.dat' index 7 using 1:2 w l lw 2 lc rgb '#888888' t '', &
       &'plot.dat' index 8 using 1:2 w l lw 2 lc rgb '#008800' t '', &
       &'plot.dat' index 9 using 1:2 w l lw 2 lc rgb '#888888' t '', &
       &'plot.dat' index 10 using 1:2 w l lw 2 lc rgb '#008800' t '', &
       &'plot.dat' index 1 using 1:3:6 w filledcurve &
       &  lt 1 lc rgb '#FF0000' fs transparent solid 0.35 t '', &
       &'plot.dat' index 1 using 1:4:5 w filledcurve &
       &  lt 1 lc rgb '#FF0000' fs transparent solid 0.35 t '', &
       &'plot.dat' index 1 using 1:2 w line lw 3 lc rgb '#FF0000' t '', &
       &'plot.dat' index 0 using 1:2:3 w errorbars &
       &  lw 3 lc rgb '#000000' pt 7 ps 0.8 t ''"
    close(io)
    write(6,'(a)')'Made gnuplot '//trim(pname)//'.gpi file and plain-text '//&
       &trim(pname)//'.dat data file.'
    write(6,'()')

  END SUBROUTINE main


    ! INPUT FILE HANDLING ROUTINES.


  SUBROUTINE check_file(fname,nline,ncolumn)
    !----------------------------------------------------!
    ! Check file contains data, and return the number of !
    ! data lines and the number of data columns in it.   !
    !----------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(out) :: nline,ncolumn
    CHARACTER(8192) line
    INTEGER ipos,ierr,ncol
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Initialize.
    nline=-1 ! flag non-existing file
    ncolumn=0

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)return
    nline=0

    ! Loop over lines.
    do
      read(io,'(a)',iostat=ierr)line
      if(ierr<0)exit
      if(ierr>0)then
        nline=-2 ! flag reading error
        exit
      endif
      line=adjustl(line)
      ! Skip comments.
      ipos=scan(line,'#!')
      if(ipos==1)cycle
      if(ipos>1)line=line(1:ipos-1)
      ! Skip empty lines.
      if(len_trim(line)==0)cycle
      ! Find how many elements there are in this line.
      ncol=nfield(line)
      if(ncolumn==0)then
        ncolumn=ncol
      else
        ncolumn=min(ncolumn,ncol)
      endif
      if(ncolumn<1)then
        nline=-3 ! flag column count problem
        exit
      endif
      nline=nline+1
    enddo

    ! Close file.
    close(io)

  END SUBROUTINE check_file


  SUBROUTINE read_file(fname,nxy,x,y,dy,ierr)
    !-------------------------------------!
    ! Read in the data in the input file. !
    !-------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: fname
    INTEGER,INTENT(in) :: nxy
    DOUBLE PRECISION,INTENT(inout) :: x(nxy),y(nxy),dy(nxy)
    INTEGER,INTENT(inout) :: ierr
    CHARACTER(8192) line
    INTEGER i,ipos
    ! Constants.
    INTEGER, PARAMETER :: io=10

    ! Open file.
    open(unit=io,file=trim(fname),status='old',iostat=ierr)
    if(ierr/=0)call quit('Problem opening "'//trim(fname)//'".')

    ! Loop over lines.
    i=0
    do
      read(io,'(a)',iostat=ierr)line
      if(ierr<0)then
        ierr=0
        exit
      endif
      if(ierr>0)call quit('Problem getting line from "'//trim(fname)//'".')
      line=adjustl(line)
      ! Skip comments.
      ipos=scan(line,'#!')
      if(ipos==1)cycle
      if(ipos>1)line=line(1:ipos-1)
      ! Skip empty lines.
      if(len_trim(line)==0)cycle
      i=i+1
      ! Read data point from string.
      x(i)=dble_field(1,line,ierr)
      if(ierr/=0)call quit('Failed to parse value of x in "'//trim(fname)//&
         &'".')
      y(i)=dble_field(2,line,ierr)
      if(ierr/=0)call quit('Failed to parse value of y in "'//trim(fname)//&
         &'".')
      dy(i)=dble_field(3,line,ierr)
      if(ierr/=0)call quit('Failed to parse value of dy in "'//trim(fname)//&
         &'".')
      if(lt_dble(dy(i),0.d0))call quit('Found negative dy in "'//trim(fname)//&
         &'".')
    enddo ! i

    ! Close file.
    close(io)

  END SUBROUTINE read_file


  SUBROUTINE obtain_ucentre_omega (rmass, norder_v, vcoeff_nat, norder, &
     &ucentre, omega)
    !----------------------------------------------------------!
    ! Given the natural-polynomial coefficients of a potential !
    ! VCOEFF_NAT(0:NORDER_V), obtain the values of ucentre and !
    ! omega that minimizes the variational energy or virial    !
    ! ratio error for a trial wave function of order NORDER.   !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder_v, norder
    DOUBLE PRECISION, INTENT(in) :: rmass, vcoeff_nat(0:norder_v)
    DOUBLE PRECISION, INTENT(inout) :: ucentre, omega
    DOUBLE PRECISION, PARAMETER :: OMEGA_TOL = 1.d-9
    LOGICAL, PARAMETER :: SET_OMEGA_BY_X2 = .false.
    LOGICAL, PARAMETER :: MINIMIZE_E = .true.
    INTEGER, PARAMETER :: MAX_ITER_CENTRE = 10
    DOUBLE PRECISION vcoeff(0:norder_v), orbcoeff(0:norder), e0, vratio, &
       &ucentre_prev, omega_prev, omega_init
    DOUBLE PRECISION xu, xv, xw, fu, fv, fw, x, f
    INTEGER iter, iorder_v
    LOGICAL rejected

    ! Obtain initial guess for omega so that:
    ! * Scaling a potential V(u) by a multiplicative constant results in the
    !   same dimensionless potential v(x).
    ! * omega is the frequency for a harmonic potential.
    omega_init = sqrt(2.d0)*rmass
    do iorder_v = norder_v, 2, -1
      if (mod(iorder_v,2)/=0) cycle
      if (le_dble(vcoeff_nat(iorder_v),0.d0)) cycle
      omega_init = sqrt(2.d0)*rmass*&
         &vcoeff_nat(iorder_v)**(2.d0/dble(iorder_v+2))
    enddo ! iorder_v

    ! Initialize ucentre.
    ucentre = 0.d0

    if (SET_OMEGA_BY_X2) then

      ! Loop over self-consistence cycles.
      omega = omega_init
      do iter = 1, MAX_ITER_CENTRE

        ! Evaluate <x> at this omega.
        call transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
        call get_ground_state (rmass, norder, norder_v, vcoeff, e0, orbcoeff, &
           &vratio)
        ucentre_prev = ucentre
        omega_prev = omega
        ucentre = ucentre_prev - &
           &eval_xpower_expval (norder, orbcoeff, 1)/sqrt(omega_prev)
        omega = 0.5d0 * omega_prev/eval_xpower_expval (norder, orbcoeff, 2)
        if (abs(ucentre-ucentre_prev)<1.d-10.and.abs(omega-omega_prev)<1.d-10) &
           &exit

      enddo

    else ! .not. SET_OMEGA_BY_X2

      ! Loop over self-consistence cycles.
      do iter = 1, MAX_ITER_CENTRE

        ! Obtain initial guess for omega so that:
        ! * Scaling a potential V(u) by a multiplicative constant results in the
        !   same dimensionless potential v(x).
        ! * omega is the frequency for a harmonic potential.
        xv = omega_init
        call transform_potential (norder_v, vcoeff_nat, ucentre, xv, vcoeff)
        call get_ground_state (rmass, norder, norder_v, vcoeff, e0, orbcoeff, &
           &vratio)
        if (MINIMIZE_E) then
          fv = e0*xv
        else
          fv = abs(vratio-1.d0)
        endif

        ! Perform line minimization for omega.

        ! Bracket to the left.
        xw = 0.d0
        fw = fv-1.d0
        do
          xu = 0.9d0*xv
          call transform_potential (norder_v, vcoeff_nat, ucentre, xu, vcoeff)
          call get_ground_state (rmass, norder, norder_v, vcoeff, e0, &
             &orbcoeff, vratio)
          if (MINIMIZE_E) then
            fu = e0*xu
          else
            fu = abs(vratio-1.d0)
          endif
          if (fu>fv) exit
          xw = xv
          fw = fv
          xv = xu
          fv = fu
        enddo

        ! Bracket to the right.
        if (fw<fu) then
          do
            xw = 1.1d0*xv
            call transform_potential (norder_v, vcoeff_nat, ucentre, xw, &
               &vcoeff)
            call get_ground_state (rmass, norder, norder_v, vcoeff, e0, &
               &orbcoeff, vratio)
            if (MINIMIZE_E) then
              fw = e0*xw
            else
              fw = abs(vratio-1.d0)
            endif
            if (fw>fv) exit
            xu = xv
            fu = fv
            xv = xw
            fv = fw
          enddo
        endif

        ! Zone in on minimum.
        do
          call parabolic_min (xu, xv, xw, fu, fv, fw, x, f, rejected)
          if (rejected) exit
          if (x<=xu.or.x>=xw) exit
          call transform_potential (norder_v, vcoeff_nat, ucentre, x, vcoeff)
          call get_ground_state (rmass, norder, norder_v, vcoeff, e0, &
             &orbcoeff, vratio)
          if (MINIMIZE_E) then
            f = e0*x
          else
            f = abs(vratio-1.d0)
          endif
          if (f<fv) then
            if (x<xv) then
              xw = xv
              fw = fv
            elseif (x>xv) then
              xu = xv
              fu = fv
            else
              exit
            endif
            xv = x
            fv = f
          elseif (f>fv) then
            if (x<xv) then
              xu = x
              fu = f
            elseif (x>xv) then
              xw = x
              fw = f
            else
              exit
            endif
          else
            exit
          endif
          if (xw-xu<OMEGA_TOL) exit
        enddo

        ! Return position of minimum.
        omega = xv

        ! Evaluate <x> at this omega and shift ucentre by it.
        call transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
        call get_ground_state (rmass, norder, norder_v, vcoeff, e0, orbcoeff, &
           &vratio)
        ucentre_prev = ucentre
        ucentre = ucentre_prev - &
           &eval_xpower_expval (norder, orbcoeff, 1)/sqrt(omega)
        if (abs(ucentre-ucentre_prev)<1.d-10) exit

      enddo

    endif ! SET_OMEGA_BY_X2 or not

  END SUBROUTINE obtain_ucentre_omega


  SUBROUTINE transform_potential (norder_v, vcoeff_nat, ucentre, omega, vcoeff)
    !---------------------------------------------------------!
    ! Given a function represented as a natural polynomial of !
    ! coefficients VCOEFF_NAT(0:NORDER_V), apply a horizontal !
    ! shift UCENTRE and rescale factor sqrt(OMEGA), and       !
    ! re-represent it in normalized Hermite polynomials of    !
    ! coefficients VCOEFF(0:NORDER_V).                        !
    !---------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder_v
    DOUBLE PRECISION, INTENT(in) :: vcoeff_nat(0:norder_v), ucentre, omega
    DOUBLE PRECISION, INTENT(inout) :: vcoeff(0:norder_v)
    DOUBLE PRECISION lu_hmatrix(0:norder_v,0:norder_v), &
       &vcoeff_nat1(0:norder_v), vcoeff_nat2(0:norder_v)
    INTEGER i, j, piv_hmatrix(0:norder_v)

    ! Prepare Hermite transformation matrix.
    call lu_decom_hermite_matrix (norder_v, lu_hmatrix, piv_hmatrix)

    ! Apply shift.
    do i = 0, norder_v
      vcoeff_nat1(i) = sum( (/ ( rchoose(i+j,i)*(-ucentre)**j*vcoeff_nat(i+j), &
          &                      j=0,norder_v-i ) /) )
    enddo ! i

    ! Apply change of variable.
    vcoeff_nat2 = (/ ( vcoeff_nat1(i)*omega**(-0.5d0*dble(i+2)), &
       &               i=0,norder_v ) /)

    ! Perform conversion.
    call convert_natpoly_to_hermite (norder_v, lu_hmatrix, piv_hmatrix, &
       &vcoeff_nat2, vcoeff)

  END SUBROUTINE transform_potential


  SUBROUTINE get_ground_state (rmass, norder, norder_v, vcoeff, e0, orbcoeff, &
     &vratio, all_eigval, all_eigvec)
    !---------------------------------------------------------!
    ! Given a one-dimensional (anharmonic) potential,         !
    !                                                         !
    !   v(x) = Sum_i VCOEFF(i)*N_i(x) ,                       !
    !                                                         !
    ! where H_i is the i-th Hermite polynomial, construct the !
    ! matrix elements of the Hamiltonian                      !
    !                                                         !
    !   H(x) = -1/2 d/dx^2 + v(x) ,                           !
    !                                                         !
    ! in the basis of the eigenfunctions of the harmonic      !
    ! oscillator of unit frequency,                           !
    !                                                         !
    !   phi_i(x) = exp(-x^2/2) * N_i(x) ,                     !
    !                                                         !
    ! and solve for the coefficients of the normalized trial  !
    ! wave function,                                          !
    !                                                         !
    !   Psi(x) = Sum_i ORBCOEFF(i)*phi_i(x) ,                 !
    !                                                         !
    ! and for the associated variational energy E0.           !
    !---------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder, norder_v
    DOUBLE PRECISION, INTENT(in) :: rmass, vcoeff(0:norder_v)
    DOUBLE PRECISION, INTENT(inout) :: e0, orbcoeff(0:norder), vratio
    DOUBLE PRECISION, INTENT(inout), OPTIONAL :: all_eigval(0:norder), &
       &all_eigvec(0:norder,0:norder)
    ! Eigenproblem arrays.
    DOUBLE PRECISION alpha(0:norder), hmatrix(0:norder,0:norder), &
       &cmatrix(0:norder,0:norder)
    ! Buffer for logarithms of factorials.
    DOUBLE PRECISION log_fact(0:norder)
    ! LAPACK work arrays.
    DOUBLE PRECISION, ALLOCATABLE :: lapack_work(:)
    INTEGER, ALLOCATABLE :: lapack_iwork(:)
    INTEGER lapack_lwork, lapack_liwork
    ! Parameters.
    DOUBLE PRECISION, PARAMETER :: TOL_ZERO = 1.d3*epsilon(1.d0)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_root_pi_over_sqrt8 = &
       &                           pi**0.25d0*sqrt(0.125d0)
    ! Virial ratio evaluation.
    DOUBLE PRECISION vircoeff(0:norder_v), xdv_expval, t_expval
    ! Misc local variables.
    INTEGER i, j, k, ierr
    DOUBLE PRECISION t1, t2, inv_rmass

    ! Get numerical constants to speed up operations.
    log_fact(0:norder) = eval_log_fact( (/ (i, i=0,norder) /) )
    inv_rmass = 1.d0/rmass

    ! Populate Hamiltonian matrix in the basis of harmonic-oscillator
    ! eigenfunctions.
    do i = 0, norder
      hmatrix(i,i) = eval_comb_Gamma (norder_v, i, i, vcoeff, log_fact) + &
         &inv_rmass * ( dble(i)+0.25d0 - &
         &fourth_root_pi_over_sqrt8*eval_Gamma(i,i,2,log_fact) )
      do j = i+1, norder
        hmatrix(i,j) = eval_comb_Gamma (norder_v, i, j, vcoeff, log_fact) - &
           &fourth_root_pi_over_sqrt8 * eval_Gamma(i,j,2,log_fact) * inv_rmass
        hmatrix(j,i) = hmatrix(i,j)
      enddo ! j
    enddo ! i
    cmatrix = hmatrix

    ! Diagonalize Hamiltonian.
    lapack_lwork = 1
    lapack_liwork = 1
    allocate(lapack_work(lapack_lwork), lapack_iwork(lapack_liwork))
    lapack_lwork = -1
    lapack_liwork = -1
    call dsyevd ('V', 'U', norder+1, cmatrix, norder+1, alpha, &
       &lapack_work, lapack_lwork, lapack_iwork, lapack_liwork, ierr)
    if (ierr/=0) call quit ('DSYEVD error '//trim(i2s(ierr))//'.')
    lapack_lwork = nint(lapack_work(1))
    lapack_liwork = lapack_iwork(1)
    deallocate(lapack_work, lapack_iwork)
    allocate(lapack_work(lapack_lwork), lapack_iwork(lapack_liwork))
    call dsyevd ('V', 'U', norder+1, cmatrix, norder+1, alpha, &
       &lapack_work, lapack_lwork, lapack_iwork, lapack_liwork, ierr)
    if (ierr/=0) call quit ('DSYEVD error '//trim(i2s(ierr))//'.')
    deallocate(lapack_work, lapack_iwork)

    ! Loop over all eigenstates to tidy them up.
    do i = 0, norder
      ! Normalize the wave function.
      t1 = 1.d0/sqrt(sum(cmatrix(0:norder,i)**2))
      cmatrix(0:norder,i) = t1*cmatrix(0:norder,i)
      ! Flush small coefficients to zero.
      where (abs(cmatrix(0:norder,i)) < TOL_ZERO) cmatrix(0:norder,i) = 0.d0
      ! Normalize the wave function again, and make first coefficient positive
      ! while at it.
      t1 = 1.d0/sqrt(sum(cmatrix(0:norder,i)**2))
      do j = 0, norder
        if (abs(cmatrix(j,i)) > TOL_ZERO) then
          if (cmatrix(j,i)<0.d0) t1 = -t1
          exit
        endif
      enddo ! j
      cmatrix(0:norder,i) = t1*cmatrix(0:norder,i)
      ! Recalculate the variational energy.
      alpha(i) = 0.d0
      do j = 0, norder
        t1 = 0.d0
        do k = j+1, norder
          t1 = t1 + cmatrix(k,i)*hmatrix(k,j)
        enddo ! k
        t1 = cmatrix(j,i)*(cmatrix(j,i)*hmatrix(j,j) + 2.d0*t1)
        alpha(i) = alpha(i)+t1
      enddo ! j
    enddo ! i

    ! Evaluate the virial ratio for the ground state.
    vircoeff = 0.d0
    do i = 0, norder_v
      if (i<norder_v-1) vircoeff(i) = vircoeff(i) + &
         &vcoeff(i+2)*sqrt(dble((i+1)*(i+2)))
      if (i>1) vircoeff(i) = vircoeff(i) + vcoeff(i)*dble(i)
    enddo ! i
    xdv_expval = 0.d0
    t_expval = 0.d0
    do i = 0, norder
      t1 = eval_comb_Gamma (norder_v, i, i, vircoeff, log_fact)
      t2 = dble(i)+0.25d0 - fourth_root_pi_over_sqrt8*eval_Gamma(i,i,2,log_fact)
      xdv_expval = xdv_expval + t1*cmatrix(i,0)**2
      t_expval = t_expval + t2*cmatrix(i,0)**2
      do j = i+1, norder
        t1 = eval_comb_Gamma (norder_v, i, j, vircoeff, log_fact)
        t2 = -fourth_root_pi_over_sqrt8*eval_Gamma(i,j,2,log_fact)
        xdv_expval = xdv_expval + 2.d0*t1*cmatrix(i,0)*cmatrix(j,0)
        t_expval = t_expval + 2.d0*t2*cmatrix(i,0)*cmatrix(j,0)
      enddo ! j
    enddo ! i
    vratio = 2.d0*t_expval*inv_rmass/xdv_expval

    ! Return ground-state components.
    ! FIXME - if ground state is degenerate, how do we choose which
    ! eigenstate to return?
    e0 = alpha(0)
    orbcoeff(0:norder) = cmatrix(0:norder,0)

    ! Copy data for all eigenstates if requested.
    if (present(all_eigval)) all_eigval(0:norder) = alpha(0:norder)
    if (present(all_eigvec)) all_eigvec(0:norder,0:norder) = &
       &cmatrix(0:norder,0:norder)

  END SUBROUTINE get_ground_state


  DOUBLE PRECISION FUNCTION eval_xpower_expval (norder, orbcoeff, iexp)
    !---------------------------------------------------------------------!
    ! Given the coefficients ORBCOEFF(0:NORDER) of a trial wave function, !
    ! evaluate the expectation value of the natural power <x^IEXP> .      !
    !---------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder, iexp
    DOUBLE PRECISION, INTENT(in) :: orbcoeff(0:norder)
    ! Misc local variables.
    INTEGER i, j, piv_hmatrix(0:iexp)
    DOUBLE PRECISION lu_hmatrix(0:iexp,0:iexp), xcoeff_nat(0:iexp), &
       &xcoeff(0:iexp), log_fact(0:max(iexp,norder))

    ! Re-represent x^iexp in the basis of Hermite polynomials.
    xcoeff_nat(0:iexp) = 0.d0
    xcoeff_nat(iexp) = 1.d0
    call lu_decom_hermite_matrix (iexp, lu_hmatrix, piv_hmatrix)
    call convert_natpoly_to_hermite (iexp, lu_hmatrix, piv_hmatrix, &
       &xcoeff_nat, xcoeff)

    ! Get numerical constants to speed up operations.
    log_fact(0:max(iexp,norder)) = &
       &eval_log_fact( (/ (i, i=0,max(iexp,norder)) /) )

    ! Evaluate expectation value.
    eval_xpower_expval = 0.d0
    do i = 0, norder
      eval_xpower_expval = eval_xpower_expval + orbcoeff(i) * orbcoeff(i) * &
         &eval_comb_Gamma (iexp, i, i, xcoeff, log_fact)
      do j = i+1, norder
        eval_xpower_expval = eval_xpower_expval + &
           &2.d0 * orbcoeff(j) * orbcoeff(i) * eval_comb_Gamma &
           &(iexp, i, j, xcoeff, log_fact)
      enddo ! j
    enddo ! i
    if (abs(eval_xpower_expval)<1.d3*epsilon(1.d0)) eval_xpower_expval = 0.d0

  END FUNCTION eval_xpower_expval


  DOUBLE PRECISION FUNCTION locate_quantile (norder, orbcoeff, p)
    !----------------------------------------------------------!
    ! Locate value of x at which the cumulative distribution   !
    ! function associated with the wave function of normalized !
    ! Hermite coefficients ORBCOEFF(0:NORDER) is P.            !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder
    DOUBLE PRECISION, INTENT(in) :: orbcoeff(0:norder), p
    DOUBLE PRECISION, PARAMETER :: ABS_TOL = 1.d-8
    DOUBLE PRECISION, PARAMETER :: REL_TOL = 1.d-3
    DOUBLE PRECISION x, f, x1, f1, x2, f2, dx

    if (p<=0.d0 .or. p>=1.d0) call quit('Asked for impossible quantile.')
    locate_quantile = 0.d0

    ! Get initial pair of points.
    x1 = 0.d0
    x2 = 0.d0
    f1 = eval_cdf (norder, orbcoeff, x1)
    f2 = f1
    dx = 1.d0
    if (f1>p) then
      do while (f1>p)
        x2 = x1
        f2 = f1
        x1 = x1 - dx
        f1 = eval_cdf (norder, orbcoeff, x1)
        dx = dx*1.5d0
      enddo
    elseif (f2<p) then
      do while (f2<p)
        x1 = x2
        f1 = f2
        x2 = x2 + dx
        f2 = eval_cdf (norder, orbcoeff, x2)
        dx = dx*1.5d0
      enddo
    endif

    ! Loop until convergence.
    do
      !x = x1 + 0.5d0*(p-f1)*((x2-x1)/(f2-f1))
      !if (x<=x1.or.x>=x2) x = 0.5d0*(x1+x2) ! guard against poor numerics
      x = 0.5d0*(x1+x2) ! pure bisection
      if (x<=x1.or.x>=x2) exit ! regard as converged
      f = eval_cdf (norder, orbcoeff, x)
      if (abs(f-p)<min(ABS_TOL,REL_TOL*min(p,1.d0-p))) exit
      if (f>p) then
        x2 = x
        f2 = f
      else
        x1 = x
        f1 = f
      endif
    enddo
    locate_quantile = x

  END FUNCTION locate_quantile


  DOUBLE PRECISION FUNCTION eval_cdf (norder, orbcoeff, x)
    !-----------------------------------------------------------!
    ! Evaluate the cumulative distribution function associated  !
    ! with the wave function of normalized Hermite coefficients !
    ! ORBCOEFF(0:NORDER) at X.                                  !
    !-----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder
    DOUBLE PRECISION, INTENT(in) :: orbcoeff(0:norder), x
    DOUBLE PRECISION exp_x2, hbasis(0:2*norder), log_fact(0:2*norder), &
       &coeff(0:2*norder), t1
    INTEGER i, j, k

    ! Compute required quantities.
    exp_x2 = exp(-x*x)
    call eval_hermite_poly_norm (2*norder, x, hbasis)
    log_fact(0:2*norder) = eval_log_fact( (/ (i, i=0,2*norder) /) )
    coeff(0) = 0.d0
    coeff(1:2*norder) = (/ ( -exp_x2 * hbasis(k-1) / sqrt(dble(2*k)), &
         &                   k=1,2*norder ) /)

    ! Main contribution.
    eval_cdf = 0.5d0*(1.d0+erf(x))

    ! Loop over corrections.
    do i = 0, norder
      t1 = eval_comb_Gamma (2*i, i, i, coeff(0:2*i), log_fact(0:2*i))
      eval_cdf = eval_cdf + orbcoeff(i)*orbcoeff(i)*t1
      do j = i+1, norder
        t1 = eval_comb_Gamma (i+j, i, j, coeff(0:i+j), log_fact(0:i+j))
        eval_cdf = eval_cdf + 2.d0*orbcoeff(i)*orbcoeff(j)*t1
      enddo ! j
    enddo ! i

  END FUNCTION eval_cdf


  ! Hermite polynomial tools.  NB, we use "normalized" Hermite polynomials,
  ! N_n(x) = pi^-1/4 2^-n/2 (n!)^-1/2 H_n(x).


  DOUBLE PRECISION FUNCTION eval_Gamma (i, j, k, log_fact)
    !--------------------------------------------------------------!
    ! Evaluate the Gamma_i,j,k symbol,                             !
    !                                                              !
    !   Gamma_i,j,k = integral exp(-x^2) N_i(x) N_j(x) N_k(x) dx . !
    !--------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: i, j, k
    DOUBLE PRECISION, INTENT(in) :: log_fact(0:)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_log_pi = 0.25d0*log(pi)
    ! Local variables.
    DOUBLE PRECISION t1
    eval_Gamma = 0.d0
    if (k<abs(i-j) .or. k>i+j .or. mod(i+j+k,2)/=0) return
    t1 = -fourth_log_pi + 0.5d0*(log_fact(i) + log_fact(j) + log_fact(k)) -&
      &log_fact((j+k-i)/2) - log_fact((k+i-j)/2) - log_fact((i+j-k)/2)
    eval_Gamma = exp(t1)
  END FUNCTION eval_Gamma


  DOUBLE PRECISION FUNCTION eval_comb_Gamma (n, i, j, coeff, log_fact)
    !--------------------------------------------------!
    ! Evaluate the linear combination of Gamma symbols !
    !                                                  !
    !   Sum_k=0^n coeff_k * Gamma_i,j,k .              !
    !                                                  !
    ! NB, Gamma evaluator inlined to skip the logic.   !
    !--------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n, i, j
    DOUBLE PRECISION, INTENT(in) :: coeff(0:n), log_fact(0:)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_log_pi = 0.25d0*log(pi)
    ! Local variables.
    DOUBLE PRECISION t1
    INTEGER k
    eval_comb_Gamma = 0.d0
    do k = i+j, abs(i-j), -2
      if (k>n) cycle
      t1 = -fourth_log_pi + 0.5d0*(log_fact(i) + log_fact(j) + log_fact(k)) -&
        &log_fact((j+k-i)/2) - log_fact((k+i-j)/2) - log_fact((i+j-k)/2)
      eval_comb_Gamma = eval_comb_Gamma + coeff(k)*exp(t1)
    enddo ! k
  END FUNCTION eval_comb_Gamma


  SUBROUTINE eval_hermite_poly_norm (n, x, h)
    !---------------------------------------------------!
    ! Evaluate N_n(x) for n = 0:N at x=X, returning the !
    ! values in H(0:N).                                 !
    !---------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    DOUBLE PRECISION, INTENT(in) :: x
    DOUBLE PRECISION, INTENT(inout) :: h(0:n)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: inv_pi_one_fourth = pi**(-0.25d0)
    DOUBLE PRECISION, PARAMETER :: sqrt2 = sqrt(2.d0)
    ! Local variables.
    DOUBLE PRECISION t0, t1, t2, sqrti, sqrti_1, sqrt2_x
    INTEGER i
    if (n<0) return
    t1 = inv_pi_one_fourth
    h(0) = t1
    if (n<1) return
    sqrt2_x = sqrt2*x
    t0 = inv_pi_one_fourth*sqrt2_x
    sqrti = 1.d0
    h(1) = t0
    do i = 2, n
      t2 = t1
      t1 = t0
      sqrti_1 = sqrti
      sqrti = sqrt(dble(i))
      t0 = (sqrt2_x*t1 - sqrti_1*t2)/sqrti
      h(i) = t0
    enddo ! i
  END SUBROUTINE eval_hermite_poly_norm


  SUBROUTINE lu_decom_hermite_matrix (norder, lu_hmatrix, piv_hmatrix)
    !-----------------------------------------------------------------!
    ! Returns the LU decomposition of the matrix of coefficients of   !
    ! natural powers in Hermite polynomials up to order NORDER, which !
    ! can then be used to express natural polynomials in Hermite      !
    ! polynomials.                                                    !
    !-----------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder
    DOUBLE PRECISION, INTENT(inout) :: lu_hmatrix(0:norder,0:norder)
    INTEGER, INTENT(inout) :: piv_hmatrix(0:norder)
    ! Numerical constants.
    DOUBLE PRECISION, PARAMETER :: pi = 4.d0*atan(1.d0)
    DOUBLE PRECISION, PARAMETER :: fourth_log_pi = 0.25d0*log(pi)
    DOUBLE PRECISION, PARAMETER :: log2 = log(2.d0)
    ! Misc local variables.
    DOUBLE PRECISION t1
    INTEGER n, i, isgn, ierr

    ! Build matrix of coefficients of N_n as a linear combination of natural
    ! powers.
    lu_hmatrix = 0.d0
    do n = 0, norder
      isgn = 1
      do i = n, 0, -2
        t1 = 0.5d0*(eval_log_fact(n) + dble(2*i-n)*log2) - &
           &eval_log_fact(i) - eval_log_fact(n/2-i/2) - fourth_log_pi
        lu_hmatrix(i,n) = dble(isgn) * exp(t1)
        isgn = -isgn
      enddo ! i
    enddo ! n

    ! LU-decompose the matrix.
    call dgetrf (norder+1, norder+1, lu_hmatrix, norder+1, piv_hmatrix, ierr)
    if (ierr/=0) call quit ('DGETRF error '//trim(i2s(ierr))//'.')

  END SUBROUTINE lu_decom_hermite_matrix


  SUBROUTINE convert_natpoly_to_hermite (norder, lu_hmatrix, piv_hmatrix, &
     &pcoeff, hcoeff)
    !----------------------------------------------------------!
    ! Convert a polynomial of coefficients PCOEFF(0:NORDER) to !
    ! a Hermite polynomial of coefficients HCOEFF(0:NORDER).   !
    !----------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: norder, piv_hmatrix(0:norder)
    DOUBLE PRECISION, INTENT(in) :: lu_hmatrix(0:norder,0:norder), &
       &pcoeff(0:norder)
    DOUBLE PRECISION, INTENT(inout) :: hcoeff(0:norder)
    INTEGER ierr
    hcoeff(0:norder) = pcoeff(0:norder)
    call dgetrs ('N', norder+1, 1, lu_hmatrix, norder+1, piv_hmatrix, hcoeff, &
       &norder+1, ierr)
    if (ierr/=0) call quit ('DGETRS error '//trim(i2s(ierr))//'.')
  END SUBROUTINE convert_natpoly_to_hermite


  ! Numerical tools.


    SUBROUTINE perform_fit (nxy, npoly, x, y, a, ierr)
    !-----------------------------------------!
    ! Perform least-squares fit of xy data to !
    ! polynomial of order npoly-1.            !
    !-----------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: nxy, npoly
    DOUBLE PRECISION, INTENT(in) :: x(nxy), y(nxy)
    DOUBLE PRECISION, INTENT(inout) :: a(npoly)
    INTEGER, INTENT(inout) :: ierr
    ! Local variables.
    INTEGER i, j, ipiv(npoly), lwork
    DOUBLE PRECISION M(npoly,npoly), Minv(npoly,npoly), c(npoly)
    DOUBLE PRECISION, ALLOCATABLE :: work(:)

    ! Construct c vector and M matrix.
    do i=1,npoly
      c(i)=sum(y(:)*x(:)**dble(i-1))
      M(i,i)=sum(x(:)**(dble(i+i-2)))
      do j=i+1,npoly
        M(j,i)=sum(x(:)**(dble(i+j-2)))
        M(i,j)=M(j,i)
      enddo ! j
    enddo ! i

    ! Invert M.
    Minv=M
    allocate(work(1))
    lwork=-1
    call dsytrf('L',npoly,Minv,npoly,ipiv,work,lwork,ierr)
    if(ierr/=0)return
    lwork=nint(work(1))
    deallocate(work)
    allocate(work(lwork),stat=ierr)
    if(ierr/=0)return
    call dsytrf('L',npoly,Minv,npoly,ipiv,work,lwork,ierr)
    if(ierr/=0)return
    deallocate(work)
    allocate(work(npoly),stat=ierr)
    if(ierr/=0)return
    call dsytri('L',npoly,Minv,npoly,ipiv,work,ierr)
    if(ierr/=0)return
    deallocate(work)

    ! Complete Minv and evaluate coefficients.
    do i=1,npoly
      do j=i+1,npoly
        Minv(i,j)=Minv(j,i)
      enddo ! j
    enddo ! i
    a=matmul(Minv,c)

  END SUBROUTINE perform_fit


  DOUBLE PRECISION FUNCTION eval_poly (npoly, a, x_target)
    !----------------------------------------------------!
    ! Evaluate the polynomial of coefficients A(1:NPOLY) !
    ! at X_TARGET.                                       !
    !----------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: npoly
    DOUBLE PRECISION, INTENT(in) :: a(npoly), x_target
    INTEGER ipoly
    DOUBLE PRECISION xx
    eval_poly = 0.d0
    if (npoly<1) return
    eval_poly = a(1)
    xx = 1.d0
    do ipoly = 2, npoly
      xx = xx*x_target
      eval_poly = eval_poly + a(ipoly)*xx
    enddo ! ipoly
  END FUNCTION eval_poly


  SUBROUTINE find_fit_minimum (npoly, a, x1, x2, x3, xmin, ymin)
    !-------------------------------------------------------!
    ! Given a polynomial with parameters A(1:NPOLY), find   !
    ! the location of its minimum (XMIN,YMIN) from initial  !
    ! bracket X1 > X2 < X3, where X1 and X3 are potentially !
    ! unreliable.                                           !
    !-------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: npoly
    DOUBLE PRECISION, INTENT(in) :: a(npoly), x1, x2, x3
    DOUBLE PRECISION, INTENT(inout) :: xmin, ymin
    LOGICAL rejected
    DOUBLE PRECISION xu, xv, xw, fu, fv, fw, xt, ft
    DOUBLE PRECISION, PARAMETER :: XTOL = 1.d-7

    ! Bracket by data range.
    xu = x1
    fu = eval_poly (npoly, a, xu)
    xv = x2
    fv = eval_poly (npoly, a, xv)
    xw = x3
    fw = eval_poly (npoly, a, xw)

    ! Tighten U.
    do while (le_dble(fu,fv))
      xu = 0.5d0*(xu+xv)
      fu = eval_poly (npoly, a, xu)
    enddo

    ! Tighten W.
    do while (le_dble(fw,fv))
      xw = 0.5d0*(xw+xv)
      fw = eval_poly (npoly, a, xw)
    enddo

    ! Zone in on minimum.
    do
      call parabolic_min (xu, xv, xw, fu, fv, fw, xt, ft, rejected)
      if (rejected) exit
      if (xt<=xu.or.xt>=xw) exit
      ft = eval_poly (npoly, a, xt)
      if (ft<fv) then
        if (xt<xv) then
          xw = xv
          fw = fv
        elseif (xt>xv) then
          xu = xv
          fu = fv
        else
          exit
        endif
        xv = xt
        fv = ft
      elseif (ft>fv) then
        if (xt<xv) then
          xu = xt
          fu = ft
        elseif (xt>xv) then
          xw = xt
          fw = ft
        else
          exit
        endif
      else
        exit
      endif
      if (xw-xu < XTOL) exit
    enddo

    ! Output values.
    xmin = xv
    ymin = fv

  END SUBROUTINE find_fit_minimum


  ELEMENTAL DOUBLE PRECISION FUNCTION eval_log_fact(k)
    !----------------------------------------------!
    ! Returns the logarithm of the factorial of k. !
    !----------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: k
    ! Local variables.
    INTEGER i
    eval_log_fact = 0.d0
    do i = 2, k
      eval_log_fact = eval_log_fact + log(dble(i))
    enddo ! i
  END FUNCTION eval_log_fact


  SUBROUTINE parabolic_min (x1, x2, x3, y1, y2, y3, x0, y0, rejected)
    !-----------------------------------------------------------------!
    ! Fit three points to a parabola and return (x,y) of the min/max. !
    !-----------------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: x1, x2, x3, y1, y2, y3
    DOUBLE PRECISION, INTENT(out) :: x0, y0
    LOGICAL, INTENT(out) :: rejected
    DOUBLE PRECISION a, b, c, numa, numb, numc, den, x1_sq, x2_sq, x3_sq,&
       &x21, x32, x31, z1, z2, z3, invden

    ! Initialize.
    x0 = x2
    y0 = y2
    rejected = .false.

    ! Check that x and y values are distinguishable and in correct order.
    if (x1>=x2 .or. x2>=x3 .or. (y2>=y1.eqv.y3>=y2))then
      rejected = .true.
      return
    endif

    ! Compute squares.
    x1_sq = x1*x1
    x2_sq = x2*x2
    x3_sq = x3*x3

    ! Renormalize for better numerics.
    x31 = x3-x1
    x21 = (x2-x1)/x31
    x32 = (x3-x2)/x31
    z1 = y1*x32
    z2 = y2
    z3 = y3*x21

    ! Solve linear system.
    den = -x1_sq*x32 + x2_sq - x3_sq*x21
    numa = -z1 + z2 - z3
    numb = z1*(x2+x3) - z2*(x1+x3) + z3*(x1+x2)
    numc = -z1*x2*x3 + z2*x3*x1 - z3*x1*x2

    ! Find x0 and y0.
    invden = 1.d0/den
    a = numa*invden
    b = numb*invden
    c = numc*invden
    x0 = -0.5d0*numb/numa
    y0 = (a*x0+b)*x0 + c

  END SUBROUTINE parabolic_min


  DOUBLE PRECISION FUNCTION rchoose(a,b)
    !-------------------------------------!
    ! This function returns a choose b as !
    ! a floating-point real number.       !
    !-------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: a,b
    INTEGER i
    rchoose=1.d0
    do i=1,b
      rchoose=rchoose*(dble(a+1-i)/dble(i))
    enddo ! i
  END FUNCTION rchoose


  LOGICAL ELEMENTAL FUNCTION eq_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X == Y within a reasonable tolerance.             !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    DOUBLE PRECISION abs_x,abs_y,big,small
    ! Parameters.
    DOUBLE PRECISION,PARAMETER :: tol_zero=1.d-12
    DOUBLE PRECISION,PARAMETER :: tol_rel=1.d-9
    if(present(tol))then
      eq_dble=abs(x-y)<=tol
    else
      abs_x=abs(x)
      abs_y=abs(y)
      if(abs_x<=tol_zero.and.abs_y<=tol_zero)then
        eq_dble=.true.
      elseif(x>0.d0.eqv.y>0.d0)then
        big=max(abs_x,abs_y)
        small=min(abs_x,abs_y)
        eq_dble=big-small<=big*tol_rel
      else
        eq_dble=.false.
      endif
    endif
  END FUNCTION eq_dble


  LOGICAL ELEMENTAL FUNCTION neq_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X /= Y within a reasonable tolerance.             !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    neq_dble=.not.eq_dble(x,y,tol)
  END FUNCTION neq_dble


  LOGICAL ELEMENTAL FUNCTION lt_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X < Y within a reasonable tolerance.              !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    lt_dble=x<y.and.neq_dble(x,y,tol)
  END FUNCTION lt_dble


  LOGICAL ELEMENTAL FUNCTION le_dble(x,y,tol)
    !---------------------------------------------------!
    ! Check if two floating-point numbers are such that !
    ! X <= Y within a reasonable tolerance.             !
    !---------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(in) :: x,y
    DOUBLE PRECISION,OPTIONAL,INTENT(in) :: tol
    le_dble=x<y.or.eq_dble(x,y,tol)
  END FUNCTION le_dble


  FUNCTION gaussian_random_number(w) RESULT(gauss_vector)
    !--------------------------------------------------------------!
    ! Given w(1:n), returns n random numbers distributed according !
    ! to a Gaussian of variance w**2 centred at zero.              !
    !--------------------------------------------------------------!
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(in) :: w(:)
    DOUBLE PRECISION :: gauss_vector(size(w))
    DOUBLE PRECISION v(2),rad,fac
    INTEGER i,n
    n=size(w)
    do i=1,n,2
      do
        call random_number(v)
        v(:)=2.d0*v(:)-1.d0
        rad=sum(v**2)
        if(rad<=1.d0.and.rad>0.d0)exit
      enddo
      fac=log(rad)/rad
      fac=sqrt(-fac-fac)
      gauss_vector(i)=v(2)*fac*w(i)
      if(i<n)gauss_vector(i+1)=v(1)*fac*w(i+1)
    enddo ! i
  END FUNCTION gaussian_random_number


  SUBROUTINE characterize_dist(M,A,w,mean,stderr,err_stderr,var,skew,kurt)
    !-------------------------------------------------!
    ! Evaluate the variance, skewness and kurtosis of !
    ! a data set.                                     !
    !-------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M),w(M)
    DOUBLE PRECISION, INTENT(out), OPTIONAL :: mean,stderr,err_stderr,var,&
       &skew,kurt
    DOUBLE PRECISION V1,V2,V3,V4,m1,m2,m3,m4,K2,K3,K4,num1,num2,denom
    ! Initialize.
    if(present(mean))mean=0.d0
    if(present(stderr))stderr=0.d0
    if(present(var))var=0.d0
    if(present(skew))skew=0.d0
    if(present(kurt))kurt=0.d0
    ! Compute mean.
    if(present(mean).or.present(stderr).or.present(err_stderr).or.&
       &present(var).or.present(skew).or.present(kurt))then
      V1=sum(w)
      m1=sum(w*A)/V1
      if(present(mean))mean=m1
    endif
    if(M<2)return
    ! Compute variance.
    if(present(stderr).or.present(err_stderr).or.present(var).or.&
       &present(skew).or.present(kurt))then
      V2=sum(w**2)
      m2=sum(w*(A-m1)**2)/V1
      K2=m2*(V1*V1/(V1*V1-V2))
      if(present(stderr))stderr=sqrt(max(0.d0,K2/dble(M)))
      if(present(err_stderr))err_stderr=sqrt(max(0.d0,&
         &0.5d0*K2/(dble(M)*dble(M-1))))
      if(present(var))var=K2
    endif
    if(K2<=0.d0)return
    if(M<3)return
    ! Compute skewness.
    if(present(skew).or.present(kurt))then
      V3=sum(w**3)
      m3=sum(w*(A-m1)**3)/V1
      K3=m3*(V1*V1*V1/(V1*V1*V1-3.d0*V1*V2+2.d0*V3))
      if(present(skew))skew=K3/K2**1.5d0
    endif
    if(M<4)return
    ! Compute kurtosis.
    if(present(kurt))then
      V4=sum(w**4)
      m4=sum(w*(A-m1)**4)/V1
      num1=V1*V1*(V1**4-4.d0*V1*V3+3.d0*V2*V2)
      num2=3.d0*V1*V1*(V1**4-2.d0*V1*V1*V2+4.d0*V1*V3-3.d0*V2*V2)
      denom=(V1*V1-V2)*(V1**4-6.d0*V1*V1*V2+8.d0*V1*V3+3.d0*V2*V2-6.d0*V4)
      K4=(m4*(num1/denom)-m2*m2*(num2/denom))
      kurt=K4/(K2*K2)
    endif
  END SUBROUTINE characterize_dist


  SUBROUTINE get_conf_intervals(M,A,mean,CI)
    !------------------------------------------------------------------!
    ! Given A(1:M), return the values of A from the dataset closest to !
    ! the 2.3% | 15.9% | 84.1% | 97.7% quantiles.  These correspond to !
    ! the one-sigma and two-sigma confidence intervals for a normal    !
    ! distribution (hence the significance of these particular         !
    ! quantiles).                                                      !
    !------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: M
    DOUBLE PRECISION, INTENT(in) :: A(M)
    DOUBLE PRECISION, INTENT(inout) :: mean, CI(4)
    INTEGER indx(M)
    mean=sum(A)/dble(M)
    call quicksort(M,A,indx)
    CI(1)=A(indx(floor(1.d0+0.022750131948179d0*dble(M))))
    CI(2)=A(indx(floor(1.d0+0.158655253931457d0*dble(M))))
    CI(3)=A(indx(floor(1.d0+0.841344746068543d0*dble(M))))
    CI(4)=A(indx(floor(1.d0+0.977249868051821d0*dble(M))))
  END SUBROUTINE get_conf_intervals


  SUBROUTINE quicksort(n,x,indx,preinitialized)
    !------------------------------------------------------------------!
    ! Apply the quicksort algorithm to create an integer index vector  !
    ! INDX such that the values X(INDX(I)) increase with increasing I. !
    !                                                                  !
    ! This is Leonard J. Moss' 1986 SORTX routine, see                 !
    !  http://www.fortran.com/quick_sort2.f                            !
    ! de-goto-ed by PLR, 06.2009.                                      !
    !                                                                  !
    ! In the comments the square-bracket notation 'X[I]' is used as    !
    ! shorthand for 'X(INDX(I))'.                                      !
    !                                                                  !
    ! For information about quicksort, see:                            !
    !  http://en.wikipedia.org/wiki/Quicksort                          !
    !------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    INTEGER, INTENT(inout) :: indx(n)
    DOUBLE PRECISION, INTENT(in) :: x(n)
    LOGICAL, INTENT(in), OPTIONAL :: preinitialized
    ! Size of subsequence below which we stop QuickSort and move to the final
    ! straight-insertion sorting stage.  Optimal value by 1986 standards is
    ! m =~ 9.  A simple test for array sizes ranging from 500 to 100000 on an
    ! Intel Core 2 processor with the ifort compiler seem to point to m =~ 30
    ! as the optimal value in 2009 -- the difference in timing is small,
    ! though.
    INTEGER, PARAMETER :: m=30
    ! Maximum stack size
    INTEGER, PARAMETER :: stack_size=128
    INTEGER l,r,i,j,p,lstk(stack_size),rstk(stack_size),istk
    LOGICAL already_init
    DOUBLE PRECISION xp

    already_init=.false.
    if(present(preinitialized))already_init=preinitialized

    if(.not.already_init)then
      ! Initialize INDX.
      do i=1,n
        indx(i)=i
      enddo ! i
    endif

    if(n>m)then ! Skip QuickSort for small vectors.

      ! Initialize left/right region boundaries and boundary stack, and start
      ! main loop.
      istk=0
      l=1
      r=n
      do

        ! Sort the subsequence X[L]..X[R].
        ! At this point, X[KL] <= X[KM] <= X[KR] for all KL < L, KR > R, and
        ! L <= KM <= R.  This is not applicable to the first iteration where
        ! there is no data for KL < L or KR > R.
        i=l
        j=r

        ! Let the pivot, P, be the midpoint of this subsequence, P=(L+R)/2.
        ! Then rearrange INDX(L), INDX(P), and INDX(R) so the corresponding X
        ! values are in increasing order.  The pivot key is then X[P].
        p=(l+r)/2
        if(x(indx(l))>x(indx(p)))call swap1(indx(p),indx(l))
        if(x(indx(p))>x(indx(r)))call swap1(indx(p),indx(r))
        if(x(indx(l))>x(indx(p)))call swap1(indx(p),indx(l))
        xp=x(indx(p))

        ! Inner loop over elements to sort in L..R.  We want to swap values
        ! between the right and left sides and/or move X[P] until all smaller
        ! values are left of P and all larger values are right of P.  At the
        ! end of this process neither the left or right side will be
        ! internally ordered yet, but X[P] will be in its final position.
        do
          ! Search for datum on left >= X[P] and datum on right <= X[P].
          ! At this point X[L] <= X[P] <= X[R], therefore we can start
          ! scanning up from L and down from R to find the required elements.
          i=i+1
          do while(x(indx(i))<xp)
            i=i+1
          enddo
          j=j-1
          do while(x(indx(j))>xp)
            j=j-1
          enddo
          if(i>=j)exit ! exit when the two scans collide
          call swap1(indx(i),indx(j))
        enddo ! loop over elements to sort in L..R.

        ! Select next subsequence to sort.  At this point, I >= J.  The
        ! elements in the left subsequence {X[KL], L <= KL < I} and right
        ! subsequence {X[KR], J < KR <= R} verify that
        ! X[KL] <= X[I] == X[P] <= X[KR].
        if(r-j>=i-l.and.i-l>m)then
          ! Both subsequences are more than M elements long. Push longer (left)
          ! on stack and QuickSort the shorter (right).
          istk=istk+1
          if(istk>stack_size)then
            write(6,*)'Quicksort: stack size too small <1>.'
            stop
          endif
          lstk(istk)=j+1
          rstk(istk)=r
          r=i-1
        elseif(i-l>r-j.and.r-j>m)then
          ! Both subsequences are more than M elements long. Push longer
          ! (right) on stack and QuickSort the shorter (left).
          istk=istk+1
          if(istk>stack_size)then
            write(6,*)'Quicksort: stack size too small <2>.'
            stop
          endif
          lstk(istk)=l
          rstk(istk)=i-1
          l=j+1
        elseif(r-j>m)then
          ! Only right subsequence is more than M elements long. QuickSort it.
          l=j+1
        elseif(i-l>m)then
          ! Only left subsequence is more than M elements long. QuickSort it.
          r=i-1
        else
          ! Both subsequences are less than M elements long. Pop the stack, or
          ! terminate QuickSort if empty
          if(istk<1)exit
          l=lstk(istk)
          r=rstk(istk)
          istk=istk-1
        endif

      enddo ! main loop

    endif ! n>m

    ! Final straight-insertion sorting stage.
    do i=2,n
      if(x(indx(i-1))<=x(indx(i)))cycle
      call swap1(indx(i-1),indx(i))
      do j=i-1,2,-1
        if(x(indx(j-1))<=x(indx(j)))exit
        call swap1(indx(j-1),indx(j))
      enddo ! j
    enddo ! i

  END SUBROUTINE quicksort


  SUBROUTINE swap1(x,y)
    !------------------------!
    ! Swap integers X and Y. !
    !------------------------!
    IMPLICIT NONE
    INTEGER, INTENT(inout) :: x,y
    INTEGER z
    z=x
    x=y
    y=z
  END SUBROUTINE swap1


  ! String utilities.


  FUNCTION field(n,line)
    !--------------------------------------------------------!
    ! Return the N-th field in string LINE, where the fields !
    ! are separated by one or more spaces.  An empty string  !
    ! is returned if N<1.                                    !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    CHARACTER(*),INTENT(in) :: line
    CHARACTER(len(line)) :: field
    CHARACTER(len(line)) remainder
    INTEGER i,k
    ! Initialize.
    field=''
    if(n<1)return
    remainder=adjustl(line)
    ! Loop over fields.
    i=0
    do
      i=i+1
      if(remainder(1:1)=='"')then
        ! Potentially start of a multi-word field.
        ! Locate end of field (quote followed by <space> or <EOL>).
        k=index(trim(remainder(2:)),'" ')+1
        if(k==1)then
          ! quote-space not found, so see if there is a quote at EOL.
          k=len_trim(remainder)
          if(remainder(k:k)/='"')k=1
        endif
        if(k>1)then
          ! Found end of field.
          if(i==n)then
            field=trim(remainder(2:k-1))
          else
            remainder=adjustl(remainder(k+1:))
          endif
          cycle
        endif
      endif
      ! Single-word field.
      ! Locate end of field.
      k=scan(trim(remainder),' ')
      if(k==0)then
        ! End is EOL.
        if(i==n)field=trim(remainder)
        return
      elseif(i==n)then
        field=trim(remainder(1:k-1))
        return
      else
        remainder=adjustl(remainder(k:))
      endif
    enddo
  END FUNCTION field


  INTEGER FUNCTION nfield(line)
    !--------------------------------------!
    ! Return the number of fields in LINE. !
    !--------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: line
    nfield=0
    do
      if(len_trim(field(nfield+1,line))==0)exit
      nfield=nfield+1
    enddo
  END FUNCTION nfield


  DOUBLE PRECISION FUNCTION dble_field(ifield,command,ierr)
    !--------------------------------------------------------!
    ! Like field, but returning the value as an real number. !
    ! ierr is set to a non-zero value if the requested field !
    ! could not be parsed as a real number.                  !
    !--------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: ifield
    CHARACTER(*),INTENT(in) :: command
    INTEGER,INTENT(inout) :: ierr
    dble_field=parse_dble(field(ifield,command),ierr)
  END FUNCTION dble_field


  DOUBLE PRECISION FUNCTION parse_dble(string,ierr)
    !-----------------------------------------------------!
    ! Parse a string to obtain a double-precision number. !
    ! Supports fractions, e.g., string="1/3".             !
    !-----------------------------------------------------!
    IMPLICIT NONE
    CHARACTER(*),INTENT(in) :: string
    INTEGER,INTENT(inout) :: ierr
    INTEGER ipos
    DOUBLE PRECISION t1,t2
    ipos=scan(string,'/')
    if(ipos==0)then
      read(string,*,iostat=ierr)parse_dble
    else
      ierr=-1
      if(ipos<2)return
      if(ipos>=len_trim(string))return
      read(string(1:ipos-1),*,iostat=ierr)t1
      if(ierr/=0)return
      read(string(ipos+1:),*,iostat=ierr)t2
      if(ierr/=0)return
      if(eq_dble(t2,0.d0))then
        ierr=-1
        return
      endif
      parse_dble=t1/t2
    endif
  END FUNCTION parse_dble


  CHARACTER(12) FUNCTION i2s(n)
    !-----------------------------------------------------------------------!
    ! Convert integers to left justified strings that can be printed in the !
    ! middle of a sentence without introducing large amounts of white space.!
    !-----------------------------------------------------------------------!
    IMPLICIT NONE
    INTEGER,INTENT(in) :: n
    INTEGER i,j
    INTEGER,PARAMETER :: ichar0=ichar('0')
    i2s=''
    i=abs(n)
    do j=len(i2s),1,-1
      i2s(j:j)=achar(ichar0+mod(i,10))
      i=i/10
      if(i==0)exit
    enddo ! j
    if(n<0)then
      i2s='-'//adjustl(i2s(2:12))
    else
      i2s=adjustl(i2s)
    endif ! n<0
  END FUNCTION i2s


  ! Generic utilities.


  SUBROUTINE quit (msg)
    !---------------------!
    ! Quit with an error. !
    !---------------------!
    IMPLICIT NONE
    CHARACTER(*), INTENT(in), OPTIONAL :: msg
    if (present(msg)) then
      write(6,*)'ERROR : '//msg
    else
      write(6,*)'Quitting.'
    endif
    stop
  END SUBROUTINE quit


END PROGRAM vib_levels
