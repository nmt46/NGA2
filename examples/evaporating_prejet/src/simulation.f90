!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg1,cfg2,amGrp1,amGrp2,grp1,grp2
   use lpt_class,         only: lpt
   use lowmach_class,     only: lowmach
   use incomp_class,      only: incomp
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use fluidTable_class,  only: fluidTable
   use sgsmodel_class,    only: sgsmodel
   use coupler_class,     only: coupler
   implicit none
   private
   
   !> Get solvers for spray evaporation domain
   type(lpt),         public :: lp
   type(lowmach),     public :: fs1
   type(vdscalar),    public :: T_sc, Yf_sc
   type(sgsmodel),    public :: sgs1
   type(timetracker), public :: time1

   !> Get solvers for jet development domain
   type(incomp),     public :: fs2
   type(sgsmodel),    public :: sgs2
   type(timetracker), public :: time2

   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out1
   type(ensight)  :: ens_out2
   type(event)    :: ens_evt1
   type(event)    :: ens_evt2

   !> Coupler
   type(coupler)  :: cpl
   
   !> Simulation monitor file
   type(monitor) :: mfile,pfile,cflFile,jetFile

   !> Monitoring quantities
   
   public :: simulation_init,simulation_run,simulation_final
   real(WP) :: meanU
   !> Case quantities

   real(WP) :: p0 ! ambient pressure
   real(WP) :: diff_T,diff_Yf ! non-turbulent diffusivities for scalar solvers

   real(WP) :: dp,dp_sig,rp,rp_sig,u_inj,u_inj_sig,T_d,inj_rate,r_jet,u_jet,T_jet,Yf_jet,u_coflow,T_coflow,Yf_coflow
   integer :: nInj

   !> Quantities for flow & scalar solvers:
   real(WP), dimension(:,:,:),   allocatable :: resU1,resV1,resW1,resT,resYf
   real(WP), dimension(:,:,:),   allocatable :: Ui1,Vi1,Wi1
   real(WP), dimension(:,:,:,:), allocatable :: SR1

   real(WP), dimension(:,:,:),   allocatable :: resU2,resV2,resW2
   real(WP), dimension(:,:,:),   allocatable :: Ui2,Vi2,Wi2
   real(WP), dimension(:,:,:,:), allocatable :: SR2

   real(WP), dimension(:,:),     allocatable :: passUVW

   type(fluidTable) :: lTab,gTab


   
contains
   
   ! function xm_locator(pg,i,j,k) result(isIn)
   !    use pgrid_class, only: pgrid
   !    implicit none
   !    class(pgrid), intent(in) :: pg
   !    integer, intent(in) :: i,j,k
   !    logical :: isIn
   !    isIn = .false.
   !    if (i.eq.pg%imin) isIn = .true.
   ! end function xm_locator

   function xp_locator2(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if ((i.eq.pg%imax+1).and.(sqrt(pg%ym(j)**2+pg%zm(k)**2)).le.r_jet) isIn = .true.
   end function xp_locator2

   function xm_locator2(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imin) isIn = .true.
   end function xm_locator2

   function gas_inj_locator2(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imin.and.(sqrt(pg%ym(j)**2+pg%zm(k)**2)).le.r_jet) isIn = .true.
   end function gas_inj_locator2


   function xp_locator1(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imax+1) isIn = .true.
   end function xp_locator1

   function xp_sc_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imax) isIn = .true.
   end function xp_sc_locator

   function ym_sc_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j.eq.pg%jmin) isIn = .true.
   end function ym_sc_locator

   function yp_sc_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (j.eq.pg%jmax) isIn = .true.
   end function yp_sc_locator

   function zm_sc_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k.eq.pg%kmin) isIn = .true.
   end function zm_sc_locator

   function zp_sc_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (k.eq.pg%kmax) isIn = .true.
   end function zp_sc_locator
   

   function gas_inj_locator1(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imin.and.(sqrt(pg%ym(j)**2+pg%zm(k)**2)).le.r_jet) isIn = .true.
   end function gas_inj_locator1

   function coflow_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imin.and.(sqrt(pg%ym(j)**2+pg%zm(k)**2)).gt.r_jet) isIn = .true.
   end function coflow_locator

      !> Define here our equation of state - rho(T,mass)
   subroutine get_rho()
      implicit none
      ! real(WP), intent(in) :: mass
      integer :: i,j,k
      real(WP) :: W_g,W_l,R_cst
      R_cst = 8.314472_WP ! Gas constant [J/(kg*K)]
      W_g = gTab%MW ! Carrier gas molar weight [kg/mol]
      W_l = lTab%MW ! Drop molar weight [kg/mol]

      ! Calculate density
      do k=T_sc%cfg%kmino_,T_sc%cfg%kmaxo_
         do j=T_sc%cfg%jmino_,T_sc%cfg%jmaxo_
            do i=T_sc%cfg%imino_,T_sc%cfg%imaxo_
               if (Yf_sc%SC(i,j,k).ge.0.0_WP) then
                  T_sc%rho(i,j,k)=p0/(R_cst*(Yf_sc%SC(i,j,k)/W_l+(1.0_WP-Yf_sc%SC(i,j,k))/W_g)*T_sc%SC(i,j,k))! + lp%srcM(i,j,k)
               else
                  T_sc%rho(i,j,k)=p0*W_g/(R_cst*T_sc%SC(i,j,k))! + lp%srcM(i,j,k)
               end if
               Yf_sc%rho(i,j,k)=T_sc%rho(i,j,k)
            end do
         end do
      end do
   end subroutine get_rho

   !> Initialization of problem solver
   subroutine simulation_init
      use coolprop
      use param, only: param_read
      implicit none

      read_global_vars : block
         call param_read('Gas jet velocity',             u_jet); meanU = u_jet
         call param_read('Gas coflow velocity',          u_coflow)
         call param_read('Gas jet radius',               r_jet)
         call param_read('Pressure',                     p0)
         call param_read('Gas jet temperature',          T_jet)
         call param_read('Gas coflow temperature',       T_coflow)
         call param_read('Gas jet fuel mass fraction',   Yf_jet)
         call param_read('Gas coflow fuel mass fraction',Yf_coflow)

         call param_read('Droplet mean diameter',              dp)
         call param_read('Droplet diameter standard deviation',dp_sig)
         call param_read('Droplet mean radius',                rp)
         call param_read('Droplet radius standard deviation',  rp_sig)
         call param_read('Droplet temperature',                T_d)
         call param_read('Droplet mean velocity',              u_inj)
         call param_read('Droplet velocity standard deviation',u_inj_sig)
         call param_read('Droplet injection rate',             inj_rate)

         call param_read('Thermal diffusivity',diff_T)
         call param_read('Fuel diffusivity',   diff_Yf)
      end block read_global_vars

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!! Initialize fluid properties !!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      initialize_fluid_properties: block
         ! use messager, only : die
         use string,    only: str_medium
         use fluidTable_class, only: Cp_ID,Lv_ID,Tb_ID,rho_ID,MW_ID,mu_ID
         character(len=str_medium) :: name
         integer :: i,j,nP,nT
         real(WP) :: Tmin,Tmax,Pmin,Pmax,dP,dT
         real(WP) :: testVal

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Liquid !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call param_read('Fuel',name)
         nP = 9; nT = 9 ! points on interpolation mesh

         lTab=fluidTable(name=name)

         Tmin = maxval([250.0_WP,1.01_WP*cprop(output='Tmin'//char(0),name1='P'//char(0),prop1=p0,name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0))]); 
         Tmax = 0.999_WP*cprop(output='T'//char(0),name1='P'//char(0),prop1=p0,name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0))
         Pmin = p0/1.1_WP; Pmax = p0*1.1_WP

         call lTab%initialize_fluidTable(Pmin=Pmin,Pmax=Pmax,Tmin=Tmin,Tmax=Tmax,nP=nP,nT=nT)

         call lTab%addProp(propID=Cp_ID)
         call lTab%addProp(propID=Lv_ID)
         call lTab%addProp(propID=Tb_ID)
         call lTab%addProp(propID=rho_ID)
         call lTab%addProp(propID=MW_ID)
         call lTab%addProp(propID=mu_ID)


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


         call param_read('Gas',name)

         gTab=fluidTable(name=name)

         Pmin = p0/1.1_WP; Pmax = p0*1.1_WP
         dP = (Pmax-Pmin)/(real(gTab%nP,WP)-1.0_WP)

         Tmin = maxval([250.0_WP,1.01_WP*cprop(output='Tmin'//char(0),name1='P'//char(0),prop1=p0,name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0))])
         Tmax = T_jet*1.5

         call gTab%initialize_fluidTable(Pmin=Pmin,Pmax=Pmax,Tmin=Tmin,Tmax=Tmax,nP=nP,nT=nT)

         call gTab%addProp(propID=Cp_ID)
         call gTab%addProp(propID=Lv_ID)
         call gTab%addProp(propID=Tb_ID)
         call gTab%addProp(propID=rho_ID)
         call gTab%addProp(propID=MW_ID)
         call gTab%addProp(propID=mu_ID)

         call param_read('Gas Prandtl number',gTab%Pr)
         call param_read('Fuel-gas Schmidt number',gTab%Sc)

      end block initialize_fluid_properties

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Create Coupler !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! Both groups prepare the coupler
      if (amGrp1.or.amGrp2) then
         ! Create the coupler
         cpl=coupler(src_grp=grp2,dst_grp=grp1,name='jet_coupler')
         ! Set the grids
         if (amGrp2) call cpl%set_src(cfg2,'x')
         if (amGrp1) call cpl%set_dst(cfg1,'x')
         ! Initialize the metrics
         call cpl%initialize()
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!! Initialize jet development solvers !!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (amGrp2) then
         initialize_timetracker2: block
            time1=timetracker(amRoot=cfg2%amRoot)
            call param_read('Max timestep size',time2%dtmax)
            call param_read('Max cfl number',time2%cflmax)
            call param_read('Max time',time2%tmax)
            time2%dt=min(time2%dtmax,0.1_WP*cfg2%dx(1)/u_jet)
            time2%itmax=2
         end block initialize_timetracker2



         initialize_res2: block
            allocate(resU2 (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); resU2=0.0_WP
            allocate(resV2 (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); resV2=0.0_WP
            allocate(resW2 (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); resV2=0.0_WP
            allocate(Ui2   (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); Ui2=0.0_WP
            allocate(Vi2   (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); Vi2=0.0_WP
            allocate(Wi2   (cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); Wi2=0.0_WP
            allocate(SR2 (6,cfg2%imino_:cfg2%imaxo_,cfg2%jmino_:cfg2%jmaxo_,cfg2%kmino_:cfg2%kmaxo_)); SR2=0.0_WP
         end block initialize_res2

         initialize_fs2: block
            use mathtools, only: twoPi
            use random, only : random_normal
            use incomp_class, only: clipped_neumann,dirichlet,bcond
            use ils_class, only: pcg_amg,pcg_pfmg,gmres_amg
            use fluidTable_class, only: mu_ID
            type(bcond),pointer :: mybc
            integer :: i,j,k,n
            real(WP) :: rhof,viscf,R_cst,r
            fs2=incomp(cfg=cfg2,name='Incompressible N-S')

            ! Setup bconds
            ! call fs2%add_bcond(name='bc_xp',  type=clipped_neumann,face='x',dir=+1,canCorrect=.true., locator=xp_locator2)
            ! call fs2%add_bcond(name='gas_inj',type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=gas_inj_locator2)

            ! Set initial density, viscosity
            R_cst = 8.314472_WP ! Gas constant [J/(kg*K)]
            fs2%rho = p0*gTab%MW/(R_cst*T_jet) ! Ideal gas law
            call gTab%evalProps(propOut=viscf,T_q=T_jet,P_q=p0,propID=mu_ID); fs2%visc=viscf

            ! Configure pressure solver
            call param_read('Pressure iteration',fs2%psolv%maxit)
            call param_read('Pressure tolerance',fs2%psolv%rcvg)
            ! Configure implicit velocity solver
            call param_read('Implicit iteration',fs2%implicit%maxit)
            call param_read('Implicit tolerance',fs2%implicit%rcvg)
            ! Setup the solver
            fs2%psolv%maxlevel=16
            call fs2%setup(pressure_ils=pcg_amg,implicit_ils=gmres_amg) ! had pressure as pcg_amg

            ! Initialize velocity field
            fs2%U = 0.0_WP
            fs2%V = 0.0_WP
            fs2%W = 0.0_WP
            do k=cfg2%kmino_,cfg2%kmaxo_
               do j=cfg2%jmino_,cfg2%jmaxo_
                  do i=cfg2%imino_,cfg2%imaxo_
                     if (cfg2%VF(i,j,k).eq.1.0_WP) fs2%U(i,j,k)=u_jet!*(1.0_WP+0.1_WP*cos(twoPi*3*sqrt(cfg2%ym(j)**2+cfg2%zm(k)**2)/r_jet))!*(1+random_normal(0.0_WP,u_jet/100.0_WP))
                     ! if (cfg2%VF(i,j,k).eq.1.0_WP) then
                     !    r = sqrt(r_jet*sin((cfg2%xm(i)-cfg2%xm(cfg2%nx/2))/r_jet)**2+cfg2%ym(j)**2+cfg2%zm(k)**2)/r_jet
                     !    fs2%U(i,j,k) = u_jet*(1.0_WP+2.0_WP*exp(-r*5.0_WP)*cfg2%ym(j)/r_jet)
                     !    fs2%V(i,j,k) = u_jet*(      -2.0_WP*exp(-r*5.0_WP)*sin((cfg2%xm(i)-cfg2%xm(cfg2%nx/2))/r_jet))
                     ! end if
                  end do
               end do
            end do

            ! call fs2%get_bcond('gas_inj',mybc)
            ! do n=1,mybc%itr%no_
            !    i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            !    fs2%U(i,j,k) = u_jet
            ! end do
            call fs2%apply_bcond(time2%t,time2%dt)
            ! Compute MFR through all boundary conditions
            call fs2%get_mfr()
            call fs2%correct_mfr()
            call fs2%interp_vel(Ui2,Vi2,Wi2)
            call fs2%get_div()
         end block initialize_fs2

         ! Create second LES model
         create_sgs2: block
            sgs2=sgsmodel(cfg=fs2%cfg,umask=fs2%umask,vmask=fs2%vmask,wmask=fs2%wmask)
         end block create_sgs2
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!! Initialize spray evaporation solvers !!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (amGrp1) then
         ! Initialize time tracker with 1 subiterations
         initialize_timetracker1: block
            time1=timetracker(amRoot=cfg1%amRoot)
            call param_read('Max timestep size',time1%dtmax)
            call param_read('Max cfl number',time1%cflmax)
            call param_read('Max time',time1%tmax)
            time1%dt=0.1_WP*cfg1%dx(1)/u_jet
            ! time1%dt=time1%dtmax
            time1%itmax=2
         end block initialize_timetracker1
      
         ! Initialize our LPT
         initialize_lpt: block
            use fluidTable_class, only: rho_ID
            use mathtools, only: twoPi
            use random, only: random_normal,random_uniform
            real(WP) :: theta
            integer :: i,np
            ! Create solver
            lp=lpt(cfg=cfg1,name='LPT')         
            ! Get droplet density from the fluid table
            call lTab%evalProps(propOut=lp%rho,T_q=T_d,P_q=p0,propID=rho_ID)
            ! Set filter scale to 3.5*dx
            lp%filter_width=3.5_WP*cfg1%min_meshsize
            np = 0!int(floor(inj_rate*time1%dt))
            nInj = np
            ! Root process initializes 1000 particles randomly
            if (lp%cfg%amRoot) then
               call lp%resize(np)
               do i=1,np
                  ! Give id
                  lp%p(i)%id=int(i,8)
                  lp%p(i)%d=random_normal(dp,dp_sig)
                  
                  ! Set the temperature 
                  lp%p(i)%T_d = T_d

                  !!!!!!!!! Assign position !!!!!!!!!!
                  ! For ring distribution:
                  if (cfg1%ny.eq.1) then
                     theta = 2.0_WP*floor(random_uniform(0.0_WP,2.0_WP))-1.0_WP ! Use theta as either -1 or 1
                     lp%p(np)%pos = random_normal(rp,rp_sig)*[0.0_WP,0.0_WP,theta]
                  elseif (cfg1%nz.eq.1) then
                     theta = 2.0_WP*floor(random_uniform(0.0_WP,2.0_WP))-1.0_WP ! Use theta as either -1 or 1
                     lp%p(np)%pos = random_normal(rp,rp_sig)*[0.0_WP,theta,0.0_WP]
                  else
                     theta = random_uniform(0.0_WP,twoPi)
                     lp%p(np)%pos = random_normal(rp,rp_sig)*[0.0_WP,cos(theta),sin(theta)]
                  end if



                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Give zero velocity
                  lp%p(i)%vel=[random_normal(u_inj,u_inj_sig),0.0_WP,0.0_WP]
                  ! Give zero collision force
                  lp%p(i)%col=0.0_WP
                  ! Give zero dt
                  lp%p(i)%dt=0.0_WP
                  ! Locate the particle on the mesh
                  lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                  ! Activate the particle
                  lp%p(i)%flag=0
                  ! print*,'ID:',lp%p(i)%id,'d',lp%p(i)%d,'pos',lp%p(i)%pos,'vel',lp%p(i)%vel
               end do
            end if
            ! Distribute particles
            call lp%sync()
            ! Get initial particle volume fraction
            call lp%update_VF()
            ! Set collision timescale
            lp%Tcol=5.0_WP*time1%dt
         end block initialize_lpt
         

         initialize_res1: block
            allocate(resU1 (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); resU1=0.0_WP
            allocate(resV1 (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); resV1=0.0_WP
            allocate(resW1 (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); resW1=0.0_WP
            allocate(resT  (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); resT=0.0_WP
            allocate(resYf (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); resYf=0.0_WP
            allocate(Ui1   (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); Ui1=0.0_WP
            allocate(Vi1   (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); Vi1=0.0_WP
            allocate(Wi1   (cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); Wi1=0.0_WP
            allocate(SR1 (6,cfg1%imino_:cfg1%imaxo_,cfg1%jmino_:cfg1%jmaxo_,cfg1%kmino_:cfg1%kmaxo_)); SR1=0.0_WP
         end block initialize_res1
         
         ! Initialize scalar solvers
         initialize_sc: block
            use ils_class,      only: gmres_amg
            use vdscalar_class, only: neumann,dirichlet,bcond
            type(bcond), pointer :: mybc
            integer :: i,j,k,n
            ! create the scalar solvers for temperature and fuel mass fraction
            T_sc=vdscalar(cfg=cfg1,scheme=1,name='Temperature solver')
            Yf_sc=vdscalar(cfg=cfg1,scheme=1,name='Fuel mass fraction solver')
            ! Configure implicit scalar solver
            call param_read('Implicit iteration',T_sc%implicit%maxit)
            call param_read('Implicit tolerance',T_sc%implicit%rcvg)
            Yf_sc%implicit%maxit=T_sc%implicit%maxit; YF_sc%implicit%rcvg=T_sc%implicit%rcvg
            ! Create bconds
            ! call T_sc%add_bcond(name='bc_xm',type=neumann,dir='xm',locator=xm_locator)
            ! call T_sc%add_bcond(name='chop_out',type=dirichlet,dir='xp',locator=xp_sc_locator)
            call T_sc%add_bcond(name='T_out',   type=neumann,  dir='xp',locator=xp_sc_locator)
            call T_sc%add_bcond(name='T_ym',    type=neumann,  dir='ym',locator=ym_sc_locator)
            call T_sc%add_bcond(name='T_yp',    type=neumann,  dir='yp',locator=yp_sc_locator)
            call T_sc%add_bcond(name='T_zm',    type=neumann,  dir='zm',locator=zm_sc_locator)
            call T_sc%add_bcond(name='T_zp',    type=neumann,  dir='zp',locator=zp_sc_locator)
            call T_sc%add_bcond(name='T_jet',   type=dirichlet,dir='xm',locator=gas_inj_locator1)
            call T_sc%add_bcond(name='T_coflow',type=dirichlet,dir='xm',locator=coflow_locator)

            call Yf_sc%add_bcond(name='Yf_out',   type=neumann, dir='xp',locator=xp_sc_locator)
            call Yf_sc%add_bcond(name='Yf_ym',    type=neumann,  dir='ym',locator=ym_sc_locator)
            call Yf_sc%add_bcond(name='Yf_yp',    type=neumann,  dir='yp',locator=yp_sc_locator)
            call Yf_sc%add_bcond(name='Yf_zm',    type=neumann,  dir='zm',locator=zm_sc_locator)
            call Yf_sc%add_bcond(name='Yf_zp',    type=neumann,  dir='zp',locator=zp_sc_locator)
            call Yf_sc%add_bcond(name='Yf_jet',   type=dirichlet,dir='xm',locator=gas_inj_locator1)
            call Yf_sc%add_bcond(name='Yf_coflow',type=dirichlet,dir='xm',locator=coflow_locator)
            ! call Yf_sc%add_bcond(name='chop_out',type=dirichlet,dir='xp',locator=xp_sc_locator)
            ! call Yf_sc%add_bcond(name='bc_xm', type=neumann,  dir='xm',locator=xm_locator)

            ! Setup the solver
            call T_sc%setup(implicit_ils=gmres_amg)
            call Yf_sc%setup(implicit_ils=gmres_amg)
            ! Set initial field values
            T_sc%SC=T_coflow
            Yf_sc%SC=Yf_coflow
            call T_sc%get_bcond('T_jet',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               T_sc%SC(i,j,k)=T_jet
            end do
            call T_sc%get_bcond('T_coflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               T_sc%SC(i,j,k)=T_coflow
            end do
            call Yf_sc%get_bcond('Yf_jet',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               Yf_sc%SC(i,j,k)=Yf_jet
            end do
            call Yf_sc%get_bcond('Yf_coflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               Yf_sc%SC(i,j,k)=Yf_coflow
            end do
            call get_rho()
            ! Set diffusivities
            T_sc%diff = diff_T
            Yf_sc%diff = diff_Yf
         end block initialize_sc

         ! Initialize the flow solver
         initialize_fs1: block
            use mathtools, only: twoPi
            use lowmach_class, only: clipped_neumann,dirichlet,bcond
            use ils_class, only: pcg_amg,pfmg,gmres_amg
            use fluidTable_class, only: mu_ID
            type(bcond),pointer :: mybc
            integer :: i,j,k,n
            real(WP) :: rhof,viscf,viscf2
            fs1=lowmach(cfg=cfg1,name='Variable density low Mach NS')
            ! Setup bconds
            ! call fs1%add_bcond(name='bc_xm',type=neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
            call fs1%add_bcond(name='bc_xp',    type=clipped_neumann,face='x',dir=+1,canCorrect=.true., locator=xp_locator1)
            call fs1%add_bcond(name='gas_inj',  type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=gas_inj_locator1)
            call fs1%add_bcond(name='gas_coflow',type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=coflow_locator)

            ! Set initial density, viscosity
            fs1%rho = Yf_sc%rho ! Take density from scalar solver
            call gTab%evalProps(propOut=viscf,T_q=T_coflow,P_q=p0,propID=mu_ID); fs1%visc=viscf
            ! Configure pressure solver
            call param_read('Pressure iteration',fs1%psolv%maxit)
            call param_read('Pressure tolerance',fs1%psolv%rcvg)
            ! Configure implicit velocity solver
            call param_read('Implicit iteration',fs1%implicit%maxit)
            call param_read('Implicit tolerance',fs1%implicit%rcvg)
            ! Setup the solver
            call fs1%setup(pressure_ils=pcg_amg,implicit_ils=gmres_amg)

            ! Initialize velocity field
            do k=lp%cfg%kmino_,lp%cfg%kmaxo_
               do j=lp%cfg%jmino_,lp%cfg%jmaxo_
                  do i=lp%cfg%imino_,lp%cfg%imaxo_
                     fs1%U(i,j,k)=0.0_WP;fs1%rhoU(i,j,k)=0.0_WP
                     fs1%V(i,j,k)=0.0_WP;fs1%rhoV(i,j,k)=0.0_WP
                     fs1%W(i,j,k)=0.0_WP;fs1%rhoW(i,j,k)=0.0_WP
                  end do
               end do
            end do
            call fs1%get_bcond('gas_inj',mybc)
            !!!!!!!!!!!!!!!!! Allocate storage of passed velocties !!!!!!!!!!!!!!!!!
            allocate(passUVW(3,mybc%itr%no_))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs1%rhoU(i,j,k)=u_jet*fs1%rho(i,j,k)
               fs1%U(i,j,k)   =u_jet
            end do
            call fs1%get_bcond('gas_coflow',mybc)
            do n=1,mybc%itr%no_
               i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
               fs1%rhoU(i,j,k)=u_coflow*fs1%rho(i,j,k)
               fs1%U(i,j,k)   =u_coflow
            end do
            call fs1%apply_bcond(time1%t,time1%dt)
            call T_sc%get_drhodt(drhodt=resT,dt=time1%dt) ! use resT to hold our drhodt :)
            call fs1%get_div(drhodt=resT)
            ! Compute MFR through all boundary conditions
            call fs1%correct_mfr(drhodt=resT)
            call fs1%rho_divide()
            call fs1%interp_vel(Ui1,Vi1,Wi1)
         end block initialize_fs1

         ! Create an LES model
         create_sgs1: block
            sgs1=sgsmodel(cfg=fs1%cfg,umask=fs1%umask,vmask=fs1%vmask,wmask=fs1%wmask)
         end block create_sgs1
      end if
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!! Initialize outputs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      ! Create partmesh object for Lagrangian particle output
      create_pmesh: block
         integer :: i
         pmesh=partmesh(nvar=2,name='lpt')
         pmesh%varname(1)='radius'
         pmesh%varname(2)='temperature'
         call lp%update_partmesh(pmesh)
         do i = 1,lp%np_
            pmesh%var(1,i)=0.5_WP*lp%p(i)%d ! Droplet Radius
            pmesh%var(2,i)=lp%p(i)%T_d ! Droplet Temperature
         end do
      end block create_pmesh     
      if (amGrp1) then
         ! Add Ensight output
         create_ensight1: block
            ! Create Ensight output from cfg1
            ens_out1=ensight(cfg=lp%cfg,name='spray')
            ! Create event for Ensight output
            ens_evt1=event(time=time1,name='Ensight output')
            call param_read('Ensight output period',ens_evt1%tper)
            ! Add variables to output
            call ens_out1%add_particle('particles',pmesh)
            call ens_out1%add_vector('velocity',fs1%U,fs1%V,fs1%W)
            call ens_out1%add_vector('srcUVW',lp%srcU,lp%srcV,lp%srcW)
            ! call ens_out1%add_scalar('epsp',lp%VF)
            call ens_out1%add_scalar('yf',Yf_sc%SC)
            call ens_out1%add_scalar('temperature',T_sc%SC)
            call ens_out1%add_scalar('density',T_sc%rho)
            call ens_out1%add_scalar('srcM',lp%srcM)
            call ens_out1%add_scalar('srcE',lp%srcE)

            ! Output to ensight
            if (ens_evt1%occurs()) call ens_out1%write_data(time1%t)
         end block create_ensight1
      end if
      if (amGrp2) then
         create_ensight2: block
            ! Create Ensight output from cfg1
            ens_out2=ensight(cfg=cfg2,name='jet_dev')
            ! Create event for Ensight output
            ens_evt2=event(time=time2,name='Ensight output')
            call param_read('Ensight output period',ens_evt2%tper)
            ! Add variables to output
            call ens_out2%add_vector('velocity',fs2%U,fs2%V,fs2%W)
            ! Output to ensight
            if (ens_evt2%occurs()) call ens_out2%write_data(time1%t)
         end block create_ensight2
      end if
      if (amGrp1) then
         ! Create a monitor file
         create_monitor1: block
            ! Prepare some info about fields
            call lp%get_max()
            call fs1%get_max()
            call Yf_sc%get_max()
            call T_sc%get_max()
            call fs1%get_cfl(time1%dt,time1%cfl)
            ! Create simulation monitor
            mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
            call mfile%add_column(time1%n,'N_time')
            call mfile%add_column(time1%wt,'wall_t')
            call mfile%add_column(time1%t,'time')
            call mfile%add_column(time1%dt,'dt')
            call mfile%add_column(time1%cfl,'CFL')
            call mfile%add_column(fs1%Umax,'FS_Max_U')
            call mfile%add_column(fs1%Vmax,'FS_Max_V')
            ! call mfile%add_column(fs1%Wmax,'FS_Max_W')
            call mfile%add_column(Yf_sc%SCmin,'Yf_min')
            call mfile%add_column(Yf_sc%SCmax,'Yf_max')
            call mfile%add_column(T_sc%SCmin,'T_min')
            call mfile%add_column(T_sc%SCmax,'T_max')
            call mfile%add_column(lp%np,'N_part')
            call mfile%add_column(nInj,'N_inj')
            call mfile%add_column(lp%nEvap,'N_evap')
            ! call mfile%add_column(lp%VFmean,'Mean_VF')
            ! call mfile%add_column(lp%Umin,'P_Umin')
            call mfile%add_column(lp%Umean,'P_Umean')
            ! call mfile%add_column(lp%Xmin,'P_Xmin')
            ! call mfile%add_column(lp%Xmax,'P_Xmax')
            call mfile%add_column(lp%Umax,'P_Umax')
            call mfile%add_column(lp%dmean,'P_dmean')
            call mfile%add_column(lp%Tmean,'P_Tmean')
            call mfile%add_column(lp%Tmax,'P_Tmax')
            call mfile%write()

            cflFile = monitor(amroot=cfg1%amRoot,name='cfl')
            call cflFile%add_column(time1%n,'N_time')
            call cflFile%add_column(time1%t,'time')
            call cflFile%add_column(max(fs1%CFLc_x,fs1%CFLc_y,fs1%CFLc_z),'max_CFLc')
            call cflFile%add_column(max(fs1%CFLv_x,fs1%CFLv_y,fs1%CFLv_z),'max_CFLv')
            call cflFile%write()
            ! pfile=monitor(amroot=lp%cfg%amRoot,name='particle')
            ! call pfile%add_column(time1%n,'T_step')
            ! call pfile%add_column(time1%t,'Time')
            ! call pfile%add_column(time1%dt,'dt')
            ! if (lp%cfg%amRoot) then
            !    call pfile%add_column(lp%p(1)%vel(1),'Part_U')
            !    call pfile%add_column(lp%p(1)%vel(2),'Part_V')
            !    call pfile%add_column(lp%p(1)%vel(3),'Part_W')
            !    call pfile%add_column(lp%p(1)%d,'Part_d')
            !    call pfile%add_column(lp%p(1)%T_d,'Part_T')
            ! end if
            ! call pfile%write()
         end block create_monitor1
      end if
      if (amGrp2) then
         create_monitor2 : block
            call fs2%get_max()
            call fs2%get_cfl(time2%dt,time2%cfl)
            jetFile=monitor(amroot=cfg2%amRoot,name='jetDev')
            call jetfile%add_column(time2%n,'N_time')
            call jetfile%add_column(time2%wt,'wall_t')
            call jetfile%add_column(time2%t,'time')
            call jetfile%add_column(time2%dt,'dt')
            call jetfile%add_column(time2%cfl,'CFL')
            call jetfile%add_column(fs2%Umax,'FS_Max_U')
            call jetfile%add_column(fs2%Vmax,'FS_Max_V')
            call jetfile%add_column(fs2%Wmax,'FS_Max_W')
            call jetfile%add_column(meanU,'meanU')
            call jetFile%write()

         end block create_monitor2
      end if

      create_coupler : block


      end block create_coupler

   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      logical :: amDone = .false.
      ! print*,'max rho1',maxval(Yf_sc%rho)
      ! Perform time integration
      do while (.not.time1%done().and.(.not.amDone))
         
         ! Increment time         
         call fs1%get_cfl(time1%dt,time1%cfl)
         call time1%adjust_dt()
         ! if (time1%dt.gt.T_dt_cur) time1%dt = T_dt_cur ! modify dt given temperature change on mesh
         call time1%increment()
         ! print*,'New dt:',time1%dt


         ! Update velocity in development region (domain 2)
         if (amGrp2) then
            do while (time2%t.lt.time1%t)
               ! print*,'Jet Development Solver, t=',time2%t
               ! Increment time
               call fs2%get_cfl(time2%dt,time2%cfl)
               call time2%adjust_dt()
               call time2%increment()
               
               ! Remember old velocity
               fs2%Uold=fs2%U
               fs2%Vold=fs2%V
               fs2%Wold=fs2%W

               ! Apply time-varying Dirichlet conditions
               ! This is where time-dpt Dirichlet would be enforced
               
               update_visc2: block
                  use fluidTable_class, only: mu_ID
                  real(WP) :: mu_normal
                  integer :: i,j,k
                  ! Get the turbulent viscosity
                  call fs2%interp_vel(Ui2,Vi2,Wi2)
                  call fs2%get_strainrate(Ui=Ui2,Vi=Vi2,Wi=Wi2,SR=SR2)
                  resU2=fs2%rho ! Abuse resU2 :O
                  call sgs2%get_visc(dt=time2%dtold,rho=resU2,Ui=Ui2,Vi=Vi2,Wi=Wi2,SR=SR2)
                  ! Get the normal dynamic viscosity
                  call gTab%evalProps(propOut=mu_normal,T_q=T_jet,P_q=p0,propID=mu_ID)
                  ! Loop some things
                  do k=fs2%cfg%kmin_,fs2%cfg%kmax_
                     do j=fs2%cfg%jmin_,fs2%cfg%jmax_
                        do i=fs2%cfg%imin_,fs2%cfg%imax_
                           ! Apply the turbulent viscosity
                           fs2%visc(i,j,k) = mu_normal + sgs2%visc(i,j,k)
                        end do
                     end do
                  end do
               end block update_visc2
               ! print*,'here0'
               ! Perform sub-iterations
               do while (time2%it.le.time2%itmax)
                  
                  ! Build mid-time velocity
                  fs2%U=0.5_WP*(fs2%U+fs2%Uold)
                  ! print*,'maxU',maxval(fs2%U),'minU',minval(fs2%U)
                  fs2%V=0.5_WP*(fs2%V+fs2%Vold)
                  ! print*,'maxV',maxval(fs2%V),'minV',minval(fs2%V)
                  fs2%W=0.5_WP*(fs2%W+fs2%Wold)
                  ! print*,'maxW',maxval(fs2%W),'minW',minval(fs2%W)
                  ! print*,'here1',time2%it
                  ! Explicit calculation of drho*u/dt from NS
                  call fs2%get_dmomdt(resU2,resV2,resW2)
                  ! print*,'here2',time2%it
                  ! Assemble explicit residual
                  ! print*,resW2
                  ! print*,'maxval(U-Uold)',maxval(fs2%V-fs2%Vold),'dt',time2%dt,'max,minval(resU2)',maxval(resV2)!,minval(resV2)
                  ! print*,'max,minval(fs2%rho*fs2%V-fs2%rho*fs2%Vold)',maxval(fs2%rho*fs2%V-fs2%rho*fs2%Vold),minval(fs2%rho*fs2%V-fs2%rho*fs2%Vold),'dt',time2%dt,'max,minval(resV2)',maxval(resV2),minval(resV2)
                  ! print*,'one more!'
                  resU2=-2.0_WP*(fs2%rho*fs2%U-fs2%rho*fs2%Uold)+time2%dt*resU2
                  resV2=-2.0_WP*(fs2%rho*fs2%V-fs2%rho*fs2%Vold)+time2%dt*resV2
                  resW2=-2.0_WP*(fs2%rho*fs2%W-fs2%rho*fs2%Wold)+time2%dt*resW2

                  forcing: block
                     use mpi_f08,  only: MPI_SUM,MPI_ALLREDUCE
                     use parallel, only: MPI_REAL_WP
                     integer :: i,j,k,ierr
                     real(WP) :: myU,myUvol,myW,myWvol,Uvol,Wvol
                     myU=0.0_WP; myUvol=0.0_WP; myW=0.0_WP; myWvol=0.0_WP
                     do k=fs2%cfg%kmin_,fs2%cfg%kmax_
                        do j=fs2%cfg%jmin_,fs2%cfg%jmax_
                           do i=fs2%cfg%imin_,fs2%cfg%imax_
                              if (fs2%umask(i,j,k).eq.0) then
                                 myU   =myU   +fs2%cfg%dxm(i)*fs2%cfg%dy(j)*fs2%cfg%dz(k)*(2.0_WP*fs2%U(i,j,k)-fs2%Uold(i,j,k))
                                 myUvol=myUvol+fs2%cfg%dxm(i)*fs2%cfg%dy(j)*fs2%cfg%dz(k)
                              end if
                           end do
                        end do
                     end do
                     call MPI_ALLREDUCE(myUvol,Uvol ,1,MPI_REAL_WP,MPI_SUM,fs2%cfg%comm,ierr)
                     call MPI_ALLREDUCE(myU   ,meanU,1,MPI_REAL_WP,MPI_SUM,fs2%cfg%comm,ierr); meanU=meanU/Uvol
                     ! print*,'meanU',meanU
                     where (fs2%umask.eq.0) resU2=resU2+0.5_WP*(u_jet-meanU)
                  end block forcing   

                  ! print*,'here3',time2%it
                  ! Form implicit residuals
                  call fs2%solve_implicit(time2%dt,resU2,resV2,resW2)
                  ! print*,'here4',time2%it
                  ! Apply these residuals
                  fs2%U=2.0_WP*fs2%U-fs2%Uold+resU2
                  fs2%V=2.0_WP*fs2%V-fs2%Vold+resV2
                  fs2%W=2.0_WP*fs2%W-fs2%Wold+resW2
                  ! print*,'here5',time2%it
                  ! Apply other boundary conditions on the resulting fields
                  call fs2%apply_bcond(time2%t,time2%dt)
                  ! print*,'here6',time2%it
                  ! ! Apply dirichlet bconds
                  ! dirichlet2 : block
                  !    use incomp_class, only: bcond
                  !    type(bcond),pointer::mybc
                  !    integer :: i,j,k,n
                  !    call fs2%get_bcond('gas_inj',mybc)
                  !    do n=1,mybc%itr%no_
                  !       i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  !       fs2%U(i,j,k) = u_jet
                  !       fs2%V(i,j,k) = 0.0_WP
                  !       fs2%W(i,j,k) = 0.0_WP
                  !    end do
                  ! end block dirichlet2
                  ! print*,'here7',time2%it
                  ! Solve Poisson equation
                  call fs2%correct_mfr()
                  call fs2%get_div()
                  fs2%psolv%rhs=-fs2%cfg%vol*fs2%div*fs2%rho/time2%dt
                  fs2%psolv%sol=0.0_WP
                  call fs2%psolv%solve()
                  ! print*,'here8',time2%it
                  ! Correct velocity
                  call fs2%get_pgrad(fs2%psolv%sol,resU2,resV2,resW2)
                  fs2%P=fs2%P+fs2%psolv%sol
                  fs2%U=fs2%U-time2%dt*resU2/fs2%rho
                  fs2%V=fs2%V-time2%dt*resV2/fs2%rho
                  fs2%W=fs2%W-time2%dt*resW2/fs2%rho
                  ! print*,'here9',time2%it
                  ! Increment sub-iteration counter
                  time2%it=time2%it+1
                  
               end do
               ! Recompute interpolated velocity and divergence
               call fs2%interp_vel(Ui2,Vi2,Wi2)
               call fs2%get_div()
               call fs2%get_max()
               call fs2%get_cfl(time2%dt,time2%cfl)
               call jetFile%write()
               ! Output to ensight
               if (ens_evt2%occurs()) call ens_out2%write_data(time2%t)
               
               ! ! Perform and output monitoring
               ! call fs2%get_max()
               ! call mfile2%write()
               ! call cflfile2%write()
               
            end do
         end if


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!! Pass velocitys between meshes !!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         pass_velocities : block
            use lowmach_class, only: bcond
            type(bcond),pointer :: mybc
            integer :: i,j,k,n

            ! Pass U
            if (amGrp2) resU2 = fs2%U
            if (amGrp2) call cpl%push(resU2)
            if (amGrp1.or.amGrp2) call cpl%transfer()
            if (amGrp1) then
               call cpl%pull(resU1)
               call fs1%get_bcond('gas_inj',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  passUVW(1,n)=resU1(i,j,k)
               end do
            end if

            ! Pass V
            if (amGrp2) resU2 = fs2%V
            if (amGrp2) call cpl%push(resU2)
            if (amGrp1.or.amGrp2) call cpl%transfer()
            if (amGrp1) then
               call cpl%pull(resU1)
               call fs1%get_bcond('gas_inj',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  passUVW(2,n)=resU1(i,j,k)
               end do
            end if

            ! Pass W
            if (amGrp2) resU2 = fs2%W
            if (amGrp2) call cpl%push(resU2)
            if (amGrp1.or.amGrp2) call cpl%transfer()
            if (amGrp1) then
               call cpl%pull(resU1)
               call fs1%get_bcond('gas_inj',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  passUVW(3,n)=resU1(i,j,k)
               end do
            end if

         end block pass_velocities

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SPRAY SOLVER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (amGrp1) then

            ! Update scalars & flow field:
            ! Remember old scalar
            T_sc%rhoold =T_sc%rho
            T_sc%SCold  =T_sc%SC
            Yf_sc%rhoold=Yf_sc%rho
            Yf_sc%SCold =Yf_sc%SC

            
            ! Remember old velocity and momentum
            fs1%rhoold=fs1%rho
            fs1%Uold=fs1%U; fs1%rhoUold=fs1%rhoU
            fs1%Vold=fs1%V; fs1%rhoVold=fs1%rhoV
            fs1%Wold=fs1%W; fs1%rhoWold=fs1%rhoW
            
            ! Apply time-varying Dirichlet conditions
            ! This is where time-dpt Dirichlet would be enforced
            ! print*,'here0'
            ! Perform sub-iterations
            do while (time1%it.le.time1%itmax)
               
               ! ============= SCALAR SOLVER =======================
               ! Build mid-time scalar
               T_sc%SC=0.5_WP*(T_sc%SC+T_sc%SCold)
               Yf_sc%SC=0.5_WP*(Yf_sc%SC+Yf_sc%SCold)
               
               ! Explicit calculation of drhoSC/dt from scalar equation
               call T_sc%get_drhoSCdt(resT,fs1%rhoU,fs1%rhoV,fs1%rhoW)
               call Yf_sc%get_drhoSCdt(resYf,fs1%rhoU,fs1%rhoV,fs1%rhoW)
               ! print*,'here0A',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))

               ! Assemble explicit residual            
               resT=time1%dt*resT-(2.0_WP*T_sc%rho*T_sc%SC-(T_sc%rho+T_sc%rhoold)*T_sc%SCold)
               resYf=time1%dt*resYf-(2.0_WP*Yf_sc%rho*Yf_sc%SC-(Yf_sc%rho+Yf_sc%rhoold)*Yf_sc%SCold)
               ! print_res : block
               !    use messager, only : die
               !    ! print*,resT
               !    ! print*,''
               !    ! print*,resYf 
               !    call die('end')  
               ! end block print_res
               ! Add mass, energy source terms
               add_mass_energy_src: block
                  use fluidTable_class, only: Cp_ID
                  integer :: i,j,k
                  real(WP) :: Cp_g
                  do k=fs1%cfg%kmin_,fs1%cfg%kmax_
                     do j=fs1%cfg%jmin_,fs1%cfg%jmax_
                        do i=fs1%cfg%imin_,fs1%cfg%imax_
                           ! Need in SC*kg/m^3
                           call gTab%evalProps(propOut=Cp_g,T_q=T_sc%SC(i,j,k),P_q=p0,propID=Cp_ID)
                           resT(i,j,k)=resT(i,j,k)+lp%srcE(i,j,k)/Cp_g!/cfg1%vol(i,j,k)
                           resYf(i,j,k)=resYf(i,j,k)+lp%srcM(i,j,k)!/cfg1%vol(i,j,k) ! *Yf_sc%rho(i,j,k)
                        end do
                     end do
                  end do
               end block add_mass_energy_src
               ! print*,'here0B'
               ! Form implicit residual
               call T_sc%solve_implicit(time1%dt,resT,fs1%rhoU,fs1%rhoV,fs1%rhoW)
               call Yf_sc%solve_implicit(time1%dt,resYf,fs1%rhoU,fs1%rhoV,fs1%rhoW)
               ! print*,'here0C',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))
               ! Apply this residual
               T_sc%SC=2.0_WP*T_sc%SC-T_sc%SCold+resT
               Yf_sc%SC=2.0_WP*Yf_sc%SC-Yf_sc%SCold+resYf

               ! clip_Yf: block
               !    integer :: i,j,k
               !    do k=fs1%cfg%kmin_,fs1%cfg%kmax_
               !       do j=fs1%cfg%jmin_,fs1%cfg%jmax_
               !          do i=fs1%cfg%imin_,fs1%cfg%imax_
               !             if (Yf_sc%SC(i,j,k).lt.0.0_WP) Yf_sc%SC(i,j,k) = 0.0_WP
               !          end do
               !       end do
               !    end do
               ! end block clip_Yf

               ! print*,'max rho2',maxval(Yf_sc%rho)
               
               ! Apply other boundary conditions on the resulting field
               scalar_dirichlet : block
                  use vdscalar_class, only: bcond
                  type(bcond),pointer::mybc
                  integer :: i,j,k,n
                  real(WP) :: W_g,W_l,R_cst
                  R_cst = 8.314472_WP ! Gas constant [J/(kg*K)]
                  W_g = gTab%MW ! Carrier gas molar weight [kg/mol]
                  W_l = lTab%MW ! Drop molar weight [kg/mol]
                  
                  call T_sc%get_bcond('T_jet',mybc)
                  do n=1,mybc%itr%no_
                     i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                     T_sc%SC(i,j,k) =T_jet
                     Yf_sc%SC(i,j,k)=Yf_jet
                     T_sc%rho(i,j,k)=p0/(R_cst*(Yf_sc%SC(i,j,k)/W_l+(1.0_WP-Yf_sc%SC(i,j,k))/W_g)*T_sc%SC(i,j,k))
                     Yf_sc%rho(i,j,k)=T_sc%rho(i,j,k)
                     T_sc%rhoSC(i,j,k) =T_jet*T_sc%rho(i,j,k)
                     Yf_sc%rhoSC(i,j,k)=Yf_jet*Yf_sc%rho(i,j,k)
                  end do  
                  
                  call T_sc%get_bcond('T_coflow',mybc)
                  do n=1,mybc%itr%no_
                     i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                     T_sc%SC(i,j,k) =T_coflow
                     Yf_sc%SC(i,j,k)=Yf_coflow
                     T_sc%rho(i,j,k)=p0/(R_cst*(Yf_sc%SC(i,j,k)/W_l+(1.0_WP-Yf_sc%SC(i,j,k))/W_g)*T_sc%SC(i,j,k))
                     Yf_sc%rho(i,j,k)=T_sc%rho(i,j,k)
                     T_sc%rhoSC(i,j,k) =T_jet*T_sc%rho(i,j,k)
                     Yf_sc%rhoSC(i,j,k)=Yf_jet*Yf_sc%rho(i,j,k)
                  end do  

               end block scalar_dirichlet

               call T_sc%apply_bcond(time1%t,time1%dt)
               call Yf_sc%apply_bcond(time1%t,time1%dt)
               ! print*,'here0D',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))
               ! print*,'max rho3',maxval(Yf_sc%rho)
               ! ===================================================
               
               ! ============ UPDATE PROPERTIES ====================
               ! Backup rhoSC
               resT=T_sc%rho*T_sc%SC
               resYf=Yf_sc%rho*Yf_sc%SC

               ! Update density
               call get_rho()
               T_sc%rho = T_sc%rho+lp%srcM
               Yf_sc%rho = Yf_sc%rho+lp%srcM
               ! print*,'here0E',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))
               ! Rescale scalar for conservation
               ! T_sc%SC=resT/T_sc%rho
               ! Yf_sc%SC=resYf/Yf_sc%rho

               ! UPDATE THE VISCOSITY & APPLY SUBGRID MODELING
               update_visc_sgs: block
                  use fluidTable_class, only: mu_ID
                  integer :: i,j,k
                  real(WP) :: Sc_t ! Turbulent schmidt number
                  ! Get the turbulent viscosity
                  call fs1%interp_vel(Ui1,Vi1,Wi1)
                  call fs1%get_strainrate(Ui=Ui1,Vi=Vi1,Wi=Wi1,SR=SR1)
                  call sgs1%get_visc(dt=time1%dtold,rho=T_sc%rho,Ui=Ui1,Vi=Vi1,Wi=Wi1,SR=SR1)

                  Sc_t = 1.2 ! from Li, et al. "Numerical and experimental investigation of turbulent n-heptane jet-in-hot-coflow flames." Fuel 283 (2021)

                  ! Loop some things
                  do k=fs1%cfg%kmin_,fs1%cfg%kmax_
                     do j=fs1%cfg%jmin_,fs1%cfg%jmax_
                        do i=fs1%cfg%imin_,fs1%cfg%imax_
                           ! Get the normal dynamic viscosity
                           call gTab%evalProps(propOut=fs1%visc(i,j,k),T_q=T_sc%SC(i,j,k),P_q=p0,propID=mu_ID)
                           ! Apply the turbulent viscosity
                           fs1%visc(i,j,k) = fs1%visc(i,j,k) + sgs1%visc(i,j,k)
               ! UPDATE THE DIFFUSIVITY
                           T_sc%diff(i,j,k) = diff_T +sgs1%visc(i,j,k)/ T_sc%rho(i,j,k)
                           Yf_sc%diff(i,j,k)= diff_Yf+sgs1%visc(i,j,k)/(Yf_sc%rho(i,j,k)*Sc_t)
                        end do
                     end do
                  end do
               end block update_visc_sgs
               ! print*,'max turbulent visc:',maxval(sgs1%visc)

               ! ===================================================
               ! ============ VELOCITY SOLVER ======================
               ! print*,'max rho4',maxval(Yf_sc%rho)
               ! Build n+1 density
               fs1%rho=0.5_WP*(Yf_sc%rho+Yf_sc%rhoold)
               
               ! Build mid-time velocity and momentum
               fs1%U=0.5_WP*(fs1%U+fs1%Uold); fs1%rhoU=0.5_WP*(fs1%rhoU+fs1%rhoUold)
               fs1%V=0.5_WP*(fs1%V+fs1%Vold); fs1%rhoV=0.5_WP*(fs1%rhoV+fs1%rhoVold)
               fs1%W=0.5_WP*(fs1%W+fs1%Wold); fs1%rhoW=0.5_WP*(fs1%rhoW+fs1%rhoWold)
               
               ! Explicit calculation of drho*u/dt from NS
               call fs1%get_dmomdt(resU1,resV1,resW1)
               ! print*,'here0F',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))
               
               ! Assemble explicit residual
               resU1=time1%dtmid*resU1-(2.0_WP*fs1%rhoU-2.0_WP*fs1%rhoUold)
               resV1=time1%dtmid*resV1-(2.0_WP*fs1%rhoV-2.0_WP*fs1%rhoVold)
               resW1=time1%dtmid*resW1-(2.0_WP*fs1%rhoW-2.0_WP*fs1%rhoWold)

               add_lpt_src: block
                  integer :: i,j,k
                  do k=fs1%cfg%kmin_,fs1%cfg%kmax_
                     do j=fs1%cfg%jmin_,fs1%cfg%jmax_
                        do i=fs1%cfg%imin_,fs1%cfg%imax_
                           resU1(i,j,k)=resU1(i,j,k)+sum(fs1%itpr_x(:,i,j,k)*lp%srcU(i-1:i,j,k)) 
                           resV1(i,j,k)=resV1(i,j,k)+sum(fs1%itpr_y(:,i,j,k)*lp%srcV(i,j-1:j,k))
                           resW1(i,j,k)=resW1(i,j,k)+sum(fs1%itpr_z(:,i,j,k)*lp%srcW(i,j,k-1:k))
                        end do
                     end do
                  end do
               end block add_lpt_src
               ! print*,'here0G',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))
               ! Form implicit residuals
               call fs1%solve_implicit(time1%dtmid,resU1,resV1,resW1)
               ! print*,'here0H',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W)),'rhoU',maxval(fs1%rhoU),'resU1',maxval(resU1)
               ! print*,'rhoU',maxval(fs1%rhoU),'resU1',maxval(resU1)
               ! Apply these residuals
               fs1%U=2.0_WP*fs1%U-fs1%Uold+resU1
               fs1%V=2.0_WP*fs1%V-fs1%Vold+resV1
               fs1%W=2.0_WP*fs1%W-fs1%Wold+resW1
               
               ! Apply other boundary conditions and update momentum
               call fs1%apply_bcond(time1%tmid,time1%dtmid)
               call fs1%rho_multiply()
               call fs1%apply_bcond(time1%tmid,time1%dtmid)
               ! print*,'here0I',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W)),'rhoU',maxval(fs1%rhoU),'resU1',maxval(resU1)
               ! This is where dirichlet BCs would go
               apply_dirichlet : block
                  use lowmach_class, only: bcond
                  type(bcond),pointer::mybc
                  integer :: i,j,k,n
                  call fs1%get_bcond('gas_inj',mybc)
                  do n=1,mybc%itr%no_
                     i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                     fs1%rhoU(i,j,k)=passUVW(1,n)*fs1%rho(i,j,k)
                     fs1%U(i,j,k)   =passUVW(1,n)
                     fs1%rhoV(i,j,k)=passUVW(2,n)*fs1%rho(i,j,k)
                     fs1%V(i,j,k)   =passUVW(2,n)
                     fs1%rhoW(i,j,k)=passUVW(3,n)*fs1%rho(i,j,k)
                     fs1%W(i,j,k)   =passUVW(3,n)
                  end do               
                  call fs1%get_bcond('gas_coflow',mybc)
                  do n=1,mybc%itr%no_
                     i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                     fs1%rhoU(i,j,k)=u_coflow*fs1%rho(i,j,k)
                     fs1%U(i,j,k)   =u_coflow
                     fs1%V(i,j,k)   =0.0_WP; fs1%W(i,j,k)=0.0_WP
                     fs1%rhoV(i,j,k)=0.0_WP; fs1%rhoW(i,j,k)=0.0_WP
                  end do  

               end block apply_dirichlet

               ! Solve Poisson equation
               call Yf_sc%get_drhodt(dt=time1%dt,drhodt=resYf)
               call fs1%correct_mfr(drhodt=(resYf + fs1%cfg%vol*lp%srcM/time1%dtmid))
               call fs1%get_div(drhodt=(resYf + fs1%cfg%vol*lp%srcM/time1%dtmid))
               fs1%psolv%rhs=-fs1%cfg%vol*fs1%div/time1%dtmid
               fs1%psolv%sol=0.0_WP
               call fs1%psolv%solve()
               call fs1%shift_p(fs1%psolv%sol)
               ! print*,'here0J',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W)),'rhoU',maxval(fs1%rhoU),'resU1',maxval(resU1)
               ! Correct momentum and rebuild velocity
               call fs1%get_pgrad(fs1%psolv%sol,resU1,resV1,resW1)
               ! print*,'here0K',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W)),'rhoU',maxval(fs1%rhoU),'resU1',maxval(resU1)
               fs1%P=fs1%P+fs1%psolv%sol
               ! print*,'here0L',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W)),'rhoU',maxval(fs1%rhoU),'resU1',maxval(resU1)
               ! print*,'rhoU',maxval(fs1%rhoU),'resU1',maxval(resU1),'dt*res',maxval(-time1%dtmid*resU1)
               fs1%rhoU=fs1%rhoU-time1%dtmid*resU1
               fs1%rhoV=fs1%rhoV-time1%dtmid*resV1
               fs1%rhoW=fs1%rhoW-time1%dtmid*resW1
               
               ! print*,'here0M',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))
               call fs1%rho_divide
               ! ===================================================
               ! print*,'here0N',time1%it,max(maxval(fs1%U),maxval(fs1%V),maxval(fs1%W))
               ! Increment sub-iteration counter
               time1%it=time1%it+1
            end do

            ! Spawn particles
            spawn_particles: block
               use random, only: random_normal,random_uniform
               use mathtools, only: twoPi
               integer :: i,np,nAdd
               real(WP) :: theta
               if (lp%cfg%amRoot) then
                  nAdd = int(floor(inj_rate*time1%t))-nInj
                  do i=1,nAdd
                     np=lp%np_+1; call lp%resize(np)
                     ! Give id
                     lp%p(i)%id=nInj+i
                     ! Set the diameter
                     lp%p(np)%d=random_normal(dp,dp_sig)
                     ! Set the temperature 
                     lp%p(np)%T_d = T_d
                     ! Position with uniform distribution
                     if (cfg1%ny.eq.1) then
                        lp%p(np)%pos = [0.0_WP,0.0_WP,random_uniform(-rp,rp)]
                     elseif (cfg1%nz.eq.1) then
                        lp%p(np)%pos = [0.0_WP,random_uniform(-rp,rp),0.0_WP]
                     else
                        theta = random_uniform(0.0_WP,twoPi)
                        lp%p(np)%pos = sqrt(random_uniform(0.0_WP,rp_sig))*[0.0_WP,cos(theta),sin(theta)]
                     end if
                     ! Position with normal distribution
                        ! if (cfg1%ny.eq.1) then
                        !    theta = 2.0_WP*floor(random_uniform(0.0_WP,2.0_WP))-1.0_WP ! Use theta as either -1 or 1
                        !    lp%p(np)%pos = random_normal(rp,rp_sig)*[0.0_WP,0.0_WP,theta]
                        ! elseif (cfg1%nz.eq.1) then
                        !    theta = 2.0_WP*floor(random_uniform(0.0_WP,2.0_WP))-1.0_WP ! Use theta as either -1 or 1
                        !    lp%p(np)%pos = random_normal(rp,rp_sig)*[0.0_WP,theta,0.0_WP]
                        ! else
                        !    theta = random_uniform(0.0_WP,twoPi)
                        !    lp%p(np)%pos = random_normal(rp,rp_sig)*[0.0_WP,cos(theta),sin(theta)]
                        ! end if

                     ! Give zero velocity
                     lp%p(np)%vel=[random_normal(u_inj,u_inj_sig),0.0_WP,0.0_WP]
                     ! Give zero collision force
                     lp%p(np)%col=0.0_WP
                     ! Give zero dt
                     lp%p(np)%dt=0.0_WP
                     ! Locate the particle on the mesh
                     lp%p(np)%ind=lp%cfg%get_ijk_global(lp%p(np)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
                     ! Activate the particle
                     lp%p(np)%flag=0
                     lp%np_=np
                  end do
                  nInj=nInj+nAdd
               end if
               ! Distribute particles
               call lp%sync()
                        ! Get initial particle volume fraction
               call lp%update_VF()

            end block spawn_particles

            ! Advance particles by dt
            ! print*,'U=',maxval(fs1%U),'V=',maxval(fs1%V),'W=',maxval(fs1%W),'rho=',maxval(fs1%rho),'visc=',maxval(fs1%visc),'T=',maxval(T_sc%SC),'Yf=',maxval(Yf_sc%SC)
            ! print*,'rank',cfg1%rank,'entering lp%advance with Yf_max,min:',maxval(Yf_sc%SC),minval(Yf_sc%SC)
            call lp%advance(dt=time1%dt,U=fs1%U,V=fs1%V,W=fs1%W,rho=fs1%rho,visc=fs1%visc,T=T_sc%SC,Yf=Yf_sc%SC,lTab=lTab,gTab=gTab,p_therm=p0)
            ! print*,'rank',cfg1%rank,'max(srcM)',maxval(lp%srcM),'min(srcM)',minval(lp%srcM),'max(srcE)',maxval(lp%srcE),'min(srcE)',minval(lp%srcE)
            ! print*,'here2'
            if ((lp%np.eq.0).and.(time1%t.gt.3.0_WP)) amDone=.true.
            ! Output to ensight
            if (ens_evt1%occurs()) then
               ensight_output1: block
                  integer :: i
                  call lp%update_partmesh(pmesh)
                  do i = 1,lp%np_
                     pmesh%var(1,i)=0.5_WP*lp%p(i)%d ! Droplet Radius
                     pmesh%var(2,i)=lp%p(i)%T_d ! Droplet Temperature
                  end do
         
                  call ens_out1%write_data(time1%t)
                  end block ensight_output1
            end if

            ! print*,'here3'
            ! Perform and output monitoring
            call lp%get_max()
            call fs1%get_max()
            call Yf_sc%get_max()
            call T_sc%get_max()
            call cflFile%write()
            call mfile%write()
            call pfile%write()
            ! print*,'here4'
         end if
      end do
   end subroutine simulation_run
   
   
   !> Finalize the NGA2 simulation
   subroutine simulation_final
      implicit none
      
      ! Get rid of all objects - need destructors
      ! monitor
      ! ensight
      ! bcond
      ! timetracker
      
      ! Deallocate work arrays
      ! deallocate(rho,visc,U,V,W)
   end subroutine simulation_final
   
end module simulation
