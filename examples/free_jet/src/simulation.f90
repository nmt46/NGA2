!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use fluidTable_class,  only: fluidTable
   use sgsmodel_class,    only: sgsmodel
   implicit none
   private
   
   !> Get a LPT solver, lowmach solver, scalar solver, LES model, and corresponding time tracker
   type(lowmach),     public :: fs
   type(vdscalar),    public :: T_sc, Yf_sc
   type(sgsmodel),    public :: sgs
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,cflFile
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Fluid phase arrays
   ! real(WP), dimension(:,:,:), allocatable :: U,V,W,T,Yf
   ! real(WP), dimension(:,:,:), allocatable :: rho,visc
   real(WP) :: p0 ! ambient pressure
   real(WP) :: diff_T,diff_Yf ! non-turbulent diffusivities for scalar solvers

   !> Case quantities
   real(WP) :: r_jet,v_jet,T_jet,Yf_jet,v_coflow,T_coflow,Yf_coflow
   integer :: nInj

   !> Quantities for flow & scalar solvers:
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW,resT,resYf
   real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   real(WP), dimension(:,:,:,:), allocatable :: SR

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

   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imax+1) isIn = .true.
   end function xp_locator

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
   
   function gas_inj_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imin.and.(sqrt(pg%ym(j)**2+pg%zm(k)**2)).le.r_jet) isIn = .true.
   end function gas_inj_locator

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
                  T_sc%rho(i,j,k)=p0/(R_cst*(Yf_sc%SC(i,j,k)/W_l+(1.0_WP-Yf_sc%SC(i,j,k))/W_g)*T_sc%SC(i,j,k))
               else
                  T_sc%rho(i,j,k)=p0*W_g/(R_cst*T_sc%SC(i,j,k))
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
         call param_read('Gas jet velocity',v_jet)
         call param_read('Gas coflow velocity',v_coflow)
         call param_read('Gas jet radius',r_jet)
         call param_read('Pressure',p0)
         call param_read('Gas jet temperature',T_jet)
         call param_read('Gas coflow temperature',T_coflow)
         call param_read('Gas jet fuel mass fraction',Yf_jet)
         call param_read('Gas coflow fuel mass fraction',Yf_coflow)

      end block read_global_vars
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=min(0.1_WP*cfg%dx(1)/v_jet,time%dtmax)
         ! time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      ! print*,'init0'
      initialize_fluid_properties: block
         ! use messager, only : die
         use string,    only: str_medium
         use fluidTable_class, only: Cp_ID,Lv_ID,Tb_ID,rho_ID,MW_ID,mu_ID
         character(len=str_medium) :: name
         integer :: nP,nT
         real(WP) :: Tmin,Tmax,Pmin,Pmax,dP,dT

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

         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         call param_read('Pressure',p0)
         call param_read('Gas',name)
         nP = 9; nT = 9 ! points on interpolation mesh

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
      ! print*,'init1'
      initialize_res: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resU=0.0_WP
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resV=0.0_WP
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resW=0.0_WP
         allocate(resT (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resT=0.0_WP
         allocate(resYf(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resYf=0.0_WP
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Ui=0.0_WP
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Vi=0.0_WP
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); Wi=0.0_WP
         allocate(SR (6,cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); SR=0.0_WP
      end block initialize_res
      ! print*,'init2'

      ! Initialize scalar solvers
      initialize_sc: block
         use ils_class,      only: gmres_amg,gmres
         use vdscalar_class, only: neumann,dirichlet,bcond
         type(bcond),pointer::mybc
         integer :: i,j,k,n
         ! create the scalar solvers for temperature and fuel mass fraction
         T_sc=vdscalar(cfg=cfg,scheme=1,name='Temperature solver')
         Yf_sc=vdscalar(cfg=cfg,scheme=1,name='Fuel mass fraction solver')
         ! print*,'init2a'
         ! Configure implicit scalar solver
         call param_read('Implicit iteration',T_sc%implicit%maxit)
         call param_read('Implicit tolerance',T_sc%implicit%rcvg)
         ! print*,'init2b'
         Yf_sc%implicit%maxit=T_sc%implicit%maxit; YF_sc%implicit%rcvg=T_sc%implicit%rcvg
         ! print*,'init2b1'
         ! Create bconds
         call T_sc%add_bcond(name='T_out',   type=neumann,  dir='xp',locator=xp_sc_locator)
         ! print*,'init2b2'
         call T_sc%add_bcond(name='T_ym',    type=neumann,  dir='ym',locator=ym_sc_locator)
         ! print*,'init2b3'
         call T_sc%add_bcond(name='T_yp',    type=neumann,  dir='yp',locator=yp_sc_locator)
         ! print*,'init2b4'
         call T_sc%add_bcond(name='T_zm',    type=neumann,  dir='zm',locator=zm_sc_locator)
         ! print*,'init2b5'
         call T_sc%add_bcond(name='T_zp',    type=neumann,  dir='zp',locator=zp_sc_locator)
         ! print*,'init2b6'
         call T_sc%add_bcond(name='T_jet',   type=dirichlet,dir='xm',locator=gas_inj_locator)
         ! print*,'init2b7'
         call T_sc%add_bcond(name='T_coflow',type=dirichlet,dir='xm',locator=coflow_locator)
         ! print*,'init2b8'
         
         call Yf_sc%add_bcond(name='Yf_out',   type=neumann, dir='xp',locator=xp_sc_locator)
         call Yf_sc%add_bcond(name='Yf_ym',    type=neumann,  dir='ym',locator=ym_sc_locator)
         call Yf_sc%add_bcond(name='Yf_yp',    type=neumann,  dir='yp',locator=yp_sc_locator)
         call Yf_sc%add_bcond(name='Yf_zm',    type=neumann,  dir='zm',locator=zm_sc_locator)
         call Yf_sc%add_bcond(name='Yf_zp',    type=neumann,  dir='zp',locator=zp_sc_locator)
         call Yf_sc%add_bcond(name='Yf_jet',   type=dirichlet,dir='xm',locator=gas_inj_locator)
         call Yf_sc%add_bcond(name='Yf_coflow',type=dirichlet,dir='xm',locator=coflow_locator)
         ! print*,'init2c'
         ! Setup the solver
         call T_sc%setup(implicit_ils=gmres_amg)
         call Yf_sc%setup(implicit_ils=gmres_amg)

         ! Set initial field values
         T_sc%SC=T_coflow
         Yf_sc%SC=Yf_coflow
         ! print*,'init2d'
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
         ! print*,'init2e'
         call get_rho()
         ! print*,'init2f'
         ! Set diffusivities
         call param_read('Thermal diffusivity',diff_T)
         call param_read('Fuel diffusivity',diff_Yf)
         T_sc%diff = diff_T
         Yf_sc%diff = diff_Yf
      end block initialize_sc
      ! print*,'init3'

      ! Initialize the flow solver
      initialize_fs: block
         use mathtools, only: twoPi
         use lowmach_class, only: clipped_neumann,dirichlet,bcond
         use ils_class, only: pcg_amg,pfmg,pcg_pfmg
         use fluidTable_class, only: mu_ID
         type(bcond),pointer :: mybc
         integer :: i,j,k,n
         real(WP) :: rhof,viscf,viscf2
         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Setup bconds
         ! call fs%add_bcond(name='bc_xm',type=neumann,face='x',dir=-1,canCorrect=.true.,locator=xm_locator)
         call fs%add_bcond(name='bc_xp',     type=clipped_neumann,face='x',dir=+1,canCorrect=.true., locator=xp_locator)
         call fs%add_bcond(name='gas_inj',   type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=gas_inj_locator)
         call fs%add_bcond(name='gas_coflow',type=dirichlet,      face='x',dir=-1,canCorrect=.false.,locator=coflow_locator)

         ! Set initial density, viscosity
         fs%rho = Yf_sc%rho ! Take density from scalar solver
         call gTab%evalProps(propOut=viscf,T_q=T_jet,P_q=p0,propID=mu_ID); fs%visc=viscf
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_amg,implicit_ils=pcg_pfmg)!pfmg)

         ! Initialize velocity field
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  fs%U(i,j,k)=0.0_WP;fs%rhoU(i,j,k)=0.0_WP
                  fs%V(i,j,k)=0.0_WP;fs%rhoV(i,j,k)=0.0_WP
                  fs%W(i,j,k)=0.0_WP;fs%rhoW(i,j,k)=0.0_WP
               end do
            end do
         end do
         call fs%get_bcond('gas_inj',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%rhoU(i,j,k)=v_jet*fs%rho(i,j,k)
            fs%U(i,j,k)   =v_jet
         end do
         call fs%get_bcond('gas_coflow',mybc)
         do n=1,mybc%itr%no_
            i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
            fs%rhoU(i,j,k)=v_coflow*fs%rho(i,j,k)
            fs%U(i,j,k)   =v_coflow
         end do
         call fs%apply_bcond(time%t,time%dt)
         call T_sc%get_drhodt(drhodt=resT,dt=time%dt) ! use resT to hold our drhodt :)
         call fs%get_div(drhodt=resT)
         ! Compute MFR through all boundary conditions
         call fs%correct_mfr(drhodt=resT)
         call fs%rho_divide()
         call fs%interp_vel(Ui,Vi,Wi)
      end block initialize_fs
      ! print*,'init4'

      ! Create an LES model
      create_sgs: block
         sgs=sgsmodel(cfg=fs%cfg,umask=fs%umask,vmask=fs%vmask,wmask=fs%wmask)
      end block create_sgs
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=cfg,name='free_jet')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_vector('velocity',fs%U,fs%V,fs%W)
         call ens_out%add_scalar('yf',Yf_sc%SC)
         call ens_out%add_scalar('temperature',T_sc%SC)
         call ens_out%add_scalar('density',T_sc%rho)

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      ! print*,'init5'

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call fs%get_max()
         call Yf_sc%get_max()
         call T_sc%get_max()
         call fs%get_cfl(time%dt,time%cfl)
         ! Create simulation monitor
         mfile=monitor(amroot=cfg%amRoot,name='simulation')
         call mfile%add_column(time%n,'N_time')
         call mfile%add_column(time%wt,'wall_t')
         call mfile%add_column(time%t,'time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(time%cfl,'CFL')
         call mfile%add_column(fs%Umax,'FS_Max_U')
         call mfile%add_column(fs%Vmax,'FS_Max_V')
         call mfile%add_column(fs%Wmax,'FS_Max_W')
         call mfile%add_column(Yf_sc%SCmin,'Yf_min')
         call mfile%add_column(Yf_sc%SCmax,'Yf_max')
         call mfile%add_column(T_sc%SCmin,'T_min')
         call mfile%add_column(T_sc%SCmax,'T_max')
         call mfile%write()

         cflFile = monitor(amroot=cfg%amRoot,name='cfl')
         call cflFile%add_column(time%n,'N_time')
         call cflFile%add_column(time%t,'time')
         call cflFile%add_column(max(fs%CFLc_x,fs%CFLc_y,fs%CFLc_z),'max_CFLc')
         call cflFile%add_column(max(fs%CFLv_x,fs%CFLv_y,fs%CFLv_z),'max_CFLv')
         call cflFile%write()

      end block create_monitor
      ! print*,'init6'
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      logical :: amDone = .false.
      ! Perform time integration
      do while (.not.time%done().and.(.not.amDone))
         
         ! Increment time         
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()

         ! Update scalars & flow field:


         ! Remember old scalar
         T_sc%rhoold =T_sc%rho
         T_sc%SCold  =T_sc%SC
         Yf_sc%rhoold=Yf_sc%rho
         Yf_sc%SCold =Yf_sc%SC

         
         ! Remember old velocity and momentum
         fs%rhoold=fs%rho
         fs%Uold=fs%U; fs%rhoUold=fs%rhoU
         fs%Vold=fs%V; fs%rhoVold=fs%rhoV
         fs%Wold=fs%W; fs%rhoWold=fs%rhoW
         
         ! Apply time-varying Dirichlet conditions
         ! This is where time-dpt Dirichlet would be enforced
         ! print*,'here0'
         ! Perform sub-iterations
         do while (time%it.le.time%itmax)
            
            ! ============= SCALAR SOLVER =======================
            ! Build mid-time scalar
            T_sc%SC=0.5_WP*(T_sc%SC+T_sc%SCold)
            Yf_sc%SC=0.5_WP*(Yf_sc%SC+Yf_sc%SCold)
            
            ! Explicit calculation of drhoSC/dt from scalar equation
            call T_sc%get_drhoSCdt(resT,fs%rhoU,fs%rhoV,fs%rhoW)
            call Yf_sc%get_drhoSCdt(resYf,fs%rhoU,fs%rhoV,fs%rhoW)
            ! print*,'here0A',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))

            ! Assemble explicit residual            
            resT=time%dt*resT-(2.0_WP*T_sc%rho*T_sc%SC-(T_sc%rho+T_sc%rhoold)*T_sc%SCold)
            resYf=time%dt*resYf-(2.0_WP*Yf_sc%rho*Yf_sc%SC-(Yf_sc%rho+Yf_sc%rhoold)*Yf_sc%SCold)
            ! print_res : block
            !    use messager, only : die
            !    ! print*,resT
            !    ! print*,''
            !    ! print*,resYf 
            !    call die('end')  
            ! end block print_res
            ! print*,'here0B'
            ! Form implicit residual
            call T_sc%solve_implicit(time%dt,resT,fs%rhoU,fs%rhoV,fs%rhoW)
            call Yf_sc%solve_implicit(time%dt,resYf,fs%rhoU,fs%rhoV,fs%rhoW)
            ! print*,'here0C',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! Apply this residual
            T_sc%SC=2.0_WP*T_sc%SC-T_sc%SCold+resT
            Yf_sc%SC=2.0_WP*Yf_sc%SC-Yf_sc%SCold+resYf
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
            call T_sc%apply_bcond(time%t,time%dt)
            call Yf_sc%apply_bcond(time%t,time%dt)
            ! print*,'here0D',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! print*,'max rho3',maxval(Yf_sc%rho)
            ! ===================================================
            
            ! ============ UPDATE PROPERTIES ====================
            ! Backup rhoSC
            ! resT=T_sc%rho*T_sc%SC
            ! resYf=Yf_sc%rho*Yf_sc%SC

            ! Update density
            call get_rho()

            ! print*,'here0E',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! Rescale scalar for conservation
            ! T_sc%SC=resT/T_sc%rho
            ! Yf_sc%SC=resYf/Yf_sc%rho

            ! UPDATE THE VISCOSITY & APPLY SUBGRID MODELING
            update_visc_sgs: block
               use fluidTable_class, only: mu_ID
               integer :: i,j,k
               real(WP) :: Sc_t ! Turbulent schmidt number
               ! Get the turbulent viscosity
               call fs%interp_vel(Ui,Vi,Wi)
               call fs%get_strainrate(Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)
               call sgs%get_visc(dt=time%dtold,rho=T_sc%rho,Ui=Ui,Vi=Vi,Wi=Wi,SR=SR)

               Sc_t = 1.2 ! from Li, et al. "Numerical and experimental investigation of turbulent n-heptane jet-in-hot-coflow flames." Fuel 283 (2021)

               ! Loop some things
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        ! Get the normal dynamic viscosity
                        call gTab%evalProps(propOut=fs%visc(i,j,k),T_q=T_sc%SC(i,j,k),P_q=p0,propID=mu_ID)
                        ! Apply the turbulent viscosity
                        fs%visc(i,j,k) = fs%visc(i,j,k) + sgs%visc(i,j,k)
            ! UPDATE THE DIFFUSIVITY
                        T_sc%diff(i,j,k) = diff_T +sgs%visc(i,j,k)/ T_sc%rho(i,j,k)
                        Yf_sc%diff(i,j,k)= diff_Yf+sgs%visc(i,j,k)/(Yf_sc%rho(i,j,k)*Sc_t)
                     end do
                  end do
               end do
            end block update_visc_sgs


            ! ============ VELOCITY SOLVER ======================
            ! print*,'max rho4',maxval(Yf_sc%rho)
            ! Build n+1 density
            ! print*,'fsSolve0'
            fs%rho=0.5_WP*(Yf_sc%rho+Yf_sc%rhoold)
            ! print*,minval(fs%rho(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_))
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            ! print*,'fsSolve1'
            ! print*,resU(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)
            ! print*,resV(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)
            ! print*,resW(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)
            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

            ! print*,'maxW',maxval(resW(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)),&
                  !  'minW',minval(resW(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)),&
                  !  'max(abs(W))',maxval(abs(resW(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_))),&
                  !  'min(abs(W))',minval(abs(resW(fs%cfg%imin_:fs%cfg%imax_,fs%cfg%jmin_:fs%cfg%jmax_,fs%cfg%kmin_:fs%cfg%kmax_)))
            ! Form implicit residuals
            ! print*,'fsSolve2'
            call fs%solve_implicit(time%dtmid,resU,resV,resW)
            ! print*,'fsSolve3'
            ! print*,'here0H',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! print*,'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            ! print*,'fsSolve4'
            ! Apply other boundary conditions and update momentum
            call fs%apply_bcond(time%tmid,time%dtmid)
            call fs%rho_multiply()
            call fs%apply_bcond(time%tmid,time%dtmid)
            ! print*,'fsSolve5'
            ! print*,'here0I',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! This is where dirichlet BCs would go
            apply_dirichlet : block
               use lowmach_class, only: bcond
               type(bcond),pointer::mybc
               integer :: i,j,k,n
               call fs%get_bcond('gas_inj',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%rhoU(i,j,k)=v_jet*fs%rho(i,j,k)
                  fs%U(i,j,k)   =v_jet
                  fs%V(i,j,k)   =0.0_WP; fs%W(i,j,k)=0.0_WP
                  fs%rhoV(i,j,k)=0.0_WP; fs%rhoW(i,j,k)=0.0_WP
               end do               
               call fs%get_bcond('gas_coflow',mybc)
               do n=1,mybc%itr%no_
                  i=mybc%itr%map(1,n); j=mybc%itr%map(2,n); k=mybc%itr%map(3,n)
                  fs%rhoU(i,j,k)=v_coflow*fs%rho(i,j,k)
                  fs%U(i,j,k)   =v_coflow
                  fs%V(i,j,k)   =0.0_WP; fs%W(i,j,k)=0.0_WP
                  fs%rhoV(i,j,k)=0.0_WP; fs%rhoW(i,j,k)=0.0_WP
               end do  

            end block apply_dirichlet
            ! print*,'fsSolve6'
            ! Solve Poisson equation
            call Yf_sc%get_drhodt(dt=time%dt,drhodt=resYf)
            call fs%correct_mfr(drhodt=resYf)
            call fs%get_div(drhodt=resYf)
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            ! print*,'fsSolve7'
            call fs%shift_p(fs%psolv%sol)
            ! print*,'here0J',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            ! print*,'here0K',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            fs%P=fs%P+fs%psolv%sol
            ! print*,'here0L',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! print*,'rhoU',maxval(fs%rhoU),'resU',maxval(resU),'dt*res',maxval(-time%dtmid*resU)
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            ! print*,'fsSolve8'
            ! print*,'here0M',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            call fs%rho_divide
            ! ===================================================
            ! print*,'here0N',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! Increment sub-iteration counter
            time%it=time%it+1
            ! print*,'fsSolve9'
         end do

         ! Output to ensight
         if (ens_evt%occurs()) then
            ensight_output: block
               call ens_out%write_data(time%t)
               end block ensight_output
         end if
         ! print*,'here3'
         ! Perform and output monitoring
         call fs%get_max()
         call Yf_sc%get_max()
         call T_sc%get_max()
         call cflFile%write()
         call mfile%write()
         ! print*,'here4'
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
      if (cfg%rank.eq.0) print*,'Clean termination'
   end subroutine simulation_final
   
end module simulation
