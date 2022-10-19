!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lpt_class,         only: lpt
   use lowmach_class,     only: lowmach
   use vdscalar_class,    only: vdscalar
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use fluidTable_class,  only: fluidTable
   implicit none
   private
   
   !> Get a LPT solver, lowmach solver, and corresponding time tracker
   type(lpt),         public :: lp
   type(lowmach),     public :: fs
   type(vdscalar),    public :: T_sc, Yf_sc
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,pfile,cflFile

   !> Monitoring quantities
   real(WP) :: Bm_debug = 0.0_WP ! The Spalding number for debugging
   
   public :: simulation_init,simulation_run,simulation_final
   
   !> Fluid phase arrays
   real(WP), dimension(:,:,:), allocatable :: U,V,W,T,Yf
   ! real(WP), dimension(:,:,:), allocatable :: rho,visc
   real(WP) :: p0 ! ambient pressure

   !> Case quantities
   real(WP) :: dp,dp_sig,rp,rp_sig,u_inj,u_inj_sig,T_d

   !> Quantities for flow & scalar solvers:
   real(WP), dimension(:,:,:), allocatable :: resU,resV,resW,resT,resYf
   real(WP), dimension(:,:,:), allocatable :: Ui,Vi,Wi
   real(WP) :: maxTemp_dt ! maximum allowable temperature change per time step [K/timeStep]
   real(WP) :: T_dt_cur ! the current temperature-based max allowable timestep [sec]

   type(fluidTable) :: lTab,gTab


   
contains
   
   function xn_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imin) isIn = .true.
   end function xn_locator

   function xp_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (i.eq.pg%imax+1) isIn = .true.
   end function xp_locator

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
                  T_sc%rho(i,j,k)=p0/(R_cst*(Yf_sc%SC(i,j,k)/W_l+(1.0_WP-Yf_sc%SC(i,j,k))/W_g)*T_sc%SC(i,j,k)) + lp%srcM(i,j,k)
               else
                  T_sc%rho(i,j,k)=p0*W_g/(R_cst*T_sc%SC(i,j,k)) + lp%srcM(i,j,k)
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
      
      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      
      initialize_fluid_properties: block
         ! use messager, only : die
         ! use fluidTable_class, only: mu_ID
         integer :: i,j
         real(WP) :: Tmin,Tmax,Pmin,Pmax,dP,dT
         real(WP) :: testVal

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Liquid !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         lTab%nP = 21
         lTab%nT = 21
         allocate(lTab%Cp(lTab%nT)); lTab%Cp = 0.0_WP
         allocate(lTab%rho(lTab%nP,lTab%nT)); lTab%rho = 0.0_WP
         allocate(lTab%L_v(lTab%nP)); lTab%L_v = 0.0_WP
         allocate(lTab%T_b(lTab%nP)); lTab%T_b = 0.0_WP
         allocate(lTab%T(lTab%nT)); lTab%T = 0.0_WP
         allocate(lTab%P(lTab%nP)); lTab%P = 0.0_WP
         allocate(lTab%mu(lTab%nP,lTab%nT)); lTab%mu = 0.0_WP

         call param_read('Pressure',p0)
         call param_read('Fuel',lTab%name)
         call param_read('Gas inlet temperature',Tmax)

         Pmin = p0/1.1_WP; Pmax = p0*1.1_WP
         dP = (Pmax-Pmin)/(real(lTab%nP,WP)-1.0_WP)

         Tmin = maxval([250.0_WP,1.01_WP*cprop(output='Tmin'//char(0),name1='P'//char(0),prop1=p0,name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0))]); 
         Tmax = 0.99_WP*cprop(output='T'//char(0),name1='P'//char(0),prop1=p0,name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0))
         dT = (Tmax-Tmin)/(real(lTab%nT,WP)-1.0_WP)

         lTab%MW = cprop(output='M'//char(0),name1='T'//char(0),prop1=Tmin,name2='P'//char(0),prop2=p0,fluidname=trim(lTab%name)//char(0))

         do i=1,lTab%nP
            lTab%P(i) = Pmin+dP*(real(i,WP)-1.0_WP)
            do j=1,lTab%nT
               if (i.eq.1) then
                  lTab%T(j) = Tmin+dT*(real(j,WP)-1.0_WP)
                  lTab%Cp(j) = cprop(output='CPMASS'//char(0),name1='T'//char(0),prop1=lTab%T(j),name2='P'//char(0),prop2=p0,fluidname=trim(lTab%name)//char(0))
               end if
               lTab%rho(i,j) = cprop(output='D'//char(0),name1='T'//char(0),prop1=lTab%T(j),name2='P'//char(0),prop2=lTab%P(i),fluidname=trim(lTab%name)//char(0))
               lTab%mu(i,j) = cprop(output='V'//char(0),name1='T'//char(0),prop1=lTab%T(j),name2='P'//char(0),prop2=lTab%P(i),fluidname=trim(lTab%name)//char(0))

            end do
            lTab%T_b(i) = cprop(output='T'//char(0),name1='P'//char(0),prop1=lTab%P(i),name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0))
            lTab%L_v(i) = cprop(output='H'//char(0),name1='P'//char(0),prop1=lTab%P(i),name2='Q'//char(0),prop2=1.0_WP,fluidname=trim(lTab%name)//char(0)) - &
                        cprop(output='H'//char(0),name1='P'//char(0),prop1=lTab%P(i),name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0)) ! latent heat of vaporization
         end do

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         gTab%nP = 11
         gTab%nT = 11
         allocate(gTab%Cp(gTab%nT)); gTab%Cp = 0.0_WP
         allocate(gTab%rho(gTab%nP,gTab%nT)); gTab%rho = 0.0_WP
         allocate(gTab%L_v(gTab%nP)); gTab%L_v = 0.0_WP
         allocate(gTab%T_b(gTab%nP)); gTab%T_b = 0.0_WP
         allocate(gTab%T(gTab%nT)); gTab%T = 0.0_WP
         allocate(gTab%P(gTab%nP)); gTab%P = 0.0_WP
         allocate(gTab%mu(gTab%nP,gTab%nT)); gTab%mu = 0.0_WP

         call param_read('Pressure',p0)
         call param_read('Gas',gTab%name)
         call param_read('Gas inlet temperature',Tmax)
         call param_read('Gas Prandtl number',gTab%Pr)
         call param_read('Fuel-gas Schmidt number',gTab%Sc)
         Pmin = p0/1.1_WP; Pmax = p0*1.1_WP
         dP = (Pmax-Pmin)/(real(gTab%nP,WP)-1.0_WP)

         Tmin = maxval([250.0_WP,1.01_WP*cprop(output='Tmin'//char(0),name1='P'//char(0),prop1=p0,name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(lTab%name)//char(0))])
         Tmax = Tmax*1.5
         dT = (Tmax-Tmin)/(real(gTab%nT,WP)-1.0_WP)

         gTab%MW = cprop(output='M'//char(0),name1='T'//char(0),prop1=Tmin,name2='P'//char(0),prop2=p0,fluidname=trim(gTab%name)//char(0))

         do i=1,gTab%nP
            gTab%P(i) = Pmin+dP*(real(i,WP)-1.0_WP)
            do j=1,gTab%nT
               if (i.eq.1) then
                  gTab%T(j) = Tmin+dT*(real(j,WP)-1.0_WP)
                  gTab%Cp(j)= cprop(output='CPMASS'//char(0),name1='T'//char(0),prop1=gTab%T(j),name2='P'//char(0),prop2=p0,fluidname=trim(gTab%name)//char(0))
               end if
               gTab%rho(i,j)= cprop(output='D'//char(0),name1='T'//char(0),prop1=gTab%T(j),name2='P'//char(0),prop2=gTab%P(i),fluidname=trim(gTab%name)//char(0))
               gTab%mu(i,j) = cprop(output='V'//char(0),name1='T'//char(0),prop1=gTab%T(j),name2='P'//char(0),prop2=gTab%P(i),fluidname=trim(gTab%name)//char(0))

            end do
            gTab%T_b(i) = cprop(output='T'//char(0),name1='P'//char(0),prop1=gTab%P(i),name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(gTab%name)//char(0))
            gTab%L_v(i) = cprop(output='H'//char(0),name1='P'//char(0),prop1=gTab%P(i),name2='Q'//char(0),prop2=1.0_WP,fluidname=trim(gTab%name)//char(0)) - &
                          cprop(output='H'//char(0),name1='P'//char(0),prop1=gTab%P(i),name2='Q'//char(0),prop2=0.0_WP,fluidname=trim(gTab%name)//char(0)) ! latent heat of vaporization
         end do
         ! print*,'Tmin',Tmin,'Tmax',Tmax,'nT',lTab%nT,'Pmin',Pmin,'Pmax',Pmax,'nP',lTab%nP
         ! print*,'max(mu)',maxval(gTab%mu),'min(mu)',minval(gTab%mu)
         ! call gTab%evalProps(propOut=testVal,T_q=280.0_WP,P_q=p0,propID=mu_ID)
         ! print*,'interped mu',testVal
         ! call die('ending here')

      end block initialize_fluid_properties
   
      ! Initialize our LPT
      initialize_lpt: block
         use fluidTable_class, only: rho_ID
         ! use random, only: random_uniform
         real(WP),dimension(3) :: v_drop
         integer :: i
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')         
         ! Get droplet diameter from the input
         call param_read('Droplet mean diameter',dp)
         call param_read('Droplet diameter standard deviation',dp_sig)
         call param_read('Injection mean radius',rp)
         call param_read('Injection radius standard deviation',rp_sig)
         call param_read('Injection mean velocity',u_inj)
         call param_read('Injection velocity standard deviation',u_inj_sig)
         ! Get droplet initial temperature
         call param_read('Droplet temperature',T_d)
         ! Get droplet density from the fluid table
         call lTab%evalProps(propOut=lp%rho,T_q=T_d,P_q=p0,propID=rho_ID)
         call param_read('Droplet velocity',v_drop)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         ! Root process initializes 1000 particles randomly
         if (lp%cfg%amRoot) then
            call lp%resize(1)
            do i=1,1
               ! Give id
               lp%p(i)%id=int(i,8)
               ! Set the diameter
               lp%p(i)%d=dp
               ! Set the temperature 
               lp%p(i)%T_d = T_d
               !! Assign random position in the domain
               ! lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
               ! &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
               ! &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
               ! lp%p(i)%pos(3)=lp%cfg%zm(lp%cfg%kmin)

               ! Assign position in center of domain
               lp%p(i)%pos = 0.0_WP
               ! Give zero velocity
               lp%p(i)%vel=v_drop![10.0_WP,0.0_WP,0.0_WP]
               ! Give zero collision force
               lp%p(i)%col=0.0_WP
               ! Give zero dt
               lp%p(i)%dt=0.0_WP
               ! Locate the particle on the mesh
               lp%p(i)%ind=lp%cfg%get_ijk_global(lp%p(i)%pos,[lp%cfg%imin,lp%cfg%jmin,lp%cfg%kmin])
               ! Activate the particle
               lp%p(i)%flag=0
            end do
         end if
         ! Distribute particles
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
         ! Set collision timescale
         lp%Tcol=5.0_WP*time%dt
      end block initialize_lpt
      

      initialize_res: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resT (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(resYf(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Ui   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Vi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
         allocate(Wi   (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_))
      end block initialize_res
      
      ! Initialize scalar solvers
      initialize_sc: block
         use ils_class,      only: gmres
         real(WP) :: diff_T,diff_Yf,T_in
         ! create the scalar solvers for temperature and fuel mass fraction
         T_sc=vdscalar(cfg=cfg,scheme=1,name='Temperature solver')
         Yf_sc=vdscalar(cfg=cfg,scheme=1,name='Fuel mass fraction solver')
         ! Configure implicit scalar solver
         call param_read('Implicit iteration',T_sc%implicit%maxit)
         call param_read('Implicit tolerance',T_sc%implicit%rcvg)
         Yf_sc%implicit%maxit=T_sc%implicit%maxit; YF_sc%implicit%rcvg=T_sc%implicit%rcvg
         ! Setup the solver
         call T_sc%setup(implicit_ils=gmres)
         call Yf_sc%setup(implicit_ils=gmres)
         ! Set initial field values
         call param_read('Gas inlet temperature',T_in)
         T_sc%SC=T_in
         Yf_sc%SC=0
         call get_rho()
         ! Set diffusivities
         call param_read('Thermal diffusivity',diff_T)
         call param_read('Fuel diffusivity',diff_Yf)
         T_sc%diff = diff_T
         Yf_sc%diff = diff_Yf
         call param_read('Max temperature change',maxTemp_dt)
      end block initialize_sc

      ! Initialize the flow solver
      initialize_fs: block
         use mathtools, only: twoPi
         use lowmach_class, only: clipped_neumann
         use ils_class, only: pcg_amg,pfmg
         use fluidTable_class, only: mu_ID
         integer :: i,j,k
         real(WP) :: rhof,viscf,T_in,viscf2
         ! Allocate arrays
         ! allocate(rho (lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
         ! allocate(visc(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
         ! allocate(U(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
         ! allocate(V(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))
         ! allocate(W(lp%cfg%imino_:lp%cfg%imaxo_,lp%cfg%jmino_:lp%cfg%jmaxo_,lp%cfg%kmino_:lp%cfg%kmaxo_))

         fs=lowmach(cfg=cfg,name='Variable density low Mach NS')
         ! Setup bconds
         call fs%add_bcond(name='bc_xn',type=clipped_neumann,face='x',dir=-1,canCorrect=.true.,locator=xn_locator)
         call fs%add_bcond(name='bc_xp',type=clipped_neumann,face='x',dir=+1,canCorrect=.true.,locator=xp_locator)

         ! Set initial density, viscosity
         fs%rho = Yf_sc%rho ! Take density from scalar solver
         call param_read('Gas inlet temperature',T_in)
         call gTab%evalProps(propOut=viscf,T_q=T_in,P_q=p0,propID=mu_ID); fs%visc=viscf
         ! Configure pressure solver
         call param_read('Pressure iteration',fs%psolv%maxit)
         call param_read('Pressure tolerance',fs%psolv%rcvg)
         ! Configure implicit velocity solver
         call param_read('Implicit iteration',fs%implicit%maxit)
         call param_read('Implicit tolerance',fs%implicit%rcvg)
         ! Setup the solver
         call fs%setup(pressure_ils=pcg_amg,implicit_ils=pfmg)

         ! Initialize velocity field
         do k=lp%cfg%kmino_,lp%cfg%kmaxo_
            do j=lp%cfg%jmino_,lp%cfg%jmaxo_
               do i=lp%cfg%imino_,lp%cfg%imaxo_
                  fs%U(i,j,k)=0.0_WP
                  fs%V(i,j,k)=0.0_WP
                  fs%W(i,j,k)=0.0_WP
               end do
            end do
         end do
         call fs%interp_vel(Ui,Vi,Wi)
      end block initialize_fs

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
      ! create_pmesh: block
      !    integer :: i
      !    pmesh=partmesh(nvar=1,name='lpt')
      !    pmesh%varname(1)='radius'
      !    call lp%update_partmesh(pmesh)
      !    do i=1,lp%np_
      !       pmesh%var(1,i)=0.5_WP*lp%p(i)%d
      !    end do
      ! end block create_pmesh

      
      
      ! Add Ensight output
      create_ensight: block
         ! Create Ensight output from cfg
         ens_out=ensight(cfg=lp%cfg,name='evaporation')
         ! Create event for Ensight output
         ens_evt=event(time=time,name='Ensight output')
         call param_read('Ensight output period',ens_evt%tper)
         ! Add variables to output
         call ens_out%add_particle('particles',pmesh)
         call ens_out%add_vector('velocity',fs%U,fs%V,fs%W)
         call ens_out%add_vector('source',lp%srcU,lp%srcV,lp%srcW)
         call ens_out%add_scalar('epsp',lp%VF)
         call ens_out%add_scalar('yf',Yf_sc%SC)
         call ens_out%add_scalar('temperature',T_sc%SC)

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call lp%get_max()
         call fs%get_max()
         call Yf_sc%get_max()
         call T_sc%get_max()
         call fs%get_cfl(time%dt,time%cfl)
         ! Create simulation monitor
         mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
         call mfile%add_column(time%n,'N_time')
         call mfile%add_column(time%wt,'wall_t')
         call mfile%add_column(time%t,'time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(fs%Umax,'FS_Max_U')
         call mfile%add_column(fs%Vmax,'FS_Max_V')
         call mfile%add_column(fs%Wmax,'FS_Max_W')
         call mfile%add_column(Yf_sc%SCmin,'Yf_min')
         call mfile%add_column(Yf_sc%SCmax,'Yf_max')
         call mfile%add_column(T_sc%SCmin,'T_min')
         call mfile%add_column(T_sc%SCmax,'T_max')
         call mfile%add_column(lp%np,'N_part')
         ! call mfile%add_column(lp%VFmean,'Mean_VF')
         if (lp%cfg%amRoot) then
            call mfile%add_column(lp%Umin,'P_Umin')
            call mfile%add_column(lp%Umean,'P_Umean')
            call mfile%add_column(lp%Umax,'P_Umax')
            call mfile%add_column(lp%dmean,'P_dmean')
            call mfile%add_column(lp%Tmean,'P_Tmean')
         end if
         ! call mfile%add_column(Bm_debug,'Spalding')
         call mfile%write()

         cflFile = monitor(amroot=lp%cfg%amRoot,name='cfl')
         call cflFile%add_column(time%n,'N_time')
         call cflFile%add_column(time%t,'time')
         call cflFile%add_column(max(fs%CFLc_x,fs%CFLc_y,fs%CFLc_z),'max_CFLc')
         call cflFile%add_column(max(fs%CFLv_x,fs%CFLv_y,fs%CFLv_z),'max_CFLv')
         call cflFile%write()

         pfile=monitor(amroot=lp%cfg%amRoot,name='particle')
         call pfile%add_column(time%n,'T_step')
         call pfile%add_column(time%t,'Time')
         call pfile%add_column(time%dt,'dt')
         if (lp%cfg%amRoot) then
            call pfile%add_column(lp%p(1)%vel(1),'Part_U')
            call pfile%add_column(lp%p(1)%vel(2),'Part_V')
            call pfile%add_column(lp%p(1)%vel(3),'Part_W')
            call pfile%add_column(lp%p(1)%d,'Part_d')
            call pfile%add_column(lp%p(1)%T_d,'Part_T')
         end if
         call pfile%write()
      end block create_monitor

   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      logical :: amDone = .false.
      T_dt_cur = time%dt
      ! print*,'max rho1',maxval(Yf_sc%rho)
      ! Perform time integration
      do while (.not.time%done().and.(.not.amDone))
         
         ! Increment time         
         call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         ! if (time%dt.gt.T_dt_cur) time%dt = T_dt_cur ! modify dt given temperature change on mesh
         call time%increment()
         ! print*,'New dt:',time%dt

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
            
            ! Add mass, energy source terms
            add_mass_energy_src: block
               use fluidTable_class, only: Cp_ID
               integer :: i,j,k
               real(WP) :: Cp_g
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        ! Need in SC*kg/m^3
                        call gTab%evalProps(propOut=Cp_g,T_q=T_sc%SC(i,j,k),P_q=p0,propID=Cp_ID)
                        resT(i,j,k)=resT(i,j,k)+lp%srcE(i,j,k)/Cp_g!/cfg%vol(i,j,k)
                        resYf(i,j,k)=resYf(i,j,k)+lp%srcM(i,j,k)!/cfg%vol(i,j,k) ! *Yf_sc%rho(i,j,k)
                     end do
                  end do
               end do
            end block add_mass_energy_src

            ! Form implicit residual
            call T_sc%solve_implicit(time%dt,resT,fs%rhoU,fs%rhoV,fs%rhoW)
            call Yf_sc%solve_implicit(time%dt,resYf,fs%rhoU,fs%rhoV,fs%rhoW)
            ! print*,'here0B',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! Apply this residual
            T_sc%SC=2.0_WP*T_sc%SC-T_sc%SCold+resT
            Yf_sc%SC=2.0_WP*Yf_sc%SC-Yf_sc%SCold+resYf

            ! clip_Yf: block
            !    integer :: i,j,k
            !    do k=fs%cfg%kmin_,fs%cfg%kmax_
            !       do j=fs%cfg%jmin_,fs%cfg%jmax_
            !          do i=fs%cfg%imin_,fs%cfg%imax_
            !             if (Yf_sc%SC(i,j,k).lt.0.0_WP) Yf_sc%SC(i,j,k) = 0.0_WP
            !          end do
            !       end do
            !    end do
            ! end block clip_Yf

            ! print*,'max rho2',maxval(Yf_sc%rho)
            
            ! Apply other boundary conditions on the resulting field
            call T_sc%apply_bcond(time%t,time%dt)
            call Yf_sc%apply_bcond(time%t,time%dt)
            ! print*,'here0C',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
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
            ! print*,'here0D',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! Rescale scalar for conservation
            ! T_sc%SC=resT/T_sc%rho
            ! Yf_sc%SC=resYf/Yf_sc%rho

            ! UPDATE THE VISCOSITY
            update_visc: block
               use fluidTable_class, only: mu_ID
               integer :: i,j,k

               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        call gTab%evalProps(propOut=fs%visc(i,j,k),T_q=T_sc%SC(i,j,k),P_q=p0,propID=mu_ID)
                        
                     end do
                  end do
               end do
            end block update_visc
            ! UPDATE THE DIFFUSIVITY
            ! ===================================================
            
            ! ============ VELOCITY SOLVER ======================
            ! print*,'max rho4',maxval(Yf_sc%rho)
            ! Build n+1 density
            fs%rho=0.5_WP*(Yf_sc%rho+Yf_sc%rhoold)
            
            ! Build mid-time velocity and momentum
            fs%U=0.5_WP*(fs%U+fs%Uold); fs%rhoU=0.5_WP*(fs%rhoU+fs%rhoUold)
            fs%V=0.5_WP*(fs%V+fs%Vold); fs%rhoV=0.5_WP*(fs%rhoV+fs%rhoVold)
            fs%W=0.5_WP*(fs%W+fs%Wold); fs%rhoW=0.5_WP*(fs%rhoW+fs%rhoWold)
            
            ! Explicit calculation of drho*u/dt from NS
            call fs%get_dmomdt(resU,resV,resW)
            ! print*,'here0E',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            
            ! Assemble explicit residual
            resU=time%dtmid*resU-(2.0_WP*fs%rhoU-2.0_WP*fs%rhoUold)
            resV=time%dtmid*resV-(2.0_WP*fs%rhoV-2.0_WP*fs%rhoVold)
            resW=time%dtmid*resW-(2.0_WP*fs%rhoW-2.0_WP*fs%rhoWold)

            add_lpt_src: block
               integer :: i,j,k
               do k=fs%cfg%kmin_,fs%cfg%kmax_
                  do j=fs%cfg%jmin_,fs%cfg%jmax_
                     do i=fs%cfg%imin_,fs%cfg%imax_
                        resU(i,j,k)=resU(i,j,k)+sum(fs%itpr_x(:,i,j,k)*lp%srcU(i-1:i,j,k)) 
                        resV(i,j,k)=resV(i,j,k)+sum(fs%itpr_y(:,i,j,k)*lp%srcV(i,j-1:j,k))
                        resW(i,j,k)=resW(i,j,k)+sum(fs%itpr_z(:,i,j,k)*lp%srcW(i,j,k-1:k))
                     end do
                  end do
               end do
            end block add_lpt_src
            ! print*,'here0F',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! Form implicit residuals
            call fs%solve_implicit(time%dtmid,resU,resV,resW)
            ! print*,'here0G',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! print*,'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! Apply these residuals
            fs%U=2.0_WP*fs%U-fs%Uold+resU
            fs%V=2.0_WP*fs%V-fs%Vold+resV
            fs%W=2.0_WP*fs%W-fs%Wold+resW
            
            ! Apply other boundary conditions and update momentum
            call fs%apply_bcond(time%tmid,time%dtmid)
            call fs%rho_multiply()
            call fs%apply_bcond(time%tmid,time%dtmid)
            ! print*,'here0H',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! This is where dirichlet BCs would go

            ! Solve Poisson equation
            call Yf_sc%get_drhodt(dt=time%dt,drhodt=resYf)
            call fs%correct_mfr(drhodt=(resYf + fs%cfg%vol*lp%srcM/time%dtmid))
            call fs%get_div(drhodt=(resYf + fs%cfg%vol*lp%srcM/time%dtmid))
            fs%psolv%rhs=-fs%cfg%vol*fs%div/time%dtmid
            fs%psolv%sol=0.0_WP
            call fs%psolv%solve()
            call fs%shift_p(fs%psolv%sol)
            ! print*,'here0I',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! Correct momentum and rebuild velocity
            call fs%get_pgrad(fs%psolv%sol,resU,resV,resW)
            ! print*,'here0J',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            fs%P=fs%P+fs%psolv%sol
            ! print*,'here0K',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W)),'rhoU',maxval(fs%rhoU),'resU',maxval(resU)
            ! print*,'rhoU',maxval(fs%rhoU),'resU',maxval(resU),'dt*res',maxval(-time%dtmid*resU)
            fs%rhoU=fs%rhoU-time%dtmid*resU
            fs%rhoV=fs%rhoV-time%dtmid*resV
            fs%rhoW=fs%rhoW-time%dtmid*resW
            
            ! print*,'here0L',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            call fs%rho_divide
            ! ===================================================
            ! print*,'here0M',time%it,max(maxval(fs%U),maxval(fs%V),maxval(fs%W))
            ! Increment sub-iteration counter
            time%it=time%it+1
         end do

         ! Spawn particles
         spawn_particles: block
            use random, only: random_normal,random_uniform
            use mathtools, only: twoPi
            integer :: i,newP,np_old
            real(WP) :: theta
            !dp,dp_sig,rp,rp_sig,u_inj,u_inj_sig
            newP = 5
            np_old = lp%np
            if (lp%cfg%amRoot) then
               call lp%resize(np_old+newP)
               do i=(np_old+1),np_old+newP
                  ! Give id
                  ! lp%p(i)%id=int(i,8)
                  ! Set the diameter
                  lp%p(i)%d=random_normal(dp,dp_sig)
                  ! Set the temperature 
                  lp%p(i)%T_d = T_d
                  !! Assign random position in the domain
                  ! lp%p(i)%pos=[random_uniform(lp%cfg%x(lp%cfg%imin),lp%cfg%x(lp%cfg%imax+1)),&
                  ! &            random_uniform(lp%cfg%y(lp%cfg%jmin),lp%cfg%y(lp%cfg%jmax+1)),&
                  ! &            random_uniform(lp%cfg%z(lp%cfg%kmin),lp%cfg%z(lp%cfg%kmax+1))]
                  ! lp%p(i)%pos(3)=lp%cfg%zm(lp%cfg%kmin)
   
                  ! Assign position in center of domain
                  theta = random_uniform(0.0_WP,twoPi)
                  lp%p(i)%pos = random_normal(rp,rp_sig)*[0.0_WP,cos(theta),sin(theta)]
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
               end do
            end if
            ! Distribute particles
            call lp%sync()
                     ! Get initial particle volume fraction
            call lp%update_VF()

         end block spawn_particles
         ! Advance particles by dt
         ! print*,'U=',maxval(fs%U),'V=',maxval(fs%V),'W=',maxval(fs%W),'rho=',maxval(fs%rho),'visc=',maxval(fs%visc),'T=',maxval(T_sc%SC),'Yf=',maxval(Yf_sc%SC)
         call lp%advance(dt=time%dt,U=fs%U,V=fs%V,W=fs%W,rho=fs%rho,visc=fs%visc,T=T_sc%SC,Yf=Yf_sc%SC,Bm_debug=Bm_debug,lTab=lTab,gTab=gTab,p_therm=p0)
         ! print*,'max(srcM)',maxval(lp%srcM),'max(srcE)',maxval(lp%srcE)
         ! Use new source terms to determine if need to reduce next time step
         ! if (maxval(abs(lp%srcE)).gt.(0.0_WP)) then
         !    temp_dt_control: block
         !       real(WP) :: Cp_g = 1200.0_WP
         !       real(WP) :: dT_cur ! maximum temperature increase per time step across mesh
         !       dT_cur = maxval(T_sc%SC-T_sc%SCold)
         !       ! dT_cur = maxval(abs(lp%srcE/(fs%rho*Cp_g))) ! [K/timeStep]
         !       if (abs(dT_cur).gt.abs(maxTemp_dt)) T_dt_cur = time%dt*abs(maxTemp_dt/dT_cur)
         !       ! print*,'dT_cur [K]',abs(dT_cur),'T_dt_cur [s]',T_dt_cur
         !    end block temp_dt_control
         ! end if
         ! print*,'here2'
         if ((lp%np.eq.0)) amDone=.true.
         ! Output to ensight
         if (ens_evt%occurs()) then
            ensight_output: block
               integer :: i
               call lp%update_partmesh(pmesh)
               do i = 1,lp%np_
                  pmesh%var(1,i)=0.5_WP*lp%p(i)%d ! Droplet Radius
                  pmesh%var(2,i)=lp%p(i)%T_d ! Droplet Temperature
               end do
      
               call ens_out%write_data(time%t)
               end block ensight_output
         end if
         ! print*,'here3'
         ! Perform and output monitoring
         call lp%get_max()
         call fs%get_max()
         call Yf_sc%get_max()
         call T_sc%get_max()
         call cflFile%write()
         call mfile%write()
         call pfile%write()
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
      
   end subroutine simulation_final
   
end module simulation
