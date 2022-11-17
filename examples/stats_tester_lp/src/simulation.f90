!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use geometry,          only: cfg
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use fluidTable_class,  only: fluidTable
   use stats_class,       only: stats_parent
   implicit none
   private
   
   !> Get a LPT solver, stats, and corresponding time tracker
   type(lpt),         public :: lp
   type(stats_parent), public :: st
   type(timetracker), public :: time
   
   !> Ensight postprocessing
   type(partmesh) :: pmesh
   type(ensight)  :: ens_out
   type(event)    :: ens_evt
   
   !> Simulation monitor file
   type(monitor) :: mfile,pfile,cflFile
      
   public :: simulation_init,simulation_run,simulation_final
   
   real(WP) :: p0 ! ambient pressure
   ! real(WP) :: diff_T,diff_Yf ! non-turbulent diffusivities for scalar solvers

   ! !> Case quantities
   real(WP) :: dp,dp_sig,r_inj,r_inj_sig,u_inj,u_inj_sig,T_d,inj_rate!,r_jet,v_coflow,T_coflow,Yf_coflow
   real(WP) :: T_jet,Yf_jet,v_jet
   integer :: nInj

   !> Quantities for flow & scalar solvers:
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW,resT,resYf,visc,rho
   ! real(WP), dimension(:,:,:),   allocatable :: Ui,Vi,Wi
   ! real(WP), dimension(:,:,:,:), allocatable :: SR

   type(fluidTable) :: lTab,gTab


   
contains

   function station_1_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      logical :: isIn
      isIn = .false.
      if (pg%xm(i).gt.0.and.pg%xm(i).lt.3.0_WP*pg%dx(1)) isIn = .true.
   end function station_1_locator

   function station_2_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: xLoc
      logical :: isIn
      isIn = .false.
      xLoc = 0.15_WP
      if (pg%xm(i).gt.xLoc.and.pg%xm(i).lt.xLoc+10.0_WP*pg%dx(1)) isIn = .true.
   end function station_2_locator

   function station_3_locator(pg,i,j,k) result(isIn)
      use pgrid_class, only: pgrid
      implicit none
      class(pgrid), intent(in) :: pg
      integer, intent(in) :: i,j,k
      real(WP) :: xLoc
      logical :: isIn
      isIn = .false.
      xLoc = 0.35_WP
      if (pg%xm(i).gt.xLoc.and.pg%xm(i).lt.xLoc+1.0_WP*pg%dx(1)) isIn = .true.
   end function station_3_locator

   ! !> Define here our equation of state - rho(T,mass)
   ! subroutine get_rho()
   !    implicit none
   !    ! real(WP), intent(in) :: mass
   !    integer :: i,j,k
   !    real(WP) :: W_g,W_l,R_cst
   !    R_cst = 8.314472_WP ! Gas constant [J/(kg*K)]
   !    W_g = gTab%MW ! Carrier gas molar weight [kg/mol]
   !    W_l = lTab%MW ! Drop molar weight [kg/mol]

   !    ! Calculate density
   !    do k=T_sc%cfg%kmino_,T_sc%cfg%kmaxo_
   !       do j=T_sc%cfg%jmino_,T_sc%cfg%jmaxo_
   !          do i=T_sc%cfg%imino_,T_sc%cfg%imaxo_
   !             if (Yf_sc%SC(i,j,k).ge.0.0_WP) then
   !                T_sc%rho(i,j,k)=p0/(R_cst*(Yf_sc%SC(i,j,k)/W_l+(1.0_WP-Yf_sc%SC(i,j,k))/W_g)*T_sc%SC(i,j,k))! + lp%srcM(i,j,k)
   !             else
   !                T_sc%rho(i,j,k)=p0*W_g/(R_cst*T_sc%SC(i,j,k))! + lp%srcM(i,j,k)
   !             end if
   !             Yf_sc%rho(i,j,k)=T_sc%rho(i,j,k)
   !          end do
   !       end do
   !    end do
   ! end subroutine get_rho

   !> Initialization of problem solver
   subroutine simulation_init
      use coolprop
      use param, only: param_read
      implicit none

      read_global_vars : block
         call param_read('Gas jet velocity',             v_jet)
         ! call param_read('Gas coflow velocity',          v_coflow)
         ! call param_read('Gas jet radius',               r_jet)
         call param_read('Pressure',                     p0)
         call param_read('Gas jet temperature',          T_jet)
         ! call param_read('Gas coflow temperature',       T_coflow)
         call param_read('Gas jet fuel mass fraction',   Yf_jet)
         ! call param_read('Gas coflow fuel mass fraction',Yf_coflow)

         call param_read('Droplet mean diameter',              dp)
         call param_read('Droplet diameter standard deviation',dp_sig)
         call param_read('Droplet mean radius',                r_inj)
         call param_read('Droplet radius standard deviation',  r_inj_sig)
         call param_read('Droplet temperature',                T_d)
         call param_read('Droplet mean velocity',              u_inj)
         call param_read('Droplet velocity standard deviation',u_inj_sig)
         call param_read('Droplet injection rate',             inj_rate)

         ! call param_read('Thermal diffusivity',diff_T)
         ! call param_read('Fuel diffusivity',   diff_Yf)


      end block read_global_vars

      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         ! time%dt=min(time%dtmax,0.1_WP*cfg%dx(1)/v_jet)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker
      print*,'initializing fluid tables'
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


         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Gas !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
   
      print*,'initializing lpt'
      ! Initialize our LPT
      initialize_lpt: block
         use fluidTable_class, only: rho_ID
         use mathtools, only: twoPi
         use random, only: random_normal,random_uniform
         real(WP) :: theta
         integer :: i,np
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')         
         ! Get droplet density from the fluid table
         call lTab%evalProps(propOut=lp%rho,T_q=T_d,P_q=p0,propID=rho_ID)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         np = 0!int(floor(inj_rate*time%dt))
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
               ! Assign position in center of domain
               theta = random_uniform(0.0_WP,twoPi)
               lp%p(i)%pos = random_normal(r_inj,r_inj_sig)*[0.0_WP,cos(theta),sin(theta)]
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
         lp%Tcol=5.0_WP*time%dt
      end block initialize_lpt

      initialize_res: block
         allocate(resU (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resU=0.0_WP
         allocate(resV (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resV=0.0_WP
         allocate(resW (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resW=0.0_WP
         allocate(resT (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resT=0.0_WP
         allocate(resYf(cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); resYf=0.0_WP
         allocate(visc (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); visc=0.0_WP
         allocate(rho  (cfg%imino_:cfg%imaxo_,cfg%jmino_:cfg%jmaxo_,cfg%kmino_:cfg%kmaxo_)); rho=0.0_WP
         
      end block initialize_res
      
      ! Initialize scalar solvers
      initialize_sc: block
      print*,'initializing scalars'
     
         resT = T_jet
         resYf= Yf_jet
      
      end block initialize_sc

      ! Initialize the velocity field
      initialize_fs: block
         
         ! use mathtools, only: twoPi
         ! use lowmach_class, only: clipped_neumann,dirichlet,bcond
         ! use ils_class, only: pcg_amg,pfmg,gmres_amg,pcg_pfmg
         use mathtools, only : twoPi
         use fluidTable_class, only: mu_ID,rho_ID
         ! type(bcond),pointer :: mybc
         integer :: i,j,k,n
         real(WP) :: viscf,rhof
         
         call gTab%evalProps(propOut=viscf,T_q=T_jet,P_q=p0,propID=mu_ID);  visc=viscf
         call gTab%evalProps(propOut=rhof, T_q=T_jet,P_q=p0,propID=rho_ID); rho=rhof

         print*,'initializing velocity field'
         do k=cfg%kmino_,cfg%kmaxo_
            do j=cfg%jmino_,cfg%jmaxo_
               do i=cfg%imino_,cfg%imaxo_
                  resU(i,j,k)=v_jet
                  resV(i,j,k)=0.1_WP*v_jet*sin(twoPi*cfg%y(j)/(cfg%y(cfg%jmax)-cfg%y(cfg%jmin)))
                  resW(i,j,k)=0.1_WP*v_jet*sin(twoPi*cfg%z(k)/(cfg%z(cfg%kmax)-cfg%z(cfg%kmin)))
               end do
            end do
         end do
         
      end block initialize_fs

      print*,'entering stats initialization'
      ! Create stats collector
      initialize_stats : block
         st = stats_parent(cfg=cfg,name='Particle Stats')

         ! Add our stations
         call st%add_station(name='x=0',dim='yz',locator=station_1_locator)
         call st%add_station(name='x=0.2',dim='yz',locator=station_2_locator)
         call st%add_station(name='x=0.39',dim='yz',locator=station_3_locator)

         ! Initialize stats
         call st%init_stats(lp=lp)

      end block initialize_stats


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
         call ens_out%add_vector('velocity',resU,resV,resW)
         call ens_out%add_vector('srcUVW',lp%srcU,lp%srcV,lp%srcW)
         ! call ens_out%add_scalar('epsp',lp%VF)
         ! call ens_out%add_scalar('yf',Yf_sc%SC)
         ! call ens_out%add_scalar('temperature',T_sc%SC)
         ! call ens_out%add_scalar('density',T_sc%rho)
         ! call ens_out%add_scalar('srcM',lp%srcM)
         ! call ens_out%add_scalar('srcE',lp%srcE)

         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight
      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call lp%get_max()
         ! call fs%get_max()
         ! call Yf_sc%get_max()
         ! call T_sc%get_max()
         ! call fs%get_cfl(time%dt,time%cfl)
         ! Create simulation monitor
         mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
         call mfile%add_column(time%n,'N_time')
         call mfile%add_column(time%wt,'wall_t')
         call mfile%add_column(time%t,'time')
         call mfile%add_column(time%dt,'dt')
         ! call mfile%add_column(time%cfl,'CFL')
         ! call mfile%add_column(fs%Umax,'FS_Max_U')
         ! call mfile%add_column(fs%Vmax,'FS_Max_V')
         ! call mfile%add_column(fs%Wmax,'FS_Max_W')
         ! call mfile%add_column(Yf_sc%SCmin,'Yf_min')
         ! call mfile%add_column(Yf_sc%SCmax,'Yf_max')
         ! call mfile%add_column(T_sc%SCmin,'T_min')
         ! call mfile%add_column(T_sc%SCmax,'T_max')
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

         cflFile = monitor(amroot=lp%cfg%amRoot,name='cfl')
         call cflFile%add_column(time%n,'N_time')
         call cflFile%add_column(time%t,'time')
         ! call cflFile%add_column(max(fs%CFLc_x,fs%CFLc_y,fs%CFLc_z),'max_CFLc')
         ! call cflFile%add_column(max(fs%CFLv_x,fs%CFLv_y,fs%CFLv_z),'max_CFLv')
         call cflFile%write()

         ! pfile=monitor(amroot=lp%cfg%amRoot,name='particle')
         ! call pfile%add_column(time%n,'T_step')
         ! call pfile%add_column(time%t,'Time')
         ! call pfile%add_column(time%dt,'dt')
         ! if (lp%cfg%amRoot) then
         !    call pfile%add_column(lp%p(1)%vel(1),'Part_U')
         !    call pfile%add_column(lp%p(1)%vel(2),'Part_V')
         !    call pfile%add_column(lp%p(1)%vel(3),'Part_W')
         !    call pfile%add_column(lp%p(1)%d,'Part_d')
         !    call pfile%add_column(lp%p(1)%T_d,'Part_T')
         ! end if
         ! call pfile%write()
      end block create_monitor
      print*,'good init'
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      logical :: amDone = .false.
      ! print*,'max rho1',maxval(Yf_sc%rho)
      ! Perform time integration
      do while (.not.time%done().and.(.not.amDone))
         
         ! Increment time         
         ! call fs%get_cfl(time%dt,time%cfl)
         call time%adjust_dt()
         call time%increment()
         
         ! Spawn particles
         spawn_particles: block
            use random, only: random_normal,random_uniform,random_lognormal
            use mathtools, only: twoPi
            integer :: i,np,nAdd
            real(WP) :: theta,r
            if (lp%cfg%amRoot) then
               nAdd = int(floor(inj_rate*time%t))-nInj
               do i=1,nAdd
                  np=lp%np_+1; call lp%resize(np)
                  ! Give id
                  lp%p(i)%id=nInj+i
                  ! Set the diameter
                  lp%p(np)%d=random_lognormal(dp,dp_sig)
                  ! Set the temperature 
                  lp%p(np)%T_d = T_d
                  ! Assign position in center of domain
                  if (cfg%ny.eq.1) then
                     lp%p(np)%pos = [random_uniform(0.0_WP,lp%cfg%dx(1)),0.0_WP,random_uniform(-r_inj,r_inj)]
                  elseif (cfg%nz.eq.1) then
                     lp%p(np)%pos = [random_uniform(0.0_WP,lp%cfg%dx(1)),random_uniform(-r_inj,r_inj),0.0_WP]
                  else
                     theta = random_uniform(0.0_WP,twoPi)
                     r = r_inj*sqrt(random_uniform(0.0_WP,1.0_WP))
                     lp%p(np)%pos = [random_uniform(0.0_WP,lp%cfg%dx(1)),r*cos(theta),r*sin(theta)]
                  end if
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
         ! print*,'U=',maxval(fs%U),'V=',maxval(fs%V),'W=',maxval(fs%W),'rho=',maxval(fs%rho),'visc=',maxval(fs%visc),'T=',maxval(T_sc%SC),'Yf=',maxval(Yf_sc%SC)
         ! print*,'rank',cfg%rank,'entering lp%advance with Yf_max,min:',maxval(Yf_sc%SC),minval(Yf_sc%SC)
         call lp%advance(dt=time%dt,U=resU,V=resV,W=resW,rho=rho,visc=visc,T=resT,Yf=resYf,lTab=lTab,gTab=gTab,p_therm=p0)
         ! print*,'rank',cfg%rank,'max(srcM)',maxval(lp%srcM),'min(srcM)',minval(lp%srcM),'max(srcE)',maxval(lp%srcE),'min(srcE)',minval(lp%srcE)
         ! print*,'here2'
         if ((lp%np.eq.0).and.(time%t.gt.3.0_WP)) amDone=.true.
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
         print*,'enter sample'
         call st%sample_stats(lp=lp,dt=time%dt)
         print*,'exit sample'
         ! call fs%get_max()
         ! call Yf_sc%get_max()
         ! call T_sc%get_max()
         call cflFile%write()
         call mfile%write()
         call pfile%write()
         ! print*,'here4'
         ! imDying : block
         !    use messager, only : die
         !    call die('Dying now for forensic reasons')
         ! end block imDying
      end do
      ! if (cfg%rank.eq.1) then
      !    slow_down : block
      !       integer :: i,j,k
      !       real(WP) :: a,b,c
      !       a = 1.0_WP
      !       b = 2.0_WP
      !       do i=1,1000
      !          do j=1,1000
      !             do k=1,1000
      !                c = a*sqrt(b*i**2+j**2)
      !             end do
      !          end do
      !       end do

      !    end block slow_down
      ! end if
      ! print_stats : block
      !    if (st%stats(1)%amIn) then
      !       print*,'Number, station 1'
      !       print*,st%stats(1)%vals(1,:,:,:)
      !       print*,'Mean diameter [micron], station 1'
      !       print*,st%stats(1)%vals(2,:,:,:)*1.0e6_WP
      !       print*,'diameter varience [mm^2], station 1'
      !       print*,st%stats(1)%vals(3,:,:,:)*1.0e6_WP
      !    end if
      !    if (st%stats(2)%amIn) then
      !       print*,'Mean diameter [micron], station 2'
      !       print*,st%stats(2)%vals(2,:,:,:)*1.0e6_WP
      !       print*,'diameter varience [mm^2], station 2'
      !       print*,st%stats(2)%vals(3,:,:,:)*1.0e6_WP
      !    end if
      !    if (st%stats(3)%amIn) then
      !       print*,'Mean diameter [micron], station 3'
      !       print*,st%stats(3)%vals(2,:,:,:)*1.0e6_WP
      !       print*,'diameter varience [mm^2], station 3'
      !       print*,st%stats(3)%vals(3,:,:,:)*1.0e6_WP
      !    end if
      ! end block print_stats
      call st%post_process()
      if (cfg%rank.eq.1) then
         slow_down2 : block
            integer :: i,j,k
            real(WP) :: a,b,c
            a = 1.0_WP
            b = 2.0_WP
            do i=1,2000
               do j=1,2000
                  do k=1,2000
                     c = a*sqrt(b*i**2+j**2)
                  end do
               end do
            end do

         end block slow_down2
      end if
      print_stats2 : block
      if (st%stats(1)%amIn) then
         print*,'Number, station 1'
         print*,st%stats(1)%vals(1,:,:,:)
         print*,'Mean diameter [micron], station 1'
         print*,st%stats(1)%vals(2,:,:,:)*1.0e6_WP
         print*,'diameter varience [mm^2], station 1'
         print*,st%stats(1)%vals(3,:,:,:)*1.0e6_WP
      end if
      if (st%stats(2)%amIn) then
         print*,'Mean diameter [micron], station 2'
         print*,st%stats(2)%vals(2,:,:,:)*1.0e6_WP
         print*,'diameter varience [mm^2], station 2'
         print*,st%stats(2)%vals(3,:,:,:)*1.0e6_WP
      end if
      if (st%stats(3)%amIn) then
         print*,'Mean diameter [micron], station 3'
         print*,st%stats(3)%vals(2,:,:,:)*1.0e6_WP
         print*,'diameter varience [mm^2], station 3'
         print*,st%stats(3)%vals(3,:,:,:)*1.0e6_WP
      end if
      end block print_stats2

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
