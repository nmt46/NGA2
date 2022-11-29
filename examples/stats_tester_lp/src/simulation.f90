!> Various definitions and tools for running an NGA2 simulation
module simulation
   use precision,         only: WP
   use string,            only: str_medium
   use geometry,          only: cfg
   use lpt_class,         only: lpt
   use timetracker_class, only: timetracker
   use ensight_class,     only: ensight
   use partmesh_class,    only: partmesh
   use event_class,       only: event
   use monitor_class,     only: monitor
   use fluidTable_class,  only: fluidTable
   use stats_class,       only: stats_parent
   use event_class,       only: event
   use datafile_class,    only: datafile
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

   !> Various things for saving restart/stats data
   type(event)    :: save_rst_evt,save_stat_evt
   type(datafile) :: df
   logical :: restarted
   character(len=str_medium) :: stats_name
   character(len=str_medium) :: lpt_name
      
   public :: simulation_init,simulation_run,simulation_final
   
   real(WP) :: p0 ! ambient pressure

   ! !> Case quantities
   real(WP) :: dp,dp_sig,r_inj,r_inj_sig,u_inj,u_inj_sig,T_d,inj_rate
   real(WP) :: T_jet,Yf_jet,v_jet
   integer :: nInj

   !> Quantities for flow & scalar "solvers":
   real(WP), dimension(:,:,:),   allocatable :: resU,resV,resW,resT,resYf,visc,rho

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

   function print_wall_time() result(wt)
      use parallel, only: wtinit,parallel_time
      real(WP) :: wt
      wt = parallel_time()-wtinit
   end function print_wall_time

   !> Initialization of problem solver
   subroutine simulation_init
      use coolprop
      use param, only: param_read
      implicit none

      read_global_vars : block
         call param_read('Gas jet velocity',             v_jet)
         call param_read('Pressure',                     p0)
         call param_read('Gas jet temperature',          T_jet)
         call param_read('Gas jet fuel mass fraction',   Yf_jet)

         call param_read('Droplet mean diameter',              dp)
         call param_read('Droplet diameter standard deviation',dp_sig)
         call param_read('Droplet mean radius',                r_inj)
         call param_read('Droplet radius standard deviation',  r_inj_sig)
         call param_read('Droplet temperature',                T_d)
         call param_read('Droplet mean velocity',              u_inj)
         call param_read('Droplet velocity standard deviation',u_inj_sig)
         call param_read('Droplet injection rate',             inj_rate)

      end block read_global_vars

      
      ! Initialize time tracker with 1 subiterations
      initialize_timetracker: block
         time=timetracker(amRoot=cfg%amRoot)
         call param_read('Max timestep size',time%dtmax)
         call param_read('Max cfl number',time%cflmax)
         call param_read('Max time',time%tmax)
         time%dt=time%dtmax
         time%itmax=2
      end block initialize_timetracker

      ! Restart if needed
      restart_sim : block
         character(len=str_medium) :: dir_restart
         ! Create event for restarting based on our timetracker
         save_rst_evt=event(time,'Restart output')
         call param_read('Restart output period',save_rst_evt%tper)
         ! Are we restarting?
         call param_read(tag='Restart from',val=dir_restart,short='r',default='')
         restarted=.false.; if (len_trim(dir_restart).gt.0) restarted=.true.
         if (restarted) then
            ! Pull datafile
            df = datafile(pg=cfg,fdata=trim(adjustl(dir_restart))//'/'//'data.f')
            ! Set names
            stats_name = trim(adjustl(dir_restart))//'/'//'data.stats'
            lpt_name=trim(adjustl(dir_restart))//'/'//'data.lpt'
            ! Set timetracker t, dt
            call df%pullval(name='t', val=time%t)
            call df%pullval(name='dt',val=time%dt)
         else ! Otherwise, let's get the datafile set up
            df = datafile(pg=cfg,filename=trim(cfg%name),nval=3,nvar=3)
            df%valname(1)='t'
            df%valname(2)='dt'
            df%valname(3)='nInj'
            df%varname(1)='U'
            df%varname(2)='V'
            df%varname(3)='W'
         end if
      end block restart_sim

      ! Store stats periodically independently of restarts
      init_stats_cache : block
         save_stat_evt = event(time,'Statistics output')
         call param_read('Statistics output period',save_stat_evt%tper)
      end block init_stats_cache

      initialize_fluid_properties: block
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
         real(WP) :: buf
         ! Create solver
         lp=lpt(cfg=cfg,name='LPT')         
         ! Get droplet density from the fluid table
         call lTab%evalProps(propOut=lp%rho,T_q=T_d,P_q=p0,propID=rho_ID)
         ! Set filter scale to 3.5*dx
         lp%filter_width=3.5_WP*cfg%min_meshsize
         if (restarted) then
            call lp%read(filename=trim(lpt_name))
            call df%pullval(name='nInj',val=buf); nInj=int(buf)
         else
            nInj = 0
         end if
         ! Distribute particles
         call lp%sync()
         ! Get initial particle volume fraction
         call lp%update_VF()
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
      initialize_vel: block
         use mathtools, only : twoPi
         use fluidTable_class, only: mu_ID,rho_ID
         use random, only : random_normal
         integer :: i,j,k
         real(WP) :: viscf,rhof
         
         call gTab%evalProps(propOut=viscf,T_q=T_jet,P_q=p0,propID=mu_ID);  visc=viscf
         call gTab%evalProps(propOut=rhof, T_q=T_jet,P_q=p0,propID=rho_ID); rho=rhof
         if (restarted) then
            call df%pullvar(name='U',var=resU)
            call df%pullvar(name='V',var=resV)
            call df%pullvar(name='W',var=resW)
         else
            do k=cfg%kmino_,cfg%kmaxo_
               do j=cfg%jmino_,cfg%jmaxo_
                  do i=cfg%imino_,cfg%imaxo_
                     resU(i,j,k)=random_normal(v_jet,v_jet*0.1_WP) ! Set random x velocity with absolutely zero regard for it solving Navier-Stokes
                     if (cfg%ny.ne.1) then
                        resV(i,j,k)=0.1_WP*v_jet*sin(twoPi*cfg%y(j)/(cfg%y(cfg%jmax)-cfg%y(cfg%jmin)))
                     else
                        resV(i,j,k)=0.0_WP
                     end if
                     if (cfg%nz.ne.1) then
                        resW(i,j,k)=0.1_WP*v_jet*sin(twoPi*cfg%z(k)/(cfg%z(cfg%kmax)-cfg%z(cfg%kmin)))
                     else
                        resW(i,j,k)=0.0_WP
                     end if
                  end do
               end do
            end do
         end if
      end block initialize_vel

      ! Create stats collector
      initialize_stats : block
         use messager, only: die
         st = stats_parent(cfg=cfg,name='Particle Stats')
         if (restarted) call st%restart_from_file(stats_name)
         ! Add our stations
         call st%add_station(name='x=0',dim='yz',locator=station_1_locator)
         call st%add_station(name='x=0.2',dim='yz',locator=station_2_locator)
         call st%add_station(name='x=0.39',dim='yz',locator=station_3_locator)
         ! Add pointers to eulerian arrays
         call st%add_array(name='U',array=resU)
         call st%add_array(name='V',array=resV)
         call st%add_array(name='W',array=resW)
         ! Add pointer to lp solver
         call st%set_lp(lp=lp)
         ! Add definitions for stats in terms of arrays and lp properties
         call st%add_definition(name='U',     def='U')
         call st%add_definition(name='U^2',   def='U*U')
         call st%add_definition(name='p_d',   def='p_d')
         call st%add_definition(name='p_d^2', def='p_d*p_d')
         call st%add_definition(name='p_T',   def='p_T')
         call st%add_definition(name='p_U',   def='p_U')
         call st%add_definition(name='p_U^2', def='p_U*p_U')
         call st%add_definition(name='p_U*U', def='p_U*U')
         call st%add_definition(name='p_rad', def='p_rad')
         call st%add_definition(name='p_rad^2',def='p_rad*p_rad')
         ! Finalize initialization (really just allocating a couple arrays)
         call st%init_stats()
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
         ! Output to ensight
         if (ens_evt%occurs()) call ens_out%write_data(time%t)
      end block create_ensight

      ! Create a monitor file
      create_monitor: block
         ! Prepare some info about fields
         call lp%get_max()
         ! Create simulation monitor
         mfile=monitor(amroot=lp%cfg%amRoot,name='simulation')
         call mfile%add_column(time%n,'N_time')
         call mfile%add_column(time%wt,'wall_t')
         call mfile%add_column(time%t,'time')
         call mfile%add_column(time%dt,'dt')
         call mfile%add_column(lp%np,'N_part')
         call mfile%add_column(nInj,'N_inj')
         call mfile%add_column(lp%nEvap,'N_evap')
         call mfile%add_column(lp%Umean,'P_Umean')
         call mfile%add_column(lp%Umax,'P_Umax')
         call mfile%add_column(lp%dmean,'P_dmean')
         call mfile%add_column(lp%Tmean,'P_Tmean')
         call mfile%add_column(lp%Tmax,'P_Tmax')
         call mfile%write()
      end block create_monitor
      print*,'good init'
   end subroutine simulation_init
   
   
   !> Perform an NGA2 simulation
   subroutine simulation_run
      implicit none
      logical :: amDone = .false.

      ! Perform time integration
      do while (.not.time%done().and.(.not.amDone))
         ! Increment time         
         call time%adjust_dt()
         call time%increment()

         ! Mess with the velocity
         if (mod(time%n,10).eq.0) then
            mess_vel : block
            use random, only :random_normal
            integer :: i,j,k
               do k=cfg%kmino_,cfg%kmaxo_
                  do j=cfg%jmino_,cfg%jmaxo_
                     do i=cfg%imino_,cfg%imaxo_
                        resU(i,j,k)=random_normal(v_jet,v_jet*0.1_WP)
                     end do
                  end do
               end do
            end block mess_vel
         end if

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
                  ! Assign uniformly distributed random position with r_inj
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
         call lp%advance(dt=time%dt,U=resU,V=resV,W=resW,rho=rho,visc=visc,T=resT,Yf=resYf,lTab=lTab,gTab=gTab,p_therm=p0)
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

         ! Perform and output monitoring
         call lp%get_max()
         ! if (cfg%amRoot) print*,'enter sample, wt=',print_wall_time()
         call st%sample_stats(dt=time%dt)
         ! if (cfg%amRoot) print*,'exit sample,  wt=',print_wall_time()
         call mfile%write()
         call pfile%write()

         ! Save a restart if needed
         save_restart : block
            character(len=str_medium) :: dirname,timestamp
            if (save_rst_evt%occurs()) then
               ! Prefix for files
               dirname='restart_'; write(timestamp,'(es12.5)') time%t
               ! Prepare a new directory
               if (cfg%amRoot) call execute_command_line('mkdir -p '//trim(adjustl(dirname))//trim(adjustl(timestamp)))
               call st%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.stats')
               call lp%write(filename=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.lpt')
               call df%pushval(name='t'   ,val=time%t)
               call df%pushval(name='dt'  ,val=time%dt)
               call df%pushval(name='nInj',val=real(nInj,WP))
               call df%pushvar(name='U'   ,var=resU)
               call df%pushvar(name='V'   ,var=resV)
               call df%pushvar(name='W'   ,var=resW)
               call df%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.f')
            end if
         end block save_restart

         ! Cache statistics data if needed
         cache_stats : block
            character(len=str_medium) :: dirname,timestamp
            if (save_stat_evt%occurs()) then
               ! File pathnames
               dirname='stats_cache/stat_'; write(timestamp,'(es12.5)') time%t
               ! Prepare the directory
               if (cfg%amRoot) call execute_command_line('mkdir -p '//trim(adjustl(dirname))//trim(adjustl(timestamp)))
               call st%write(fdata=trim(adjustl(dirname))//trim(adjustl(timestamp))//'/'//'data.stats')
            end if
         end block cache_stats
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
      deallocate(rho,visc,resU,resV,resW)
      if (cfg%rank.eq.0) print*,'Clean termination'
   end subroutine simulation_final
   
end module simulation
