! Collect specialized statistics for sprays
module stats_class
    use precision, only: WP
    use string,    only: str_medium
    use config_class, only: config
    use mpi_f08
    implicit none
    private


    public :: stats_parent,station ! Show stats_parent to the compiler so it stops yelling at me

    ! Contains information on sampling regions for stats gathering
    ! Stations are assumed to be rectangular prisms in the cartesian mesh, though
    !   arbitrary volumes fixed to the lattice are allowed
    type :: station
        character(len=str_medium) :: name = 'Unknown_station'               ! Name of our station
        logical :: amIn = .false.                                           ! Is this processor in this station?
        integer :: dim                                                      ! Dimension of station, x=1,y=2,z=3,xy=4,yz=5,zx=6,xyz=7
        integer :: imin ,imax ,jmin ,jmax ,kmin ,kmax                       ! Global min/maxes in x,y,z. Indices correspond to spray_parent%cfg
        integer :: imin_,imax_,jmin_,jmax_,kmin_,kmax_                      ! Local min/maxes in x,y,z. Indices correspond to spray_parent%cfg
        integer :: n_stats                                                  ! Number of statistics being gathered
        character(len=str_medium), dimension(:), allocatable :: names       ! Names associated with gather statistics
        real(WP), dimension(:,:,:,:), allocatable :: vals                   ! Values for each statistic (before post-processing)
        type(MPI_Comm) :: comm = MPI_COMM_NULL                              ! Communicator for processors on this station
        logical :: amRoot
    end type station

    
    type :: stats_parent
        type(config), pointer :: cfg                                ! The config around which spray_parent is built
        character(len=str_medium) :: name='Unknown_stats_parent'    ! Name of our stat parent
        integer :: nLoc                                             ! Number of stations for gathering statistics
        type(station), dimension(:),allocatable :: stats            ! Array of our statistics
        real(WP) :: sum_time                                        ! Amount of time over which we have averaged
        logical, dimension(5) :: my_solvers                         ! List of whether we use each solver by logical value in each index
                                                                        ! 1=lowmach,2=T_sc,3=Yf_sc,4=incomp,5=lp
    contains
        procedure :: add_station            ! Add a station for stats measuring
        procedure :: init_stats             ! Initialize statistics. To be called after adding all stations
        procedure :: sample_stats           ! Sample stats over all stations
        procedure :: cache_data             ! Store data to the disk
        procedure :: read_data              ! Read raw data from disk
        procedure :: post_process           ! Process data from binary file into usable forms
        procedure :: reset_stats            ! Reset statistics, maintaining stations
        

    end type stats_parent

    interface stats_parent
        procedure constructor
    end interface stats_parent
    
contains

    function constructor(cfg,name) result(self)
        implicit none
        type(stats_parent) :: self
        class(config), target, intent(in) :: cfg
        character(len=*), optional :: name

        ! Set the stats_parent's name
        if (present(name)) self%name=trim(adjustl(name))

        ! Point to our config object
        self%cfg=>cfg

        ! Initialize stations as null
        self%nLoc=0
        
    end function constructor

    ! Subroutine for adding stations to the stats sampling
    subroutine add_station(this,name,dim,locator)
        use string, only : lowercase
        use messager, only : die
        use mpi_f08,  only: MPI_ALLREDUCE,MPI_MAX,MPI_MIN,MPI_SUM,MPI_GROUP_SIZE,MPI_Group_rank,MPI_Comm_split,MPI_UNDEFINED
        use parallel, only: MPI_REAL_WP 
        implicit none
        class(stats_parent),intent(inout) :: this
        character(len=*),intent(in) :: dim
        character(len=*),intent(in),optional :: name
        integer :: i,j,k,ierr,n
        integer :: nproc,myRank,key,color
        real(WP) :: buf,nproc_stat
        integer, dimension(:), allocatable :: new_ranks
        type(station), dimension(:), allocatable :: tmp
        interface
            logical function locator(pargrid,ind1,ind2,ind3)
                use pgrid_class, only: pgrid
                class(pgrid), intent(in) :: pargrid
                integer, intent(in) :: ind1,ind2,ind3
            end function locator
        end interface
        print*,'add0, i=',this%nLoc+1
        ! Reallocate station array
        if (allocated(this%stats)) then
            allocate(tmp(this%nLoc+1))
            tmp(1:this%nLoc)=this%stats
            call move_alloc(tmp,this%stats)
        else
            allocate(this%stats(1))
        end if

        this%nLoc=this%nLoc+1

        

        if (present(name)) this%stats(this%nLoc)%name=trim(adjustl(name))

        ! Parse dimension of station
        select case (lowercase(dim))
        case ('0D');      this%stats(this%nLoc)%dim=0
        case ('x');       this%stats(this%nLoc)%dim=1
        case ('y');       this%stats(this%nLoc)%dim=2
        case ('z');       this%stats(this%nLoc)%dim=3
        case ('xy','yx'); this%stats(this%nLoc)%dim=4
        case ('yz','zy'); this%stats(this%nLoc)%dim=5
        case ('zx','xz'); this%stats(this%nLoc)%dim=6
        case ('xyz','xzy','yxz','yzx','zxy','zyx'); this%stats(this%nLoc)%dim=7
        case default; call die('[stats add_station] Unrecognized dimension of station.')
        end select
        print*,'add1, i=',this%nLoc
        ! Loop over mesh to find local bounding indices for the station
        this%stats(this%nLoc)%imin_ = HUGE(1)
        this%stats(this%nLoc)%jmin_ = HUGE(1)
        this%stats(this%nLoc)%kmin_ = HUGE(1)
        this%stats(this%nLoc)%imax_ =-HUGE(1)
        this%stats(this%nLoc)%jmax_ =-HUGE(1)
        this%stats(this%nLoc)%kmax_ =-HUGE(1)
        do i=this%cfg%imin_,this%cfg%imax_
            do j=this%cfg%jmin_,this%cfg%jmax_
                do k=this%cfg%kmin_,this%cfg%kmax_
                    if (locator(this%cfg,i,j,k)) then
                        this%stats(this%nLoc)%amIn = .true.
                        this%stats(this%nLoc)%imin_ = min(this%stats(this%nLoc)%imin_,i)
                        this%stats(this%nLoc)%jmin_ = min(this%stats(this%nLoc)%jmin_,j)
                        this%stats(this%nLoc)%kmin_ = min(this%stats(this%nLoc)%kmin_,k)
                        this%stats(this%nLoc)%imax_ = max(this%stats(this%nLoc)%imax_,i)
                        this%stats(this%nLoc)%jmax_ = max(this%stats(this%nLoc)%jmax_,j)
                        this%stats(this%nLoc)%kmax_ = max(this%stats(this%nLoc)%kmax_,k)
                    end if
                end do
            end do
        end do
        print*,'add2, i=',this%nLoc
        nproc_stat = 0.0_WP
        if (this%stats(this%nLoc)%amIn) nproc_stat = 1.0_WP
        ! Get the global bounding indices
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%imin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(this%nLoc)%imin=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%jmin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(this%nLoc)%jmin=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%kmin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(this%nLoc)%kmin=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%imax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(this%nLoc)%imax=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%jmax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(this%nLoc)%jmax=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%kmax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(this%nLoc)%kmax=int(buf)
        print*,'add3, i=',this%nLoc
        ! Figure out how many processors are within this station
        call MPI_ALLREDUCE(nproc_stat,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); nproc_stat=buf
        print*,'add3a, i=',this%nLoc
        allocate(new_ranks(nproc_stat))
        print*,'add3b, i=',this%nLoc

        ! Figure out our group and communicator
        call MPI_GROUP_SIZE(this%cfg%group,nproc,ierr)
        print*,'add3c, i=',this%nLoc

        call MPI_Group_rank(this%cfg%group,myRank,ierr)
        print*,'add4, i=',this%nLoc
        do i=0,nproc-1
            if (myRank.eq.i) then
                if (this%stats(this%nLoc)%amIn) then
                    key = i
                    color = 0
                else
                    key = 0 
                    color = MPI_UNDEFINED
                end if
            end if
        end do
        call MPI_Comm_split(this%cfg%comm,color,key,this%stats(this%nLoc)%comm)
        ! Get new rank and group
        this%stats(this%nLoc)%amRoot = .false.
        if (this%stats(this%nLoc)%amIn) then
            call MPI_Comm_rank(this%stats(this%nLoc)%comm,myRank,ierr)
            if (myRank.eq.0) this%stats(this%nLoc)%amRoot = .true.
        end if
        
        deallocate(new_ranks)
        print*,'add5, i=',this%nLoc
    end subroutine add_station

    ! Subroutine for initializing stats arrays based on the stations which have already been initialized
    subroutine init_stats(this,lm,spns,lp,T_sc,Yf_sc)
        use lowmach_class, only :  lowmach
        use incomp_class, only :   incomp
        use lpt_class, only :      lpt
        use vdscalar_class, only : vdscalar
        implicit none
        class(stats_parent),intent(inout) :: this
        type(lowmach), optional, intent(in) :: lm
        type(incomp), optional, intent(in) :: spns
        type(lpt), optional, intent(in) :: lp
        type(vdscalar), optional, intent(in) :: T_sc,Yf_sc
        integer :: i,j,k,nVar,n,iLoc

        ! Assume all stations gather the same statistics!!!
        nVar = 0
        this%sum_time = 0.0_WP
        this%my_solvers = .false.
        if (present(lm)) then
            this%my_solvers(1)=.true.
            nVar=nVar+14
        end if
        if (present(T_sc)) then
            this%my_solvers(2)=.true.
            nVar=nVar+4
        end if
        if (present(Yf_sc)) then
            this%my_solvers(3)=.true.
            nVar=nVar+4
        end if
        if (present(spns)) then
            this%my_solvers(4)=.true.
            nVar=nVar+6
        end if
        if (present(lp)) then
            this%my_solvers(5)=.true.
            nVar=nVar+12
        end if

        ! Loop over all stations, only entering if the processor is in the station, then allocate all the things
        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                this%stats(iLoc)%n_stats = nVar
                allocate(this%stats(iLoc)%names(nVar))
                ! Allocate arrays based on dim
                select case (this%stats(iLoc)%dim)
                case (0); allocate(this%stats(iLoc)%vals(nVar,1,1,1))
                case (1); allocate(this%stats(iLoc)%vals(nVar,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,1,1))
                case (2); allocate(this%stats(iLoc)%vals(nVar,1,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,1))
                case (3); allocate(this%stats(iLoc)%vals(nVar,1,1,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                case (4); allocate(this%stats(iLoc)%vals(nVar,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,1))
                case (5); allocate(this%stats(iLoc)%vals(nVar,1,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                case (6); allocate(this%stats(iLoc)%vals(nVar,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,1,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                case (7); allocate(this%stats(iLoc)%vals(nVar,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                end select
                this%stats(iLoc)%vals = 0.0_WP

                n = 0
                if (present(lm)) then
                    this%stats(iLoc)%names(n+1) = 'rho'
                    this%stats(iLoc)%names(n+2) = 'rho^2'
                    this%stats(iLoc)%names(n+3) = 'U'
                    this%stats(iLoc)%names(n+4) = 'V'
                    this%stats(iLoc)%names(n+5) = 'W'
                    this%stats(iLoc)%names(n+6) = 'U^2'
                    this%stats(iLoc)%names(n+7) = 'V^2'
                    this%stats(iLoc)%names(n+8) = 'W^2'
                    this%stats(iLoc)%names(n+9) = 'rhoU'
                    this%stats(iLoc)%names(n+10) = 'rhoV'
                    this%stats(iLoc)%names(n+11) = 'rhoW'
                    this%stats(iLoc)%names(n+12) = 'rhoU^2'
                    this%stats(iLoc)%names(n+13) = 'rhoV^2'
                    this%stats(iLoc)%names(n+14) = 'rhoW^2'
                    n = n+14
                end if

                if (present(T_sc)) then
                    this%stats(iLoc)%names(n+1) = 'T'
                    this%stats(iLoc)%names(n+2) = 'T^2'
                    this%stats(iLoc)%names(n+3) = 'rhoT'
                    this%stats(iLoc)%names(n+4) = 'rhoT^2'
                    n = n+4
                end if

                if (present(Yf_sc)) then
                    this%stats(iLoc)%names(n+1) = 'Yf'
                    this%stats(iLoc)%names(n+2) = 'Yf^2'
                    this%stats(iLoc)%names(n+3) = 'rhoYf'
                    this%stats(iLoc)%names(n+4) = 'rhoYf^2'
                    n = n+4
                end if

                if (present(spns)) then
                    this%stats(iLoc)%names(n+1) = 'U'
                    this%stats(iLoc)%names(n+2) = 'V'
                    this%stats(iLoc)%names(n+3) = 'W'
                    this%stats(iLoc)%names(n+4) = 'U^2'
                    this%stats(iLoc)%names(n+5) = 'V^2'
                    this%stats(iLoc)%names(n+6) = 'W^2'
                    n = n+6
                end if

                if (present(lp)) then
                    this%stats(iLoc)%names(n+1) = 'lp_number'
                    this%stats(iLoc)%names(n+2) = 'lp_D'
                    this%stats(iLoc)%names(n+3) = 'lp_D^2'
                    this%stats(iLoc)%names(n+4) = 'lp_D^3'
                    this%stats(iLoc)%names(n+5) = 'lp_T'
                    this%stats(iLoc)%names(n+6) = 'lp_T^2'
                    this%stats(iLoc)%names(n+7) = 'lp_U'
                    this%stats(iLoc)%names(n+8) = 'lp_V'
                    this%stats(iLoc)%names(n+9) = 'lp_W'
                    this%stats(iLoc)%names(n+10) = 'lp_U^2'
                    this%stats(iLoc)%names(n+11) = 'lp_V^2'
                    this%stats(iLoc)%names(n+12) = 'lp_W^2'
                    n = n+12
                end if

                
            end if
        end do
        print*,'Successfully Initialized',this%nLoc,'Stations for Statistics'
    end subroutine init_stats


    ! Subroutine for collecting statistics
    subroutine sample_stats(this,lm,spns,lp,T_sc,Yf_sc,dt)
        use lowmach_class, only :  lowmach
        use incomp_class, only :   incomp
        use lpt_class, only :      lpt
        use vdscalar_class, only : vdscalar
        use messager, only :       die
        implicit none
        class(stats_parent),intent(inout) :: this
        type(lowmach), optional, intent(in) :: lm
        type(incomp), optional, intent(in) :: spns
        type(lpt), optional, intent(in) :: lp
        type(vdscalar), optional, intent(in) :: T_sc,Yf_sc
        integer :: i,j,k,n,iLoc,iS,jS,kS
        real(WP), intent(in) :: dt
        real(WP) :: corr ! Corrects for sampling volume when stat cell is larger than mesh cell
        integer, dimension(3,2) :: indMod ! Modify indices for different station dimensions


        this%sum_time = this%sum_time+dt
        ! Loop over all stations, only entering if the processor is in the station, then gather stats
        ! print*,'sample0'
        do iLoc=1,this%nLoc
            ! print*,'sample1'
            if (this%stats(iLoc)%amIn) then
                n = 0
                ! Set our modifiers by the dimensionality of this station
                select case (this%stats(iLoc)%dim)
                case(0) ! 0D
                    corr = real((1+this%stats(iLoc)%imax_-this%stats(iLoc)%imin_) * &
                                (1+this%stats(iLoc)%jmax_-this%stats(iLoc)%jmin_) * &
                                (1+this%stats(iLoc)%kmax_-this%stats(iLoc)%kmin_),WP)
                    indMod(:,1) = [1,1,1] ! First column is a static offset
                    indMod(:,2) = [0,0,0] ! Second column multiplies into ijk
                case(1) ! X
                    corr = real((1+this%stats(iLoc)%jmax_-this%stats(iLoc)%jmin_) * &
                                (1+this%stats(iLoc)%kmax_-this%stats(iLoc)%kmin_),WP)
                    indMod(:,1) = [0,1,1] ! First column is a static offset
                    indMod(:,2) = [1,0,0] ! Second column multiplies into ijk
                case(2) ! Y
                    corr = real((1+this%stats(iLoc)%imax_-this%stats(iLoc)%imin_) * &
                                (1+this%stats(iLoc)%kmax_-this%stats(iLoc)%kmin_),WP)
                    indMod(:,1) = [1,0,1] ! First column is a static offset
                    indMod(:,2) = [0,1,0] ! Second column multiplies into ijk
                case(3) ! Z
                    corr = real((1+this%stats(iLoc)%imax_-this%stats(iLoc)%imin_) * &
                                (1+this%stats(iLoc)%jmax_-this%stats(iLoc)%jmin_),WP)
                    indMod(:,1) = [1,1,0] ! First column is a static offset
                    indMod(:,2) = [0,0,1] ! Second column multiplies into ijk
                case(4) ! XY
                    corr = real(1+this%stats(iLoc)%kmax_-this%stats(iLoc)%kmin_,WP)
                    indMod(:,1) = [0,0,1] ! First column is a static offset
                    indMod(:,2) = [1,1,0] ! Second column multiplies into ijk
                case(5) ! YZ
                    corr = real(1+this%stats(iLoc)%imax_-this%stats(iLoc)%imin_,WP)
                    indMod(:,1) = [1,0,0] ! First column is a static offset
                    indMod(:,2) = [0,1,1] ! Second column multiplies into ijk
                case(6) ! ZX
                    corr = real(1+this%stats(iLoc)%jmax_-this%stats(iLoc)%jmin_,WP)
                    indMod(:,1) = [0,1,0] ! First column is a static offset
                    indMod(:,2) = [1,0,1] ! Second column multiplies into ijk
                case(7) ! XYZ
                    corr = 1.0_WP
                    indMod(:,1) = [0,0,0] ! First column is a static offset
                    indMod(:,2) = [1,1,1] ! Second column multiplies into ijk
                end select


                do i=this%stats(iLoc)%imin_,this%stats(iLoc)%imax_
                    do j=this%stats(iLoc)%jmin_,this%stats(iLoc)%jmax_
                        do k=this%stats(iLoc)%kmin_,this%stats(iLoc)%kmax_
                            iS = sum(indMod(1,:)*[1,i])
                            jS = sum(indMod(2,:)*[1,j])
                            kS = sum(indMod(3,:)*[1,k])
                            corr =  sum([this%cfg%dx(i),this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin)]*real(  indMod(1,:),WP))/ &
                                    sum([this%cfg%dx(i),this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin)]*real(1-indMod(1,:),WP))* &
                                    sum([this%cfg%dy(j),this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin)]*real(  indMod(2,:),WP))/ &
                                    sum([this%cfg%dy(j),this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin)]*real(1-indMod(2,:),WP))* &
                                    sum([this%cfg%dz(k),this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin)]*real(  indMod(3,:),WP))/ &
                                    sum([this%cfg%dz(k),this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin)]*real(1-indMod(3,:),WP))
                            ! corr =  (this%cfg%dx(i)*indMod(1,1)+,real(this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin),WP)]*real(  indMod(1,:),WP))/
                            !         sum([this%cfg%dx(i),real(this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin),WP)]*real(1-indMod(1,:),WP))*
                            !         sum([this%cfg%dy(j),real(this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin),WP)]*real(  indMod(2,:),WP))/
                            !         sum([this%cfg%dy(j),real(this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin),WP)]*real(1-indMod(2,:),WP))*
                            !         sum([this%cfg%dz(k),real(this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin),WP)]*real(  indMod(3,:),WP))/
                            !         sum([this%cfg%dz(k),real(this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin),WP)]*real(1-indMod(3,:),WP))*
                            if (present(lm)) then
                                this%stats(iLoc)%vals(n+1,iS,jS,kS) =  this%stats(iLoc)%vals(n+1,iS,jS,kS)+ dt*corr*lm%rho(i,j,k)     !'rho'
                                this%stats(iLoc)%vals(n+2,iS,jS,kS) =  this%stats(iLoc)%vals(n+2,iS,jS,kS)+ dt*corr*lm%rho(i,j,k)**2  !'rho^2'
                                this%stats(iLoc)%vals(n+3,iS,jS,kS) =  this%stats(iLoc)%vals(n+3,iS,jS,kS)+ dt*corr*lm%U(i,j,k)       !'U'
                                this%stats(iLoc)%vals(n+4,iS,jS,kS) =  this%stats(iLoc)%vals(n+4,iS,jS,kS)+ dt*corr*lm%V(i,j,k)       !'V'
                                this%stats(iLoc)%vals(n+5,iS,jS,kS) =  this%stats(iLoc)%vals(n+5,iS,jS,kS)+ dt*corr*lm%W(i,j,k)       !'W'
                                this%stats(iLoc)%vals(n+6,iS,jS,kS) =  this%stats(iLoc)%vals(n+6,iS,jS,kS)+ dt*corr*lm%U(i,j,k)**2    !'U^2'
                                this%stats(iLoc)%vals(n+7,iS,jS,kS) =  this%stats(iLoc)%vals(n+7,iS,jS,kS)+ dt*corr*lm%V(i,j,k)**2    !'V^2'
                                this%stats(iLoc)%vals(n+8,iS,jS,kS) =  this%stats(iLoc)%vals(n+8,iS,jS,kS)+ dt*corr*lm%W(i,j,k)**2    !'W^2'
                                this%stats(iLoc)%vals(n+9,iS,jS,kS) =  this%stats(iLoc)%vals(n+9,iS,jS,kS)+ dt*corr*lm%rhoU(i,j,k)    !'rhoU'
                                this%stats(iLoc)%vals(n+10,iS,jS,kS) = this%stats(iLoc)%vals(n+10,iS,jS,kS)+dt*corr*lm%rhoV(i,j,k)    !'rhoV'
                                this%stats(iLoc)%vals(n+11,iS,jS,kS) = this%stats(iLoc)%vals(n+11,iS,jS,kS)+dt*corr*lm%rhoW(i,j,k)    !'rhoW'
                                this%stats(iLoc)%vals(n+12,iS,jS,kS) = this%stats(iLoc)%vals(n+12,iS,jS,kS)+dt*corr*lm%rhoU(i,j,k)**2 !'rhoU^2'
                                this%stats(iLoc)%vals(n+13,iS,jS,kS) = this%stats(iLoc)%vals(n+13,iS,jS,kS)+dt*corr*lm%rhoV(i,j,k)**2 !'rhoV^2'
                                this%stats(iLoc)%vals(n+14,iS,jS,kS) = this%stats(iLoc)%vals(n+14,iS,jS,kS)+dt*corr*lm%rhoW(i,j,k)**2 !'rhoW^2'
                                n = n+14
                            end if

                            if (present(T_sc)) then
                                this%stats(iLoc)%vals(n+1,iS,jS,kS) = this%stats(iLoc)%vals(n+1,iS,jS,kS)+dt*corr*T_sc%SC(i,j,k)       !'T'
                                this%stats(iLoc)%vals(n+2,iS,jS,kS) = this%stats(iLoc)%vals(n+2,iS,jS,kS)+dt*corr*T_sc%SC(i,j,k)**2    !'T^2'
                                this%stats(iLoc)%vals(n+3,iS,jS,kS) = this%stats(iLoc)%vals(n+3,iS,jS,kS)+dt*corr*T_sc%rhoSC(i,j,k)    !'rhoT'
                                this%stats(iLoc)%vals(n+4,iS,jS,kS) = this%stats(iLoc)%vals(n+4,iS,jS,kS)+dt*corr*T_sc%rhoSC(i,j,k)**2 !'rhoT^2'
                                n = n+4
                            end if

                            if (present(Yf_sc)) then
                                this%stats(iLoc)%vals(n+1,iS,jS,kS) = this%stats(iLoc)%vals(n+1,iS,jS,kS)+dt*corr*Yf_sc%SC(i,j,k)       !'Yf'
                                this%stats(iLoc)%vals(n+2,iS,jS,kS) = this%stats(iLoc)%vals(n+2,iS,jS,kS)+dt*corr*Yf_sc%SC(i,j,k)**2    !'Yf^2'
                                this%stats(iLoc)%vals(n+3,iS,jS,kS) = this%stats(iLoc)%vals(n+3,iS,jS,kS)+dt*corr*Yf_sc%rhoSC(i,j,k)    !'rhoYf'
                                this%stats(iLoc)%vals(n+4,iS,jS,kS) = this%stats(iLoc)%vals(n+4,iS,jS,kS)+dt*corr*Yf_sc%rhoSC(i,j,k)**2 !'rhoYf^2'
                                n = n+4
                            end if

                            if (present(spns)) then
                                this%stats(iLoc)%vals(n+1,iS,jS,kS) = this%stats(iLoc)%vals(n+1,iS,jS,kS)+dt*corr*spns%U(i,j,k)    !'U'
                                this%stats(iLoc)%vals(n+2,iS,jS,kS) = this%stats(iLoc)%vals(n+2,iS,jS,kS)+dt*corr*spns%V(i,j,k)    !'V'
                                this%stats(iLoc)%vals(n+3,iS,jS,kS) = this%stats(iLoc)%vals(n+3,iS,jS,kS)+dt*corr*spns%W(i,j,k)    !'W'
                                this%stats(iLoc)%vals(n+4,iS,jS,kS) = this%stats(iLoc)%vals(n+4,iS,jS,kS)+dt*corr*spns%U(i,j,k)**2 !'U^2'
                                this%stats(iLoc)%vals(n+5,iS,jS,kS) = this%stats(iLoc)%vals(n+5,iS,jS,kS)+dt*corr*spns%V(i,j,k)**2 !'V^2'
                                this%stats(iLoc)%vals(n+6,iS,jS,kS) = this%stats(iLoc)%vals(n+6,iS,jS,kS)+dt*corr*spns%W(i,j,k)**2 !'W^2'
                                n = n+6
                            end if
                        end do
                    end do
                end do
                ! print*,'sample4'
                if (present(lp)) then
                    do i=1,lp%np_
                        if ((lp%p(i)%ind(1).ge.this%stats(iLoc)%imin_.and.lp%p(i)%ind(1).le.this%stats(iLoc)%imax_).and. &
                            (lp%p(i)%ind(2).ge.this%stats(iLoc)%jmin_.and.lp%p(i)%ind(2).le.this%stats(iLoc)%jmax_).and. &
                            (lp%p(i)%ind(3).ge.this%stats(iLoc)%kmin_.and.lp%p(i)%ind(3).le.this%stats(iLoc)%kmax_).and. &
                                lp%p(i)%flag.eq.0) then ! If particle is inside our station
                            iS = sum(indMod(1,:)*[1,lp%p(i)%ind(1)])
                            jS = sum(indMod(2,:)*[1,lp%p(i)%ind(2)])
                            kS = sum(indMod(3,:)*[1,lp%p(i)%ind(3)])
                                ! print*,'made it into saving stats, i=',i
                            this%stats(iLoc)%vals(n+1,iS,jS,kS) =  this%stats(iLoc)%vals(n+1,iS,jS,kS)+ dt                   !'lp_number'
                            this%stats(iLoc)%vals(n+2,iS,jS,kS) =  this%stats(iLoc)%vals(n+2,iS,jS,kS)+ dt*lp%p(i)%d         !'lp_D'
                            this%stats(iLoc)%vals(n+3,iS,jS,kS) =  this%stats(iLoc)%vals(n+3,iS,jS,kS)+ dt*lp%p(i)%d**2      !'lp_D^2'
                            this%stats(iLoc)%vals(n+4,iS,jS,kS) =  this%stats(iLoc)%vals(n+4,iS,jS,kS)+ dt*lp%p(i)%d**3      !'lp_D^3'
                            this%stats(iLoc)%vals(n+5,iS,jS,kS) =  this%stats(iLoc)%vals(n+5,iS,jS,kS)+ dt*lp%p(i)%T_d       !'lp_T'
                            this%stats(iLoc)%vals(n+6,iS,jS,kS) =  this%stats(iLoc)%vals(n+6,iS,jS,kS)+ dt*lp%p(i)%T_d**2    !'lp_T^2'
                            this%stats(iLoc)%vals(n+7,iS,jS,kS) =  this%stats(iLoc)%vals(n+7,iS,jS,kS)+ dt*lp%p(i)%vel(1)    !'lp_U'
                            this%stats(iLoc)%vals(n+8,iS,jS,kS) =  this%stats(iLoc)%vals(n+8,iS,jS,kS)+ dt*lp%p(i)%vel(2)    !'lp_V'
                            this%stats(iLoc)%vals(n+9,iS,jS,kS) =  this%stats(iLoc)%vals(n+9,iS,jS,kS)+ dt*lp%p(i)%vel(3)    !'lp_W'
                            this%stats(iLoc)%vals(n+10,iS,jS,kS) = this%stats(iLoc)%vals(n+10,iS,jS,kS)+dt*lp%p(i)%vel(1)**2 !'lp_U^2'
                            this%stats(iLoc)%vals(n+11,iS,jS,kS) = this%stats(iLoc)%vals(n+11,iS,jS,kS)+dt*lp%p(i)%vel(2)**2 !'lp_V^2'
                            this%stats(iLoc)%vals(n+12,iS,jS,kS) = this%stats(iLoc)%vals(n+12,iS,jS,kS)+dt*lp%p(i)%vel(3)**2 !'lp_W^2'
                        end if
                    end do
                    n = n+12
                end if

            end if
        end do
        if (this%cfg%amRoot) print*,'Statistics collected :)'
    end subroutine sample_stats

    ! Normalize collected statistics over time
    subroutine post_process(this)
        use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
        use parallel, only: MPI_REAL_WP
        class(stats_parent), intent(inout) :: this
        integer :: n,i,iLoc,ierr,numVals
        real(WP), dimension(:,:,:,:), allocatable :: buf,globVals

        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                select case (this%stats(iLoc)%dim)
                case (0) ! 0D
                    allocate(buf(this%stats(iLoc)%n_stats,1,1,1))
                    allocate(globVals(this%stats(iLoc)%n_stats,1,1,1)); globVals = 0.0_WP
                    globVals = this%stats(iLoc)%vals
                case (1) ! X
                    allocate(buf(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,1))
                    allocate(globVals(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,1)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,:,:) = this%stats(iLoc)%vals
                case (2) ! Y
                    allocate(buf(this%stats(iLoc)%n_stats,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1))
                    allocate(globVals(this%stats(iLoc)%n_stats,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1)); globVals = 0.0_WP
                    globVals(:,:,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,:) = this%stats(iLoc)%vals
                case (3) ! Z
                    allocate(buf(this%stats(iLoc)%n_stats,1,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%stats(iLoc)%n_stats,1,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,:,:,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                case (4) ! XY
                    allocate(buf(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1))
                    allocate(globVals(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,:) = this%stats(iLoc)%vals
                case (5) ! YZ
                    allocate(buf(this%stats(iLoc)%n_stats,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%stats(iLoc)%n_stats,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,:,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                case (6) ! ZX
                    allocate(buf(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,:,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                case (7) ! XYZ
                    allocate(buf(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%stats(iLoc)%n_stats,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                end select
                numVals = size(globVals,1)*size(globVals,2)*size(globVals,3)*size(globVals,4)
                call MPI_ALLREDUCE(globVals,buf,numVals,MPI_REAL_WP,MPI_SUM,this%stats(iLoc)%comm,ierr); globVals = buf
                deallocate(buf)
                n = 0



                if (this%my_solvers(1)) then
                    globVals(n+1:n+14,:,:,:) = globVals(n+1:n+14,:,:,:)/this%sum_time
                    n = n+14
                end if

                if (this%my_solvers(2)) then
                    globVals(n+1:n+4,:,:,:) = globVals(n+1:n+4,:,:,:)/this%sum_time
                    n = n+4
                end if

                if (this%my_solvers(3)) then
                    globVals(n+1:n+4,:,:,:) = globVals(n+1:n+4,:,:,:)/this%sum_time
                    n = n+4
                end if

                if (this%my_solvers(4)) then
                    globVals(n+1:n+6,:,:,:) = globVals(n+1:n+6,:,:,:)/this%sum_time
                    n = n+6
                end if
                ! Add all reduces
                if (this%my_solvers(5)) then
                    do i=n+2,n+12
                        where (globVals(n+1,:,:,:).gt.0.0_WP) ! Only normalize where we found particles
                            globVals(i,:,:,:) = globVals(i,:,:,:)/globVals(n+1,:,:,:)
                        elsewhere
                            globVals(i,:,:,:) = 0.0_WP
                        end where
                    end do
                    globVals(n+1,:,:,:) = globVals(n+1,:,:,:)/this%sum_time
                    n = n+12
                end if
                select case (this%stats(iLoc)%dim)
                case (0) ! 0D
                    this%stats(iLoc)%vals = globVals
                case (1) ! X
                    this%stats(iLoc)%vals = globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,:,:)
                case (2) ! Y
                    this%stats(iLoc)%vals = globVals(:,:,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,:)
                case (3) ! Z
                    this%stats(iLoc)%vals = globVals(:,:,:,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_)
                case (4) ! XY
                    this%stats(iLoc)%vals = globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,:)
                case (5) ! YZ
                    this%stats(iLoc)%vals = globVals(:,:,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_)
                case (6) ! ZX
                    this%stats(iLoc)%vals = globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,:,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_)
                case (7) ! XYZ
                    this%stats(iLoc)%vals = globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_)
                end select
                deallocate(globVals)
            end if
        end do
    end subroutine post_process

    ! Store current state of statistics to the disk for retrieval later. Doesn't affect current arrays
    subroutine cache_data(this)
        class(stats_parent), intent(inout) :: this

    end subroutine cache_data

    ! Read the data contained in a prior cache
    subroutine read_data(this)
        class(stats_parent), intent(inout) :: this

    end subroutine read_data


    subroutine reset_stats(this)
        class(stats_parent), intent(inout) :: this
        integer :: iLoc

        this%sum_time = 0.0_WP
        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                this%stats(iLoc)%vals=0.0_WP
            end if
        end do

    end subroutine reset_stats



end module stats_class