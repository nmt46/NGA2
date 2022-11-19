! Collect specialized statistics for sprays
module stats_class
    use precision, only: WP
    use string,    only: str_medium
    use config_class, only: config
    use mpi_f08
    use lpt_class, only: lpt
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
        ! character(len=str_medium), dimension(:), allocatable :: names       ! Names associated with gather statistics
        real(WP), dimension(:,:,:,:), allocatable :: vals                   ! Values for each statistic (before post-processing)
        type(MPI_Comm) :: comm = MPI_COMM_NULL                              ! Communicator for processors on this station
        logical :: amRoot
    end type station

    ! Array of pointers to the data fed into the statistics and associated names
    type :: arr_point
        real(WP), dimension(:,:,:), pointer :: val                  ! Pointer to this array
        character(len=str_medium) :: name = 'Unknown_arr_point'     ! Name of this array
    end type arr_point

    ! Contains the parameters for creating one statistical measure 
    type :: stat_def
        integer :: nTA                                  ! Number of arrays in arp used for this statistic
        integer,  dimension(:), allocatable :: iAs      ! Ind(ex)(ices) for each array within arp
        integer,  dimension(6) :: expP                  ! Exponents corresponding to number,lp%p%d,T_d,U,V,W
        character(len=str_medium) :: name               ! Name of this definition

    end type stat_def

    ! Parent class for all the statistics related things
    type :: stats_parent
        type(config), pointer :: cfg                                ! The config around which spray_parent is built
        character(len=str_medium) :: name='Unknown_stats_parent'    ! Name of our stat parent
        integer :: nLoc                                             ! Number of stations for gathering statistics
        type(station), dimension(:),allocatable :: stats            ! Array of our stations
        real(WP) :: sum_time                                        ! Amount of time over which we have averaged
        type(arr_point), dimension(:), allocatable :: arp           ! Array of pointers to the arrays from whch statistics are gathered
        integer :: nAr                                              ! Number of array pointers we have
        type(stat_def), dimension(:), allocatable :: def            ! Array of the definitions for each statistical measure
        integer :: nDef                                             ! Number of stats definitions
        integer :: iNum = -1                                        ! Index in def which corresponds to particle number. If -1, none are
        type(lpt), pointer :: lp                                    ! Pointer to LPT solver we use for particles
        logical, dimension(5) :: my_solvers                         ! List of whether we use each solver by logical value in each index
                                                                        ! 1=lowmach,2=T_sc,3=Yf_sc,4=incomp,5=lp

    contains
        procedure :: add_station            ! Add a station for stats measuring
        procedure :: add_array              ! Add a pointer to a solver array
        procedure :: add_definition         ! Add the definition for one statistical measure. To be called after adding all arrays
        procedure :: set_lp                 ! Set up the lp pointer
        procedure :: init_stats             ! Initialize statistics. To be called after adding all stations, arrays, and defitions
        procedure :: sample_stats           ! Sample stats over all stations
        procedure :: cache_data             ! Store data to the disk
        procedure :: read_data              ! Read raw data from disk
        procedure :: post_process           ! Process data from binary file into usable forms
        procedure :: reset_stats            ! Reset statistics, maintaining stations
        procedure :: get_def                ! Get the index of a stat definition
        

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

        ! Initialize various counters
        self%nLoc = 0
        self%nAr  = 0
        self%nDef = 0

        ! Default lp pointer to null
        self%lp => NULL()
    end function constructor

    ! Function for splitting an statistic definition with * as delimiter
    function parse_def(def) result(parsed_def)
        character(len=*) :: def
        integer :: i,i1,i2,n,nT
        character(len=str_medium), dimension(:), allocatable :: parsed_def
        n = len(trim(def))
        nT = 0
        do i=1,n
            if (def(i:i).eq.'*') nT = nT+1
        end do
        allocate(parsed_def(nT+1))
        i1 = 1
        i2 = 1
        do i=1,n
            if (def(i:i).eq.'*') then
                parsed_def(i2) = def(i1:i-1)
                i1 = i+1
                i2 = i2+1
            end if
        end do
        parsed_def(i2) = def(i1:len(trim(def)))
    end function parse_def

    ! Convert integer to character array
    function int2char(myInt) result(myChar)
        integer, intent(in) :: myInt
        character(len=64) :: myChar
        integer :: i,n,buf
        i = 0
        n = len(myChar)
        buf = myInt
        do while (buf.gt.0)
            myChar(n-i:n-i) = char(48+mod(buf,10))
            buf = buf/10
            i = i+1
        end do
        myChar = trim(myChar)
    end function int2char
    
    ! Given a stats definition name, return the index of that definition
    function get_def(this,name) result(iDef)
        use messager, only: die
        use string, only: lowercase
        class(stats_parent),intent(in) :: this
        character(len=*), intent(in) :: name
        integer :: iDef
        integer :: i
        iDef = -1
        do i=1,this%nDef
            if (lowercase(trim(this%def(i)%name)).eq.lowercase(trim(name))) iDef = i
        end do
        if (iDef.eq.-1) call die('[stats_class get_def] Unknown statistics definition was passed.')
    end function get_def

    ! Set lp, initialize statistic definition for particle number (required for all other particle stats)
    subroutine set_lp(this,lp)
        implicit none
        class(stats_parent), intent(inout) :: this
        type(lpt), intent(in), target :: lp

        this%lp => lp
        call this%add_definition(name='p_n',def='p_n')
    end subroutine set_lp

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
        nproc_stat = 0.0_WP
        if (this%stats(this%nLoc)%amIn) nproc_stat = 1.0_WP
        ! Get the global bounding indices
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%imin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(this%nLoc)%imin=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%jmin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(this%nLoc)%jmin=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%kmin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(this%nLoc)%kmin=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%imax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(this%nLoc)%imax=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%jmax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(this%nLoc)%jmax=int(buf)
        call MPI_ALLREDUCE(real(this%stats(this%nLoc)%kmax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(this%nLoc)%kmax=int(buf)
        ! Figure out how many processors are within this station
        call MPI_ALLREDUCE(nproc_stat,buf,1,MPI_REAL_WP,MPI_SUM,this%cfg%comm,ierr); nproc_stat=buf
        allocate(new_ranks(int(nproc_stat)))

        ! Figure out our group and communicator
        call MPI_GROUP_SIZE(this%cfg%group,nproc,ierr)

        call MPI_Group_rank(this%cfg%group,myRank,ierr)
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
    end subroutine add_station

    ! Add array to use when collecting statistics. Must use overlapping pgrid dimensions
    subroutine add_array(this,array,name)
        implicit none
        class(stats_parent),intent(inout) :: this
        type(arr_point), dimension(:), allocatable :: tmp
        real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: array
        character(len=*),intent(in) :: name

        if (allocated(this%arp)) then
            allocate(tmp(this%nAr+1))
            tmp(1:this%nAr)=this%arp
            call move_alloc(tmp,this%arp)
        else
            allocate(this%arp(1))
        end if
        this%nAr = this%nAr+1

        this%arp(this%nAr)%name = trim(adjustl(name))
        this%arp(this%nAr)%val => array

    end subroutine add_array

    subroutine add_definition(this,name,def)
        use messager, only : die
        use string, only : lowercase
        implicit none
        class(stats_parent), intent(inout) :: this
        character(len=*),intent(in) :: name
        type(stat_def), dimension(:), allocatable :: tmp
        character(len=*), intent(in) :: def
        character(len=str_medium),dimension(:),allocatable :: def_parsed
        character(len=str_medium),dimension(6) :: p_names
        integer :: i,j,n,iAr,iP

        ! Allocate the next definition
        if (allocated(this%def)) then
            allocate(tmp(this%nDef+1))
            tmp(1:this%nDef)=this%def
            call move_alloc(tmp,this%def)
        else
            allocate(this%def(1))
        end if
        ! Increment definition counter
        this%nDef = this%nDef+1
        
        ! Assign a name
        this%def(this%nDef)%name = name

        ! Zero-out some things
        this%def(this%nDef)%nTA = 0
        this%def(this%nDef)%expP = 0

        ! Make string array of the names of terms for particles
        p_names(1) = 'p_n'; p_names(2) = 'p_d'; p_names(3) = 'p_T'
        p_names(4) = 'p_U'; p_names(5) = 'p_V'; p_names(6) = 'p_W'

        ! Get a string array of the components of the definition
        def_parsed = parse_def(def)
        do i = 1,size(def_parsed)
            do iAr = 1,this%nAr
                ! Need to count number of arrays to allocate iAs
                if (lowercase(trim(def_parsed(i))).eq.lowercase(trim(this%arp(iAr)%name))) this%def(this%nDef)%nTA = this%def(this%nDef)%nTA+1
            end do
            do iP = 1,6
                if (lowercase(trim(def_parsed(i))).eq.lowercase(trim(p_names(iP)))) then
                    this%def(this%nDef)%expP(iP) = this%def(this%nDef)%expP(iP)+1
                    if (iP.eq.1) then
                        this%iNum = this%nDef
                    end if
                end if
            end do
        end do
        if (this%def(this%nDef)%nTA.gt.0) then
            allocate(this%def(this%nDef)%iAs(this%def(this%nDef)%nTA))
            j = 1
            do i = 1,size(def_parsed)
                do iAr = 1,this%nAr
                    ! Set this value of iAr to iAs
                    if (lowercase(trim(def_parsed(i))).eq.lowercase(trim(this%arp(iAr)%name))) then
                        this%def(this%nDef)%iAs(j) = iAr
                        j = j+1
                    end if
                end do
            end do
        else
            allocate(this%def(this%nDef)%iAs(1))
            this%def(this%nDef)%iAs = 0
        end if
        
        ! Check that user has given us a particle solver if they want particle statistics
        if(sum(abs(this%def(this%nDef)%expP)).ne.0.and..not.associated(this%lp)) call die('[stats_class add_definition] No LP solver is associated with stat_parent but requested statistic uses LP')
        ! if(lowercase(trim(name)).eq.'u') print*,this%def(this%nDef)%expP
    end subroutine add_definition

    

    ! Subroutine for initializing stats arrays based on the stations which have already been initialized
    subroutine init_stats(this)
        implicit none
        class(stats_parent),intent(inout) :: this
        integer :: i,j,k,n,iLoc

        this%sum_time = 0.0_WP
        
        ! Loop over all stations, only entering if the processor is in the station, then allocate all the things
        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                ! Allocate arrays based on dim
                select case (this%stats(iLoc)%dim)
                case (0); allocate(this%stats(iLoc)%vals(this%nDef,1,1,1))
                case (1); allocate(this%stats(iLoc)%vals(this%nDef,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,1,1))
                case (2); allocate(this%stats(iLoc)%vals(this%nDef,1,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,1))
                case (3); allocate(this%stats(iLoc)%vals(this%nDef,1,1,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                case (4); allocate(this%stats(iLoc)%vals(this%nDef,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,1))
                case (5); allocate(this%stats(iLoc)%vals(this%nDef,1,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                case (6); allocate(this%stats(iLoc)%vals(this%nDef,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,1,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                case (7); allocate(this%stats(iLoc)%vals(this%nDef,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                end select
                this%stats(iLoc)%vals = 0.0_WP
            end if
        end do
        print*,'Successfully Initialized',this%nLoc,'Stations for Statistics'
    end subroutine init_stats


    ! Subroutine for collecting statistics
    subroutine sample_stats(this,dt)
        use messager, only :       die
        implicit none
        class(stats_parent),intent(inout) :: this
        real(WP), intent(in) :: dt
        integer :: i,j,k,iDef,iLoc,iS,jS,kS,iA
        real(WP) :: corr ! Corrects for sampling volume when stat cell is larger than mesh cell
        integer, dimension(3,2) :: indMod ! Modify indices for different station dimensions
        real(WP), dimension(:,:,:), allocatable:: buf,buf2 ! Buffer for building stats


        this%sum_time = this%sum_time+dt
        ! Loop over all stations, only entering if the processor is in the station, then gather stats
        ! print*,'sample0'
        do iLoc=1,this%nLoc
            ! print*,'sample1'
            if (this%stats(iLoc)%amIn) then
                ! Set our modifiers by the dimensionality of this station
                select case (this%stats(iLoc)%dim)
                case(0) ! 0D
                    allocate(buf(1,1,1))
                    indMod(:,1) = [1,1,1] ! First column is a static offset
                    indMod(:,2) = [0,0,0] ! Second column multiplies into ijk
                case(1) ! X
                    allocate(buf(this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,1,1))
                    indMod(:,1) = [0,1,1] ! First column is a static offset
                    indMod(:,2) = [1,0,0] ! Second column multiplies into ijk
                case(2) ! Y
                    allocate(buf(1,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,1))
                    indMod(:,1) = [1,0,1] ! First column is a static offset
                    indMod(:,2) = [0,1,0] ! Second column multiplies into ijk
                case(3) ! Z
                    allocate(buf(1,1,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                    indMod(:,1) = [1,1,0] ! First column is a static offset
                    indMod(:,2) = [0,0,1] ! Second column multiplies into ijk
                case(4) ! XY
                    allocate(buf(this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,1))
                    indMod(:,1) = [0,0,1] ! First column is a static offset
                    indMod(:,2) = [1,1,0] ! Second column multiplies into ijk
                case(5) ! YZ
                    allocate(buf(1,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                    indMod(:,1) = [1,0,0] ! First column is a static offset
                    indMod(:,2) = [0,1,1] ! Second column multiplies into ijk
                case(6) ! ZX
                    allocate(buf(this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,1,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                    indMod(:,1) = [0,1,0] ! First column is a static offset
                    indMod(:,2) = [1,0,1] ! Second column multiplies into ijk
                case(7) ! XYZ
                    allocate(buf(this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                    indMod(:,1) = [0,0,0] ! First column is a static offset
                    indMod(:,2) = [1,1,1] ! Second column multiplies into ijk
                end select

                do iDef=1,this%nDef
                    buf = 1.0_WP
                    if (this%def(iDef)%nTA.gt.0) then
                        allocate(buf2(this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))
                        do i=this%stats(iLoc)%imin_,this%stats(iLoc)%imax_
                            do j=this%stats(iLoc)%jmin_,this%stats(iLoc)%jmax_
                                do k=this%stats(iLoc)%kmin_,this%stats(iLoc)%kmax_
                                    ! Get indices which correspond to stat arrays
                                    ! iS = sum(indMod(1,:)*[1,i])
                                    ! jS = sum(indMod(2,:)*[1,j])
                                    ! kS = sum(indMod(3,:)*[1,k])
                                    ! Correct for volume of cfg vs stat cell differences
                                    corr =  sum([this%cfg%dx(i),this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin)]*real(indMod(1,:),WP))/ &
                                                            (this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin))                      * &
                                            sum([this%cfg%dy(j),this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin)]*real(indMod(2,:),WP))/ &
                                                            (this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin))                      * &
                                            sum([this%cfg%dz(k),this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin)]*real(indMod(3,:),WP))/ &
                                                            (this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin))
                                    ! Apply the volume correction to the buffer
                                    buf2(i,j,k) = corr
                                    ! Multiply the buffer by the values in the original arrays
                                    do iA=1,this%def(iDef)%nTA
                                        ! print*,'iA',iA,'ijk',i,j,k
                                        buf2(i,j,k) = buf2(i,j,k) * this%arp(this%def(iDef)%iAs(iA))%val(i,j,k)
                                    end do
                                end do
                            end do
                        end do
                        select case(this%stats(iLoc)%dim)
                        case (0) ! 0D
                            buf(1,1,1) = sum(buf2)
                        case (1) ! X
                            buf(:,1,1) = sum(sum(buf2,3),2)
                        case (2) ! Y
                            buf(1,:,1) = sum(sum(buf2,3),1)
                        case (3) ! Z
                            buf(1,1,:) = sum(sum(buf2,2),1)
                        case (4) ! XY
                            buf(:,:,1) = sum(buf2,3)
                        case (5) ! YZ
                            buf(1,:,:) = sum(buf2,1)
                        case (6) ! ZX
                            buf(:,1,:) = sum(buf2,2)
                        case (7) ! XYZ
                            buf = buf2
                        end select 
                        deallocate(buf2)
                    end if
                    ! At this point the buffer is either a) all 1.0_WP (only particle stats) or 
                    !   b) the factor to multiply into each particle when collecting particle statistics
                    ! Loop over all the particles, if we are doing particles
                    if (sum(abs(this%def(iDef)%expP)).ne.0) then
                        do i=1,this%lp%np_
                            ! If particle is inside our station
                            if ((this%lp%p(i)%ind(1).ge.this%stats(iLoc)%imin_.and.this%lp%p(i)%ind(1).le.this%stats(iLoc)%imax_).and. &
                                (this%lp%p(i)%ind(2).ge.this%stats(iLoc)%jmin_.and.this%lp%p(i)%ind(2).le.this%stats(iLoc)%jmax_).and. &
                                (this%lp%p(i)%ind(3).ge.this%stats(iLoc)%kmin_.and.this%lp%p(i)%ind(3).le.this%stats(iLoc)%kmax_).and. &
                                 this%lp%p(i)%flag.eq.0) then 

                                ! Get the indices that correspond to the stats arrays
                                iS = sum(indMod(1,:)*[1,this%lp%p(i)%ind(1)])
                                jS = sum(indMod(2,:)*[1,this%lp%p(i)%ind(2)])
                                kS = sum(indMod(3,:)*[1,this%lp%p(i)%ind(3)])
                                ! Add to our stats
                                this%stats(iLoc)%vals(iDef,iS,jS,kS) = this%stats(iLoc)%vals(iDef,iS,jS,kS)+dt*buf(iS,jS,kS)* &
                                    product([1.0_WP,this%lp%p(i)%d,this%lp%p(i)%T_d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]**this%def(iDef)%expP)
                            end if
                        end do
                    else
                        this%stats(iLoc)%vals(iDef,:,:,:) = this%stats(iLoc)%vals(iDef,:,:,:) + buf*dt
                    end if
                end do
                deallocate(buf)
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
                ! Allocate array to hold the global (for station) statistics
                select case (this%stats(iLoc)%dim)
                case (0) ! 0D
                    allocate(buf(     this%nDef,1,1,1))
                    allocate(globVals(this%nDef,1,1,1)); globVals = 0.0_WP
                    globVals = this%stats(iLoc)%vals
                case (1) ! X
                    allocate(buf(     this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,1))
                    allocate(globVals(this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,1)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,:,:) = this%stats(iLoc)%vals
                case (2) ! Y
                    allocate(buf(     this%nDef,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1))
                    allocate(globVals(this%nDef,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1)); globVals = 0.0_WP
                    globVals(:,:,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,:) = this%stats(iLoc)%vals
                case (3) ! Z
                    allocate(buf(     this%nDef,1,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%nDef,1,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,:,:,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                case (4) ! XY
                    allocate(buf(     this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1))
                    allocate(globVals(this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,1)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,:) = this%stats(iLoc)%vals
                case (5) ! YZ
                    allocate(buf(     this%nDef,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%nDef,1,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,:,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                case (6) ! ZX
                    allocate(buf(     this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,1,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,:,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                case (7) ! XYZ
                    allocate(buf(     this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))
                    allocate(globVals(this%nDef,this%stats(iLoc)%imin:this%stats(iLoc)%imax,this%stats(iLoc)%jmin:this%stats(iLoc)%jmax,this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)); globVals = 0.0_WP
                    globVals(:,this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_) = this%stats(iLoc)%vals
                end select

                ! Sum over all partitions within the station
                numVals = size(globVals,1)*size(globVals,2)*size(globVals,3)*size(globVals,4)
                call MPI_ALLREDUCE(globVals,buf,numVals,MPI_REAL_WP,MPI_SUM,this%stats(iLoc)%comm,ierr); globVals = buf
                deallocate(buf)


                do i = 1,this%nDef
                    ! if (this%stats(iLoc)%amRoot) then
                    !     print*,'i',i
                    !     print*,globVals(i,:,:,:)
                    ! end if

                    ! If no particles, this is simple
                    if (sum(abs(this%def(i)%expP)).eq.0) then
                        globVals(i,:,:,:) = globVals(i,:,:,:)/this%sum_time
                    elseif (i.ne.this%iNum) then ! Don't divide out particle number yet... that'd mess everything else up
                        ! If there are particles, need to be careful to not divide by zero
                        where (globVals(this%iNum,:,:,:).gt.0.0_WP) ! Only normalize where we found particles
                            globVals(i,:,:,:) = globVals(i,:,:,:)/globVals(this%iNum,:,:,:)
                        elsewhere
                            globVals(i,:,:,:) = 0.0_WP
                        end where
                    end if
                end do
                if (this%iNum.ne.-1) then ! Normalize the particle number, but only if it actually exists
                    globVals(this%iNum,:,:,:) = globVals(this%iNum,:,:,:)/this%sum_time
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

        ! Required assumption: 
        ! Save glob_vals

        ! Save all definition information
        ! Save arr_point info
            ! Perhaps can have "repoint" subroutine for when reading this data?
        ! 

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