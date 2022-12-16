! Collect specialized statistics for sprays
module stats_class
    use precision, only: WP
    use string,    only: str_medium
    use config_class, only: config
    use mpi_f08
    use lpt_class, only: lpt
    implicit none
    private


    public :: stats_parent ! Show stats_parent to the compiler so it stops yelling at me

    ! Contains information on sampling regions for stats gathering
    ! Stations are to be rectangular prisms in the cartesian mesh, though
    !   arbitrary volumes fixed to the lattice are allowed (i.e. does not have to span the mesh in any direction)
    type :: station
        character(len=str_medium) :: name = 'Unknown_station'               ! Name of our station. Must be unique (case-insensitive).
        logical :: amIn = .false.                                           ! Is this processor in this station?
        integer :: dim                                                      ! Dimension of station, x=1,y=2,z=3,xy=4,yz=5,zx=6,xyz=7
        integer, dimension(3,2) :: msk                                      ! Mask for modifying mesh indices to station%vals indices
        integer :: imin ,imax ,jmin ,jmax ,kmin ,kmax                       ! Global min/maxes in x,y,z. Indices correspond to spray_parent%cfg
        integer :: imin_,imax_,jmin_,jmax_,kmin_,kmax_                      ! Local min/maxes in x,y,z. Indices correspond to spray_parent%cfg
        real(WP), dimension(:,:,:,:), allocatable :: vals                   ! Values for each statistic (before post-processing)
        real(WP), dimension(:,:,:,:,:), allocatable :: binned_vals          ! Values for each diameter-binned statistic (before post-processing)
        type(MPI_Comm) :: all_comm = MPI_COMM_NULL                          ! Communicator for processors on this station
        type(MPI_Comm) :: sum_comm = MPI_COMM_NULL                          ! Communicator for summing between processors on this station
        type(MPI_Comm) :: write_comm = MPI_COMM_NULL                        ! Communicator for processors on this station which write stats data
        type(MPI_Datatype) :: view  = MPI_DATATYPE_NULL                     ! Local to station-global array mapping info - real(WP)
        logical :: amRoot                                                   ! Is this process the root of all_comm?
        logical :: should_write = .false.                                   ! Is this process to be involved in writing data
    end type station

    ! Array of pointers to the data fed into the statistics and associated names
    type :: arr_point
        real(WP), dimension(:,:,:), pointer :: val                  ! Pointer to this array
        character(len=str_medium) :: name = 'Unknown_arr_point'     ! Name of this array. Must be unique (case-insensitive).
    end type arr_point

    ! Contains the parameters for creating one statistical measure 
    type :: stat_def
        integer :: nTA                                  ! Number of arrays in arp used for this statistic
        integer,  dimension(:), allocatable :: iAs      ! Ind(ex)(ices) for each array within arp
        integer,  dimension(7) :: expP                  ! Exponents corresponding to number,lp%p%d,T_d,U,V,W,radial velocity
        character(len=str_medium) :: name               ! Name of this definition. Must be unique (case-insensitive).
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

        !!!!!!!! Binning: All particle stats are binned if use_bins is true !!!!!!!!
        logical :: use_bins = .false.                               ! Are we binning particle data by diameter?
        integer :: n_bins = 0                                       ! Number of bins for particle diameters
        integer, dimension(:), allocatable :: i_b2def               ! Map of indices from binned_vals to the definition array
        integer :: iNum_b                                           ! Index in def_b which corresponds to particle number
        integer :: nDef_b = 0                                       ! Number of definitions which include particle statistics
        real(WP), dimension(:), allocatable :: d_bins               ! Diameters for the lower cutoff of each bin

    contains
        procedure :: restart_from_file      ! Load all the data from the restart file
        procedure :: add_station            ! Add a station for stats measuring
        procedure :: add_array              ! Add a pointer to a solver array
        procedure :: add_definition         ! Add the definition for one statistical measure. To be called after adding all arrays
        procedure :: add_bins               ! Add bins for particle statistics based on diameter
        procedure :: set_lp                 ! Set up the lp pointer
        procedure :: init_stats             ! Initialize statistics. To be called after adding all stations, arrays, bins, and definitions
        procedure :: sample_stats           ! Sample stats over all stations
        procedure :: write                  ! Store data to the disk
        procedure :: post_process           ! Process data from binary file into usable forms
        procedure :: reset_stats            ! Reset statistics, maintaining stations
        procedure :: get_def                ! Get the index of a stat definition
        procedure :: get_reduced            ! Return the reduced statistics for a given station (normalizes results by time and where applicable, particle number)
        

    end type stats_parent

    interface stats_parent
        procedure construct
    end interface stats_parent
    
contains

    ! Constructor for a brand new statistic object
    function construct(cfg,name) result(self)
        implicit none
        type(stats_parent) :: self
        class(config), target, intent(in) :: cfg
        character(len=*), intent(in), optional :: name

        ! Set the stats_parent's name
        if (present(name)) self%name=trim(adjustl(name))

        ! Point to our config object
        self%cfg=>cfg

        ! Initialize various counters
        self%nLoc = 0
        self%nAr  = 0
        self%nDef = 0
        self%sum_time = 0.0_WP

        ! Default lp pointer to null
        self%lp => NULL()
    end function construct

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

    ! Given a stats definition name, return the index of that definition
    function get_bin_def(this,name) result(iDef)
        use messager, only: die
        use string, only: lowercase
        class(stats_parent),intent(in) :: this
        character(len=*), intent(in) :: name
        integer :: iDef
        integer :: i
        iDef = -1
        do i=1,this%nDef_b
            if (lowercase(trim(this%def(this%i_b2def(i))%name)).eq.lowercase(trim(name))) iDef = this%i_b2def(i)
        end do
        if (iDef.eq.-1) call die('[stats_class get_bin_def] Unknown statistics definition was passed.')
    end function get_bin_def

    ! Given a station name, return the index of that station
    function get_loc(this,name) result(iLoc)
        use messager, only: die
        use string, only: lowercase
        class(stats_parent),intent(in) :: this
        character(len=*), intent(in) :: name
        integer :: iLoc,i
        iLoc = -1
        do i=1,this%nLoc
            if (lowercase(trim(this%stats(i)%name)).eq.lowercase(trim(name))) iLoc = i
        end do
        if (iLoc.eq.-1) call die('[stats_class get_loc] Unknown station was passed.')
    end function get_loc


    ! Load in information to the stats_parent from a restart file
    ! Not performed here:
    !   - pointing to array pointers
    !   - pointing to lpt
    ! It is safe to call the normal add_(station,definition,array),init_stats assuming the names are not changed from the initial run
    ! The restarted statistics must use the same sgrid as the restart file, but the partitioning may change
    subroutine restart_from_file(this,filename)
        use messager, only: die
        use parallel, only: info_mpiio,MPI_REAL_WP
        use mpi_f08
        implicit none
        class(stats_parent) :: this
        character(len=*), intent(in) :: filename
        character(len=18) :: the_date_time
        integer :: i,j,k,ierr,n,iLoc,i_b
        integer :: key,color,myRank
        real(WP) :: corr
        real(WP), dimension(:), allocatable :: x,y,z
        real(WP) :: buf
        integer, dimension(3) :: dims
        type(MPI_File) :: ifile
        type(MPI_Status) :: status
        integer(kind=MPI_OFFSET_KIND) :: pos
        integer, dimension(3) :: gsizes,lsizes,lstart
        integer :: nx_,ny_,nz_,my_sum_rank

        ! Open the file
        call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
        if (ierr.ne.0) call die('[stats constructor] Problem encountered while reading data file: '//trim(filename))
        
        ! Christen the stats object
        call MPI_FILE_READ_ALL(ifile,this%name,str_medium,MPI_CHARACTER,status,ierr)

        ! Dumping date and time into dummy variable
        call MPI_FILE_READ_ALL(ifile,the_date_time,18,MPI_CHARACTER,status,ierr) 

        ! Get various counters
        call MPI_FILE_READ_ALL(ifile,this%nLoc,1,MPI_INTEGER,status,ierr)
        call MPI_FILE_READ_ALL(ifile,this%nAr ,1,MPI_INTEGER,status,ierr)
        call MPI_FILE_READ_ALL(ifile,this%nDef,1,MPI_INTEGER,status,ierr)

        ! Allocate the array pointers, definitions, and locations
        allocate(this%def(  this%nDef))
        allocate(this%arp(  this%nAr))
        allocate(this%stats(this%nLoc))

        ! Get summation time
        call MPI_FILE_READ_ALL(ifile,this%sum_time,1,MPI_REAL_WP,status,ierr)

        ! Get bin info
        call MPI_FILE_READ_ALL(ifile,this%n_bins,1,MPI_INTEGER,status,ierr)
        call MPI_FILE_READ_ALL(ifile,this%nDef_b,1,MPI_INTEGER,status,ierr)
        if (this%n_bins.gt.0) then
            allocate(this%i_b2def(this%nDef_b))
            allocate(this%d_bins(this%n_bins))
            call MPI_FILE_READ_ALL(ifile,this%i_b2def,this%nDef_b,MPI_INTEGER,status,ierr)
            call MPI_FILE_READ_ALL(ifile,this%d_bins,this%n_bins,MPI_REAL_WP,status,ierr)
            this%use_bins = .true.
        end if


        ! Loop over definitions
        do n=1,this%nDef
            ! Get the name
            call MPI_FILE_READ_ALL(ifile,this%def(n)%name,str_medium,MPI_CHARACTER,status,ierr)
            ! Read expP
            call MPI_FILE_READ_ALL(ifile,this%def(n)%expP,7,MPI_INTEGER,status,ierr)
            ! Read nTA
            call MPI_FILE_READ_ALL(ifile,this%def(n)%nTA,1,MPI_INTEGER,status,ierr)
            ! Read iAs
            if (this%def(n)%nTA.gt.0) then
                allocate(this%def(n)%iAs(this%def(n)%nTA))
                call MPI_FILE_READ_ALL(ifile,this%def(n)%iAs,this%def(n)%nTA,MPI_INTEGER,status,ierr)
            end if
            ! If this is particle number, note as such
            if (this%def(n)%nTA.eq.0.and.sum(abs(this%def(n)%expP)).eq.1.and.this%def(n)%expP(1).eq.1) this%iNum = n
        end do
        do n=1,this%nDef_b ! only runs if nDef_b>0
            ! Write the map from binned definition index to global definition index
            if (this%i_b2def(n).eq.this%iNum) this%iNum_b = n
        end do

        ! Loop over array pointers
        do n=1,this%nAr
            call MPI_FILE_READ_ALL(ifile,this%arp(n)%name,str_medium,MPI_CHARACTER,status,ierr)
        end do

        ! To facilitate looping, close the file now, saving our position
        call MPI_FILE_GET_POSITION(ifile,pos,ierr)
        call MPI_FILE_CLOSE(ifile,ierr)

        ! Loop over stations
        do iLoc=1,this%nLoc
            ! Reopen the file on all processes and reset our position
            call MPI_FILE_OPEN(this%cfg%comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
            call MPI_FILE_SEEK(ifile,pos,MPI_SEEK_SET,ierr)

            ! Get the name
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%name,str_medium,MPI_CHARACTER,status,ierr)
            ! Get the dimension
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%dim,1,MPI_INTEGER,status,ierr)
            ! Get the mesh-based indices
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%imin,1,MPI_INTEGER,status,ierr)
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%imax,1,MPI_INTEGER,status,ierr)
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%jmin,1,MPI_INTEGER,status,ierr)
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%jmax,1,MPI_INTEGER,status,ierr)
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%kmin,1,MPI_INTEGER,status,ierr)
            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%kmax,1,MPI_INTEGER,status,ierr)

            ! Figure out if the config coordinates match. This should also catch if mesh stretching changed between runs
            allocate(x(this%stats(iLoc)%imin:this%stats(iLoc)%imax))
            allocate(y(this%stats(iLoc)%jmin:this%stats(iLoc)%jmax))
            allocate(z(this%stats(iLoc)%kmin:this%stats(iLoc)%kmax))

            do i=this%stats(iLoc)%imin,this%stats(iLoc)%imax
                call MPI_FILE_READ_ALL(ifile,x(i),1,MPI_REAL_WP,status,ierr)
            end do
            do j=this%stats(iLoc)%jmin,this%stats(iLoc)%jmax
                call MPI_FILE_READ_ALL(ifile,y(j),1,MPI_REAL_WP,status,ierr)
            end do
            do k=this%stats(iLoc)%kmin,this%stats(iLoc)%kmax
                call MPI_FILE_READ_ALL(ifile,z(k),1,MPI_REAL_WP,status,ierr)
            end do
            buf = sum(abs(this%cfg%x(this%stats(iLoc)%imin:this%stats(iLoc)%imax)-x))+&
                  sum(abs(this%cfg%y(this%stats(iLoc)%jmin:this%stats(iLoc)%jmax)-y))+&
                  sum(abs(this%cfg%z(this%stats(iLoc)%kmin:this%stats(iLoc)%kmax)-z))
            if ((this%stats(iLoc)%imin.lt.this%cfg%imin).or.(this%stats(iLoc)%imax.gt.this%cfg%imax).or. &
                (this%stats(iLoc)%jmin.lt.this%cfg%jmin).or.(this%stats(iLoc)%jmax.gt.this%cfg%jmax).or. &
                (this%stats(iLoc)%kmin.lt.this%cfg%kmin).or.(this%stats(iLoc)%kmax.gt.this%cfg%kmax).or. &
                buf.gt.1.0e-6_WP*(abs(max(maxval(x)-minval(x),maxval(y)-minval(y),maxval(z)-minval(z)))))  then
                if (this%cfg%amRoot) print*,'Residual of pgrid coordinates minus coordinates from file:',buf
                if (this%cfg%amRoot) call die('[stats constructor] Coordinates from file are incompatible with passed pgrid.')
            end if
            deallocate(x)
            deallocate(y)
            deallocate(z)
            ! Dimensions of the information to be read... more useful for external reading of the data
            call MPI_FILE_READ_ALL(ifile,dims,3,MPI_INTEGER,status,ierr)
            call MPI_FILE_GET_POSITION(ifile,pos,ierr)
            ! Close the file now so that only processes within this station open it later
            call MPI_FILE_CLOSE(ifile,ierr)
            ! Set the index mask for this station
            select case (this%stats(iLoc)%dim)
            case(0) ! 0D
                this%stats(iLoc)%msk(:,1) = [1,1,1] ! First column is a static offset
            case(1) ! X
                this%stats(iLoc)%msk(:,1) = [0,1,1] ! First column is a static offset
            case(2) ! Y
                this%stats(iLoc)%msk(:,1) = [1,0,1] ! First column is a static offset
            case(3) ! Z
                this%stats(iLoc)%msk(:,1) = [1,1,0] ! First column is a static offset
            case(4) ! XY
                this%stats(iLoc)%msk(:,1) = [0,0,1] ! First column is a static offset
            case(5) ! YZ
                this%stats(iLoc)%msk(:,1) = [1,0,0] ! First column is a static offset
            case(6) ! ZX
                this%stats(iLoc)%msk(:,1) = [0,1,0] ! First column is a static offset
            case(7) ! XYZ
                this%stats(iLoc)%msk(:,1) = [0,0,0] ! First column is a static offset
            end select
            this%stats(iLoc)%msk(:,2) = 1-this%stats(iLoc)%msk(:,1) ! Second column multiplies into ijk

            ! Create communicator for all the processors in this station

            
            ! Figure out if we're even in this station
            if (.not.((this%cfg%imin_.gt.this%stats(iLoc)%imax).or.(this%cfg%imax_.lt.this%stats(iLoc)%imin).or. &
                (this%cfg%jmin_.gt.this%stats(iLoc)%jmax).or.(this%cfg%jmax_.lt.this%stats(iLoc)%jmin).or. &
                (this%cfg%kmin_.gt.this%stats(iLoc)%kmax).or.(this%cfg%kmax_.lt.this%stats(iLoc)%kmin))) then

                    this%stats(iLoc)%amIn = .true.
                ! Set the local array indices
                ! X
                if (this%cfg%imin_.gt.this%stats(iLoc)%imin) then
                    this%stats(iLoc)%imin_ = this%cfg%imin_
                else
                    this%stats(iLoc)%imin_ = this%stats(iLoc)%imin
                end if
                if (this%cfg%imax_.lt.this%stats(iLoc)%imax) then
                    this%stats(iLoc)%imax_ = this%cfg%imax_
                else
                    this%stats(iLoc)%imax_ = this%stats(iLoc)%imax
                end if
                ! Y
                if (this%cfg%jmin_.gt.this%stats(iLoc)%jmin) then
                    this%stats(iLoc)%jmin_ = this%cfg%jmin_
                else
                    this%stats(iLoc)%jmin_ = this%stats(iLoc)%jmin
                end if
                if (this%cfg%jmax_.lt.this%stats(iLoc)%jmax) then
                    this%stats(iLoc)%jmax_ = this%cfg%jmax_
                else
                    this%stats(iLoc)%jmax_ = this%stats(iLoc)%jmax
                end if
                ! Z
                if (this%cfg%kmin_.gt.this%stats(iLoc)%kmin) then
                    this%stats(iLoc)%kmin_ = this%cfg%kmin_
                else
                    this%stats(iLoc)%kmin_ = this%stats(iLoc)%kmin
                end if
                if (this%cfg%kmax_.lt.this%stats(iLoc)%kmax) then
                    this%stats(iLoc)%kmax_ = this%cfg%kmax_
                else
                    this%stats(iLoc)%kmax_ = this%stats(iLoc)%kmax
                end if
                
                key = this%cfg%rank
                color = 0
            else
                this%stats(iLoc)%amIn = .false.
                this%stats(iLoc)%amRoot = .false.
                key = 0
                color = MPI_UNDEFINED
            end if            
            ! Make new communicator for this station
            call MPI_Comm_split(this%cfg%comm,color,key,this%stats(iLoc)%all_comm)
            ! Get new rank and group
            this%stats(iLoc)%amRoot = .false.
            myRank = -1
            if (this%stats(iLoc)%amIn) then
                call MPI_Comm_rank(this%stats(iLoc)%all_comm,myRank,ierr)
                if (myRank.eq.0) this%stats(iLoc)%amRoot = .true.
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!! Setup for processors in this station !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ! From now on only processors in this station need to do things
            if (this%stats(iLoc)%amIn) then
                ! Create communicators and views for summing and writing here
                select case (this%stats(iLoc)%dim)
                case(0) ! 0D
                    ! Need:
                    !   - all processes sum with each other---need 1 sum_comm
                    !   - only the station-global root process writes
                    key = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1+this%cfg%npy*(this%cfg%kproc-1))
                    color = 1
                case(1) ! X
                    ! Need:
                    !   - processes aligned in y and z sum with each other---need npx independent sum_comm
                    !   - root process in y and z (i.e. the one in the negative y - negative z corner) writes
                    key = (this%cfg%jproc-1)+this%cfg%npy*(this%cfg%kproc-1)
                    color = this%cfg%iproc-1
                case(2) ! Y
                    ! Need:
                    !   - processes aligned in x and z sum with each other---need npy independent sum_comm
                    !   - root process in x and z (i.e. the one in the negative x - negative z corner) writes
                    key = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%kproc-1)
                    color = this%cfg%jproc-1

                case(3) ! Z
                    ! Need:
                    !   - processes aligned in x and y sum with each other---need npz independent sum_comm
                    !   - root process in x and y (i.e. the one in the negative x - negative y corner) writes
                    key = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1)
                    color = this%cfg%kproc-1

                case(4) ! XY
                    ! Need:
                    !   - processes aligned in z sum with each other---need npx*npy independent sum_comm
                    !   - root process in z for each sum_comm writes
                    key = this%cfg%kproc-1
                    color = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1)

                case(5) ! YZ
                    ! Need:
                    !   - processes aligned in x sum with each other---need npy*npz independent sum_comm
                    !   - root process in x for each sum_comm writes
                    key   = this%cfg%iproc-1
                    color = this%cfg%jproc-1+this%cfg%npy*(this%cfg%kproc-1)

                case(6) ! ZX
                    ! Need:
                    !   - processes aligned in y sum with each other---need npx*npz independent sum_comm
                    !   - root process in y for each sum_comm writes
                    key = this%cfg%jproc-1
                    color = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%kproc-1)

                case(7) ! XYZ
                    ! Need:
                    !   - Everyone gets their own sum_comm :)
                    !   - Everyone gets to write :)
                    key = 1
                    color = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1+this%cfg%npy*(this%cfg%kproc-1))
                end select
                ! Create the summing communicator
                call MPI_Comm_split(this%stats(iLoc)%all_comm,color,key,this%stats(iLoc)%sum_comm)

                ! The root process of each summing communicator needs to write---create a communicator amongst those roots
                call MPI_Comm_rank(this%stats(iLoc)%sum_comm,my_sum_rank,ierr)
                if (my_sum_rank.eq.0) then
                    this%stats(iLoc)%should_write = .true.
                    key = 0
                    color = 1
                else
                    key = 0
                    color = MPI_UNDEFINED
                end if

                ! Create communicator for writing (the root of each summing group writes)
                call MPI_Comm_split(this%stats(iLoc)%all_comm,color,key,this%stats(iLoc)%write_comm)

                ! Create view to link local and station-global coordinates
                gsizes=[sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax-this%stats(iLoc)%imin+1]),&
                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax-this%stats(iLoc)%jmin+1]),&
                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax-this%stats(iLoc)%kmin+1])]

                lsizes=[sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_-this%stats(iLoc)%imin_+1]),&
                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_-this%stats(iLoc)%jmin_+1]),&
                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_-this%stats(iLoc)%kmin_+1])]

                lstart=[sum(this%stats(iLoc)%msk(1,:)*[0,this%stats(iLoc)%imin_-this%stats(iLoc)%imin]),&
                        sum(this%stats(iLoc)%msk(2,:)*[0,this%stats(iLoc)%jmin_-this%stats(iLoc)%jmin]),&
                        sum(this%stats(iLoc)%msk(3,:)*[0,this%stats(iLoc)%kmin_-this%stats(iLoc)%kmin])]
                call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_REAL_WP,this%stats(iLoc)%view,ierr)
                call MPI_TYPE_COMMIT(this%stats(iLoc)%view,ierr)

                ! Allocate vals according to dimension of locator
                allocate(this%stats(iLoc)%vals( &
                    sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                    sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                    sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),this%nDef))

                nx_ = size(this%stats(iLoc)%vals,1)
                ny_ = size(this%stats(iLoc)%vals,2)
                nz_ = size(this%stats(iLoc)%vals,3)

                ! Reopen our file with the station's communicator:
                call MPI_FILE_OPEN(this%stats(iLoc)%all_comm,trim(filename),MPI_MODE_RDONLY,info_mpiio,ifile,ierr)
                
                ! Read the non-binned stats data
                do n=1,this%nDef
                    call MPI_FILE_SET_VIEW(ifile,pos,MPI_REAL_WP,this%stats(iLoc)%view,'native',info_mpiio,ierr)
                    call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%vals(:,:,:,n),nx_*ny_*nz_,MPI_REAL_WP,status,ierr)
                    pos=pos+int(dims(1),MPI_OFFSET_KIND)*int(dims(2),MPI_OFFSET_KIND)*int(dims(3),MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
                end do

                ! Read the binned data
                if (this%use_bins) then
                    allocate(this%stats(iLoc)%binned_vals( &
                        sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),this%n_bins,this%nDef_b))
                    do n=1,this%nDef_b
                        do i_b=1,this%n_bins
                            call MPI_FILE_SET_VIEW(ifile,pos,MPI_REAL_WP,this%stats(iLoc)%view,'native',info_mpiio,ierr)
                            call MPI_FILE_READ_ALL(ifile,this%stats(iLoc)%binned_vals(:,:,:,i_b,n),nx_*ny_*nz_,MPI_REAL_WP,status,ierr)
                            pos=pos+int(dims(1),MPI_OFFSET_KIND)*int(dims(2),MPI_OFFSET_KIND)*int(dims(3),MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
                        end do
                    end do

                end if

                ! Correct for differences in local vs global station size
                corr =  sum([(this%cfg%x(this%stats(iLoc)%imax_+1)-this%cfg%x(this%stats(iLoc)%imin_)),this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin)]*real(this%stats(iLoc)%msk(1,:),WP))/ &
                            (this%cfg%x(this%stats(iLoc)%imax +1)-this%cfg%x(this%stats(iLoc)%imin))                      * &
                        sum([(this%cfg%y(this%stats(iLoc)%jmax_+1)-this%cfg%y(this%stats(iLoc)%jmin_)),this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin)]*real(this%stats(iLoc)%msk(2,:),WP))/ &
                            (this%cfg%y(this%stats(iLoc)%jmax +1)-this%cfg%y(this%stats(iLoc)%jmin))                      * &
                        sum([(this%cfg%z(this%stats(iLoc)%kmax_+1)-this%cfg%z(this%stats(iLoc)%kmin_)),this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin)]*real(this%stats(iLoc)%msk(3,:),WP))/ &
                            (this%cfg%z(this%stats(iLoc)%kmax +1)-this%cfg%z(this%stats(iLoc)%kmin))

                ! Un-normalize the locally-stored statistics array, correcting for volume and integration time
                if (this%iNum.ne.-1) then ! Un-normalize the particle number, but only if it actually exists
                    ! We also multiply by corr such that when re-summing later we get back to the original reduced value,
                    !   avoiding a situation in which samples prior to a restart are weighted higher than those after a restart
                    this%stats(iLoc)%vals(:,:,:,this%iNum) = this%stats(iLoc)%vals(:,:,:,this%iNum)*this%sum_time*corr
                end if
                ! Deal with everything that isn't particle number
                do i=1,this%nDef
                    if (sum(abs(this%def(i)%expP)).eq.0) then ! Eulerian stats are un-normalized by sum_time
                        this%stats(iLoc)%vals(:,:,:,i) = this%stats(iLoc)%vals(:,:,:,i)*this%sum_time*corr
                    elseif (i.ne.this%iNum) then ! Multiply by (particle number*time) now
                        this%stats(iLoc)%vals(:,:,:,i) = this%stats(iLoc)%vals(:,:,:,i)*this%stats(iLoc)%vals(:,:,:,this%iNum)
                    end if
                end do

                ! Deal with bins
                if (this%use_bins) then
                    if (this%iNum_b.ne.-1) then ! Un-normalize the particle number, but only if it actually exists
                        this%stats(iLoc)%binned_vals(:,:,:,:,this%iNum_b) = this%stats(iLoc)%binned_vals(:,:,:,:,this%iNum_b)*this%sum_time*corr
                    end if
                    do i_b=1,this%nDef_b ! index in binned_vals
                        if (i_b.ne.this%iNum_b) then ! Multiply by (particle number*time) now
                            this%stats(iLoc)%binned_vals(:,:,:,:,i_b) = this%stats(iLoc)%binned_vals(:,:,:,:,i_b)*this%stats(iLoc)%binned_vals(:,:,:,:,this%iNum_b)
                        end if
                    end do
                end if
                call MPI_FILE_CLOSE(ifile,ierr)
            else ! If this processor isn't in the station, need to advance the readhead to catch up
                ! Account for non-binned data:
                pos=pos+int(dims(1),MPI_OFFSET_KIND)*int(dims(2),MPI_OFFSET_KIND)*int(dims(3),MPI_OFFSET_KIND)*&
                        int(this%nDef,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
                ! Account for binned data (0 if no bins):
                pos=pos+int(dims(1),MPI_OFFSET_KIND)*int(dims(2),MPI_OFFSET_KIND)*int(dims(3),MPI_OFFSET_KIND)*&
                        int(this%nDef_b,MPI_OFFSET_KIND)*int(this%n_bins,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)
            end if
        end do        

    end subroutine restart_from_file

    ! Set lp, initialize statistic definition for particle number (required for all other particle stats)
    subroutine set_lp(this,lp)
        implicit none
        class(stats_parent), intent(inout) :: this
        type(lpt), intent(in), target :: lp

        this%lp => lp
        call this%add_definition(name='p_n',def='p_n') ! if this is already defined, nothing happens (check happens inside add_def)
    end subroutine set_lp

    ! Subroutine for adding stations to the stats sampling
    subroutine add_station(this,name,dim,locator)
        use string, only : lowercase
        use messager, only : die
        use mpi_f08
        use parallel, only: MPI_REAL_WP 
        implicit none
        class(stats_parent),intent(inout) :: this
        character(len=*),intent(in) :: dim
        character(len=*),intent(in),optional :: name
        integer :: i,j,k,ierr,n,iLoc
        integer :: nproc,myRank,key,color
        real(WP) :: buf
        type(station), dimension(:), allocatable :: tmp
        interface
            logical function locator(pargrid,ind1,ind2,ind3)
                use pgrid_class, only: pgrid
                class(pgrid), intent(in) :: pargrid
                integer, intent(in) :: ind1,ind2,ind3
            end function locator
        end interface

        ! Check if this station already exists from a reset
        iLoc = -1
        do i=1,this%nLoc
            if (trim(lowercase(this%stats(i)%name)).eq.trim(lowercase(name))) iLoc=i
        end do
        ! If it doesn't exist, keep going
        if (iLoc.eq.-1) then
            ! Reallocate station array
            if (allocated(this%stats)) then
                allocate(tmp(this%nLoc+1))
                tmp(1:this%nLoc)=this%stats
                call move_alloc(tmp,this%stats)
            else
                allocate(this%stats(1))
            end if

            this%nLoc=this%nLoc+1
            iLoc = this%nLoc
            

            if (present(name)) this%stats(iLoc)%name=trim(adjustl(name))

            ! Parse dimension of station
            select case (lowercase(dim))
            case ('0d')
                this%stats(iLoc)%dim=0
                this%stats(iLoc)%msk(:,1) = [1,1,1]
            case ('x')
                this%stats(iLoc)%dim=1
                this%stats(iLoc)%msk(:,1) = [0,1,1]
            case ('y')
                this%stats(iLoc)%dim=2
                this%stats(iLoc)%msk(:,1) = [1,0,1]
            case ('z')
                this%stats(iLoc)%dim=3
                this%stats(iLoc)%msk(:,1) = [1,1,0]
            case ('xy','yx')
                this%stats(iLoc)%dim=4
                this%stats(iLoc)%msk(:,1) = [0,0,1]
            case ('yz','zy')
                this%stats(iLoc)%dim=5
                this%stats(iLoc)%msk(:,1) = [1,0,0]
            case ('zx','xz')
                this%stats(iLoc)%dim=6
                this%stats(iLoc)%msk(:,1) = [0,1,0]
            case ('xyz','xzy','yxz','yzx','zxy','zyx','3d')
                this%stats(iLoc)%dim=7
                this%stats(iLoc)%msk(:,1) = [0,0,0]
            case default
                unrec_stat_dim : block
                    character(len=100) :: my_message
                    my_message = '[stats add_station] Unrecognized dimension of station '//trim(name)
                    call die(my_message)
                end block unrec_stat_dim
            end select
            this%stats(iLoc)%msk(:,2) = 1-this%stats(iLoc)%msk(:,1)
            ! Loop over mesh to find local bounding indices for the station
            this%stats(iLoc)%imin_ = HUGE(1)
            this%stats(iLoc)%jmin_ = HUGE(1)
            this%stats(iLoc)%kmin_ = HUGE(1)
            this%stats(iLoc)%imax_ =-HUGE(1)
            this%stats(iLoc)%jmax_ =-HUGE(1)
            this%stats(iLoc)%kmax_ =-HUGE(1)
            do i=this%cfg%imin_,this%cfg%imax_
                do j=this%cfg%jmin_,this%cfg%jmax_
                    do k=this%cfg%kmin_,this%cfg%kmax_
                        if (locator(this%cfg,i,j,k)) then
                            this%stats(iLoc)%amIn = .true.
                            this%stats(iLoc)%imin_ = min(this%stats(iLoc)%imin_,i)
                            this%stats(iLoc)%jmin_ = min(this%stats(iLoc)%jmin_,j)
                            this%stats(iLoc)%kmin_ = min(this%stats(iLoc)%kmin_,k)
                            this%stats(iLoc)%imax_ = max(this%stats(iLoc)%imax_,i)
                            this%stats(iLoc)%jmax_ = max(this%stats(iLoc)%jmax_,j)
                            this%stats(iLoc)%kmax_ = max(this%stats(iLoc)%kmax_,k)
                        end if
                    end do
                end do
            end do
            ! Get the global bounding indices
            call MPI_ALLREDUCE(real(this%stats(iLoc)%imin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(iLoc)%imin=int(buf)
            call MPI_ALLREDUCE(real(this%stats(iLoc)%jmin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(iLoc)%jmin=int(buf)
            call MPI_ALLREDUCE(real(this%stats(iLoc)%kmin_,WP),buf,1,MPI_REAL_WP,MPI_MIN,this%cfg%comm,ierr); this%stats(iLoc)%kmin=int(buf)
            call MPI_ALLREDUCE(real(this%stats(iLoc)%imax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(iLoc)%imax=int(buf)
            call MPI_ALLREDUCE(real(this%stats(iLoc)%jmax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(iLoc)%jmax=int(buf)
            call MPI_ALLREDUCE(real(this%stats(iLoc)%kmax_,WP),buf,1,MPI_REAL_WP,MPI_MAX,this%cfg%comm,ierr); this%stats(iLoc)%kmax=int(buf)

            check_in_domain : block
                integer :: int_am_in
                integer :: my_sum
                character(len=str_medium) :: my_message
                int_am_in = 0
                if (this%stats(iLoc)%amIn) int_am_in = 1.0_WP
                call MPI_ALLREDUCE(int_am_in,my_sum,1,MPI_INTEGER,MPI_SUM,this%cfg%comm,ierr)
                if (my_sum.eq.0) then
                    my_message = '[stats add_station] '//name//' is not within the domain.'
                    call die(my_message)
                end if
            end block check_in_domain
            
            create_all_comm : block
            ! Figure out our group and communicator
            call MPI_GROUP_SIZE(this%cfg%group,nproc,ierr)
            call MPI_Group_rank(this%cfg%group,myRank,ierr)
            do i=0,nproc-1
                if (myRank.eq.i) then
                    if (this%stats(iLoc)%amIn) then
                        key = i
                        color = 0
                    else
                        key = 0 
                        color = MPI_UNDEFINED
                    end if
                end if
            end do
            call MPI_Comm_split(this%cfg%comm,color,key,this%stats(iLoc)%all_comm)
            ! Get new rank and group
            this%stats(iLoc)%amRoot = .false.
            myRank = -1
            if (this%stats(iLoc)%amIn) then
                call MPI_Comm_rank(this%stats(iLoc)%all_comm,myRank,ierr)
                if (myRank.eq.0) this%stats(iLoc)%amRoot = .true.
            end if
            end block create_all_comm

            create_sum_write_comms : block
            integer :: my_sum_rank
            integer, dimension(3) :: gsizes,lsizes,lstart
            if (this%stats(iLoc)%amIn) then
                ! Pick the key and color for the sum communicator based on a select case:
                select case (this%stats(iLoc)%dim)
                case(0) ! 0D
                    ! Need:
                    !   - all processes sum with each other---need 1 sum_comm
                    !   - only the station-global root process writes
                    key = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1+this%cfg%npy*(this%cfg%kproc-1))
                    color = 1
                case(1) ! X
                    ! Need:
                    !   - processes aligned in y and z sum with each other---need npx independent sum_comm
                    !   - root process in y and z (i.e. the one in the negative y - negative z corner) writes
                    key = (this%cfg%jproc-1)+this%cfg%npy*(this%cfg%kproc-1)
                    color = this%cfg%iproc-1
                case(2) ! Y
                    ! Need:
                    !   - processes aligned in x and z sum with each other---need npy independent sum_comm
                    !   - root process in x and z (i.e. the one in the negative x - negative z corner) writes
                    key = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%kproc-1)
                    color = this%cfg%jproc-1

                case(3) ! Z
                    ! Need:
                    !   - processes aligned in x and y sum with each other---need npz independent sum_comm
                    !   - root process in x and y (i.e. the one in the negative x - negative y corner) writes
                    key = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1)
                    color = this%cfg%kproc-1

                case(4) ! XY
                    ! Need:
                    !   - processes aligned in z sum with each other---need npx*npy independent sum_comm
                    !   - root process in z for each sum_comm writes
                    key = this%cfg%kproc-1
                    color = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1)

                case(5) ! YZ
                    ! Need:
                    !   - processes aligned in x sum with each other---need npy*npz independent sum_comm
                    !   - root process in x for each sum_comm writes
                    key   = this%cfg%iproc-1
                    color = this%cfg%jproc-1+this%cfg%npy*(this%cfg%kproc-1)

                case(6) ! ZX
                    ! Need:
                    !   - processes aligned in y sum with each other---need npx*npz independent sum_comm
                    !   - root process in y for each sum_comm writes
                    key = this%cfg%jproc-1
                    color = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%kproc-1)

                case(7) ! XYZ
                    ! Need:
                    !   - Everyone gets their own sum_comm :)
                    !   - Everyone gets to write :)
                    key = 1
                    color = (this%cfg%iproc-1)+this%cfg%npx*(this%cfg%jproc-1+this%cfg%npy*(this%cfg%kproc-1))
                end select

                call MPI_Comm_split(this%stats(iLoc)%all_comm,color,key,this%stats(iLoc)%sum_comm)

                ! The root process of each summing communicator needs to write---create a communicator amongst those roots
                call MPI_Comm_rank(this%stats(iLoc)%sum_comm,my_sum_rank,ierr)
                if (my_sum_rank.eq.0) then
                    this%stats(iLoc)%should_write = .true.
                    key = 0
                    color = 1
                else
                    key = 0
                    color = MPI_UNDEFINED
                end if
                ! Create communicator for writing (the root of each summing group writes)
                call MPI_Comm_split(this%stats(iLoc)%all_comm,color,key,this%stats(iLoc)%write_comm)

                ! Create view to link local and station-global coordinates
                gsizes=[sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax-this%stats(iLoc)%imin+1]),&
                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax-this%stats(iLoc)%jmin+1]),&
                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax-this%stats(iLoc)%kmin+1])]

                lsizes=[sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_-this%stats(iLoc)%imin_+1]),&
                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_-this%stats(iLoc)%jmin_+1]),&
                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_-this%stats(iLoc)%kmin_+1])]

                lstart=[sum(this%stats(iLoc)%msk(1,:)*[0,this%stats(iLoc)%imin_-this%stats(iLoc)%imin]),&
                        sum(this%stats(iLoc)%msk(2,:)*[0,this%stats(iLoc)%jmin_-this%stats(iLoc)%jmin]),&
                        sum(this%stats(iLoc)%msk(3,:)*[0,this%stats(iLoc)%kmin_-this%stats(iLoc)%kmin])]
                call MPI_TYPE_CREATE_SUBARRAY(3,gsizes,lsizes,lstart,MPI_ORDER_FORTRAN,MPI_REAL_WP,this%stats(iLoc)%view,ierr)
                call MPI_TYPE_COMMIT(this%stats(iLoc)%view,ierr)
            end if
            end block create_sum_write_comms
        end if
            
    end subroutine add_station

    ! Add array to use when collecting statistics. Must use overlapping pgrid dimensions
    subroutine add_array(this,array,name)
        use string, only : lowercase
        implicit none
        class(stats_parent),intent(inout) :: this
        type(arr_point), dimension(:), allocatable :: tmp
        real(WP), dimension(this%cfg%imino_:,this%cfg%jmino_:,this%cfg%kmino_:), target, intent(in) :: array
        character(len=*),intent(in) :: name
        integer :: iA,i
        
        ! Check if this array already exists from a reset
        iA = -1
        do i=1,this%nAr
            if (trim(lowercase(this%arp(i)%name)).eq.trim(lowercase(name))) iA=i
        end do
        if (iA.eq.-1) then ! Allocate things if this is a new array

            if (allocated(this%arp)) then
                allocate(tmp(this%nAr+1))
                tmp(1:this%nAr)=this%arp
                call move_alloc(tmp,this%arp)
            else
                allocate(this%arp(1))
            end if
            this%nAr = this%nAr+1
            iA = this%nAr
            this%arp(iA)%name = trim(adjustl(name))
        end if
        ! Either way, point to the user-passed array
        this%arp(iA)%val => array

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
        character(len=str_medium),dimension(7) :: p_names
        integer :: i,j,n,iAr,iP,iDef

        ! Check if this definition already exists from a reset
        iDef = -1
        do i=1,this%nDef
            if (trim(lowercase(this%def(i)%name)).eq.trim(lowercase(name))) iDef=i
        end do
        ! If this is a new definition, add all the relevant parameters. Otherwise we're done
        if (iDef.eq.-1) then
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
            iDef = this%nDef
            ! Assign a name
            this%def(iDef)%name = name

            ! Zero-out some things
            this%def(iDef)%nTA = 0
            this%def(iDef)%expP = 0

            ! Make string array of the names of terms for particles
            p_names(1) = 'p_n'; p_names(2) = 'p_d'; p_names(3) = 'p_T'
            p_names(4) = 'p_U'; p_names(5) = 'p_V'; p_names(6) = 'p_W'; p_names(7) = 'p_rad'

            ! Get a string array of the components of the definition
            def_parsed = parse_def(def)
            do i = 1,size(def_parsed)
                do iAr = 1,this%nAr
                    ! Need to count number of arrays to allocate iAs
                    if (lowercase(trim(def_parsed(i))).eq.lowercase(trim(this%arp(iAr)%name))) this%def(iDef)%nTA = this%def(iDef)%nTA+1
                end do
                do iP = 1,7
                    if (lowercase(trim(def_parsed(i))).eq.lowercase(trim(p_names(iP)))) then
                        this%def(iDef)%expP(iP) = this%def(iDef)%expP(iP)+1
                    end if
                end do
            end do
            if (this%def(iDef)%expP(1).eq.1.and.sum(abs(this%def(iDef)%expP)).eq.1.and.this%def(iDef)%nTA.eq.0) this%iNum = iDef
            if (this%def(iDef)%nTA.gt.0) then
                allocate(this%def(iDef)%iAs(this%def(iDef)%nTA))
                j = 1
                do i = 1,size(def_parsed)
                    do iAr = 1,this%nAr
                        ! Set this value of iAr to iAs
                        if (lowercase(trim(def_parsed(i))).eq.lowercase(trim(this%arp(iAr)%name))) then
                            this%def(iDef)%iAs(j) = iAr
                            j = j+1
                        end if
                    end do
                end do
            else
                allocate(this%def(iDef)%iAs(1))
                this%def(iDef)%iAs = 0
            end if
            
            ! Check that user has given us a particle solver if they want particle statistics
            if(sum(abs(this%def(iDef)%expP)).ne.0.and..not.associated(this%lp)) call die('[stats_class add_definition] No LP solver is associated with stat_parent but requested statistic uses LP')
        end if
    end subroutine add_definition

    ! Set the values for binning diameters... the more complex stuff is saved for init_stats
    subroutine add_bins(this,d_bins)
        implicit none
        class(stats_parent),intent(inout) :: this
        real(WP), dimension(:), intent(in) :: d_bins

        ! If we're here, that means we're using bins
        this%use_bins = .true.
        if (.not.allocated(this%d_bins)) then
            this%n_bins = size(d_bins)
            allocate(this%d_bins(this%n_bins))
            this%d_bins = d_bins
        end if

    end subroutine add_bins

    ! Subroutine for initializing stats arrays based on the stations which have already been initialized
    subroutine init_stats(this)
        use messager, only: die
        implicit none
        class(stats_parent),intent(inout) :: this
        real(WP), dimension(:,:,:,:), allocatable :: tmp
        integer :: i,j,k,n,iLoc,iDef
        integer :: n_use_bins
        integer, dimension(:), allocatable :: tmp_i_bins

        ! Set up the this%i_b2def array
        if (this%use_bins) then
            if (.not.allocated(this%i_b2def)) then
                allocate(tmp_i_bins(this%nDef)) ! Over-allocate i_b2def at first
                n_use_bins = 0
                do iDef=1,this%nDef
                    if (sum(abs(this%def(iDef)%expP)).gt.0) then
                        n_use_bins = n_use_bins+1
                        tmp_i_bins(n_use_bins) = iDef
                    end if
                end do
                if(n_use_bins.eq.0) call die('[stats init_stats] Binned statistics requested but no definitions use particles.')
                allocate(this%i_b2def(n_use_bins))
                this%i_b2def = tmp_i_bins(1:n_use_bins)
                this%nDef_b = n_use_bins
                deallocate(tmp_i_bins)
            end if
            ! Determine iNum_b
            this%iNum_b = -1
            do i=1,this%nDef_b
                if (this%iNum.eq.this%i_b2def(i)) this%iNum_b = i
            end do
            if (this%iNum_b.eq.-1) call die('[stats init_stats] Somehow you made it this far before it died, congrats. Associate lp to use bins.')
        end if
        ! Loop over all stations, only entering if the processor is in the station, then allocate all the things
        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                ! Allocate arrays based on dim
                if (.not.allocated(this%stats(iLoc)%vals)) then ! If not allocated yet, allocate and set to
                    allocate(this%stats(iLoc)%vals( &
                        sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),this%nDef))
                    this%stats(iLoc)%vals = 0.0_WP
                elseif(size(this%stats(iLoc)%vals,4).ne.this%nDef) then
                    ! In the rare case that the user decides to add a new definition after a restart, forcibly tell them that they can't do that.
                    !   This would require changing how summing time is implemented.
                    call die('[stats init_stats] Additional definitions have been added since restart, this is not implemented.')
                end if ! If already initialized from restart, don't overwrite
                ! Binned data
                
                if (.not.this%use_bins) cycle ! Don't do this if we aren't binning :)
                if (.not.allocated(this%stats(iLoc)%binned_vals)) then
                    allocate(this%stats(iLoc)%binned_vals( &
                        sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),&
                        this%n_bins, this%nDef_b))
                    this%stats(iLoc)%binned_vals = 0.0_WP
                elseif(size(this%stats(iLoc)%binned_vals,5).ne.this%nDef_b) then
                    call die('[stats init_stats] Additional binned definitions have been added since restart, this is not implemented.')
                end if
            end if
        end do
        if (this%cfg%amRoot) print*,'Successfully Initialized',this%nLoc,'Stations for Statistics'
    end subroutine init_stats


    ! Subroutine for collecting statistics
    subroutine sample_stats(this,dt)
        use messager, only : die
        use mpi_f08, only: MPI_ALLREDUCE,MPI_SUM
        use parallel, only: MPI_REAL_WP
        implicit none
        class(stats_parent),intent(inout) :: this
        real(WP), intent(in) :: dt
        integer :: i,j,k,iDef,iLoc,iS,jS,kS,iA,ierr,iDef_b
        real(WP) :: corr ! Corrects for sampling volume when stat cell is larger than mesh cell
        integer, dimension(3,2) :: msk ! Modify indices for different station dimensions
        real(WP), dimension(:,:,:), allocatable:: buf,buf2 ! Buffer for building stats

        ! Increment summing time
        this%sum_time = this%sum_time+dt

        ! Loop over all stations, only entering if the processor is in the station, then gather stats
        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                ! Allocate buffer
                allocate(buf(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                             sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                             sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_])))
                iDef_b = 0
                do iDef=1,this%nDef
                    buf = 1.0_WP
                    if (this%def(iDef)%nTA.gt.0) then
                        ! Use buf2 to gather across p-mesh of this station, then dump into buffer with size of stats(iLoc)%vals
                        allocate(buf2(this%stats(iLoc)%imin_:this%stats(iLoc)%imax_,this%stats(iLoc)%jmin_:this%stats(iLoc)%jmax_,this%stats(iLoc)%kmin_:this%stats(iLoc)%kmax_))

                        do i=this%stats(iLoc)%imin_,this%stats(iLoc)%imax_
                            do j=this%stats(iLoc)%jmin_,this%stats(iLoc)%jmax_
                                do k=this%stats(iLoc)%kmin_,this%stats(iLoc)%kmax_
                                    ! Correct for volume of cfg vs stat cell differences
                                    corr =  sum([this%cfg%dx(i),this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin)]*real(this%stats(iLoc)%msk(1,:),WP))/ &
                                                            (this%cfg%x(this%stats(iLoc)%imax+1)-this%cfg%x(this%stats(iLoc)%imin))                      * &
                                            sum([this%cfg%dy(j),this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin)]*real(this%stats(iLoc)%msk(2,:),WP))/ &
                                                            (this%cfg%y(this%stats(iLoc)%jmax+1)-this%cfg%y(this%stats(iLoc)%jmin))                      * &
                                            sum([this%cfg%dz(k),this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin)]*real(this%stats(iLoc)%msk(3,:),WP))/ &
                                                            (this%cfg%z(this%stats(iLoc)%kmax+1)-this%cfg%z(this%stats(iLoc)%kmin))
                                    ! Apply the volume correction to the buffer
                                    buf2(i,j,k) = corr
                                    ! Multiply the buffer by the values in the original arrays
                                    do iA=1,this%def(iDef)%nTA
                                        buf2(i,j,k) = buf2(i,j,k) * this%arp(this%def(iDef)%iAs(iA))%val(i,j,k)
                                    end do
                                end do
                            end do
                        end do
                        ! Do our sums according to the dimension of the station
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
                    !   b) the volume weighted mesh quantity to multiply into each particle when collecting particle statistics

                    ! Loop over all the particles, if we are doing particles
                    if (sum(abs(this%def(iDef)%expP)).ne.0) then
                        iDef_b = iDef_b+1 ! Increment the binned definition counter

                        ! Unlike the purely eulerian stats, we need to complete the averaging _now_ for the lagrangian-eulerian statistics
                        if (this%def(iDef)%nTA.gt.0) then
                            allocate(buf2(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                                        sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                                        sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_])))
                            i = size(buf,1)*size(buf,2)*size(buf,3)
                            ! Combine across partitions if includes eulerian info; otherwise buf remains = 1.0_WP
                            ! Weighted averaging with corr is done now.
                            call MPI_ALLREDUCE(buf,buf2,i,MPI_REAL_WP,MPI_SUM,this%stats(iLoc)%all_comm,ierr); buf=buf2
                            deallocate(buf2)
                        end if
                        do i=1,this%lp%np_
                            ! If particle is inside our station
                            if ((this%lp%p(i)%ind(1).ge.this%stats(iLoc)%imin_.and.this%lp%p(i)%ind(1).le.this%stats(iLoc)%imax_).and. &
                                (this%lp%p(i)%ind(2).ge.this%stats(iLoc)%jmin_.and.this%lp%p(i)%ind(2).le.this%stats(iLoc)%jmax_).and. &
                                (this%lp%p(i)%ind(3).ge.this%stats(iLoc)%kmin_.and.this%lp%p(i)%ind(3).le.this%stats(iLoc)%kmax_).and. &
                                 this%lp%p(i)%flag.eq.0) then 

                                ! Get the indices that correspond to the stats arrays
                                iS = sum(this%stats(iLoc)%msk(1,:)*[1,this%lp%p(i)%ind(1)])
                                jS = sum(this%stats(iLoc)%msk(2,:)*[1,this%lp%p(i)%ind(2)])
                                kS = sum(this%stats(iLoc)%msk(3,:)*[1,this%lp%p(i)%ind(3)])
                                ! Add to our stats
                                this%stats(iLoc)%vals(iS,jS,kS,iDef) = this%stats(iLoc)%vals(iS,jS,kS,iDef)+dt*buf(iS,jS,kS)* &
                                    product([1.0_WP,this%lp%p(i)%d,this%lp%p(i)%T_d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3),&
                                    sum([this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]*[this%lp%p(i)%pos(2),this%lp%p(i)%pos(3)]/ &
                                    (epsilon(1.0_WP)+sqrt(this%lp%p(i)%pos(2)**2+this%lp%p(i)%pos(3)**2)))]**this%def(iDef)%expP)
                                
                                if (.not.this%use_bins) cycle
                                ! If binning, add to the bins now
                                this%stats(iLoc)%binned_vals(iS,jS,kS,:,iDef_b) = &
                                    this%stats(iLoc)%binned_vals(iS,jS,kS,:,iDef_b)+dt*buf(iS,jS,kS)* &
                                    product([1.0_WP,this%lp%p(i)%d,this%lp%p(i)%T_d,this%lp%p(i)%vel(1),this%lp%p(i)%vel(2),this%lp%p(i)%vel(3),&
                                    sum([this%lp%p(i)%vel(2),this%lp%p(i)%vel(3)]*[this%lp%p(i)%pos(2),this%lp%p(i)%pos(3)]/ &
                                    (epsilon(1.0_WP)+sqrt(this%lp%p(i)%pos(2)**2+this%lp%p(i)%pos(3)**2)))]**this%def(iDef)%expP)*&
                                    ! Now multiply by a vector which is zero everywhere except for the bin corresponding to p%d, where it is one
                                    -0.25_WP*(sign(1.0_WP,this%lp%p(i)%d- this%d_bins(1:this%n_bins)              )+1.0_WP)* &
                                             (sign(1.0_WP,this%lp%p(i)%d-[this%d_bins(2:this%n_bins),huge(1.0_WP)])-1.0_WP)
                            end if
                        end do
                    else
                        this%stats(iLoc)%vals(:,:,:,iDef) = this%stats(iLoc)%vals(:,:,:,iDef) + buf*dt
                    end if
                end do
                deallocate(buf)
            end if
        end do
        call MPI_Barrier(this%cfg%comm)
        ! if (this%cfg%amRoot) print*,'Statistics collected :)'
    end subroutine sample_stats

    ! Normalize collected statistics over time
    subroutine post_process(this)
        use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
        use parallel, only: MPI_REAL_WP
        class(stats_parent), intent(inout) :: this
        integer :: n,i,i_b,iLoc,ierr,numVals
        real(WP), dimension(:,:,:,:), allocatable :: buf,reduceVals
        real(WP), dimension(:,:,:,:,:), allocatable :: buf_b,reduceVals_b ! Binned versions

        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                ! Allocate array to hold the global (for station) statistics
                allocate(buf(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax]),&
                             sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax]),&
                             sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax]),this%nDef))
                allocate(reduceVals(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax]),&
                                  sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax]),&
                                  sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax]),this%nDef))
                reduceVals = 0.0_WP
                reduceVals(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                         sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                         sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),:) = this%stats(iLoc)%vals
                ! Sum over all partitions within the station
                numVals = size(reduceVals,1)*size(reduceVals,2)*size(reduceVals,3)*size(reduceVals,4)
                call MPI_ALLREDUCE(reduceVals,buf,numVals,MPI_REAL_WP,MPI_SUM,this%stats(iLoc)%all_comm,ierr); reduceVals = buf
                deallocate(buf)

                do i = 1,this%nDef
                    ! If no particles, this is simple
                    if (sum(abs(this%def(i)%expP)).eq.0) then
                        reduceVals(:,:,:,i) = reduceVals(:,:,:,i)/this%sum_time
                    elseif (i.ne.this%iNum) then ! Don't divide out particle number yet... that'd mess everything else up
                        ! If there are particles, need to be careful to not divide by zero
                        where (reduceVals(:,:,:,this%iNum).gt.0.0_WP) ! Only normalize where we found particles
                            reduceVals(:,:,:,i) = reduceVals(:,:,:,i)/reduceVals(:,:,:,this%iNum)
                        elsewhere
                            reduceVals(:,:,:,i) = 0.0_WP
                        end where
                    end if
                end do
                if (this%iNum.ne.-1) then ! Normalize the particle number, but only if it actually exists
                    reduceVals(:,:,:,this%iNum) = reduceVals(:,:,:,this%iNum)/this%sum_time
                end if
                ! Send global values back to the vals array
                this%stats(iLoc)%vals = reduceVals(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                                                 sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                                                 sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),:)
                deallocate(reduceVals)
                if (.not.this%use_bins) cycle
                ! Deal with bins now
                allocate(buf_b(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax]),&
                                sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax]),&
                                sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax]),this%n_bins,this%nDef))
                allocate(reduceVals_b(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax]),&
                                    sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax]),&
                                    sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax]),this%n_bins,this%nDef_b))
                reduceVals_b = 0.0_WP
                reduceVals_b(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                            sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                            sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),:,:) = this%stats(iLoc)%binned_vals
                ! Sum over all partitions within the station
                numVals = size(reduceVals_b,1)*size(reduceVals_b,2)*size(reduceVals_b,3)*size(reduceVals_b,4)*size(reduceVals_b,5)
                call MPI_ALLREDUCE(reduceVals_b,buf_b,numVals,MPI_REAL_WP,MPI_SUM,this%stats(iLoc)%all_comm,ierr); reduceVals_b = buf_b
                deallocate(buf_b)

                do i = 1,this%nDef_b
                    i_b = this%i_b2def(i)
                    if (i_b.ne.this%iNum_b) then ! Don't divide out particle number yet... that'd mess everything else up
                        ! If there are particles, need to be careful to not divide by zero
                        where (reduceVals_b(:,:,:,:,this%iNum_b).gt.0.0_WP) ! Only normalize where we found particles
                            reduceVals_b(:,:,:,:,i) = reduceVals_b(:,:,:,:,i)/reduceVals_b(:,:,:,:,this%iNum_b)
                        elsewhere
                            reduceVals_b(:,:,:,:,i) = 0.0_WP
                        end where
                    end if
                end do
                if (this%iNum_b.ne.-1) then ! Normalize the particle number, but only if it actually exists
                    reduceVals_b(:,:,:,:,this%iNum_b) = reduceVals_b(:,:,:,:,this%iNum_b)/this%sum_time
                end if
                ! Send global values back to the vals array
                this%stats(iLoc)%binned_vals = &
                    reduceVals_b( sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                                sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                                sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),:,:)
                deallocate(reduceVals_b)
            end if
        end do
    end subroutine post_process

    ! Get reduced local values for a single locator/station
    subroutine get_reduced(this,iLoc,reduceVals,reduceVals_b)
        use mpi_f08,  only: MPI_ALLREDUCE,MPI_SUM
        use parallel, only: MPI_REAL_WP
        class(stats_parent), intent(inout) :: this
        integer :: n,i,ierr,numVals,i_b
        integer, intent(in) :: iLoc
        real(WP), dimension(:,:,:,:), allocatable :: buf
        real(WP), dimension(:,:,:,:), allocatable, intent(out) :: reduceVals
        real(WP), dimension(:,:,:,:,:), allocatable :: buf_b ! Buffer for binned global values
        real(WP), dimension(:,:,:,:,:), allocatable, optional, intent(out) :: reduceVals_b ! Binned global values

        if (this%stats(iLoc)%amIn) then
            ! Allocate array to hold the global (for station) statistics
            allocate(buf(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                         sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                         sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),this%nDef))
            allocate(reduceVals(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                              sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                              sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),this%nDef))
            reduceVals = this%stats(iLoc)%vals
            ! Sum over all partitions within the station

            numVals = size(reduceVals,1)*size(reduceVals,2)*size(reduceVals,3)*size(reduceVals,4)
            call MPI_ALLREDUCE(reduceVals,buf,numVals,MPI_REAL_WP,MPI_SUM,this%stats(iLoc)%sum_comm,ierr); reduceVals = buf
            deallocate(buf)
            do i = 1,this%nDef
                ! If no particles, this is simple
                if (sum(abs(this%def(i)%expP)).eq.0) then
                    reduceVals(:,:,:,i) = reduceVals(:,:,:,i)/this%sum_time
                elseif (i.ne.this%iNum) then ! Don't divide out particle number yet... that'd mess everything else up
                    ! If there are particles, need to be careful to not divide by zero
                    where (reduceVals(:,:,:,this%iNum).gt.0.0_WP) ! Only normalize where we found particles
                        reduceVals(:,:,:,i) = reduceVals(:,:,:,i)/reduceVals(:,:,:,this%iNum)
                    elsewhere
                        reduceVals(:,:,:,i) = 0.0_WP
                    end where
                end if
            end do
            if (this%iNum.ne.-1) then ! Normalize the particle number, but only if it actually exists
                reduceVals(:,:,:,this%iNum) = reduceVals(:,:,:,this%iNum)/this%sum_time
            end if
            ! Deal with bins now
            if (.not.this%use_bins) return
            allocate(buf_b(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                            sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                            sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),this%n_bins,this%nDef_b))
            allocate(reduceVals_b(sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imin_]):sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax_]),&
                                sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmin_]):sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax_]),&
                                sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmin_]):sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax_]),this%n_bins,this%nDef_b))
            reduceVals_b = this%stats(iLoc)%binned_vals
            ! Sum over all partitions within the station which align in the summing directions
            numVals = size(reduceVals_b,1)*size(reduceVals_b,2)*size(reduceVals_b,3)*size(reduceVals_b,4)*size(reduceVals_b,5)

            call MPI_ALLREDUCE(reduceVals_b,buf_b,numVals,MPI_REAL_WP,MPI_SUM,this%stats(iLoc)%sum_comm,ierr); reduceVals_b = buf_b
            deallocate(buf_b)

            do i = 1,this%nDef_b
                i_b = this%i_b2def(i)
                if (i_b.ne.this%iNum_b) then ! Don't divide out particle number yet... that'd mess everything else up
                    ! If there are particles, need to be careful to not divide by zero
                    where (reduceVals_b(:,:,:,:,this%iNum_b).gt.0.0_WP) ! Only normalize where we found particles
                        reduceVals_b(:,:,:,:,i) = reduceVals_b(:,:,:,:,i)/reduceVals_b(:,:,:,:,this%iNum_b)
                    elsewhere
                        reduceVals_b(:,:,:,:,i) = 0.0_WP
                    end where
                end if
            end do
            if (this%iNum_b.ne.-1) then ! Normalize the particle number, but only if it actually exists
                reduceVals_b(:,:,:,:,this%iNum_b) = reduceVals_b(:,:,:,:,this%iNum_b)/this%sum_time
            end if
        end if
    end subroutine get_reduced
    
    ! Store current state of statistics to the disk for retrieval later. Doesn't affect current arrays
    subroutine write(this,fdata)
        use messager, only: die
        use parallel, only: MPI_REAL_WP,info_mpiio
        use mpi_f08
        implicit none
        class(stats_parent), intent(inout) :: this
        character(len=*), intent(in), optional :: fdata
        character(len=str_medium) :: filename
        integer :: i,j,k,iLoc,n,ierr,iunit,i_b
        integer(kind=MPI_OFFSET_KIND) :: pos
        type(MPI_File) :: ifile
        type(MPI_Status):: status
        real(WP), dimension(:,:,:,:), allocatable :: reduceVals
        real(WP), dimension(:,:,:,:,:), allocatable :: reduceVals_b
        character(len=8) :: dateChar
        character(len=10) :: timeChar

        ! Create our filename
        if (present(fdata)) then
            filename=trim(adjustl(fdata))
        else
            filename=trim(adjustl(this%name))
        end if
        
        ! Write all non-grid things on root
        if (this%cfg%amRoot) then
            ! Open file
            open(newunit=iunit,file=trim(filename),form='unformatted',status='replace',access='stream',iostat=ierr)
            if (ierr.ne.0) call die('[stats write] Problem encountered while serial-opening data file: '//trim(filename))

            !!!!!!!!! Write parent header info !!!!!!!!!
            ! Write name
            write(iunit) this%name
            ! Write the date and time of writing
            call date_and_time(date=dateChar,time=timeChar)
            write(iunit) dateChar
            write(iunit) timeChar
            ! Write nLoc
            write(iunit) this%nLoc
            ! Write nAr
            write(iunit) this%nAr
            ! Write nDef
            write(iunit) this%nDef
            ! Write sum_time
            write(iunit) this%sum_time

            !!!!!!!!!!!!! Deal with bins !!!!!!!!!!!!
            ! Write number of bins; zero if not using bins
            write(iunit) this%n_bins
            ! Write number of definitions using bins; zero if not using bins
            write(iunit) this%nDef_b
            ! Write the indices of bins which correspond to definitions; nothing is written if there are no bins
            do n=1,this%nDef_b
                write(iunit) this%i_b2def(n)
            end do
            ! Write the lower cutoff values used for the bins
            do n=1,this%n_bins
                write(iunit) this%d_bins(n)
            end do

            !!!!!!!!! Write definition info !!!!!!!!!
            do n=1,this%nDef
                ! Write name
                write(iunit) this%def(n)%name
                ! Write expP
                do i=1,7
                    write(iunit) this%def(n)%expP(i)
                end do
                ! Write nTA
                write(iunit) this%def(n)%nTA
                ! Write iAs
                do i=1,this%def(n)%nTA
                    write(iunit) this%def(n)%iAs(i)
                end do
            end do

            !!!!!!!!! Write array pointer info !!!!!!!!!
            ! This is simple, just save the name
            do n=1,this%nAr
                write(iunit) this%arp(n)%name
            end do

            close(iunit)
        end if

        ! Loop over stations
        do iLoc=1,this%nLoc
            ! Stop the processors in a later station but not the current one from stepping on the processors in the current station
            call MPI_Barrier(this%cfg%comm)

            ! Only do stuff if we're in this station
            if (this%stats(iLoc)%amIn) then
                ! Gather together global data---every process does this
                if (this%use_bins) then
                    call this%get_reduced(iLoc=iLoc,reduceVals=reduceVals,reduceVals_b=reduceVals_b)
                else
                    call this%get_reduced(iLoc=iLoc,reduceVals=reduceVals)
                end if
                ! If we are the station root, write station info
                if (this%stats(iLoc)%amRoot) then
                    open(newunit=iunit,file=trim(filename),form='unformatted',status='old',access='stream',iostat=ierr,position='append')
                    ! Write name
                    write(iunit) this%stats(iLoc)%name
                    ! Write dim
                    write(iunit) this%stats(iLoc)%dim
                    ! Write global indices
                    write(iunit) this%stats(iLoc)%imin
                    write(iunit) this%stats(iLoc)%imax
                    write(iunit) this%stats(iLoc)%jmin
                    write(iunit) this%stats(iLoc)%jmax
                    write(iunit) this%stats(iLoc)%kmin
                    write(iunit) this%stats(iLoc)%kmax

                    ! Write coordinates corresponding to indices for use in post processing
                    do i=this%stats(iLoc)%imin,this%stats(iLoc)%imax
                        write(iunit) this%cfg%x(i)
                    end do
                    do j=this%stats(iLoc)%jmin,this%stats(iLoc)%jmax
                        write(iunit) this%cfg%y(j)
                    end do
                    do k=this%stats(iLoc)%kmin,this%stats(iLoc)%kmax
                        write(iunit) this%cfg%z(k)
                    end do
                    ! Write sizes of upcoming data... not strictly necessary given above knowledge of dimension and whatnot but makes life easier later (mostly for post-processing)
                    write(iunit) sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax-this%stats(iLoc)%imin+1]) ! X
                    write(iunit) sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax-this%stats(iLoc)%jmin+1]) ! Y
                    write(iunit) sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax-this%stats(iLoc)%kmin+1]) ! Z
                    close(iunit)
                end if

                ! All the processes which should be writing write now:
                if (this%stats(iLoc)%should_write) then
                    write_parallel : block
                    integer :: nx_,ny_,nz_,nx,ny,nz ! These are for the sizes of the arrays, NOT the stations coordinates!
                    nx_ = size(reduceVals,1)
                    ny_ = size(reduceVals,2)
                    nz_ = size(reduceVals,3)
                    nx  = sum(this%stats(iLoc)%msk(1,:)*[1,this%stats(iLoc)%imax-this%stats(iLoc)%imin+1])
                    ny  = sum(this%stats(iLoc)%msk(2,:)*[1,this%stats(iLoc)%jmax-this%stats(iLoc)%jmin+1])
                    nz  = sum(this%stats(iLoc)%msk(3,:)*[1,this%stats(iLoc)%kmax-this%stats(iLoc)%kmin+1])
                    call MPI_FILE_OPEN(this%stats(iLoc)%write_comm,trim(filename),IOR(MPI_MODE_WRONLY,MPI_MODE_APPEND),info_mpiio,ifile,ierr)
                    call MPI_FILE_GET_POSITION(ifile,pos,ierr)

                    do n=1,this%nDef
                        call MPI_FILE_SET_VIEW(ifile,pos,MPI_REAL_WP,this%stats(iLoc)%view,'native',info_mpiio,ierr)
                        call MPI_FILE_WRITE_ALL(ifile,reduceVals(:,:,:,n),nx_*ny_*nz_,MPI_REAL_WP,status,ierr)
                        pos=pos+int(nx,MPI_OFFSET_KIND)*int(ny,MPI_OFFSET_KIND)*int(nz,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)                                    
                    end do

                    ! Write binned global data (if needed)
                    do n=1,this%nDef_b ! this loop never runs if there are no bins
                        do i_b=1,this%n_bins
                            call MPI_FILE_SET_VIEW(ifile,pos,MPI_REAL_WP,this%stats(iLoc)%view,'native',info_mpiio,ierr)
                            call MPI_FILE_WRITE_ALL(ifile,reduceVals_b(:,:,:,i_b,n),nx_*ny_*nz_,MPI_REAL_WP,status,ierr)
                            pos=pos+int(nx,MPI_OFFSET_KIND)*int(ny,MPI_OFFSET_KIND)*int(nz,MPI_OFFSET_KIND)*int(WP,MPI_OFFSET_KIND)                                    
                        end do
                    end do
                    call MPI_FILE_CLOSE(ifile,ierr)
                    end block write_parallel
                end if
                ! Deallocate reduceVals
                deallocate(reduceVals)
            end if
        end do
    end subroutine write

    subroutine reset_stats(this)
        class(stats_parent), intent(inout) :: this
        integer :: iLoc

        this%sum_time = 0.0_WP
        do iLoc=1,this%nLoc
            if (this%stats(iLoc)%amIn) then
                this%stats(iLoc)%vals=0.0_WP
                if (this%use_bins) this%stats(iLoc)%binned_vals=0.0_WP
            end if
        end do

    end subroutine reset_stats

end module stats_class