!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   use mpi_f08,      only: MPI_Group
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg1 ! Config for main flow domain
   type(config), public :: cfg2 ! Config for jet development region
   type(MPI_Group), public :: grp1, grp2 ! MPI groups associated with cfg1 and cfg2


   logical, public :: amGrp1, amGrp2 ! Boolean if this processor is within group1 or grp2

   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      use parallel, only: group
      use mpi_f08, only : MPI_GROUP_SIZE,MPI_Group,MPI_Group_incl
      implicit none
      type(sgrid) :: grid1,grid2
      integer, dimension(3) :: partition1,partition2
      integer :: nproc,ierr,nDom2,i
      integer, allocatable, dimension(:) :: ranks2
      integer :: nproc2,myRank

      
      ! Determine which processors are doing what

      ! Read in partition
      call param_read('1 Partition',partition1,short='p')
      call param_read('2 Partition',partition2,short='p')

      call MPI_GROUP_SIZE(group,nproc,ierr)
      call MPI_Group_rank(group,myRank,ierr)
      nDom2 = product(partition2)

      ! Assume that grid 1 will use all processors in group.
      amGrp1=.true.
      grp1 = group
      ! Create subgroup for grid 2 if necessary
      if (nproc.gt.nDom2) then
         ! Create vector of the ranks to be included in group2
         allocate(ranks2(nDom2))
         do i=1,nDom2
            ranks2(i) = i-1
         end do
         call MPI_Group_incl(group,nDom2,ranks2,grp2,ierr)
         call MPI_GROUP_SIZE(grp2,nproc2,ierr)
         deallocate(ranks2)
      else
         grp2 = group
      end if
      

      if (myRank.le.(nDom2-1)) then
         amGrp2 = .true.
      else
         amGrp2 = .false.
      end if
      
      ! Initialize normal geometry things for grid 1
      if (amGrp1) then
         ! Create a grid from input params
         create_grid1: block
            use sgrid_class, only: cartesian
            integer :: i,j,k,nx,ny,nz
            real(WP) :: Lx,Ly,Lz
            real(WP), dimension(:), allocatable :: x,y,z
            
            ! Read in grid definition
            call param_read('Lx',Lx); call param_read('nx1',nx); allocate(x(nx+1))
            call param_read('Ly',Ly); call param_read('ny1',ny); allocate(y(ny+1))
            call param_read('Lz',Lz); call param_read('nz1',nz); allocate(z(nz+1))
            
            ! Create simple rectilinear grid
            do i=1,nx+1
               x(i)=real(i-1,WP)/real(nx,WP)*Lx
            end do
            do j=1,ny+1
               y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
            end do
            do k=1,nz+1
               z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
            end do
            
            ! General serial grid object
            grid1=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='evaporating_spray')
            ! Add walls
         end block create_grid1
         
         
         ! Create a config from that grid on our entire group
         create_cfg1: block                        
            ! Create partitioned grid
            cfg1=config(grp=group,decomp=partition1,grid=grid1)
         end block create_cfg1
         
         ! Create masks for this config
         create_walls1: block
            integer :: i,j,k
            cfg1%VF=1.0_WP
            !do k=cfg1%kmino_,cfg1%kmaxo_
            !   do j=cfg1%jmino_,cfg1%jmaxo_
            !      do i=cfg1%imino_,cfg1%imaxo_
            !         if (abs(cfg1%xm(i)).lt.0.1_WP.and.abs(cfg1%ym(j)).lt.0.2_WP) cfg1%VF(i,j,k)=0.0_WP
            !      end do
            !   end do
            !end do
         end block create_walls1
      end if

      ! Initialize normal geometry things for grid 2
      if (amGrp2) then
         ! Create a grid for jet development from input params
         create_grid2: block
            use sgrid_class, only: cartesian
            integer :: i,j,k,nx,ny,nz
            real(WP) :: Lx,Ly,Lz
            real(WP) :: r_jet,D,L_over_D
            real(WP), dimension(:), allocatable :: x,y,z
            
            ! Figure out how big to make jet domain:
            call param_read('nx2',nx); allocate(x(nx+1))
            call param_read('ny2',ny); allocate(y(ny+1))
            call param_read('nz2',nz); allocate(z(nz+1))
            call param_read('Gas jet radius',r_jet); D = r_jet*2.0_WP
            call param_read('Jet development L over D',L_over_D)
            Lx = D*L_over_D
            if (ny.eq.1) then
               Ly=Lx/nx!call param_read('Ly',Ly)
            else
               Ly = D*real(ny,WP)/(real(ny-2,WP))
            end if
            if (nz.eq.1) then
               Lz=Lx/nx!call param_read('Lz',Lz)
            else
               Lz = D*real(nz,WP)/(real(nz-2,WP))
            end if

            ! Center grid on jet; x=0 at exit of pre-jet region
            do i=1,nx+1
               x(i)=(real(i-1,WP)/real(nx,WP)-1.0_WP)*Lx
            end do
            do j=1,ny+1
               y(j)=real(j-1,WP)/real(ny,WP)*Ly-0.5_WP*Ly
            end do
            do k=1,nz+1
               z(k)=real(k-1,WP)/real(nz,WP)*Lz-0.5_WP*Lz
            end do
            ! General serial grid object
            grid2=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.true.,yper=.false.,zper=.false.,name='jet_development')
         end block create_grid2
         
         
         ! Create a config from that grid on our entire group
         create_cfg2: block
            ! Create partitioned grid
            cfg2=config(grp=grp2,decomp=partition2,grid=grid2)            
         end block create_cfg2
         
         
         ! Create masks for this config
         create_walls2: block
            ! Need walls to surround jet
            real(WP) :: r_injG
            integer :: i,j,k
            cfg2%VF=0.0_WP
            call param_read('Gas jet radius',r_injG)
            do k=cfg2%kmino_,cfg2%kmaxo_
               do j=cfg2%jmino_,cfg2%jmaxo_
                  if ((sqrt(cfg2%ym(j)**2+cfg2%zm(k)**2)).le.r_injG) cfg2%VF(:,j,k)=1.0_WP
                  ! do i=cfg2%imino_,cfg2%imaxo_
                  !    if ((sqrt(cfg2%ym(j)**2+cfg2%zm(k)**2)).le.r_injG) cfg2%VF(i,j,k)=1.0_WP
                  !    ! Add a little chunk of wall in the middle to help induce turbulence (Fake surface roughness)
                  !    if ((abs(cfg2%zm(k)).lt.cfg2%dz(k)).and.((cfg2%nx/2).eq.i).and.((cfg2%ym(j).lt.r_injG).and.(cfg2%ym(j).gt.r_injG-cfg2%dy(j)))) cfg2%VF(i,j,k)=0.0_WP
                  ! end do
               end do
            end do
         end block create_walls2
      end if
      
   end subroutine geometry_init
   
   
end module geometry
