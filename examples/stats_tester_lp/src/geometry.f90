!> Various definitions and tools for initializing NGA2 config
module geometry
   use config_class, only: config
   use precision,    only: WP
   implicit none
   private
   
   !> Single config
   type(config), public :: cfg
   
   public :: geometry_init
   
contains
   
   
   !> Initialization of problem geometry
   subroutine geometry_init
      use sgrid_class, only: sgrid
      use param,       only: param_read
      implicit none
      type(sgrid) :: grid
      
      
      ! Create a grid from input params
      create_grid: block
         use sgrid_class, only: cartesian
         integer :: i,j,k,nx,ny,nz
         real(WP) :: Lx,Ly,Lz,Sx,Sy,Sz
         real(WP), dimension(:), allocatable :: x,y,z
         
         ! Read in grid definition
         call param_read('Lx',Lx); call param_read('nx',nx); allocate(x(nx+1))
         call param_read('Ly',Ly); call param_read('ny',ny); allocate(y(ny+1))
         call param_read('Lz',Lz); call param_read('nz',nz); allocate(z(nz+1))

         call param_read('Stretching in x', Sx)
         call param_read('Stretching in y', Sy)
         call param_read('Stretching in z', Sz)
         ! Create stretched rectilinear grid
         do i = 1, nx + 1
            x(i) = Lx*(sinh(Sx*(real(i - 1, WP)/real(nx, WP)))/sinh(Sx))
            end do
            do j = 1, ny + 1
            y(j) = 0.5_WP*Ly*(sinh(Sy*(real(j - 1, WP)/real(ny, WP) - 0.5_WP))/sinh(Sy*0.5_WP))
            end do
            do k = 1, nz + 1
            z(k) = 0.5_WP*Lz*(sinh(Sz*(real(k - 1, WP)/real(nz, WP) - 0.5_WP))/sinh(Sz*0.5_WP))
         end do

         ! General serial grid object
         grid=sgrid(coord=cartesian,no=2,x=x,y=y,z=z,xper=.false.,yper=.false.,zper=.false.,name='EvaporatingSpray')
         ! Add walls
      end block create_grid
      
      
      ! Create a config from that grid on our entire group
      create_cfg: block
         use parallel, only: group
         integer, dimension(3) :: partition
         
         ! Read in partition
         call param_read('Partition',partition,short='p')
         
         ! Create partitioned grid
         cfg=config(grp=group,decomp=partition,grid=grid)
         
      end block create_cfg
      
      
      ! Create masks for this config
      create_walls: block
         integer :: i,j,k
         cfg%VF=1.0_WP
      end block create_walls
      
      
   end subroutine geometry_init
   
   
end module geometry
