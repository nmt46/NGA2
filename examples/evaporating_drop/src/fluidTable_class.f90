module fluidTable_class
    use precision, only: WP
    use string,    only: str_medium

    public :: fluidTable

    ! Parameters for calling a property from the interpolation function
    integer, public, parameter :: MW_ID  = 0 ! Molecular weight
    integer, public, parameter :: Cp_ID  = 1 ! Specific heat
    integer, public, parameter :: Lv_ID  = 2 ! Latent Heat
    integer, public, parameter :: Tb_ID  = 3 ! Boiling Temperature
    integer, public, parameter :: rho_ID = 4 ! Density
    integer, public, parameter :: mu_ID  = 5 ! Viscosity



    type :: fluidTable
        character(len=str_medium) :: name ! Valid names: 'air', 'ethanol', 'n-Heptane', 'O2[0.21]&N2[0.79]', etc
        integer :: nP, nT
        real(WP), dimension(:), allocatable :: P,T ! associated presure and temperature for properties
        real(WP), dimension(:,:), allocatable :: rho,mu ! functions of both pressure and temperature
        real(WP), dimension(:), allocatable :: L_v,T_b ! Functions of exclusively pressure
        real(WP), dimension(:), allocatable :: Cp ! functions of exclusively temperature
        real(WP) :: MW,Pr,Sc ! invariant of pressure and temperature

    contains
        procedure :: evalProps ! evaluate fluid properties

    end type fluidTable


contains 

    subroutine evalProps(this,propOut,T_q,P_q,propID)
        use messager, only: die
        implicit none
        class(fluidTable), intent(inout) :: this
        real(WP), intent(out) :: propOut ! Prop evaluated at specific point
        real(WP), intent(in) :: T_q,P_q ! Query temperature and pressure
        integer, intent(in) :: propID ! ID corresponding to type of pressure, temperature dependence
        real(WP), dimension(:,:), allocatable :: propArray ! work array of the property table
        integer :: i,j
        integer :: nI,nJ
        if (propID.eq.0) then     ! MW
            propOut = this%MW
            nI = 0; nJ = 0;

        elseif (propID.eq.1) then ! Cp
            allocate(propArray(1,this%nT)); propArray(1,:) = this%Cp
            nI = 1; nJ = this%nT

        elseif (propID.eq.2) then ! Lv
            allocate(propArray(this%nP,1)); propArray(:,1) = this%L_v
            nI = this%nP; nJ = 1

        elseif (propID.eq.3) then ! Tb
            allocate(propArray(this%nP,1)); propArray(:,1) = this%T_b
            nI = this%nP; nJ = 1

        elseif (propID.eq.4) then ! Rho
            allocate(propArray(this%nP,this%nT)); propArray = this%rho
            nI = this%nP; nJ = this%nT

        elseif (propID.eq.5) then ! Mu
            allocate(propArray(this%nP,this%nT)); propArray = this%mu
            nI = this%nP; nJ = this%nT

        else ! Unknown
            call die('Unknown fluid property called in fluidTable.')
        end if
        if (nI.eq.1) then ! f(T)
            ! No extrapolation, just take low/high value if T_q outside T
            if (T_q.lt.minval(this%T)) then 
                propOut = propArray(1,1)
            elseif (T_q.ge.maxval(this%T)) then
                propOut = propArray(1,this%nT)
            else
                ! Loop over values in T to find propOut
                do j=1,nJ-1
                    if ((T_q.ge.this%T(j)).and.(T_q.lt.this%T(j+1))) then
                        propOut = propArray(1,j)+(T_q-this%T(j))*(propArray(1,j+1)-propArray(1,j))/(this%T(j+1)-this%T(j))
                    end if
                end do
            end if


        elseif (nJ.eq.1) then ! f(P)
            ! No extrapolation, just take low/high value if P_q outside P
            if (P_q.lt.minval(this%P)) then 
                propOut = propArray(1,1)
            elseif (P_q.ge.maxval(this%P)) then
                propOut = propArray(this%nP,1)
            else
                ! Loop over values in P to find propOut
                do i=1,nI-1
                    if ((P_q.ge.this%P(i)).and.(P_q.lt.this%P(i+1))) then
                        propOut = propArray(i,1)+(P_q-this%P(i))*(propArray(i+1,1)-propArray(i,1))/(this%P(i+1)-this%P(i))
                    end if
                end do
            end if


        elseif ((nI.gt.1).and.(nJ.gt.1)) then ! f(T,P)
            if (P_q.lt.minval(this%P)) then ! if P_q < min(P), interpolate along P(1), T(x)
                ! Loop over values in T to find propOut
                if (T_q.lt.minval(this%T)) then ! No extrapolation, just take low/high value if T_q outside T
                    propOut = propArray(1,1)
                elseif (T_q.ge.maxval(this%T)) then
                    propOut = propArray(1,this%nT)
                else
                    ! Loop over values in T to find propOut
                    do j=1,nJ-1
                        if ((T_q.ge.this%T(j)).and.(T_q.lt.this%T(j+1))) then
                            propOut = propArray(1,j)+(T_q-this%T(j))*(propArray(1,j+1)-propArray(1,j))/(this%T(j+1)-this%T(j))
                        end if
                    end do
                end if
            elseif (P_q.ge.maxval(this%P)) then ! if P_q > max(P), interpolate along P(end), T(x)
                ! Loop over values in T to find propOut
                if (T_q.lt.minval(this%T)) then ! No extrapolation, just take low/high value if T_q outside T
                    propOut = propArray(nI,1)
                elseif (T_q.ge.maxval(this%T)) then
                    propOut = propArray(nI,this%nT)
                else
                    ! Loop over values in T to find propOut
                    do j=1,nJ-1
                        if ((T_q.ge.this%T(j)).and.(T_q.lt.this%T(j+1))) then
                            propOut = propArray(1,j)+(T_q-this%T(j))*(propArray(1,j+1)-propArray(1,j))/(this%T(j+1)-this%T(j))
                        end if
                    end do
                end if
            else ! P_q within P
                do i=1,nI-1 ! Loop over P
                    if ((P_q.ge.this%P(i)).and.(P_q.lt.this%P(i+1))) then
                        do j=1,nJ-1
                            if ((T_q.ge.this%T(j)).and.(T_q.lt.this%T(j+1))) then
                                propOut = ((this%T(j+1)-T_q)*(propArray(i,j)*  (this%P(i+1)-P_q)+propArray(i+1,j)*  (P_q-this%P(i)))+&
                                           (T_q-this%T(j))*  (propArray(i,j+1)*(this%P(i+1)-P_q)+propArray(i+1,j+1)*(P_q-this%P(i))))/&
                                           ((this%P(i+1)-this%P(i))*(this%T(j+1)-this%T(j)))
                            end if

                        end do
                    end if
                end do
            end if
        end if
        

    end subroutine evalProps

end module fluidTable_class