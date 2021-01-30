
module updatePseudoSpeciesChem

use iso_c_binding

real(c_double) :: temperature ! temperature [K]

contains

subroutine solveODE(Y, neq, T_init, dt) bind(C, name='solveODE')

    implicit none

    !external FEX, JEX

    ! Inputs / Outputs

    ! neq : integer
    !   lengths of Y
    integer(c_long), intent(in) :: neq ! c_long is an integer C-type
    ! Y : 1d array
    !   vector concentration [mol/m^3]
    real(c_double), intent(inout) :: Y(*)
    ! T_init : float
    !   temperature [K]
    real(c_double), intent(in) :: T_init
    ! dt : float
    !   convection time step [s]. The chemistry will be update
    !   on the convection time step (cf operator splitting procedure)
    real(c_double), intent(in) :: dt

    ! Parameters and variables for subroutine DVODE
    integer(kind=4) :: i,itol,itask,istate,iopt,mf
    integer(kind=8) :: liw,lrw
    real(kind=8) :: t, t_out
    real(kind=8) :: rtol,atol
    real(kind=8), allocatable, dimension(:) :: rwork
    integer(kind=4), allocatable, dimension(:) :: iwork
    real(kind=8) :: RPAR ! not used in my case because I don't exchange float data
    integer(kind=4) :: IPAR ! not used in my case because I don't exchange integer data
    ! real(kind=8) :: Yold(neq), YDOT(1:neq)

    ! Physics constants
    real(kind=8), parameter :: Na=6.0221409e+23 ! Avogadro

    ! Local variables
    integer(kind=8) :: j

    ! Convert vector Y from [mol/m^3] to [particules/m^3]
    do j=1,neq
      Y(j)=Y(j)*Na
      ! Yold(j)=Y(j)
    enddo

    temperature = T_init ! make it accessible to FEX and JEX
    ! neq = int(n,4)

    itol=1
    rtol=1e-9 ! relative tolerance
    atol=1.0e-10*SUM(Y(1:neq)) ! absolute tolerance
    itask=1
    istate=1
    iopt=1 ! to enable or disable the cutomization of DVODE with IWORK and RWORK
    lrw=22+9*NEQ+2*NEQ**2 ! advised by line 547 of dvode.f90
    liw=30+neq ! advised by line 565 od dvode.f90
    allocate(rwork(1:lrw))
    allocate(iwork(1:liw))
    rwork(:)=0.0
    ! RWORK(6)=1.0e-11 ! HMAX
    RWORK(7)=0.0 ! HMIN
    ! If IOPT=1, then IWORK is used to pass some options to the ODE solver
    iwork(:)=0 ! set the options to default
    iwork(6)=1000000 ! MXSTEP = maximum number of (internally defined) steps (default : 500)
    ! Warning, in the previous option, all the iworks must have the same integer kind
    ! (kind=4 or kind=8) in all routines (the current one and vode.f)
    ! otherwise iwork(6) will be understood in iwork(12) for example
    mf=22
    t=0.0

    t_out=dt
    ! print *,'dt = ', dt
    call DVODE(FEX,NEQ,Y,t,t_out,ITOL,RTOL,ATOL,ITASK,ISTATE, &
               IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)

    ! Check successful integration
    ! if (ISTATE .ne. 2) then
    !   print *,"A problem occurs in DVODE"
    !   print *,"ATOL = ",atol
    !   print *,"ISTATE = ", ISTATE
    !   print *,"temperature = ",temperature
    !   print *,"Yold = ", Yold(1:neq)
    !   print *,"Y = ", Y(1:neq)
    !   ! call FEX(NEQ, t, Y, YDOT, RPAR, IPAR)
    !   ! print *,"YDOT = ", YDOT(1:neq)
    !   print *,"t = ", t, "t_out = ", t_out
    !   !
    !   ! print *,"Internal loops"
    !   ! t_out = 0.0
    !   ! t = 0.0
    !   ! istate = 1
    !   Y(1:neq) = Yold(1:neq)
    !   ! do j=1,10000
    !   !   t_out = t_out + dt/10000
    !   !   print *,"t = ", t, "t_out = ", t_out
    !   !   call DVODE(FEX,NEQ,Y,t,t_out,ITOL,RTOL,ATOL,ITASK,ISTATE, &
    !   !              IOPT,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
    !   ! enddo
    !   stop
    ! endif

    ! Check nucleus conservation
    ! if (abs(Y(1)+2*sum(Y(2:neq))-(Yold(1)+2*sum(Yold(2:neq))))/(Yold(1)+2*sum(Yold(2:neq)))>1.0e-6) then
    !   print *,"Error in nuclei conservation"
    !   print *,Y(1)+2*sum(Y(2:neq)),Yold(1)+2*sum(Yold(2:neq))
    !   print *,abs(Y(1)+2*sum(Y(2:neq))-(Yold(1)+2*sum(Yold(2:neq))))/(Yold(1)+2*sum(Yold(2:neq)))
    !   ! stop
    ! endif

    ! Convert back vector Y from [particules/m^3] to [mol/m^3]
    do j=1,neq
      Y(j)=Y(j)/Na
    enddo

end subroutine solveODE

! RHS and Jac

subroutine FEX(NEQ, t, Y, YDOT, RPAR, IPAR)
    implicit none
    integer(kind=4) :: NEQ, IPAR
    real(kind=8) :: RPAR, t, Y(1:NEQ), YDOT(1:NEQ)
    real(kind=8) :: temp
    call func(temperature, Y, YDOT)
    return
end

subroutine JEX(NEQ, t, Y, ML, MU, PD, NRPD, RPAR, IPAR)
    implicit none
    integer(kind=4) :: NEQ, ML, MU, NRPD, IPAR
    real(kind=8) :: t, Y(NEQ), PD(NRPD,NEQ), RPAR
    ! TODO generate a Jacobian fortran subroutine with python
    stop
    return
end

end module updatePseudoSpeciesChem
