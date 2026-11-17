module problem_mod
    use precision_mod
    implicit none
    private
    public :: dim, n_adams, t_end, h, x0, f

    integer, parameter :: dim = 2
    integer, parameter :: n_adams = 4
    real(dp), parameter :: t_end = 6.283185307179586476925286766559_dp  ! 2pi
    real(dp), parameter :: h = 0.01_dp
    real(dp), parameter :: x0(dim) = [1.0_dp, 0.0_dp]

contains

    ! Гармонический осциллятор:
    ! dx/dt = v
    ! dv/dt = -omega^2 * x
    pure subroutine f(t, X, dXdt)
        real(dp), intent(in) :: t, X(dim)
        real(dp), intent(out) :: dXdt(dim)
        real(dp), parameter :: omega = 1.0_dp
        dXdt(1) = X(2)      + 0.0_dp * t
        dXdt(2) = -omega*omega*X(1) + 0.0_dp * t
    end subroutine f

end module problem_mod
