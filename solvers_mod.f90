module solvers_mod
    use precision_mod
    use problem_mod
    implicit none
    private
    public :: integrate_rk4, integrate_adams_extrap, integrate_adams_interp

contains

    !  RK4
    subroutine rk4_step(t, X, X_next)
        real(dp), intent(in) :: t, X(dim)
        real(dp), intent(out) :: X_next(dim)
        real(dp) :: k1(dim), k2(dim), k3(dim), k4(dim)

        call f(t, X, k1);                  k1 = h * k1
        call f(t + 0.5_dp*h, X + 0.5_dp*k1, k2); k2 = h * k2
        call f(t + 0.5_dp*h, X + 0.5_dp*k2, k3); k3 = h * k3
        call f(t + h,       X + k3,        k4); k4 = h * k4

        X_next = X + (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4) / 6.0_dp
    end subroutine rk4_step

    subroutine integrate_rk4(filename)
        character(len=*), intent(in) :: filename
        integer :: i, n_steps, unit
        real(dp) :: X(dim), X_old(dim)

        n_steps = nint(t_end / h)
        open(newunit=unit, file=filename, status='replace', action='write')

        X = x0
        write(unit, '(*(ES25.16,1X))') 0.0_dp, X

        do i = 1, n_steps
            X_old = X
            call rk4_step(real(i-1, dp)*h, X_old, X)
            write(unit, '(*(ES25.16,1X))') real(i, dp)*h, X
        end do

        close(unit)
    end subroutine integrate_rk4

    !  Коэффициенты Адамса
    ! A_nj = (-1)^j / (j! (n-1-j)!) * int_0^1 prod_{k=0}^{n-1}(z+k)/(z+j) dz
    subroutine adams_extrap_coeffs(n, A)
        integer, intent(in) :: n
        real(dp), intent(out) :: A(0:n-1)
        real(dp), allocatable :: poly(:), quot(:)
        real(dp) :: integral, rem, sign
        integer :: j, k, idx

        allocate(poly(0:n), quot(0:n-1))

        do j = 0, n - 1
            poly = 0.0_dp
            poly(0) = 1.0_dp
            do idx = 0, n - 1
                call poly_mul_linear(poly, idx, real(idx, dp))
            end do
            call poly_div_linear(poly, n, real(j, dp), quot, rem)
            integral = 0.0_dp
            do k = 0, n - 1
                integral = integral + quot(k) / real(k+1, dp)
            end do
            sign = merge(-1.0_dp, 1.0_dp, modulo(j, 2) == 1)
            A(j) = sign * integral / (factorial_int(j) * factorial_int(n-1-j))
        end do

        deallocate(poly, quot)
    end subroutine adams_extrap_coeffs

    ! B_nj = (-1)^{j+1} / ((j+1)! (n-2-j)!) *
    !        int_0^1 [(z-1)z(z+1)...(z+n-2)]/(z+j) dz
    subroutine adams_interp_coeffs(n, B)
        integer, intent(in) :: n
        real(dp), intent(out) :: B(0:n-1)   ! B(0)=B_{n,-1}, B(k)=B_{n,k-1}
        real(dp), allocatable :: poly(:), quot(:)
        real(dp) :: integral, rem, sign
        integer :: j, k, idx

        allocate(poly(0:n), quot(0:n-1))

        do j = -1, n - 2
            poly = 0.0_dp
            poly(0) = 1.0_dp
            do idx = 0, n - 1
                call poly_mul_linear(poly, idx, real(idx - 1, dp))
            end do
            call poly_div_linear(poly, n, real(j, dp), quot, rem)
            integral = 0.0_dp
            do k = 0, n - 1
                integral = integral + quot(k) / real(k+1, dp)
            end do
            sign = merge(1.0_dp, -1.0_dp, modulo(j, 2) == 1)
            if (j == -1) then
                B(0) = sign * integral / (factorial_int(0) * factorial_int(n-1))
            else
                B(j+1) = sign * integral / (factorial_int(j+1) * factorial_int(n-2-j))
            end if
        end do

        deallocate(poly, quot)
    end subroutine adams_interp_coeffs

    ! poly = poly * (z + a)
    subroutine poly_mul_linear(poly, deg, a)
        integer, intent(in) :: deg
        real(dp), intent(inout) :: poly(0:)
        real(dp), intent(in) :: a
        real(dp) :: tmp(0:deg+1)
        integer :: i
        tmp = 0.0_dp
        do i = 0, deg
            tmp(i)   = tmp(i)   + poly(i) * a
            tmp(i+1) = tmp(i+1) + poly(i)
        end do
        poly(0:deg+1) = tmp
    end subroutine poly_mul_linear

    subroutine poly_div_linear(poly_in, deg, a, quot, rem)
        integer, intent(in) :: deg
        real(dp), intent(in) :: poly_in(0:deg), a
        real(dp), intent(out) :: quot(0:deg-1), rem
        real(dp) :: carry
        integer :: k
        carry = poly_in(deg)
        quot(deg-1) = carry
        do k = deg - 1, 1, -1
            carry = poly_in(k) - a * carry
            quot(k-1) = carry
        end do
        rem = poly_in(0) - a * carry
    end subroutine poly_div_linear

    pure real(dp) function factorial_int(n) result(f)
        integer, intent(in) :: n
        integer :: i
        f = 1.0_dp
        do i = 2, n
            f = f * real(i, dp)
        end do
    end function factorial_int

    !  Экстраполяционный метод Адамса
    subroutine integrate_adams_extrap(filename)
        character(len=*), intent(in) :: filename
        integer :: i, j, n_steps, unit
        real(dp), allocatable :: X_hist(:,:), F_hist(:,:)
        real(dp) :: A(0:n_adams-1), rhs(dim)

        n_steps = nint(t_end / h)
        call adams_extrap_coeffs(n_adams, A)

        allocate(X_hist(0:n_steps, dim), F_hist(0:n_steps, dim))

        X_hist(0, :) = x0
        call f(0.0_dp, x0, F_hist(0, :))

        do i = 1, n_adams - 1
            call rk4_step(real(i-1, dp)*h, X_hist(i-1, :), X_hist(i, :))
            call f(real(i, dp)*h, X_hist(i, :), F_hist(i, :))
        end do

        !Adams
        do i = n_adams - 1, n_steps - 1
            rhs = 0.0_dp
            do j = 0, n_adams - 1
                rhs = rhs + A(j) * F_hist(i - j, :)
            end do
            X_hist(i+1, :) = X_hist(i, :) + h * rhs
            call f(real(i+1, dp)*h, X_hist(i+1, :), F_hist(i+1, :))
        end do

        open(newunit=unit, file=filename, status='replace', action='write')
        do i = 0, n_steps
            write(unit, '(*(ES25.16,1X))') real(i, dp)*h, X_hist(i, :)
        end do
        close(unit)

        deallocate(X_hist, F_hist)
    end subroutine integrate_adams_extrap

    !  Интерполяционный метод Адамса с методом Ньютона
    subroutine integrate_adams_interp(filename)
        character(len=*), intent(in) :: filename
        integer :: i, k, n_steps, unit
        real(dp), allocatable :: X_hist(:,:), F_hist(:,:)
        real(dp) :: B(0:n_adams-1), sum_known(dim), pred(dim), X_new(dim)
        real(dp) :: t_next
        logical :: ok

        n_steps = nint(t_end / h)
        call adams_interp_coeffs(n_adams, B)

        allocate(X_hist(0:n_steps, dim), F_hist(0:n_steps, dim))

        X_hist(0, :) = x0
        call f(0.0_dp, x0, F_hist(0, :))

        do i = 1, n_adams - 1
            call rk4_step(real(i-1, dp)*h, X_hist(i-1, :), X_hist(i, :))
            call f(real(i, dp)*h, X_hist(i, :), F_hist(i, :))
        end do

        ! Adams
        do i = n_adams - 1, n_steps - 1
            t_next = real(i+1, dp) * h

            sum_known = 0.0_dp
            do k = 1, n_adams - 1
                sum_known = sum_known + B(k) * F_hist(i - (k-1), :)
            end do

            pred = X_hist(i, :) + h * F_hist(i, :)   ! предиктор Эйлера
            call newton_implicit(t_next, pred, X_hist(i, :), sum_known, B(0), X_new, ok)
            if (.not. ok) then
                write(*,'(A,I0)') 'Внимание: Ньютон не сошёлся на шаге ', i
            end if

            X_hist(i+1, :) = X_new
            call f(t_next, X_new, F_hist(i+1, :))
        end do

        open(newunit=unit, file=filename, status='replace', action='write')
        do i = 0, n_steps
            write(unit, '(*(ES25.16,1X))') real(i, dp)*h, X_hist(i, :)
        end do
        close(unit)

        deallocate(X_hist, F_hist)
    end subroutine integrate_adams_interp

    ! X - X_prev - h*(B_minus1*f(t_next,X) + sum_known) = 0
    subroutine newton_implicit(t_next, X_guess, X_prev, sum_known, B_minus1, X_out, ok)
        real(dp), intent(in) :: t_next, X_guess(dim), X_prev(dim)
        real(dp), intent(in) :: sum_known(dim), B_minus1
        real(dp), intent(out) :: X_out(dim)
        logical, intent(out) :: ok
        real(dp) :: X(dim), G(dim), J(dim,dim), dX(dim), f_next(dim), rhs(dim)
        real(dp), parameter :: tol = 1.0e-12_dp
        integer, parameter :: max_iter = 100
        integer :: iter, i

        X = X_guess
        ok = .false.

        do iter = 1, max_iter
            call f(t_next, X, f_next)
            G = X - X_prev - h * (B_minus1 * f_next + sum_known)
            if (maxval(abs(G)) < tol * max(1.0_dp, maxval(abs(X)))) then
                ok = .true.
                exit
            end if

            call jacobian_f(t_next, X, J)
            J = -h * B_minus1 * J
            do i = 1, dim
                J(i, i) = J(i, i) + 1.0_dp
            end do

            rhs = -G
            call solve_linear(J, rhs, dX)
            X = X + dX
        end do

        X_out = X
    end subroutine newton_implicit

    subroutine jacobian_f(t, X, J)
        real(dp), intent(in) :: t, X(dim)
        real(dp), intent(out) :: J(dim, dim)
        real(dp) :: Xp(dim), f0(dim), fp(dim), delta
        integer :: k
        call f(t, X, f0)
        do k = 1, dim
            delta = sqrt(epsilon(1.0_dp)) * max(abs(X(k)), 1.0_dp)
            Xp = X
            Xp(k) = Xp(k) + delta
            call f(t, Xp, fp)
            J(:, k) = (fp - f0) / delta
        end do
    end subroutine jacobian_f

    ! Решение A*x = b методом Гаусса без выбора ведущего элемента
    subroutine solve_linear(A, b, x)
        real(dp), intent(inout) :: A(dim, dim), b(dim)
        real(dp), intent(out) :: x(dim)
        integer :: i, j, k
        real(dp) :: factor, tmp

        do k = 1, dim - 1
            do i = k + 1, dim
                factor = A(i, k) / A(k, k)
                do j = k, dim
                    A(i, j) = A(i, j) - factor * A(k, j)
                end do
                b(i) = b(i) - factor * b(k)
            end do
        end do

        do i = dim, 1, -1
            tmp = b(i)
            do j = i + 1, dim
                tmp = tmp - A(i, j) * x(j)
            end do
            x(i) = tmp / A(i, i)
        end do
    end subroutine solve_linear

end module solvers_mod
