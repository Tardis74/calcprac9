program main
    use solvers_mod
    implicit none

    call integrate_rk4('rk.dat')
    call integrate_adams_extrap('ae.dat')
    call integrate_adams_interp('ai.dat')

    print '(A)', 'Записаны файлы: rk.dat, ae.dat, ai.dat'
end program main
