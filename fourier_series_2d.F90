module constants  
  implicit none 

  double precision, parameter :: pi = 3.1415926536  
  double precision, parameter :: e = 2.7182818285 
end module constants 

subroutine print_array(a, nx, ny)
  implicit none
  integer, intent(in) :: nx, ny
  double precision, intent(in) :: a(nx, ny)
  integer :: i, j

  write(*,*) "===================="
  do i=1, nx
    do j=1, ny
      write(*,'(10F10.2)', advance='no') a(i, j)
      !write(*,*) a(i, j)
    end do
    write(*,*) ""
  end do
  write(*,*) "===================="

end subroutine print_array

subroutine read_arrays(nx, ny, fname, sig, X, Y)
  implicit none
  integer, intent(in) :: nx, ny
  character*7, intent(in) :: fname
  double precision, intent(out) :: sig(nx, ny)
  double precision, intent(out) :: X(nx, ny), Y(nx, ny)
  character*160 :: fname_full

  fname_full = trim("data/") // trim(fname)
  open(unit=22,file=fname_full,form="unformatted",access="stream",status="unknown",action="read")
  read(22) sig
  close(22)

  fname_full = trim("data/") // "X.bin"
  open(unit=23,file=fname_full,form="unformatted",access="stream",status="unknown",action="read")
  read(23) X
  close(23)

  fname_full = trim("data/") // "Y.bin"
  open(unit=24,file=fname_full,form="unformatted",access="stream",status="unknown",action="read")
  read(24) Y
  close(24)
end subroutine read_arrays

subroutine write_arrays(nx, ny, N, sig_fs, alpha, beta, gamma, delta)
  implicit none
  integer, intent(in) :: nx, ny, N
  double precision, intent(in) :: sig_fs(nx, ny)
  double precision, intent(in) :: alpha(N, N), beta(N, N), gamma(N, N), delta(N, N)
  integer :: i

  open (unit=30,form="unformatted",file="data/sig_fs.bin",action="write")
  do i=1,nx
     write(30) sig_fs(i,:) ! write one row at a time
  end do
  close(30)

  open (unit=31,form="unformatted",file="data/alpha.bin",action="write")
  do i=1,N
     write(31) alpha(i,:) ! write one row at a time
  end do
  close(31)

  open (unit=32,form="unformatted",file="data/beta.bin",action="write")
  do i=1,N
     write(32) beta(i,:) ! write one row at a time
  end do
  close(32)

  open (unit=33,form="unformatted",file="data/gamma.bin",action="write")
  do i=1,N
     write(33) gamma(i,:) ! write one row at a time
  end do
  close(33)

  open (unit=34,form="unformatted",file="data/delta.bin",action="write")
  do i=1,N
     write(34) delta(i,:) ! write one row at a time
  end do
  close(34)
end subroutine write_arrays

program main
  use constants, only: pi

  implicit none

  integer :: nx, ny, N
  double precision :: Lx, Ly

  double precision, allocatable :: sig(:, :), sig_fs(:, :)
  double precision, allocatable :: X(:, :), Y(:, :)
  double precision, allocatable :: snxx(:, :), cnxx(:, :), snyy(:, :), cnyy(:, :)
  double precision, allocatable :: snx(:), cnx(:), sny(:), cny(:)
  double precision, allocatable :: alpha(:, :), beta(:, :), gamma(:, :), delta(:, :)
  double precision :: alpha_tmp, beta_tmp, gamma_tmp, delta_tmp

  double precision :: kappa
  integer :: i, j, k, l
  double precision :: dA

  namelist /settings/ nx, ny, N, Lx, Ly

  open(unit=1, file="data/settings.nml", status="OLD")
  read(unit=1, nml=settings)

  dA = (Lx / (nx - 1)) * (Ly / (ny - 1))

  allocate(sig(nx, ny))
  allocate(sig_fs(nx, ny))
  allocate(X(nx, ny))
  allocate(Y(nx, ny))
  allocate(snxx(nx, ny))
  allocate(cnxx(nx, ny))
  allocate(snyy(nx, ny))
  allocate(cnyy(nx, ny))

  allocate(snx(nx))
  allocate(cnx(nx))
  allocate(sny(ny))
  allocate(cny(ny))

  allocate(alpha(N, N))
  allocate(beta(N, N))
  allocate(gamma(N, N))
  allocate(delta(N, N))

  write(*,*)"Running fourier_seried_2d"
  write(*,*) "dA", dA
  write(*,*) "Lx", Lx
  write(*,*) "Ly", Ly
  write(*,*) "N", N

  call read_arrays(nx, ny, "sig.bin", sig, X, Y)

  ! Sanity check. Rem. fortran using 1 based indexing.
  write(*,'(10F10.8)') sig(6, 6)

  !call print_array(sig, nx, ny)
  !write(*,*) "column: ", sig(:, 1)
  !write(*,*) "row   : ", sig(1, :)
  !write(*,*) "X(1, :): ", X(1, :)
  !write(*,*) "Y(:, 1): ", Y(:, 1)

  !call print_array(sig, nx, ny)
  !call print_array(X, nx, ny)
  !call print_array(Y, nx, ny)

  if (.FALSE.) then
    sig_fs = 0
    do i=1, N
      do j=1, N
        cnxx = cos(2 * pi * (i - 1) * X / Lx)
        snxx = sin(2 * pi * (i - 1) * X / Lx)
        cnyy = cos(2 * pi * (j - 1) * Y / Ly)
        snyy = sin(2 * pi * (j - 1) * Y / Ly)

        if ((i - 1) == 0 .and. (j - 1) == 0) then
          kappa = 1
        else if ((i - 1) == 0 .or. (j - 1) == 0) then
          kappa = 2
        else
          kappa = 4
        end if

        alpha(i, j) = sum(kappa * sig * cnxx * cnyy) * dA / (Lx * Ly)
        beta(i, j)  = sum(kappa * sig * cnxx * snyy) * dA / (Lx * Ly)
        gamma(i, j) = sum(kappa * sig * snxx * cnyy) * dA / (Lx * Ly)
        delta(i, j) = sum(kappa * sig * snxx * snyy) * dA / (Lx * Ly)

        sig_fs = sig_fs + alpha(i, j) * cnxx * cnyy
        sig_fs = sig_fs + beta(i, j)  * cnxx * snyy
        sig_fs = sig_fs + gamma(i, j) * snxx * cnyy
        sig_fs = sig_fs + delta(i, j) * snxx * snyy
      end do
    end do
  else
    sig_fs = 0
    do i=1, N
      do j=1, N
        do k=1, nx
          cnx(k) = cos(2 * pi * (i - 1) * X(1, k) / Lx)
          snx(k) = sin(2 * pi * (i - 1) * X(1, k) / Lx)
        end do

        do l=1, ny
          cny(l) = cos(2 * pi * (j - 1) * Y(l, 1) / Ly)
          sny(l) = sin(2 * pi * (j - 1) * Y(l, 1) / Ly)
        end do

        if ((i - 1) == 0 .and. (j - 1) == 0) then
          kappa = 1
        else if ((i - 1) == 0 .or. (j - 1) == 0) then
          kappa = 2
        else
          kappa = 4
        end if

        alpha_tmp = 0
        do k=1, nx
          do l=1, ny
            alpha_tmp = alpha_tmp + sig(l, k) * cnx(k) * cny(l)
          end do
        end do
        alpha_tmp = kappa * alpha_tmp * dA / (Lx * Ly)
        alpha(i, j) = alpha_tmp

        beta_tmp = 0
        do k=1, nx
          do l=1, ny
            beta_tmp = beta_tmp + sig(l, k) * cnx(k) * sny(l)
          end do
        end do
        beta_tmp = kappa * beta_tmp * dA / (Lx * Ly)
        beta(i, j) = beta_tmp

        gamma_tmp = 0
        do k=1, nx
          do l=1, ny
            gamma_tmp = gamma_tmp + sig(l, k) * snx(k) * cny(l)
          end do
        end do
        gamma_tmp = kappa * gamma_tmp * dA / (Lx * Ly)
        gamma(i, j) = gamma_tmp

        delta_tmp = 0
        do k=1, nx
          do l=1, ny
            delta_tmp = delta_tmp + sig(l, k) * snx(k) * sny(l)
          end do
        end do
        delta_tmp = kappa * delta_tmp * dA / (Lx * Ly)
        delta(i, j) = delta_tmp

        do k=1, nx
          do l=1, ny
            sig_fs(l, k) = sig_fs(l, k) + alpha_tmp * cnx(k) * cny(l) &
                                        + beta_tmp  * cnx(k) * sny(l) &
                                        + gamma_tmp * snx(k) * cny(l) &
                                        + delta_tmp * snx(k) * sny(l)
          end do
        end do
      end do
    end do
  end if

  !call print_array(alpha, N, N)
  !call print_array(beta, N, N)
  !call print_array(gamma, N, N)
  !call print_array(delta, N, N)
  !call print_array(sig_fs, nx, ny)

  call write_arrays(nx, ny, N, sig_fs, alpha, beta, gamma, delta)

  deallocate(sig)
  deallocate(sig_fs)
  deallocate(X)
  deallocate(Y)
  deallocate(snx)
  deallocate(cnx)
  deallocate(sny)
  deallocate(cny)

  deallocate(alpha)
  deallocate(beta)
  deallocate(gamma)
  deallocate(delta)
end program main
