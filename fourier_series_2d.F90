module constants  
  implicit none 

  double precision, parameter :: pi = 3.1415926536  
  double precision, parameter :: e = 2.7182818285 
end module constants 

subroutine print_arr(a, nx, ny)
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

end subroutine print_arr

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
  double precision, allocatable :: snx(:, :), cnx(:, :), sny(:, :), cny(:, :)
  double precision, allocatable :: alpha(:, :), beta(:, :), gamma(:, :), delta(:, :)

  double precision :: kappa
  integer :: i, j
  double precision :: dA

  namelist /settings/ nx, ny, N, Lx, Ly

  open(unit=1, file="data/settings.nml", status="OLD")
  read(unit=1, nml=settings)

  dA = (Lx / (nx - 1)) * (Ly / (ny - 1))

  allocate(sig(nx, ny))
  allocate(sig_fs(nx, ny))
  allocate(X(nx, ny))
  allocate(Y(nx, ny))
  allocate(snx(nx, ny))
  allocate(cnx(nx, ny))
  allocate(sny(nx, ny))
  allocate(cny(nx, ny))

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

  !call print_arr(sig, nx, ny)
  !call print_arr(X, nx, ny)
  !call print_arr(Y, nx, ny)

  sig_fs = 0
  do i=1, N
    do j=1, N
      cnx = cos(2 * pi * (i - 1) * X / Lx)
      snx = sin(2 * pi * (i - 1) * X / Lx)
      cny = cos(2 * pi * (j - 1) * Y / Ly)
      sny = sin(2 * pi * (j - 1) * Y / Ly)

      if ((i - 1) == 0 .and. (j - 1) == 0) then
        kappa = 1
      else if ((i - 1) == 0 .or. (j - 1) == 0) then
        kappa = 2
      else
        kappa = 4
      end if

      alpha(i, j) = sum(kappa * sig * cnx * cny) * dA / (Lx * Ly)
      beta(i, j)  = sum(kappa * sig * cnx * sny) * dA / (Lx * Ly)
      gamma(i, j) = sum(kappa * sig * snx * cny) * dA / (Lx * Ly)
      delta(i, j) = sum(kappa * sig * snx * sny) * dA / (Lx * Ly)

      sig_fs = sig_fs + alpha(i, j) * cnx * cny
      sig_fs = sig_fs + beta(i, j)  * cnx * sny
      sig_fs = sig_fs + gamma(i, j) * snx * cny
      sig_fs = sig_fs + delta(i, j) * snx * sny
    end do
  end do

  !call print_arr(alpha, N, N)
  !call print_arr(beta, N, N)
  !call print_arr(gamma, N, N)
  !call print_arr(delta, N, N)
  !call print_arr(sig_fs, nx, ny)

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
