module constants  
  implicit none 

  real, parameter :: pi = 3.1415926536  
  real, parameter :: e = 2.7182818285 
end module constants 

subroutine print_arr(a, nx, ny)
  implicit none
  integer, intent(in) :: nx, ny
  real, intent(in) :: a(nx, ny)
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

subroutine read_sig(nx, ny, fname, sig, X, Y)
  implicit none
  integer, intent(in) :: nx, ny
  character*7, intent(in) :: fname
  real, intent(out) :: sig(nx, ny)
  real, intent(out) :: X(nx, ny), Y(nx, ny)
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
end subroutine read_sig

subroutine write_sig_fs(nx, ny, sig_fs)
  implicit none
  integer, intent(in) :: nx, ny
  real, intent(in) :: sig_fs(nx, ny)
  integer :: i

  open (unit=30,form="unformatted",file="data/sig_fs.bin",action="write")
  do i=1,nx
     write(30) sig_fs(i,:) ! write one row at a time
  end do
  close(30)
end subroutine write_sig_fs

subroutine sin_nx(nx, ny, Lx, Ly, X, n, snx)
  use constants, only: pi
  implicit none
  integer, intent(in) :: nx, ny
  real, intent(in) :: X(nx, ny)
  real, intent(in) :: Lx, Ly, n
  real, intent(out) :: snx(nx, ny)

  snx = sin(2 * pi * n * X / Lx)
end subroutine sin_nx

subroutine cos_nx(nx, ny, Lx, Ly, X, n, cnx)
  use constants, only: pi
  implicit none
  integer, intent(in) :: nx, ny
  real, intent(in) :: X(nx, ny)
  real, intent(in) :: Lx, Ly, n
  real, intent(out) :: cnx(nx, ny)

  cnx = cos(2 * pi * n * X / Lx)
end subroutine cos_nx

subroutine sin_ny(nx, ny, Lx, Ly, Y, n, sny)
  use constants, only: pi
  implicit none
  integer, intent(in) :: nx, ny
  real, intent(in) :: Y(nx, ny)
  real, intent(in) :: Lx, Ly, n
  real, intent(out) :: sny(nx, ny)

  sny = sin(2 * pi * n * Y / Lx)
end subroutine sin_ny

subroutine cos_ny(nx, ny, Lx, Ly, Y, n, cny)
  use constants, only: pi
  implicit none
  integer, intent(in) :: nx, ny
  real, intent(in) :: Y(nx, ny)
  real, intent(in) :: Lx, Ly, n
  real, intent(out) :: cny(nx, ny)

  cny = cos(2 * pi * n * Y / Lx)
end subroutine cos_ny

program main
  use iso_fortran_env
  implicit none

  integer, parameter :: nx=100, ny=100, N=100
  integer :: i, j
  real, parameter :: Lx=20, Ly=20
  real, parameter :: dA=0.04
  real :: sig(nx, ny), sig_fs(nx, ny)
  real :: X(nx, ny), Y(nx, ny)
  real :: snx(nx, ny), cnx(nx, ny), sny(nx, ny), cny(nx, ny)
  real :: alpha(N, N), beta(N, N), gamma(N, N), delta(N, N)
  real :: kappa

  write(*,*)"Running fourier_seried_2d"

  call read_sig(nx, ny, "sig.bin", sig, X, Y)

  ! Sanity check. Rem. fortran using 1 based indexing.
  write(*,'(10F10.8)') sig(6, 6)

  !call print_arr(sig, nx, ny)
  !call print_arr(X, nx, ny)
  !call print_arr(Y, nx, ny)

  do i=1, N
    do j=1, N
      call sin_nx(nx, ny, Lx, Ly, X, 1. * i, snx)
      call cos_nx(nx, ny, Lx, Ly, X, 1. * i, cnx)
      call sin_ny(nx, ny, Lx, Ly, Y, 1. * j, sny)
      call cos_ny(nx, ny, Lx, Ly, Y, 1. * j, cny)
      if (i == 0 .and. j == 0) then
        kappa = 1
      else if (i == 0 .or. j == 0) then
        kappa = 2
      else
        kappa = 4
      end if

      alpha(i, j) = sum(kappa * sig * cnx * cny) * dA / (Lx * Ly)
      beta(i, j) = sum(kappa * sig * cnx * sny) * dA / (Lx * Ly)
      gamma(i, j) = sum(kappa * sig * snx * cny) * dA / (Lx * Ly)
      delta(i, j) = sum(kappa * sig * snx * sny) * dA / (Lx * Ly)
    end do
  end do

  call print_arr(alpha, N, N)
  call print_arr(beta, N, N)
  call print_arr(gamma, N, N)
  call print_arr(delta, N, N)

  sig_fs = 0
  do i=1, N
    do j=1, N
      call sin_nx(nx, ny, Lx, Ly, X, 1. * i, snx)
      call cos_nx(nx, ny, Lx, Ly, X, 1. * i, cnx)
      call sin_ny(nx, ny, Lx, Ly, Y, 1. * j, sny)
      call cos_ny(nx, ny, Lx, Ly, Y, 1. * j, cny)

      sig_fs = sig_fs + alpha(i, j) * cnx * cny
      sig_fs = sig_fs + beta(i, j) * cnx * sny
      sig_fs = sig_fs + gamma(i, j) * snx * cny
      sig_fs = sig_fs + delta(i, j) * snx * sny
    end do
  end do

  call print_arr(sig_fs, nx, ny)
  call write_sig_fs(nx, ny, sig_fs)
end program main
