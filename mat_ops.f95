subroutine get_pair_dist(mat, N, arr)
! generate pairs of distances from an xyz matrix mat of size (N, 3)
! 22/09/15
  integer, intent(in) :: N
  real, intent(in), dimension(N, 3) :: mat
  real, intent(out), dimension(N*(N-1)/2) :: arr
  integer :: i, j, cnt

  cnt = 1
  do i = 1, N
    do j = i+1, N
      arr(cnt) = sqrt((mat(i, 1)-mat(j, 1))**2 + (mat(i, 2)-mat(j, 2))**2 + &
                      (mat(i, 3)-mat(j, 3))**2)
      cnt = cnt + 1
    enddo
  enddo
end subroutine
