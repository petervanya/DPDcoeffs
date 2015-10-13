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


subroutine get_pair_dist2(mat, N, arr)
! generate pairs of cistances together with the particle types
! 11/10/15
    integer, intent(in) :: N
    real, intent(in), dimension(N, 4) :: mat
    real, intent(out), dimension(N*(N-1)/2, 3) :: arr
    integer :: i, j, cnt

    cnt = 1
    do i = 1, N
        do j = i+1, N
            arr(cnt, 1) = mat(i, 1)
            arr(cnt, 2) = mat(j, 1)
            arr(cnt, 3) = sqrt((mat(i, 2)-mat(j, 2))**2 + (mat(i, 3)-mat(j, 3))**2 + &
                          (mat(i, 4)-mat(j, 4))**2)
            cnt = cnt + 1
        enddo
    enddo
end subroutine

subroutine get_local_op(xyz, N, rc, phi)
!
! 12/10/15
    integer, intent(in) :: N
    real, intent(in) :: rc
    real, intent(in), dimension(N, 4) :: xyz
    real, intent(out) :: phi
    integer, dimension(N-1) :: na, nb
    integer :: i, j

    na = 0
    nb = 0
    phi = 0.0

    do i = 1, N
        do j = i+1, N
            dist = sqrt((xyz(i, 2)-xyz(j, 2))**2 + (xyz(i, 3)-xyz(j, 3))**2 + &
                          (xyz(i, 4)-xyz(j, 4))**2)
            if (dist < rc) then
                if (xyz(i, 1) == 1) then
                    na(i) = na(i) + 1
                    na(j) = na(j) + 1
                else if (xyz(i, 1) == 2) then
                    nb(i) = nb(i) + 1
                    nb(j) = nb(j) + 1
                endif
            endif
        enddo
    enddo
    
    do i = 1, N-1
        phi = phi + float(na(i) - nb(i))**2/(na(i) + nb(i))**2
    enddo
    phi = phi/N
end subroutine


