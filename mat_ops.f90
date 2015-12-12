subroutine get_pair_dist(mat, N, arr)
! generate pairs of distances from an xyz matrix mat of size (N, 3)
! 22/09/15
    integer, parameter :: k = selected_int_kind(16)
    integer(kind=k), intent(in) :: N
    real(8), intent(in), dimension(N, 3) :: mat
    real(8), intent(out), dimension(N*(N-1)/2) :: arr
    integer :: i, j, cnt
    ! f2py depend(N) mat, arr

    cnt = 1
    do i = 1, N
        do j = i+1, N
            arr(cnt) = sqrt(sum((mat(i, :) - mat(j, :))**2))
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
! Calculate and return order parameter (according to Goyal PhD thesis)
! 12/10/15
    integer, intent(in) :: N
    real, intent(in) :: rc
    real, intent(in), dimension(N, 4) :: xyz
    real, intent(out) :: phi
    integer, dimension(N) :: n1, n2
    integer :: i, j
    
    n1 = 0
    n2 = 0
    phi = 0.0

    do i = 1, N
        do j = i+1, N
            dist = sqrt((xyz(i, 2)-xyz(j, 2))**2 + (xyz(i, 3)-xyz(j, 3))**2 + &
                        (xyz(i, 4)-xyz(j, 4))**2)
            if (dist < rc) then
                if (xyz(j, 1) == 1) then   ! neighbour of i-th particle is type 1
                    n1(i) = n1(i) + 1
                endif
                if (xyz(j, 1) == 2) then   ! neighbour of i-th particle is type 2
                    n2(i) = n2(i) + 1
                endif
                if (xyz(i, 1) == 1) then   ! vice versa
                    n1(j) = n1(j) + 1
                endif
                if (xyz(i, 1) == 2) then   ! vice versa
                    n2(j) = n2(j) + 1
                endif
            endif
        enddo
    enddo
    
    do i = 1, N
        phi = phi + float(n1(i) - n2(i))**2/(n1(i) + n2(i))**2
!        if (n1(i) == 0) then 
!            print *, "n1", i
!        endif
!        if (n2(i) == 0) then
!            print *, "n2", i
!        endif 
    enddo
!    print *, phi
    phi = phi/N
end subroutine


subroutine get_local_op2(xyz, N, rc, n1, n2)
! Return two nearest neighbour occupation arrays n1, n2 to later calculate 
! the order parameter (according to Goyal PhD thesis)
! 12/10/15
    integer, intent(in) :: N
    real, intent(in) :: rc
    real, intent(in), dimension(N, 4) :: xyz
    integer, intent(out), dimension(N) :: n1, n2
    integer :: i, j
    
    n1 = 0
    n2 = 0
    phi = 0.0

    do i = 1, N
        do j = i+1, N
            dist = sqrt((xyz(i, 2)-xyz(j, 2))**2 + (xyz(i, 3)-xyz(j, 3))**2 + &
                        (xyz(i, 4)-xyz(j, 4))**2)
            if (dist < rc) then
                if (xyz(j, 1) == 1) then   ! neighbour of i-th particle is type 1
                    n1(i) = n1(i) + 1
                endif
                if (xyz(j, 1) == 2) then   ! neighbour of i-th particle is type 2
                    n2(i) = n2(i) + 1
                endif
                if (xyz(i, 1) == 1) then   ! vice versa
                    n1(j) = n1(j) + 1
                endif
                if (xyz(i, 1) == 2) then   ! vice versa
                    n2(j) = n2(j) + 1
                endif
            endif
        enddo
    enddo
end subroutine

