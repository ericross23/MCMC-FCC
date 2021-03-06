module latticesubs
contains
    !-------------------------------------------------------------------------------------



    !------------------------------------------------------------------------------------
    subroutine initializefcc(box,length,x0)

        implicit none

        real(8),allocatable,dimension(:,:) :: x0
        real(8),dimension(4,3)             :: xo
        real(8)                            :: b, box
        integer                            :: i,j,k,k1,length,ip

        ! Setting up the four atoms in the unit cell
        b = box/dfloat(length)
        xo(:,:)=b/2.0d0
        xo(1,:)=0.0d0
        xo(2,3)=0.0d0
        xo(3,2)=0.0d0
        xo(4,1)=0.0d0

        ip=0
        do i=1,length
            do j=1,length
                do k=1,length
                    do k1=1,4
                        ip = ip+1
                        x0(ip,1) = xo(k1,1) + float(i-1) * b + 0.25 * b
                        x0(ip,2) = xo(k1,2) + float(j-1) * b + 0.25 * b
                        x0(ip,3) = xo(k1,3) + float(k-1) * b + 0.25 * b
                    end do
                end do
            end do
        end do

    end subroutine
    !-------------------------------------------------------------------------------




    !----------------------------------------------------------------------------------------------------
    subroutine initialpotential(initial_U,box,N,x0,d,U_i,d_temp,U_i_temp,sigma)
        implicit none
        real(8),allocatable,dimension(:,:) :: d, x0, U_i, d_temp, U_i_temp
        real(8)                            :: initial_U, box, dcut, sigma
        integer                            :: l, m, N

        initial_U = 0
        dcut = 0.98 * 0.5 * box
        do l=1,N-1
            do m=l+1,N
                call distance(l,m,d,x0,box,U_i,sigma)
                initial_U = initial_U + U_i(l,m)
            end do
        end do
        initial_U = initial_U * 4 - ((1/(dcut**12)) - (1/(dcut**6))) * 4
        d_temp = d
        U_i_temp = U_i

    end subroutine
    !-------------------------------------------------------------------------------




    !--------------------------------------------------------------------------------
    subroutine distance(l,m,d,x,box,U_i,sigma)
        implicit none
        real(8),allocatable,dimension(:,:) :: x, d, U_i
        real(8)                            :: dx, dy, dz, box, sigma
        integer                            :: l, m

        dx = x(l,1) - x(m,1)
        dy = x(l,2) - x(m,2)
        dz = x(l,3) - x(m,3)
        if(dx > box/2.0) dx = dx - box
        if(dx <= -box/2.0) dx = dx + box
        if(dy > box/2.0) dy = dy - box
        if(dy <= -box/2.0) dy = dy + box
        if(dz > box/2.0) dz = dz - box
        if(dz <= -box/2.0) dz = dz + box
        d(l,m) = sqrt(dx**2 + dy**2 + dz**2)
        U_i(l,m) = ((sigma/(d(l,m)**12)) - (sigma/(d(l,m)**6)))

    end subroutine
    !-------------------------------------------------------------------------------




    !-------------------------------------------------------------------------------
    subroutine moveparticle(x0,xt,box,k,N,step_size)
        implicit none
        real(8),allocatable,dimension(:,:) :: x0, xt, deltax
        real(8),parameter                  :: pi = 4 * atan(1.0_8)
        real(8)                            :: box, theta, phi, randnum, step_size
        integer                            :: k, N

        allocate(deltax(N,3))
        deltax(:,:) = 0

        call random_number(randnum)
        theta = acos((randnum * 2) - 1)
        call random_number(randnum)
        phi   = randnum * 2 * pi

        deltax(k,1) = step_size*sin(theta)*cos(phi)
        deltax(k,2) = step_size*sin(theta)*sin(phi)
        deltax(k,3) = step_size*cos(theta)

        xt(k,1) = x0(k,1) + deltax(k,1)
        xt(k,2) = x0(k,2) + deltax(k,2)
        xt(k,3) = x0(k,3) + deltax(k,3)

        call pbc(xt,box,N)

    end subroutine

    subroutine pbc(xt,box,N)
        implicit none
        real(8),allocatable,dimension(:,:) :: xt
        real(8)                            :: box
        integer                            :: N, i, j

        do i = 1, N
            do j = 1, 3
                if (xt(i,j) > box) xt(i,j) = xt(i,j) - box
                if (xt(i,j) < 0.0) xt(i,j) = xt(i,j) + box
            end do
        end do

    end subroutine


    subroutine newpotential(U,d,xt,initial_U,box,j,k,N,U_i,d_temp,U_i_temp,sigma)
        implicit none
        real(8),allocatable,dimension(:,:) :: U, d, d_temp, xt, U_i, U_i_temp
        real(8)                            :: initial_U, box, dcut, sigma
        integer,intent(in)                 :: j, k, N
        integer                            :: m

        if(j==1 .and. k==1) then
            U(j,k) = initial_U
        else if(k==1) then
            U(j,k) = U(j-1,N)
        else
            U(j,k) = U(j,k-1)
        end if

        U_i_temp = U_i
        d_temp = d
        dcut = 0.98 * 0.5 * box
        U(j,k) = U(j,k)/4 + ((1/(dcut**12)) - (1/(dcut**6)))
        do m=k+1,N
            call distance(k,m,d,xt,box,U_i,sigma)
            U(j,k) = U(j,k) + U_i(k,m) - U_i_temp(k,m)
        end do
        U(j,k) = U(j,k) * 4 - ((1/(dcut**12)) - (1/(dcut**6))) * 4
    end subroutine


    subroutine keepmove(accepted_u,U,x0,xt,initial_U,kt,accepted,i,j,k,N,d,U_i,U_i_temp,d_temp)
        implicit none
        real(8),allocatable,dimension(:,:,:):: accepted_u
        real(8),allocatable,dimension(:,:)  :: U, x0, xt, d, d_temp, U_i, U_i_temp
        real(8)                             :: initial_U, test_U, ptrial, kT, randnum
        integer,allocatable,dimension(:)    :: accepted
        integer                             :: i, j, k, N, flag

        if(j==1 .and. k==1) then
            test_U = initial_U
        else if(k==1) then
            test_U = U(j-1,N)
        else
            test_U = U(j,k-1)
        end if

        if(U(j,k) < test_U) then
            x0(k,1) = xt(k,1)
            x0(k,2) = xt(k,2)
            x0(k,3) = xt(k,3)
            flag = 0
            accepted(i) = accepted(i) + 1
            if(i==3) accepted_u(1,j,k) = U(j,k)
        else
            ptrial = EXP(-(U(j,k) - test_U)/kT)
            call random_number(randnum)
            if (randnum<ptrial) then
                x0(k,1) = xt(k,1)
                x0(k,2) = xt(k,2)
                x0(k,3) = xt(k,3)
                flag = 0
                accepted(i) = accepted(i) + 1
                if(i==3) accepted_u(1,j,k) = U(j,k)
            else
                U(j,k) = test_U
                xt(k,1) = x0(k,1)
                xt(k,2) = x0(k,2)
                xt(k,3) = x0(k,3)
                if(i==3) accepted_u(1,j,k) = 0
                U_i(:,:) = U_i_temp(:,:)
                d(:,:) = d_temp(:,:)
            end if
        end if
    end subroutine


    subroutine averageE(avgE,U,j,k,N,counter,kT,nMC)
        implicit none

        real(8),allocatable,dimension(:,:) :: avgE, U
        real(8)                            :: kT
        integer                            :: j, k, counter, N, nMC

        if(j==nMC/2 .and. k==1) then
            avgE(j,k) = U(j,k)
        else if(k==1) then
            avgE(j,k) = (avgE(j-50,N) * (counter - 1) + U(j,k)) / counter
        else
            avgE(j,k) = (avgE(j,k-1) * (counter - 1) + U(j,k)) / counter
        end if

        write(10,*) avgE(j,k), kT

    end subroutine

!    !--------------------------------------------------------------------------------
!
!
!
!    !--------------------------------------------------------------------------------
    subroutine deviation(AvgE,nMC,N,variance,kT)
        implicit none

        integer                            :: i, j, nMC, N, SampleSize
        real(8),allocatable,dimension(:,:) :: AvgE, variance
        real(8)                            :: M2, delta, mean, kT

        SampleSize = 0
        mean = 0.0
        M2 = 0.0

        do i=nMC/2, nMC, 50
            do j=1, N
                SampleSize = SampleSize + 1
                delta = AvgE(i,j) - mean
                mean = mean + delta/SampleSize
                M2 = M2 + delta * (AvgE(i,j) - mean)
                if(SampleSize < 2) then
                    variance(i,j) = 0.0
                else
                    variance(i,j) = M2 / (SampleSize - 1)
                end if
                write(20,*) variance(i,j), kT
            end do
        end do

    end subroutine

    subroutine disorder(total_lambda,xt,N,j,k,box,kT)

        implicit none

        real(8),allocatable,dimension(:,:) :: xt, total_lambda
        real(8),dimension(3)               :: lambda
        real(8)                            :: box,kT
        real(8),parameter                  :: pi = 4 * atan(1.0_8)
        integer                            :: i, j, k, l, N

        lambda(:) = 0
        do i=1,3
            do l=1,N
                lambda(i) = lambda(i) + cos((4*pi*xt(l,i))/box)
            end do
        end do

        do l=1,3
            lambda(l) = lambda(l) / N
        end do

        total_lambda(j,k) = sum(lambda) / 3.0
        write(40,*) total_lambda(j,k), kT

    end subroutine

    subroutine adjuststep(step_size,accepted,i,counter,j,acceptance_ratio_new,delta_acceptance_ratio_new,kT)
        implicit none

        real(8)                          :: step_size, acceptance_ratio_old, acceptance_ratio_new
        real(8)                          :: delta_acceptance_ratio_old, delta_acceptance_ratio_new, kT
        integer,allocatable,dimension(:) :: accepted
        integer                          :: i, counter, j


        if(j==1) then
            acceptance_ratio_old = real(accepted(i)) / real(counter)
        else
            acceptance_ratio_old = acceptance_ratio_new
        end if

        acceptance_ratio_new = real(accepted(i)) / real(counter)

        if(j==1) then
            delta_acceptance_ratio_old = 0
            delta_acceptance_ratio_new = 0
        else
            delta_acceptance_ratio_old = delta_acceptance_ratio_new
        end if
        delta_acceptance_ratio_new = acceptance_ratio_new - acceptance_ratio_old

        if(acceptance_ratio_new > 0.60) then
            if(acceptance_ratio_new - acceptance_ratio_old > 0) then
                if(delta_acceptance_ratio_new - delta_acceptance_ratio_old > 0) then
                    step_size = step_size * 1.20
                end if
            end if
        else if(acceptance_ratio_new < 0.40) then
            if(acceptance_ratio_new - acceptance_ratio_old < 0) then
                if(delta_acceptance_ratio_new - delta_acceptance_ratio_old < 0) then
                    step_size = step_size * 0.80
                end if
            end if
        end if

        if(step_size < 1.0d-5) step_size = 0.02

        write(60,*) acceptance_ratio_new, step_size, kT
    end subroutine

    subroutine gr(switch,box,N,xt,density,kT,g,nhis,ngr,delg)
        implicit none

        integer                            :: switch, i, j, ngr, nhis, N, ig
        real(8),parameter                  :: pi = 4 * atan(1.0_8)
        real(8),allocatable,dimension(:)   :: g
        real(8),allocatable,dimension(:,:) :: xt
        real(8)                            :: delg, dx, dy, dz, r, density, vb, kT, box, nid

        if(switch == 0) then
            ngr = 0
            delg = box/(2*nhis)
            do i=0,nhis
                g(i) = 0
            end do
        else if(switch == 1) then
            ngr = ngr + 1
            do i=1, N-1
                do j=i+1, N
                    dx = xt(i,1) - xt(j,1)
                    dx = dx - box*nint(dx/box)
                    dy = xt(i,2) - xt(j,2)
                    dy = dy - box*nint(dy/box)
                    dz = xt(i,3) - xt(j,3)
                    dz = dz - box*nint(dz/box)
                    r = sqrt(dx**2 + dy**2 + dz**2)
                    if (r < box/2) then
                        ig = int(r/delg)
                        g(ig) = g(ig) + 2
                    end if
                end do
            end do
        else if(switch == 2) then
            do i=1, nhis
                r = delg * (i+0.5)
                vb = ((i+1)**3-i**3) * delg**3
                nid = (4./3) * pi * vb * density
                g(i) = g(i) / (ngr * N * nid)
                write(70,*) r, g(i), kT
            end do
            print*, ngr
        end if
        return

    end subroutine

end module
