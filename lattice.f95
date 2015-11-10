program lattice
    use latticesubs
    implicit none
    real(8),allocatable,dimension(:,:,:):: accepted_u
    real(8),allocatable,dimension(:,:)  :: x0, xt, d, U, U_i, d_temp, U_i_temp, AvgE, variance, total_lambda
    real(8)                             :: kT_init, kT_final, kT_delta, kT, box, lattice_parameter, initial_U, start, finish
    real(8)                             :: volume, density, sigma, step_size, acceptance_ratio_new
    real(8)                             :: d_a_r_n
    integer,allocatable,dimension(:)    :: accepted
    integer                             :: nkT, N, nMC, length, i, j, k, counter

    call random_seed()
    call cpu_time(start)

    ! Program inputs and parameters----------------------!
    ! ---------------------------------------------------!
    ! ---------------------------------------------------!
    kT_init  = 0.8
    kT_final = 2.0
    kT_delta = 1.2
    nkT = int((kT_final - kT_init) / kT_delta) + 2
    kT = kT_init
    N=32
    nMC=100000
    length = int((real(N)/4) ** (1.0/3.0))
    box = 2.0d0 ** (2.0d0/3.0d0)*dfloat(length)
    lattice_parameter = box/dfloat(length)
    volume = box ** 3.0d+0
    density = 0.8
    sigma = density ** (1.0/3.0)

    ! ---------------------------------------------------!
    ! ---------------------------------------------------!
    ! Allocation of dynamic arrays-----------------------!
    ! ---------------------------------------------------!
    ! ---------------------------------------------------!
    allocate(accepted_u(1,nMC,N))
    allocate(x0(N,3))
    allocate(xt(N,3))
    allocate(d(N,N))
    allocate(U(nMC,N))
    allocate(U_i(N,N))
    allocate(accepted(nkT))
    allocate(d_temp(N,N))
    allocate(U_i_temp(N,N))
    allocate(AvgE(nMC,N))
    allocate(total_lambda(nMC,N))
    allocate(variance(nMC,N))

    open(10,file='Potential_8.csv')
    open(20,file='Deviation_8.csv')
    open(40,file='AfterLambda_8.csv')
    open(50,file='AcceptedU_8.csv')
    open(60,file='Acceptance_Ratio_8.csv')


    ! ---------------------------------------------------!
    ! ---------------------------------------------------!
    do i=1,nkT
        step_size = 0.02
        accepted = 0
        counter = 1
        call initializefcc(box,length,x0)
        call initialpotential(initial_U,box,N,x0,d,U_i,d_temp,U_i_temp,sigma)
        xt = x0
        do j=1,nMC
            do k=1,N
                call moveparticle(x0,xt,box,k,N,step_size)
                call newpotential(U,d,xt,initial_U,box,j,k,N,U_i,d_temp,U_i_temp,sigma)
                call keepmove(accepted_u,U,x0,xt,initial_U,kt,accepted,i,j,k,N,d,U_i,U_i_temp,d_temp)
                counter = counter + 1
            end do
            if(mod(j,50) == 0 .or. j==1) call adjuststep(step_size,accepted,i,counter,j,acceptance_ratio_new, d_a_r_n,kT)
        end do
        accepted = 0
        initial_U = U(nMC,N)
        U = 0
        counter = 1
        do j=1,nMC
            do k=1,N
                call moveparticle(x0,xt,box,k,N,step_size)
                call newpotential(U,d,xt,initial_U,box,j,k,N,U_i,d_temp,U_i_temp,sigma)
                call keepmove(accepted_u,U,x0,xt,initial_U,kt,accepted,i,j,k,N,d,U_i,U_i_temp,d_temp)
                call disorder(total_lambda,xt,N,j,k,box)
                call averageE(avgE,U,j,k,N,counter)
                counter = counter + 1
            end do
            if(mod(j,50) == 0 .or. j==1) call adjuststep(step_size,accepted,i,counter,j,acceptance_ratio_new, d_a_r_n,kT)
        end do
        call deviation(AvgE,nMC,N,variance)

        ! Write results to disk
        do j=1,nMC, 50
            do k=1,N
                write(10,*) AvgE(j,k), kT
                write(20,*) variance(j,k), kT
                write(40,*) total_lambda(j,k), kT
            end do
        end do
        print*, kT

        kT = kT + kT_delta
    end do

    call cpu_time(finish)
    print*,"Time equals: ",finish-start," seconds."
end program
