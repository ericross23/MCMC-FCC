program lattice
    use latticesubs
    implicit none
    real(8),allocatable,dimension(:,:,:):: accepted_u
    real(8),allocatable,dimension(:,:)  :: x0, xt, d, U, U_i, d_temp, U_i_temp, AvgE, variance, total_lambda
    real(8),allocatable,dimension(:)    :: g
    real(8)                             :: kT_init, kT_final, kT_delta, kT, box, lattice_parameter, initial_U, start, finish
    real(8)                             :: volume, density, sigma, step_size, acceptance_ratio_new
    real(8)                             :: d_a_r_n, delg
    integer,allocatable,dimension(:)    :: accepted
    integer                             :: nkT, N, nMC, length, i, j, k, counter, nhis, switch, ngr

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
    N=256
    nMC=10000
    length = int((real(N)/4) ** (1.0/3.0))
    density = 0.8
    box = (N/density) ** (1.0/3.0)
    lattice_parameter = box/dfloat(length)
    volume = box ** 3.0d+0

    sigma = ((density * volume)/N) ** (1.0/3.0)
    nhis = 80
    print*, box

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
    allocate(g(0:nhis))

    open(10,file='AveragePotential_8.csv')
    open(20,file='Deviation_8.csv')
    open(30,file='Potential_8.csv')
    open(40,file='AfterLambda_8.csv')
    open(50,file='AcceptedU_8.csv')
    open(60,file='Acceptance_Ratio_8.csv')
    open(70,file='RadialDistribution_8.csv')

    ! ---------------------------------------------------!
    ! ---------------------------------------------------!



    do i=1,nkT

        switch = 0
        call gr(switch,box,N,xt,density,kT,g,nhis,ngr,delg)
        switch = 1

        step_size = 0.02
        accepted = 0
        counter = 1
        AvgE = 0
        U = 0
        call initializefcc(box,length,x0)
        call initialpotential(initial_U,box,N,x0,d,U_i,d_temp,U_i_temp,sigma)
        xt = x0
        do j=1,nMC
            do k=1,N
                if(j == nMC/2) then
                    counter = 1
                    accepted = 0
                end if
                call moveparticle(x0,xt,box,k,N,step_size)
                call newpotential(U,d,xt,initial_U,box,j,k,N,U_i,d_temp,U_i_temp,sigma)
                call keepmove(accepted_u,U,x0,xt,initial_U,kt,accepted,i,j,k,N,d,U_i,U_i_temp,d_temp)
                if(j >= nMC/2 .and. mod(j,50) == 0) call disorder(total_lambda,xt,N,j,k,box,kT)
                if(j >= nMC/2 .and. mod(j,50) == 0) call averageE(avgE,U,j,k,N,counter,kT,nMC)
                counter = counter + 1
            end do
            call gr(switch,box,N,xt,density,kT,g,nhis,ngr,delg)
            if(mod(j,50) == 0) call adjuststep(step_size,accepted,i,counter,j,acceptance_ratio_new,d_a_r_n,kT)
        end do

        call deviation(AvgE,nMC,N,variance,kT)
        switch = 2
        call gr(switch,box,N,xt,density,kT,g,nhis,ngr,delg)
        print*, kT

        kT = kT + kT_delta
    end do

    call cpu_time(finish)
    print*,"Time equals: ",finish-start," seconds."
end program
