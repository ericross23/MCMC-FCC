program lattice
    use latticesubs
    implicit none
    real(8),allocatable,dimension(:,:,:):: accepted_u
    real(8),allocatable,dimension(:,:)  :: x0, xt, d, U, U_i, d_temp, U_i_temp, AvgE, a_k, q_k, total_lambda
    real(8)                             :: kT_init, kT_final, kT_delta, kT, box, lattice_parameter, initial_U, start, finish
    real(8)                             :: volume, density, sigma, step_size, acceptance_ratio_new
    real(8)                             :: d_a_r_n
    integer,allocatable,dimension(:)    :: accepted
    integer                             :: nkT, N, nMC, length, i, j, k, counter, Uflag, Dflag

    call random_seed()
    call cpu_time(start)

    ! Program inputs and parameters----------------------!
    ! ---------------------------------------------------!
    ! ---------------------------------------------------!
    kT_init  = 0.1
    kT_final = 2.1
    kT_delta = 0.5
    nkT = int((kT_final - kT_init) / kT_delta) + 2
    kT = kT_init
    N=256
    nMC=50000
    length = int((real(N)/4) ** (1.0/3.0))
    box = 2.0d0 ** (2.0d0/3.0d0)*dfloat(length)
    lattice_parameter = box/dfloat(length)
    volume = box ** 3.0d+0
    density = 0.4
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
    allocate(a_k(nMC,N))
    allocate(q_k(nMC,N))
    allocate(total_lambda(nMC,N))

    open(10,file='Potential_4.csv')
    open(20,file='Deviation_4.csv')
    open(40,file='AfterLambda_4.csv')
    open(50,file='AcceptedU_4.csv')
    open(60,file='Acceptance_Ratio_4.csv')


    ! ---------------------------------------------------!
    ! ---------------------------------------------------!
    do i=1,nkT
        step_size = 0.02
        Uflag = 0
        Dflag = 0
        accepted = 0
        counter = 1
        call initializefcc(box,length,x0)
        call initialpotential(initial_U,box,N,x0,d,U_i,d_temp,U_i_temp,sigma)
        xt = x0
        do j=1,nMC
            do k=1,N
                call moveparticle(x0,xt,box,k,N,step_size)
                call newpotential(U,d,xt,initial_U,box,j,k,N,U_i,d_temp,U_i_temp,sigma)
                call keepmove(accepted_u,U,x0,xt,initial_U,kt,accepted,i,j,k,N,d,U_i,U_i_temp,d_temp,Dflag,q_k,a_k,counter)
                call averageE(avgE,U,j,k,Uflag,N,kT,counter)
                call disorder(total_lambda,xt,N,j,k,box,Dflag,kT)
                counter = counter + 1
            end do
            if(mod(j,50) == 0 .or. j==1) call adjuststep(step_size,accepted,i,counter,j,acceptance_ratio_new, d_a_r_n,kT)
        end do
        accepted = 0
        initial_U = U(nMC,N)
        U = 0
        Uflag = 1
        Dflag = 1
        counter = 1
        do j=1,nMC
            do k=1,N
                call moveparticle(x0,xt,box,k,N,step_size)
                call newpotential(U,d,xt,initial_U,box,j,k,N,U_i,d_temp,U_i_temp,sigma)
                call keepmove(accepted_u,U,x0,xt,initial_U,kt,accepted,i,j,k,N,d,U_i,U_i_temp,d_temp,Dflag,q_k,a_k,counter)
                call averageE(avgE,U,j,k,Uflag,N,kT,counter)
                call disorder(total_lambda,xt,N,j,k,box,Dflag,kT)
                counter = counter + 1
            end do
            if(mod(j,50) == 0 .or. j==1) call adjuststep(step_size,accepted,i,counter,j,acceptance_ratio_new, d_a_r_n,kT)
        end do
        AvgE = 0
        print*,"After Equilibrium ", "kT: ", kT, "Accepted: ", accepted(i), "Step Size: ", step_size
        kT = kT + kT_delta


    end do

    do j=1,nMC
        do k=1,N
            write(50,*) accepted_u(1,j,k)
        end do
    end do
    call cpu_time(finish)
    print*,"Time equals: ",finish-start," seconds."
end program
