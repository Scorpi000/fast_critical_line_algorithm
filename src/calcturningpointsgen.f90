! g95 -c calcturningpointsgen.f90


SUBROUTINE calcturningpointsgen(mu_param, sigma_param, nn, mu_tps, sig_tps, &
                        weight_tps, lambda_tps, n_tps, error_flag, max_tps, &
                        lower_param, upper_param)

  implicit none

  ! using allocatable arrays:
  integer, intent(in)                  :: nn
  integer                              :: n_tpx
  integer                              :: kx
  integer, intent(out)                 :: error_flag
  integer, intent(in)                  :: max_tps

  real(8), dimension(:,:), allocatable :: Ai

  real(8), dimension(nn) :: lower, upper
  real(8), dimension(nn), intent(in) :: lower_param, upper_param
  real(8), dimension(nn) :: weight_max, weight_min
  integer, dimension(nn) :: ind_mu

  real(8), dimension(nn), intent(in)    :: mu_param
  real(8), dimension(nn,nn), intent(in) :: sigma_param
  real(8), dimension(nn)    :: mu
  real(8), dimension(nn,nn) :: sigma


  ! the turning points
  real(8), dimension(max_tps), intent(out)     :: mu_tps, sig_tps, lambda_tps
  real(8), dimension(nn,max_tps), intent(out)  :: weight_tps
  
  integer :: print_on = 0
  
  !----------------------------------

  integer                   :: i_min, i_max
  real(8)                   :: mu_max, mu_min, sig, lambda_min, lambda_max

  ! the actual number of turning points:
  integer, intent(out)      :: n_tps


  ! these declarations are needed only for testing:
  integer :: i, j, n1, i0
  real(8) :: rnd, rnd1, rnd2, aux, aux1

  real(8) :: auxv1(1), auxv2(1)
  real(8) :: rand_shift = 0.0d0   ! default value
  !-------------------
  lower = lower_param
  upper = upper_param
  mu = mu_param
  sigma = sigma_param

  error_flag = 0

  ! for the version with allocations:
!  nn = 500
  kx = nn
  n_tpx = max_tps  ! 10*nn
  allocate( Ai(kx,kx) )
  !----------------

  if( print_on /= 0 ) then
    print *,'nn=', nn
  END IF

  !----------------
!!$  !!! for testing: set mu and sigma:
!!$  rand_shift = 0.5d0
!!$  print *,'rand_shift=', rand_shift
!!$  do i = 1, nn
!!$     call random_number(rnd)
!!$     mu(i) = rnd + 0.01d0*i
!!$  end do
!!$  call timing('')
!!$  call matgen(nn, sigma, rand_shift)
!!$  call timing('generate sigma')

  ! create lower and upper bounds
  !~ call create_bounds
  call check_bounds

  ! to have only zero lower bounds:
  !~ lower = 0.0d0
  !~ upper = 1.1d0
     

  !----------------

!!$  call write_data
!!$  stop

  !~ if( print_on /= 0 ) then
    !~ print *,'read in data'
  !~ END IF
  !~ call timing('')
  !~ call read_data
  !~ call timing('Read_data')

!  call timing('')
  call index_sort(nn, mu, ind_mu)
!  call timing('index_sort')
!  print *,ind_mu(1:5)
!  print *,mu(ind_mu(1:3))

!  call timing('')
  call index_sort0(nn, mu, ind_mu)
!  call timing('index_sort0')
!  print *,ind_mu(1:5)
!  print *,mu(ind_mu(1:3))

!  stop

  call timing('')
  call turning(n_tps, lambda_tps, mu_tps, sig_tps, weight_tps)
  call timing('turning')

!  print *,'check weight_tps:'
! uncomment the following lines to run weight checking in fortran
  !~ do i = 1, n_tps
     !~ call check_weight(weight_tps(:,i))
     !~ if( error_flag /= 0 ) then
       !~ return
     !~ end if
  !~ end do


  if( print_on /= 0 ) then
    print *,'A: lambda_min=', lambda_min
    print *,'A: lambda_max=', lambda_max
  !  print *,'sum(lower)=', sum(lower)
  !  print *,'sum(upper)=', sum(upper)

    print *,'---------------------------------'
    print *,'Number of turning points: n_tps=', n_tps

    n1 = min(10, nn)
    if( n1 < nn)  then
       do i = 1, n_tps

          if( i > 10 .and. i < n_tps-10) then
             if(i==11) print *,'.....................'
             cycle
          end if

          print '(i4,f10.2,2f7.3,2x,10f6.3," ...")', i, lambda_tps(i), &
               mu_tps(i), sig_tps(i),                                  &
  !             weight_tps(1:n1,i)
               (weight_tps(1:n1,i)-lower(1:n1)) / (upper(1:n1)-lower(1:n1))
       end do
    else
       do i = 1, n_tps
          print '(i4,f10.2,2f7.3,2x,10f6.3)', i, lambda_tps(i), &
               mu_tps(i), sig_tps(i),                          &
  !             weight_tps(:,i)
               (weight_tps(:,i)-lower) / (upper-lower)
       end do
    end if

    print *,'B: lambda_min= ', lambda_min
    print *,'B: lambda_last=', lambda_tps(n_tps)
  end if

  ! check the violations on the bounds
  aux = 0.
  do i = 1, n_tps
     auxv1 = minval(weight_tps(:,i)-lower)
     auxv2 = minval(upper - weight_tps(:,i))
     aux = min(aux, auxv1(1), auxv2(1))
  end do
  if( print_on /= 0 ) then
    print *,'max. violation of bounds:', -aux
    print *,'C: lambda_min= ', lambda_min
    print *,'C: lambda_last=', lambda_tps(n_tps)
  end if

  deallocate(Ai)

contains

  subroutine turning(n_tp, lambda_tp, mu_tp, sig_tp, weight_tp)
    implicit none

    logical, parameter :: correct_weights = .true.
    logical, parameter :: check_inversion = .false.
    
    ! total number of turning points 
    integer                    :: n_tp

    ! values at the turning points:
    real(8), dimension(n_tpx)     :: mu_tp, sig_tp, lambda_tp
    real(8), dimension(nn,n_tpx)  :: weight_tp

    real(8), dimension(nn) :: w, weight, Sinv_eS, Sinv_muS, Sinv_SRwR
    integer, dimension(nn) :: ivS, ivRu, ivRl, indv, indv_temp
    integer, dimension(nn) :: ivS_min, ivRu_min, ivRl_min
    integer :: kk_max, i, i_aux, i_new, kk, count, i1_new, i2_new
    real(8) :: lambda_current, lambda1_new, lambda2_new, lambda_new, muP
    real(8) :: aux_l, aux_r
    integer :: dir1, i_last
    real(8), parameter :: epsx = 1.0d-5

    ! for accelerating get_muP:
    real(8), dimension(nn) :: sigw0
    real(8)                :: w0sigw0

    ! ii=indRl(1:kl) gives the index of the weight at the lower bound 
    integer, dimension(nn) :: indRl, indRu
    integer                :: kl, ku
    real(8), parameter     :: lambda_aux = 1.0d20

    count = 1
    kk_max = 0

    ! determine the end point, used only to check the program
    call get_starting_weight_min(ivS_min, ivRl_min, ivRu_min, i_min, &
         weight_min)
    call get_lambda_min(i_min, weight_min, i_aux, lambda_min)

    ! determine the starting point 
    ! (this will be repeated with general routine)
    call get_starting_weight_max(ivS, ivRl, ivRu, i_max, weight_max)
    if( error_flag /= 0 ) then
      return
    end if
    call get_lambda_max(i_max, weight_max, i_aux, lambda_max)

    ! get lambda_max with the general routine:
    ! this is needed to initialize some variables (not only to check!)
    Ai(1,1) = 1.0d0/sigma(i_max,i_max)
    Sinv_eS(1) = Ai(1,1)
    Sinv_muS(1) = Ai(1,1)*mu(i_max)

    aux_l = sum(sigma(:,i_max)*ivRl*lower) 
    aux_r = sum(sigma(:,i_max)*ivRu*upper) 
    Sinv_SRwR(1) = Ai(1,1)*(aux_l + aux_r)

    indv(1) = i_max
    lambda_current = lambda_aux

    call get_indR(ivRl, ivRu, kl, indRl, ku, indRu)

    call get_lambda2(ivS, ivRl, ivRu, 1, indv(1:1),    &
         kl, indRl(1:kl), ku, indRu(1:ku),             &
         Sinv_eS(1:1), Sinv_muS(1:1), Sinv_SRwR(1:1),  &
         lambda_current, i_new, lambda_new)

    lambda_max = lambda_new
    lambda_current = lambda_max
    i_last = i_new

    kk = 1
    indv(kk+1) = i_new
    ivS(i_new) = 1
    ivRl(i_new) = 0
    ivRu(i_new) = 0

    call extend_Ainv(kk, indv(1:kk+1))
    kk = kk + 1

    if(check_inversion) then
       call check_inverse(kk, indv(1:kk))
       if( error_flag /= 0 ) then
         return
       end if
    end if

    call get_indR(ivRl, ivRu, kl, indRl, ku, indRu)

    call get_weight1(ivS, ivRl, ivRu, kk, indv(1:kk),              &
         kl, indRl(1:kl), ku, indRu(1:ku), lambda_current, weight)
    call get_muP(weight, muP, sig)

    ! set the values at the turning point
    lambda_tp(count)   = lambda_current
    mu_tp(count)       = muP
    sig_tp(count)      = sig
    weight_tp(:,count) = weight

!    print '("cnt=",i5,"  kk=",i3," l_c=",g16.10,"  sig, mu=",2f16.12)', &
!         count, kk, lambda_current, sig, muP

    sigw0 = matmul(sigma, weight)
    w0sigw0 = dot_product(weight, sigw0)


    ! The iteration starts here:
    ! ---------

    do

       call get_indR(ivRl, ivRu, kl, indRl, ku, indRu)

       ! NOTE: this also calculates the vectors Sinv_* which are
       !  reused by get_lambda2
       call get_lambda1(ivS, ivRl, ivRu, kk, indv(1:kk),   &
            kl, indRl(1:kl), ku, indRu(1:ku),              &
            lambda_current, i1_new, lambda1_new,           &
            Sinv_eS, Sinv_muS, Sinv_SRwR, dir1)
       if( error_flag /= 0 ) then
         return
       end if

       ! NOTE: one has to call get_lambda1 to calculte the vectors Sinv_*
       call get_lambda2(ivS, ivRl, ivRu, kk, indv(1:kk),      &
            kl, indRl(1:kl), ku, indRu(1:ku),                 &
            Sinv_eS(1:kk), Sinv_muS(1:kk), Sinv_SRwR(1:kk),   &
            lambda_current, i2_new, lambda2_new)


       ! Hopefully it never enters this
       if( lambda1_new == lambda_aux .and. lambda2_new == lambda_aux) then
          if( print_on /= 0 ) then
            print *,'**** finished !!!???'
            print *,'lambda1_new=   ',lambda1_new 
            print *,'lambda2_new=   ',lambda2_new 
            print *,'lambda_current=',lambda_current 
          end if
          error_flag=1
          return
          !~ stop
       end if

       if( lambda1_new > lambda2_new .and. kk > 1 ) then

          lambda_current = lambda1_new
          i_last = i1_new
          ivS(i1_new) = 0
          if(dir1 == -1) then
             ivRl(i1_new) = 1
          elseif(dir1 == 1) then
             ivRu(i1_new) = 1
          else
             print *,'*** in turning: dir1=', dir1
             error_flag = 9
             return
             ! stop
          end if

          call shrink_Ainv(kk, indv(1:kk), i1_new, indv_temp(1:kk-1))
          if( error_flag /= 0 ) then
            return
          end if

          kk = kk - 1
          indv(1:kk) = indv_temp(1:kk)

          if(check_inversion) then
             call check_inverse(kk, indv(1:kk))
             if( error_flag /= 0 ) then
               return
             end if
          end if

       else

          lambda_current = lambda2_new
          i_last = i2_new
          indv(kk+1) = i2_new
          ivS(i2_new) = 1
          ivRl(i2_new) = 0
          ivRu(i2_new) = 0

          call extend_Ainv(kk, indv)
          kk = kk + 1

          if(check_inversion) then
             call check_inverse(kk, indv(1:kk))
             if( error_flag /= 0 ) then
               return
             end if
          end if

       end if

       kk_max = max(kk, kk_max)

       call get_indR(ivRl, ivRu, kl, indRl, ku, indRu)

       call get_weight1(ivS, ivRl, ivRu, kk, indv(1:kk),  &
            kl, indRl(1:kl), ku, indRu(1:ku),             &
            lambda_current, weight)
 
       if( correct_weights) then
          call correct_wghts(ivS, ivRl, ivRu, weight)
       end if

!       call get_muP(weight, muP, sig)
       call get_muP_fast(weight_tp(:,count), sigw0, w0sigw0, mu_tp(count),  &
            weight, muP, sig)

       count = count + 1

!       print '("cnt=",i5,"  kk=",i3," l_c=",g16.10,"  sig, mu=",2f16.12)', &
!            count, kk, lambda_current, sig, muP

       ! set the values at the turning point
       if(count > n_tpx) then
          if( print_on /= 0 ) then
            print *,'**** n_tpx is too small! ', count
          end if
          error_flag = 4
          RETURN
       else
          lambda_tp(count)   = lambda_current
          mu_tp(count)       = muP
          sig_tp(count)      = sig
          weight_tp(:,count) = weight
       end if

       if(lambda_current < lambda_min + epsx) then
          n_tp = count
          if(kk==1) then
             return
          else
            if( print_on /= 0 ) then
               print *,'kk=', kk
               print *,'lambda_current < lambda_min:'
               print *,'lambda_current=', lambda_current
               print *,'D: lambda_min=    ', lambda_min
            end if
            error_flag = 11
!             stop
            return
          end if
       end if
    end do
  end subroutine turning

  subroutine get_lambda_min(ix, wx, i1_min, lambda1_min)

    implicit none
    integer, intent(in)                :: ix      ! i_min
    real(8), dimension(nn), intent(in) :: wx      ! weight_min
    integer, intent(out)               :: i1_min
    real(8), intent(out)               :: lambda1_min

    real(8) :: lambda1, msm, mse, ese, det, aux, aux_mu, aux_e
    real(8), dimension(2,2) :: sigma0
    real(8), dimension(2)   :: mu0, e, sm, se, xx, v, sv, w0
    real(8) :: aux1, aux2, auxw, esv
    real(8), parameter        :: lambda_aux = 1.0d20
    integer :: i, j

    sigma0(1,1) = sigma(ix,ix)
    mu0(1) = mu(ix)
    e = (/ 1.0d0, 1.0d0 /)
    w0(1) = wx(ix)

    lambda1_min = lambda_aux
    i1_min = 0
    
    aux_mu = wx(ix)*mu(ix)
    aux_e =  wx(ix)

    do i = 1, nn
       if(i == ix) cycle

       sigma0(1,2) = sigma(ix,i)
       sigma0(2,1) = sigma0(1,2)
       sigma0(2,2) = sigma(i,i)
       mu0(2) = mu(i)
       w0(2) = wx(i)
       auxw = w0(1) + w0(2)

       v(1) = sum(sigma(:,ix)*wx) - sum(sigma0(:,1)*w0)
       v(2) = sum(sigma(:,i)*wx) - sum(sigma0(:,2)*w0)
       
       se = solve(2, sigma0, e)
       sm = solve(2, sigma0, mu0)
       sv = solve(2, sigma0, v)

       msm = dot_product(mu0, sm)
       mse = dot_product(mu0, se)
       ese = dot_product(e, se)
       esv = dot_product(e, sv)

       lambda1 = ((auxw + esv)*se(2) - (w0(2) + sv(2))*ese) / &
            (mse*se(2) - ese*sm(2))
       if( lambda1 < lambda1_min ) then
          i1_min = i
          lambda1_min = lambda1
       end if
    end do
  end subroutine get_lambda_min


  subroutine get_lambda_max(ix, wx, i1_max, lambda1_max)
    implicit none
    integer, intent(in)                :: ix      ! i_max
    real(8), dimension(nn), intent(in) :: wx      ! weight_max
    integer, intent(out)               :: i1_max
    real(8), intent(out)               :: lambda1_max

    real(8) :: lambda1, msm, mse, ese, det, aux, aux_mu, aux_e
    real(8), dimension(2,2) :: sigma0
    real(8), dimension(2)   :: mu0, e, sm, se, xx, v, sv, w0
    real(8) :: aux1, aux2, auxw, esv
    real(8), parameter        :: lambda_aux = 1.0d20
    integer :: i, j

    sigma0(1,1) = sigma(ix,ix)
    mu0(1) = mu(ix)
    e = (/ 1.0d0, 1.0d0 /)
    w0(1) = wx(ix)

    lambda1_max = -lambda_aux
    i1_max = 0
    
    aux_mu = wx(ix)*mu(ix)
    aux_e =  wx(ix)

    do i = 1, nn
       if(i == ix) cycle

       sigma0(1,2) = sigma(ix,i)
       sigma0(2,1) = sigma0(1,2)
       sigma0(2,2) = sigma(i,i)
       mu0(2) = mu(i)
       w0(2) = wx(i)
       auxw = w0(1) + w0(2)

       v(1) = sum(sigma(:,ix)*wx) - sum(sigma0(:,1)*w0)
       v(2) = sum(sigma(:,i)*wx) - sum(sigma0(:,2)*w0)
       
       se = solve(2, sigma0, e)
       sm = solve(2, sigma0, mu0)
       sv = solve(2, sigma0, v)

       msm = dot_product(mu0, sm)
       mse = dot_product(mu0, se)
       ese = dot_product(e, se)
       esv = dot_product(e, sv)

       lambda1 = ((auxw + esv)*se(2) - (w0(2) + sv(2))*ese) / &
            (mse*se(2) - ese*sm(2))
       if( lambda1 > lambda1_max ) then
          i1_max = i
          lambda1_max = lambda1
       end if
    end do
  end subroutine get_lambda_max


  function solve(n, A, b)
    ! solves the equation A*x=b
    implicit none 
    integer :: n
    real(8) :: A(n,n), b(n), solve(n), u(n)
    real(8) :: C(n,n+1), aux
    integer :: i, j

    C(:,1:n) = A
    C(:,n+1) = b

    do i = 1, n
       aux = 1.0d0/C(i,i)
       C(i,:) = aux*C(i,:)
       do j = i+1, n
          C(j,:) = C(j,:) - C(j,i)*C(i,:)
       end do
    end do

    u = C(:,n+1)
    do i = n, 2, -1
       do j = i-1, 1, -1
          !          C(j,:) = C(j,:) - C(j,i)*C(i,:)
          u(j) = u(j) - C(j,i)*u(i)
       end do
    end do

    !    solve = C(:,n+1)
    solve = u
  end function solve

  subroutine get_lambda1(jvS, jvRl, jvRu, k, indw, kl, indRl, ku, indRu,     &
       lambda_c, i_new, lambda_new, Ai_eS, Ai_muS, Ai_SRwR, dir_new)
    ! find the largest lambda, when one of the weights in S
    ! reach the boundary

    ! iw(i) = 1, i=1..nn  for i in S, and 0 for i not in S
    ! ind(i), i=1..k is the meaning of the indices of AS_inv

    implicit none
    integer, dimension(nn),  intent(in) :: jvS, jvRl, jvRu
    integer,                 intent(in) :: k, kl, ku
    integer, dimension(k),   intent(in) :: indw
    integer, dimension(kl),  intent(in) :: indRl
    integer, dimension(ku),  intent(in) :: indRu
    real(8), intent(in)                 :: lambda_c
    integer, intent(out)                :: i_new
    real(8), intent(out)                :: lambda_new
    real(8), dimension(k), intent(out)  :: Ai_eS, Ai_muS, Ai_SRwR
    integer, intent(out)                :: dir_new

    real(8), dimension(k) :: muS, eS, a, c, va
    real(8) :: alpha, beta, ac, ceS, cmuS, aux_e, aux_mu
    real(8) :: eS_Ai_eS, eS_Ai_muS, muii, aux_ee, aux_emu, lambda_i
    real(8) :: eS_Ai_SRwR, aux1, den, bi
    integer :: i, ii, dir, i0

    real(8), dimension(k) :: vp
    real(8) :: gamma
    real(8), parameter :: eps0 = 1.0d-10
    real(8), parameter        :: lambda_aux = 1.0d20

    eS = 1.0d0
    Ai_eS = matmul( Ai(1:k,1:k), eS)
    eS_Ai_eS = sum(Ai_eS(1:k))

    do i = 1, k
       muS(i) = mu(indw(i))
    end do
    Ai_muS = matmul( Ai(1:k,1:k), muS)
    eS_Ai_muS = sum(Ai_muS(1:k))

    do i = 1, k
       ii = indw(i)
       va(i) = sum(sigma(:,ii)*jvRl*lower) + sum(sigma(:,ii)*jvRu*upper)
    end do
    Ai_SRwR = matmul( Ai(1:k,1:k), va)
    eS_Ai_SRwR = sum(Ai_SRwR(1:k))

    if( k==1 ) return

    aux1 = 1.0d0 - sum(lower(indRl(1:kl))) - &
         sum(upper(indRu(1:ku))) + eS_Ai_SRwR

    lambda_new = -lambda_aux
    i_new = 0
    i0 = 0
    do i = 1, k
       ii = indw(i)

       den = eS_Ai_muS*Ai_eS(i) - eS_Ai_eS*Ai_muS(i)
       if( den > 0.0d0 ) then
          bi = upper(ii)
          dir = 1
       elseif( den < 0.0d0) then
          bi = lower(ii)
          dir = -1
       else
          if( print_on /= 0 ) then
            print *,'*** in get_lambda1: den=', den
          end if
          error_flag = 5
          return
          !~ stop
       end if

       lambda_i = (aux1*Ai_eS(i) - eS_Ai_eS*(bi+Ai_SRwR(i)))/den

!       if( lambda_i >= lambda_c - eps0 ) cycle
       if( lambda_i > lambda_new) then
          lambda_new = lambda_i
          i_new = ii
          dir_new = dir
          i0 = i
       end if
    end do
    if(i_new == 0) then
      if( print_on /= 0 ) then
        print *,'*** in get_lambda1: i_new=', i_new
      end if
      error_flag=6
      return
       !~ stop
    end if
  end subroutine get_lambda1

  subroutine get_lambda2(jvS, jvRl, jvRu, k, indw, kl, indRl, ku, indRu, &
       Ai_eS, Ai_muS, Ai_SRwR, lambda_c, i_new, lambda_new)
    ! find the largest lambda, when one of the weights on the bounds
    ! wants to come inside

    ! jvS(i) = 1, i=1..nn  for i in S,
    ! jvRl(i), jvRu(i) = 1 if w(i) is on the lower/upper bound
    ! indw(i), i=1..k is mapping of the indices of AS_inv to 1..nn

    implicit none
    integer, dimension(nn), intent(in) :: jvS, jvRl, jvRu
    integer,                intent(in) :: k, kl, ku
    integer, dimension(k),  intent(in) :: indw
    integer, dimension(kl), intent(in) :: indRl
    integer, dimension(ku), intent(in) :: indRu
    real(8), dimension(k),  intent(in) :: Ai_eS, Ai_muS, Ai_SRwR
    real(8),                intent(in) :: lambda_c
    integer,               intent(out) :: i_new
    real(8),               intent(out) :: lambda_new

    real(8), dimension(kl) :: low
    real(8), dimension(ku) :: upp  ! NOTE: ku can be zero as well, but it works

    real(8), dimension(k) :: muS, eS, a, c
    real(8) :: alpha, beta, ac, ceS, cmuS, aux_e, aux_mu
    real(8) :: eS_Ai_eS, eS_Ai_muS, aux_ee, aux_emu, lambda_i
    real(8) :: aux, aux1, auxi, eS_Ai_SRwR, cSRwR, den
    integer :: i, ii, ip
    real(8), parameter :: eps0 = 1.0d-10
    real(8), parameter        :: lambda_aux = 1.0d20

    eS = 1.0d0
    eS_Ai_eS = sum(Ai_eS(1:k))

    eS_Ai_muS = sum(Ai_muS(1:k))

    eS_Ai_SRwR = sum(Ai_SRwR(1:k))

    low = lower(indRl(1:kl))
    upp = upper(indRu(1:ku))

    aux1 = 1.0d0 - sum(low) - sum(upp) + eS_Ai_SRwR

    lambda_new = -lambda_aux
    i_new = 0

    do ii = 1, nn

       if(jvS(ii) == 1) cycle

       a = sigma(indw(1:k),ii)
       alpha = sigma(ii,ii)

       ceS = dot_product(a, Ai_eS)
       cmuS = dot_product(a, Ai_muS)
       cSRwR = dot_product(a, Ai_SRwR)

       auxi = sum(sigma(indRl(1:kl),ii)*low) +  &
             sum(sigma(indRu(1:ku),ii)*upp)

       den = eS_Ai_muS*(1.0d0 - ceS) - eS_Ai_eS*(mu(ii) - cmuS)
       
       lambda_i = (aux1*(1.0d0 - ceS) + (cSRwR - auxi)*eS_Ai_eS) / den

       if( lambda_i >= lambda_c) cycle

       ! to avoid going back
       if( den < 0.0d0 .and. jvRl(ii) == 1 ) cycle 
       if( den > 0.0d0 .and. jvRu(ii) == 1 ) cycle 

       if( lambda_i > lambda_new) then
          lambda_new = lambda_i
          i_new = ii
       end if
    end do
  end subroutine get_lambda2

  subroutine extend_Ainv(k, indw)
    ! extend the inverse from k*k to (k+1)*(k+1)

    ! ind(i), i=1..k is the meaning of the indices of Ai
    ! ind(k+1) is the new index to be added

    implicit none
    integer,                   intent(in)  :: k
    integer, dimension(k+1),   intent(in)  :: indw

    real(8), dimension(k) :: a, c
    real(8)               :: alpha, beta, ac, aux
    integer               :: j, ii

    ii = indw(k+1)
    a = sigma(indw(1:k),ii)
    alpha = sigma(ii,ii)

    c = matmul(Ai(1:k,1:k), a)

    ac = dot_product(a,c)
    beta = 1.0d0/(alpha-ac)

    do j = 1, k
       aux = beta*c(j)
       Ai(1:k,j) = Ai(1:k,j) + aux*c(1:k)
    end do
    Ai(k+1,1:k) = -beta*c(1:k)
    Ai(1:k,k+1) = -beta*c(1:k)
    Ai(k+1,k+1) = beta
  end subroutine extend_Ainv

  subroutine shrink_Ainv(k, indw, ii0, indw_new)
    ! shrink the inverse from k*k to (k-1)*(k-1)

    ! ind(i), i=1..k is the meaning of the indices of Ai
    ! ii0 is the index to leave out

    implicit none
    integer,                   intent(in)  :: k, ii0
    integer, dimension(k),     intent(in)  :: indw
    integer, dimension(k-1),   intent(out) :: indw_new

    real(8), dimension(k-1) :: b
    integer, dimension(k-1) :: ivp
    real(8)                 :: beta, aux
    integer                 :: i0, j, i
    !-----------

    if(k <= 1) then
       if( print_on /= 0 ) then
         print *,'*** shrink_Ainv called with k=', k
       end if
       error_flag = 7
       return
       !~ stop
    end if

    i0 = 0
    j = 0
    do i = 1, k
       if(indw(i)==ii0) then
          i0 = i
       else
          j = j + 1
          indw_new(j) = indw(i)
          ivp(j) = i
       end if
    end do
    if(i0==0) then
       if( print_on /= 0 ) then
         print *,'*** ii0 missing:  ii0=', ii0,'  k=', k
         print *,'indw=', indw
       end if
       error_flag = 12
       return
       !~ stop
    end if

    b = Ai(ivp(:),i0)
    beta = Ai(i0,i0)

    do j = 1, k-1
       aux = -1.0d0/beta*b(j)
       Ai(1:k-1,j) = Ai(ivp(:),ivp(j)) + aux*b(:)
    end do
  end subroutine shrink_Ainv


  subroutine check_inverse(k, indw)
    ! check that Ai is the inverse of sigma_S

    ! ind(i), i=1..k is the meaning of the indices of Ai

    implicit none
    integer,                   intent(in)  :: k
    integer, dimension(k),     intent(in)  :: indw

    real(8), dimension(k,k) :: sig_S, res
    real(8), dimension(1)   :: auxv
    integer                 :: i

    sig_S = sigma(indw(:), indw(:))
    res = matmul( Ai(1:k,1:k), sig_S)

    do i = 1, k
       res(i,i) = res(i,i) - 1.0d0
    end do
    auxv = maxval( abs(res))
    if( auxv(1) > 1.0d-10 ) then
       if( print_on /= 0 ) then
         print *,'**** in AS_inv*sig_S: k=', k,'   error=', auxv
       end if
       error_flag = 8
       return
       !~ stop
    end if
  end subroutine check_inverse

  subroutine get_muP(w, muP, sig1)
    ! finds the value of muP and sigma for a given weight
    implicit none
    real(8), dimension(nn),  intent(in) :: w
    real(8),                intent(out) :: muP, sig1
    real(8)                             :: aux

    muP = dot_product(w, mu)
    aux = dot_product(w, matmul(sigma, w))
    sig1 = sqrt(aux)
  end subroutine get_muP

  subroutine get_muP_fast(w0, sigw0, w0sigw0, muP0, w, muP, sig1)
    ! finds the value of muP and sigma for a given weight (fast version)
    implicit none
    real(8), dimension(nn),     intent(in) :: w0
    real(8), dimension(nn),  intent(inout) :: sigw0
    real(8),                 intent(inout) :: w0sigw0
    real(8),                 intent(in)    :: muP0
    real(8), dimension(nn),  intent(in)    :: w
    real(8),                intent(out)    :: muP, sig1
    real(8)                                :: aux, sig1a, dwsigw0, dwsigdw
    real(8), dimension(nn)                 :: dw, vp
    integer, dimension(nn)                 :: iv
    integer                                :: i, j, icnt, ip, jp

    iv = 0
    icnt = 0
    dw = 0.
    do i = 1, nn
       if(w0(i) /= w(i)) then
          icnt = icnt + 1
          iv(icnt) = i
          dw(icnt) = w(i) - w0(i)
       end if
    end do

    aux = sum( dw(1:icnt)*mu(iv(1:icnt)))
    muP = muP0 + aux

    dwsigw0 = sum(sigw0(iv(1:icnt))*dw(1:icnt))
    dwsigdw = 0.
    do ip = 1, icnt
       i = iv(ip)
       aux = sum( sigma(iv(1:icnt),i)*dw(1:icnt))
       dwsigdw = dwsigdw + dw(ip)*aux
    end do

    do i = 1, nn
       aux = sum( sigma(iv(1:icnt),i)*dw(1:icnt))
       sigw0(i) = sigw0(i) + aux    ! replaced by the new value!!!
    end do

    aux = w0sigw0 + 2.0d0*dwsigw0 + dwsigdw
    w0sigw0 = aux                   ! replaced by the new value!!!
    sig1 = sqrt(aux)

  end subroutine get_muP_fast

  subroutine get_weight1(jvS, jvRl, jvRu, k, indw, kl, indRl, &
       ku, indRu, lambda, weight1)
    ! finds the weight for the given lambda

    ! jvS,...
    ! indw(i), i=1..k  --> 1..nn

    implicit none
    integer, dimension(nn),  intent(in) :: jvS, jvRl, jvRu
    integer,                 intent(in) :: k, kl, ku
    integer, dimension(k),   intent(in) :: indw
    integer, dimension(kl),  intent(in) :: indRl
    integer, dimension(ku),  intent(in) :: indRu
    real(8),                 intent(in) :: lambda
    real(8), dimension(nn), intent(out) :: weight1

    real(8), dimension(kl) :: low
    real(8), dimension(ku) :: upp
    real(8), dimension(k) :: muS, eS, Ai_eS, Ai_muS, va, Ai_va
    real(8) :: eS_Ai_eS, eS_Ai_muS, gamma, eS_Ai_va
    real(8) :: auxv(1), err, aux, auxv1(1), auxv2(1), aux1
    integer :: i, ii

    eS = 1.0d0
    Ai_eS = matmul( Ai(1:k,1:k), eS)
    eS_Ai_eS = sum(Ai_eS(1:k))

    do i = 1, k
       muS(i) = mu(indw(i))
    end do

    Ai_muS = matmul( Ai(1:k,1:k), muS)
    eS_Ai_muS = sum(Ai_muS(1:k))

    low = lower(indRl(1:kl))
    upp = upper(indRu(1:ku))

    do i = 1, k
       ii = indw(i)
       va(i) = sum(sigma(indRl(1:kl),ii)*low) + &
            sum(sigma(indRu(1:ku),ii)*upp)
    end do
    Ai_va = matmul( Ai(1:k,1:k), va)
    eS_Ai_va = sum(Ai_va(1:k))

    aux1 = 1.0d0 - sum(low) - sum(upp) + eS_Ai_va

    gamma = (aux1 - lambda*eS_Ai_muS)/eS_Ai_eS

    weight1 = 0.
    do i = 1, k 
       ii = indw(i)
       weight1(ii) = lambda*Ai_muS(i) + gamma*Ai_eS(i) - Ai_va(i)
    end do

    do i = 1, nn
       if(jvRl(i) == 1) then
          weight1(i) = lower(i)
       elseif(jvRu(i) == 1) then
          weight1(i) = upper(i)
       end if
    end do

  end subroutine get_weight1

  !-------------

  function invert(n, A)
    implicit none 
    integer :: n
    real(8) :: A(n,n), invert(n,n), u(n)
    real(8) :: C(n,2*n), aux
    integer :: i, j

    C(:,1:n) = A
    C(:,n+1:) = 0.
    do i = 1, n
       C(i,n+i) = 1.0d0
    end do

    do i = 1, n
       aux = 1.0d0/C(i,i)
       C(i,:) = aux*C(i,:)
       do j = i+1, n
          C(j,:) = C(j,:) - C(j,i)*C(i,:)
       end do
    end do

    do i = n, 2, -1
       do j = i-1, 1, -1
          !          C(j,:) = C(j,:) - C(j,i)*C(i,:)
          C(j,n+1:) = C(j,n+1:) - C(j,i)*C(i,n+1:)
       end do
    end do
    invert = C(:,n+1:2*n)
  end function invert

  !-------------------------
  ! Needed only for testing

  subroutine randvec(np, rvec, shift)
    implicit none
    integer :: np
    real(8) :: rvec(np)
    real(8) :: shift, rnd
    integer :: i

    rvec = 0.d0
    do i = 1, np
       call random_number(rnd)
       rvec(i) = rnd - 0.5d0 + shift
    end do
  end subroutine randvec

  subroutine matgen(np, aap, shift)
    ! generate a positive definite random matrix
    implicit none
    integer :: np
    real(8) :: aap(np,np), auxv(np)
    real(8) :: shift, rnd, aux
    integer :: j, na, ia

    aap = 0.0d0
    !    na = 2*np
    na = np + 5
    do ia = 1, na
       call randvec(np,auxv, shift)
       do j = 1, np
          aux = auxv(j)
          aap(:,j) = aap(:,j) + auxv(:)*aux
       end do
    end do
  end subroutine matgen
  !-------------------------------

  subroutine index_sort0(n, x, ind)
    ! indexing the array x(1:n) in increasing order
    ! a slower but simpler version
    implicit none
    integer, intent(in)  :: n
    real(8), intent(in)  :: x(n)
    integer, intent(out) :: ind(n)
    integer              :: inda(n)
    integer              :: i, iauxv(1), ip
    real(8)              :: aux0, auxv(1), xp(n)

    xp = x
    auxv = maxval(x)
    aux0 = auxv(1) + 1.0d0

    do i = 1, n
       iauxv = minloc(xp)
       ip = iauxv(1)
       ind(i) = ip
       xp(ip) = aux0
    end do
  end subroutine index_sort0

  subroutine index_sort(n, x, ind)
    ! indexing the array x(1:n) in increasing order
    implicit none
    integer, intent(in)  :: n
    real(8), intent(in)  :: x(n)
    integer, intent(out) :: ind(n)
    integer              :: inda(n)
    integer              :: i

    if( n == 1) then
       ind(1) =  1
       return
    end if
    do i = 1, n
       inda(i) = i
    end do
    call i_sort(n, x, inda, ind)
  end subroutine index_sort

  recursive subroutine i_sort(n, x, ind, indp)
    implicit none
    integer, intent(in)  :: n
    real(8), intent(in)  :: x(n)
    integer, intent(in)  :: ind(n)
    integer, intent(out) :: indp(n)
    integer              :: ind0(n), ind1(n)
    integer :: i0, i1, i
    real(8) :: x0, y0(n), y1(n)

    if( n == 1) then
       indp(1) =  ind(1)
       return
    end if

    x0 = x(1)
    i0 =  0
    i1 =  0
    do i = 2, n
       if(x(i) < x0) then
          i0 = i0 + 1 
          y0(i0) = x(i)
          ind0(i0) = ind(i)
       else
          i1 = i1 + 1
          y1(i1) = x(i)
          ind1(i1) = ind(i)
       end if
    end do

    if( i0 > 0 ) then
       call i_sort(i0, y0(1:i0), ind0(1:i0), indp(1:i0))
    end if
    if( i1 > 0 ) then
       call i_sort(i1, y1(1:i1), ind1(1:i1), indp(i0+2:n))
    end if
    indp(i0+1) = ind(1)
  end subroutine i_sort

  !--------------------------
  ! for testing

  subroutine create_bounds
    implicit none
    integer :: i
    real(8) :: rnd1, rnd2, aux, aux1

    ! set lower and upper bounds:
    do i = 1, nn
       call random_number(rnd1)
       call random_number(rnd2)
       lower(i) = min(rnd1, rnd2)
       upper(i) = 2.0d0*max(rnd1, rnd2)
    end do
    aux = sum(lower)
    lower = (0.5d0/aux)*lower
    upper = (0.5d0/aux)*upper
    
    aux = sum(upper)
    if(aux < 1.0d0) then
       aux1 = (1.0d0 - aux)*2.0d0/nn
       upper = upper + aux1
    end if
  end subroutine create_bounds


  subroutine check_bounds
    implicit none
    integer :: i, ierror
    real(8) :: aux

    ierror = 0

    do i = 1, nn
       if(lower(i) > upper(i) .or. lower(i) < 0.0d0) then ! TODO: find out why negative lower bound not accepted?
       !~ if(lower(i) > upper(i)) then
          ierror = ierror + 1
         if( print_on /= 0 ) then
            print *,'bounds: ', i, lower(i), upper(i)
         end if
       end if
    end do

    aux = sum(lower)
    if(aux > 1.0d0) then
       ierror = ierror + 1
       if( print_on /= 0 ) then
         print *,'sum(lower)=', aux
       end if
    end if

    aux = sum(upper)
    if(aux < 1.0d0) then
       ierror = ierror + 1
       if( print_on /= 0 ) then
         print *,'sum(lower)=', aux
       end if
    end if
    
    if(ierror > 0) then
       error_flag = 12
       return
       ! stop
    endif
  end subroutine check_bounds

  !-----------------------

  subroutine get_starting_weight_max(jvS, jvRl, jvRu, i_mx, weight_mx)
    ! fill in the weights to have a maximal muP (== mu_max !!!)
    ! and find i_max, the only index for which the weight is not
    ! on its bounds.
    ! it is assumed that the bounds satisfy the appropriate conditions
    ! for this one should call check_bounds (or make it sure otherwise)

    implicit none
    integer, dimension(nn), intent(out) :: jvS, jvRl, jvRu
    integer,                intent(out) :: i_mx
    real(8), dimension(nn), intent(out) :: weight_mx
    integer :: i, j, cnt
    real(8) :: remainder, aux
    real(8), parameter :: eps = 1.0d-14

    weight_mx = lower
    remainder = 1.0d0 - sum(lower)

    do i = nn, 1, -1
       j = ind_mu(i)
       aux = min( upper(j)-lower(j), remainder, 1.0d0-weight_mx(j))
       weight_mx(j) = weight_mx(j) + aux
       remainder = remainder - aux
       if(remainder < eps) exit
    end do

    ! correct for numerical error
    jvS = 0
    jvRl = 0
    jvRu = 0
    do i = 1, nn
       if( abs(weight_mx(i) - lower(i)) < eps) then
          weight_mx(i) = lower(i)
          jvRl(i) = 1
       elseif( abs(weight_mx(i) - upper(i)) < eps) then
          weight_mx(i) = upper(i)
          jvRu(i) = 1
       else
          i_mx = i
          jvS(i) = 1
       end if
    end do

    mu_max = dot_product(mu, weight_mx)

    if( sum(jvS) /= 1) then
       if( print_on /= 0 ) then
         print *,'**** in get_starting_weight_max: sum(jvS)=', sum(jvS) 
       end if
       error_flag = 10
       return
       !~ stop
    end if

  end subroutine get_starting_weight_max

  subroutine get_starting_weight_min(jvSn, jvRln, jvRun, i_mn, weight_mn)
    ! fill in the weights to have a minimal muP (== mu_min !!!)
    ! and find i_min, the only index for which the weight is not
    ! on its bounds.
    ! it is assumed that the bounds satisfy the appropriate conditions
    ! for this one should call check_bounds (or make it sure otherwise)

    implicit none

    integer, dimension(nn), intent(out) :: jvSn, jvRln, jvRun
    integer,                intent(out) :: i_mn
    real(8), dimension(nn), intent(out) :: weight_mn
    integer :: i, j, cnt
    real(8) :: remainder, aux
    real(8), parameter :: eps = 1.0d-14

    weight_mn = lower
    remainder = 1.0d0 - sum(lower)

    do i = 1, nn
       j = ind_mu(i)
       aux = min( upper(j)-lower(j), remainder, 1.0d0-weight_mn(j))
       weight_mn(j) = weight_mn(j) + aux
       remainder = remainder - aux
       if(remainder < eps) exit
    end do

    ! correct for numerical error
    jvSn = 0
    jvRln = 0
    jvRun = 0
    do i = 1, nn
       if( abs(weight_mn(i) - lower(i)) < eps) then
          weight_mn(i) = lower(i)
          jvRln(i) = 1
       elseif( abs(weight_mn(i) - upper(i)) < eps) then
          weight_mn(i) = upper(i)
          jvRun(i) = 1
       else
          i_mn = i
          jvSn(i) = 1
       end if
    end do

    mu_min = dot_product(mu, weight_mn)

    if( sum(jvSn) /= 1) then
       if( print_on /= 0 ) then
         print *,'**** in get_starting_weight_min: sum(jvSn)=', sum(jvSn) 
       end if
       error_flag = 13
       return
       !~ stop
    end if
  end subroutine get_starting_weight_min


  subroutine check_weight(w)
    implicit none
    real(8), dimension(nn) :: w
    real(8) :: auxv(1)
    !~ real(8), parameter :: eps_cw = 1.0d-12
    real(8), parameter :: eps_cw = 1.0d-10
    integer :: i

    auxv = minval( w - lower )
    if( auxv(1) < -eps_cw ) then
       if( print_on /= 0 ) then
         print *,'*** min(w-lower)=', auxv(1)
       end if
       error_flag = 14
       return
       !~ stop
    end if

    auxv = minval( upper - w )
    if( auxv(1) < -eps_cw ) then
       if( print_on /= 0 ) then
         print *,'*** min(upper-w)=', auxv(1)
         do i = 1, nn
            print *, i, lower(i), w(i), upper(i)
         end do
       end if
       error_flag = 15
       return
       !~ stop
    end if
  end subroutine check_weight

  subroutine correct_wghts(jvS, jvRl, jvRu, w)
    implicit none
    integer, dimension(nn), intent(in) :: jvS, jvRl, jvRu
    real(8), dimension(nn)             :: w
    real(8) :: aux, aux1
    integer :: i

    do i = 1, nn
       if(jvRl(i) == 1) then
          w(i) = lower(i)
       elseif(jvRu(i) == 1) then
          w(i) = upper(i)
       end if
    end do

    aux1 = sum(w*jvS)
    aux = (1.0d0 - sum(jvRl*lower) - sum(jvRu*upper))/aux1
    do i = 1, nn
       if(jvS(i) == 1) then
          w(i) = aux*w(i)
       end if
    end do
  end subroutine correct_wghts

  subroutine write_data
    implicit none
    integer :: i, j

    open(20,file='sigma0.dat')
    write(20,*) nn
    do i = 1, nn
       do j = 1, i
          write(20,*) sigma(i,j)
       end do
    end do
    close(20)
    
    open(20,file='mu_bounds0.dat')
    write(20,*) nn
    do i = 1, nn
       write(20,*) mu(i), lower(i), upper(i)
    end do
    close(20)
  end subroutine write_data

  subroutine read_data
    implicit none
    integer :: i, j, nnp
    real(8) :: aux
    character(len=100) :: sigma_file, mu_file

    sigma_file = 'sigma_1000_sh_05.dat'
    mu_file = 'mu_1000.dat'

    open(20,file=sigma_file)
    read(20,*) nnp
    if(nn > nnp) then
       if( print_on /= 0 ) then
         print *, '****  nn too large: nn=', nn, '  nnp=', nnp
       end if
       error_flag = 16
       return
       !~ stop
    end if
    do i = 1, nn
       do j = 1, i
          read(20,*) aux
          sigma(i,j) = aux
          sigma(j,i) = aux
       end do
    end do
    close(20)
    
!!$    open(20,file='mu_bounds0.dat')
!!$    read(20,*) nnp
!!$    if(nn > nnp) then
!!$       print *, '****  nn too large: nn=', nn, '  nnp=', nnp
!!$       stop
!!$    end if
!!$    do i = 1, nn
!!$       read(20,*) mu(i), lower(i), upper(i)
!!$    end do
!!$    close(20)

    open(20,file=mu_file)
    read(20,*) nnp
    if(nn > nnp) then
       if( print_on /= 0 ) then
          print *, '****  nn too large: nn=', nn, '  nnp=', nnp
       end if
       error_flag = 17
       return
       !~ stop
    end if
    do i = 1, nn
       read(20,*) mu(i)
    end do
    close(20)
  end subroutine read_data

  subroutine print_iv(jvS, jvRl, jvRu, i_last, i_new)
    implicit none
    integer, dimension(nn) :: jvS, jvRl, jvRu
    integer :: i_last, i_new, i
    character(len=3)   :: txt
    character(len=100) :: text

    text = 'ivS :'
    do i = 1, nn
       if(jvS(i)==1) then
          write(txt,'(i3)') i
          text = trim(text) // txt
          if(i==i_new) text = trim(text) // '*'
          if(i==i_last) text = trim(text) // '+'
       end if
    end do
    if( print_on /= 0 ) then
      print *, trim(text)
    end if
    
    text = 'ivRl:'
    do i = 1, nn
       if(jvRl(i)==1) then
          write(txt,'(i3)') i
          text = trim(text) // txt
          if(i==i_new) text = trim(text) // '*'
          if(i==i_last) text = trim(text) // '+'
       end if
    end do
    if( print_on /= 0 ) then
      print *, trim(text)
    end if

    text = 'ivRu:'
    do i = 1, nn
       if(jvRu(i)==1) then
          write(txt,'(i3)') i
          text = trim(text) // txt
          if(i==i_new) text = trim(text) // '*'
          if(i==i_last) text = trim(text) // '+'
       end if
    end do
    if( print_on /= 0 ) then
      print *, trim(text)
    end if

  end subroutine print_iv

  subroutine timing(text)

    ! prints out the time in seconds from the previous
    ! together with the <text>

    ! call timing('')                 ! resets t0 and t00
    ! call xxxxx
    ! call timing('some_text')        ! prints: 'some_text' t-t0 sec

    implicit none
    character(len=*)  :: text

    integer, save :: itime0
    integer :: itime, irate, imax, dtime

    if( text == '' ) then
       call system_clock(itime0)
    else
       call system_clock(itime, irate, imax)
       dtime = itime - itime0
       if(dtime < 0) dtime = dtime + imax
       if( print_on /= 0 ) then
         print '(a,f9.3," [sec]")', trim(text), float(dtime)/irate
       end if
    end if
  end subroutine timing

  subroutine get_indR(jvRl, jvRu, kl, indRl, ku, indRu)
    ! constructs indRl(1:kl), indRu(1:ku) which are indexing 
    ! the weights with values at the lower and upper bounds
    implicit none
    integer, dimension(nn) :: jvRl, jvRu, indRl, indRu
    integer :: kl, ku
    integer :: i, j

    kl = sum(jvRl)
    j = 0
    do i = 1, nn
       if(jvRl(i) == 1) then
          j = j + 1
          indRl(j) = i
       end if
    end do

    ku = sum(jvRu)
    j = 0
    do i = 1, nn
       if(jvRu(i) == 1) then
          j = j + 1
          indRu(j) = i
       end if
    end do
  end subroutine get_indR

end SUBROUTINE calcturningpointsgen
