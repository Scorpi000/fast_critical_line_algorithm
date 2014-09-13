! g95  turn.f90 -o turn
! g95 -O2 turn.f90 -o turn
! g95 -g turn.f90 -o turn
! f2py -c calcturningpoints.f90 calcturningpointsgen.f90 -m cla_fortran

SUBROUTINE calcturningpoints(mu, sigma, nn, mu_tps, sig_tps, weight_tps, lambda_tps, n_tps, error_flag, max_tps)
  implicit none
  
  ! using allocatable arrays:
  integer, intent(in)                  :: nn, max_tps
  integer                              :: kx
  integer, intent(out)                 :: error_flag
  integer                              :: n_tpx
  
  real(8), dimension(:,:), allocatable :: Ai
  !~ real(8), dimension(:),   allocatable :: mu
  !~ real(8), dimension(:,:), allocatable :: sigma
  real(8), dimension(nn), intent(in)  :: mu
  real(8), dimension(nn,nn), intent(in) :: sigma

  ! the turning points
  real(8), dimension(max_tps), intent(out)  :: mu_tps, sig_tps, lambda_tps
  real(8), dimension(nn,max_tps), intent(out)  :: weight_tps
  !----------------------------------

  integer                   :: i_min, i_max, ii_last
  real(8), parameter        :: lambda_aux = 1.0d20
  real(8)                   :: mu_max, mu_min, sig, lambda_min, lambda_max

  ! the actual number of turning points:
  integer, intent(out)                   :: n_tps


  ! these declarations are needed only for testing:
  integer :: i, n1
  real(8) :: rnd
  
  integer :: print_on = 0
  !-------------------

  ! for the version with allocations:
  !~ nn = 12
  error_flag = 0
  
  kx = nn
  n_tpx = max_tps
  !~ allocate( Ai(kx,kx), mu(nn), sigma(nn,nn) )
  allocate( Ai(kx,kx) )
  !----------------

  !----------------
  !!! for testing: set mu and sigma:
  !~ do i = 1, nn
     !~ call random_number(rnd)
     !~ mu(i) = rnd + 0.01d0*i
  !~ end do
  !~ call matgen(nn, sigma, 0.0d0)
  !----------------

  call turning(n_tps, lambda_tps, mu_tps, sig_tps, weight_tps)
  if( error_flag /= 0 ) then
    return
  endif
  
  if( print_on /= 0 ) then
    print *,'mu_min=    ', mu_min, i_min
    print *,'mu_max=    ', mu_max, i_max
    print *,'lambda_min=', lambda_min
    print *,'lambda_max=', lambda_max
    print *,'minval(weight_tps)=', minval(weight_tps)
    print *,'maxval(weight_tps)=', maxval(weight_tps)

    print *,'---------------------------------'
    print *,'Number of turning points: n_tps=', n_tps
  
    n1 = min(10, nn)
    if( n1 < nn)  then
       do i = 1, n_tps
          print '(i4,f10.2,2f7.3,3x,10f6.2," ...")', i, lambda_tps(i), &
               mu_tps(i), sig_tps(i), weight_tps(1:n1,i)
       end do
    else
       do i = 1, n_tps
          print '(i4,f10.2,2f7.3,3x,10f6.2)', i, lambda_tps(i), &
               mu_tps(i), sig_tps(i), weight_tps(:,i)
       end do
    end if
  end if

  
  deallocate(Ai)

contains

  subroutine turning(n_tp, lambda_tp, mu_tp, sig_tp, weight_tp)
    implicit none

    ! total number of turning points 
    integer                    :: n_tp

    ! values at the turning points:
    real(8), dimension(n_tpx)     :: mu_tp, sig_tp, lambda_tp
    real(8), dimension(nn,n_tpx)  :: weight_tp


    real(8), dimension(nn) :: w, weight, Sinv_eS, Sinv_muS
    integer, dimension(nn) :: iv, indv, indv_temp
    integer :: kk_max, i, i_aux, i_new, kk, count, i1_new, i2_new
    real(8) :: lambda_current, lambda1_new, lambda2_new, lambda_new, muP

    kk_max = 0
    i_min = 1
    i_max = 1
    mu_min = mu(1)
    mu_max = mu(1)
    do i = 2, nn
       if( mu(i) < mu_min ) then
          i_min = i
          mu_min = mu(i)
       elseif( mu(i) > mu_max ) then
          i_max = i
          mu_max = mu(i)
       end if
    end do

    ! starting solution
    w= 0.
    w(i_max) = 1.0d0

    iv = 0
    iv(i_max) = 1
    ii_last = i_max

    call get_lambda_min(i_aux, lambda_min)

    call get_lambda_max(i_aux, lambda_max)

    ! get lambda_max with the general routine:
    Ai(1,1) = 1.0d0/sigma(i_max,i_max)
    Sinv_eS(1) = Ai(1,1)
    Sinv_muS(1) = Ai(1,1)*mu(i_max)
    indv(1) = i_max
    lambda_current = lambda_aux
    call get_lambda2(iv, 1, indv(1:1), Sinv_eS(1:1), Sinv_muS(1:1), &
         lambda_current, i_new, lambda_new)
    lambda_max = lambda_new
    lambda_current = lambda_max
    ii_last = i_new


    kk = 1
    indv(kk+1) = i_new
    iv(i_new) = 1
    call extend_Ainv(kk, indv(1:kk+1))
    kk = kk + 1

    !   call check_inverse(kk, indv(1:kk))

    count = 1

    call get_weight1(kk, indv(1:kk), lambda_current, weight)
    call get_muP(iv, kk, indv(1:kk), lambda_current, muP, sig)

    ! set the values at the turning point
    lambda_tp(count)   = lambda_current
    mu_tp(count)       = muP
    sig_tp(count)      = sig
    weight_tp(:,count) = weight


!!$    print '("cnt=",i5,"  kk=",i3," l_c=",g16.10,"  sig, mu=",2f16.12)', &
!!$         count, kk, lambda_current, sig, muP


    ! The iteration starts here:
    ! ---------

    do 
       call get_lambda1(iv, kk, indv(1:kk), lambda_current, &
            i1_new, lambda1_new, Sinv_eS, Sinv_muS)
       ! print *,'i1_new=', i1_new,'  lambda1_new=', lambda1_new, ii_last

       call get_lambda2(iv, kk, indv(1:kk), Sinv_eS(1:kk), Sinv_muS(1:kk),&
            lambda_current, i2_new, lambda2_new)
       ! print *,'i2_new=', i2_new,'  lambda2_new=', lambda2_new, ii_last

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

       if( lambda1_new > lambda2_new ) then
          lambda_current = lambda1_new
          ii_last = i1_new
          iv(i1_new) = 0
          call shrink_Ainv(kk, indv(1:kk), i1_new, indv_temp(1:kk-1))
          if( error_flag /= 0 ) then
            return
          end if
          kk = kk - 1
          indv(1:kk) = indv_temp(1:kk)

          !         call check_inverse(kk, indv(1:kk))

       else
          lambda_current = lambda2_new
          ii_last = i2_new
          indv(kk+1) = i2_new
          iv(i2_new) = 1
          call extend_Ainv(kk, indv)
          kk = kk + 1

          !         call check_inverse(kk, indv(1:kk))

       end if

       kk_max = max(kk, kk_max)
       call get_weight1(kk, indv(1:kk), lambda_current, weight)
       call get_muP(iv, kk, indv(1:kk), lambda_current, muP, sig)

       count = count + 1
!!$       print '("cnt=",i5,"  kk=",i3," l_c=",g16.10,"  sig, mu=",2f16.12)', &
!!$            count, kk, lambda_current, sig, muP

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
          n_tp=count
       end if


       if(kk==1) then
          ! print '(10x,"lambda_min= ", g16.10,19x,"mu_min=",f16.12)', &
          !        lambda_min, mu_min
          ! print *,'finished!,  kk_max=', kk_max
          !~ n_tp = count
          return
       end if
    end do

  end subroutine turning

  subroutine get_lambda_min(i1_min, lambda1_min)
    implicit none

    ! global: i_min, mu_min 

    integer :: i1_min
    real(8) :: lambda1_min
    real(8) :: lambda1, msm, mse, ese, det, aux
    real(8), dimension(2,2) :: sigma0
    real(8), dimension(2)   :: mu0, e, sm, se, xx

    integer :: i
    !----------

    !    print *,'i_min=', i_min, '  mu_min=', mu_min

    sigma0(1,1) = sigma(i_min,i_min)
    mu0(1) = mu(i_min)
    e = 1.0d0

    lambda1_min = lambda_aux
    i1_min = 0
    do i = 1, nn
       if(i==i_min) cycle

       sigma0(1,2) = sigma(i_min,i)
       sigma0(2,1) = sigma0(1,2)
       sigma0(2,2) = sigma(i,i)

       se = solve(2, sigma0, e)

       mu0(2) = mu(i)
       sm = solve(2, sigma0, mu0)

       msm = dot_product(mu0, sm)
       mse = dot_product(mu0, se)
       ese = dot_product(e, se)

       det = msm*ese - mse**2
       lambda1 = (mu_min*ese - mse)/det
       if( lambda1 < lambda1_min ) then
          i1_min = i
          lambda1_min = lambda1
       end if
    end do

  end subroutine get_lambda_min

  subroutine get_lambda_max(i1_max, lambda1_max)
    implicit none
    real(8) :: lambda1_max
    real(8) :: lambda1, msm, mse, ese, det, aux
    real(8), dimension(2,2) :: sigma0
    real(8), dimension(2)   :: mu0, e, sm, se, xx

    integer :: i, i1_max

    sigma0(1,1) = sigma(i_min,i_min)
    mu0(1) = mu(i_min)
    e = (/ 1.0d0, 1.0d0 /)

    sigma0(1,1) = sigma(i_max,i_max)
    mu0(1) = mu(i_max)
    e = (/ 1.0d0, 1.0d0 /)

    lambda1_max = -1.0d10
    i1_max = 0
    do i = 1, nn
       if(i==i_max) cycle
       sigma0(1,2) = sigma(i_max,i)
       sigma0(2,1) = sigma0(1,2)
       mu0(2) = mu(i)
       sm = solve(2, sigma0, mu0)
       se = solve(2, sigma0, e)

       msm = dot_product(mu0, sm)
       mse = dot_product(mu0, se)
       ese = dot_product(e, se)

       det = msm*ese - mse**2
       lambda1 = (mu_max*ese - mse)/det
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

  subroutine get_lambda1(iw, k, indw, lambda_c, i_new, lambda_new, &
       Ai_eS, Ai_muS)
    ! find the largest lambda, when one of the formerly non-zero weight
    ! wants to become zero 

    ! iw(i) = 1, i=1..nn  for i in S, and 0 for i not in S
    ! ind(i), i=1..k is the meaning of the indices of AS_inv

    implicit none
    integer, dimension(nn),  intent(in) :: iw
    integer,                 intent(in) :: k
    integer, dimension(k),   intent(in) :: indw
    real(8), intent(in) :: lambda_c
    integer, intent(out) :: i_new
    real(8), intent(out) :: lambda_new
    real(8), dimension(k), intent(out) :: Ai_eS, Ai_muS

    real(8), dimension(k) :: muS, eS, a, c
    real(8) :: alpha, beta, ac, ceS, cmuS, aux_e, aux_mu
    real(8) :: eS_Ai_eS, eS_Ai_muS, muii, aux_ee, aux_emu, lambda_ii
    integer :: i, ii

    eS = 1.0d0
    Ai_eS = matmul( Ai(1:k,1:k), eS)
    eS_Ai_eS = sum(Ai_eS(1:k))

    do i = 1, k
       muS(i) = mu(indw(i))
    end do
    Ai_muS = matmul( Ai(1:k,1:k), muS)
    eS_Ai_muS = sum(Ai_muS(1:k))

    lambda_new = -lambda_aux
    i_new = 0
    do i = 1, k
       ii = indw(i)
       if(ii == ii_last) cycle
       lambda_ii = Ai_eS(i)/(eS_Ai_muS*Ai_eS(i) - eS_Ai_eS*Ai_muS(i))
       if( lambda_ii >= lambda_c) cycle
       if( lambda_ii > lambda_new) then
          lambda_new = lambda_ii
          i_new = ii
       end if
    end do
  end subroutine get_lambda1

  subroutine get_lambda2(iw, k, indw, Ai_eS, Ai_muS, &
       lambda_c, i_new, lambda_new)
    ! find the largest lambda, when one of the formerly zero weight
    ! wants to become positive 

    ! iw(i) = 1, i=1..nn  for i in S, and 0 for i not in S
    ! ind(i), i=1..k is the meaning of the indices of AS_inv

    implicit none
    integer, dimension(nn), intent(in) :: iw
    integer, intent(in) :: k
    integer, dimension(k), intent(in) :: indw
    real(8), dimension(k),   intent(in) :: Ai_eS, Ai_muS
    real(8), intent(in) :: lambda_c
    integer, intent(out) :: i_new
    real(8), intent(out) :: lambda_new

    real(8), dimension(k) :: muS, eS, a, c
    real(8) :: alpha, beta, ac, ceS, cmuS, aux_e, aux_mu
    real(8) :: eS_Ai_eS, eS_Ai_muS, muii, aux_ee, aux_emu, lambda_ii
    real(8) :: aux
    integer :: i, ii

    eS = 1.0d0
    eS_Ai_eS = sum(Ai_eS(1:k))

    eS_Ai_muS = sum(Ai_muS(1:k))

    lambda_new = -lambda_aux
    i_new = 0
    do ii = 1, nn
       if(ii == ii_last) cycle
       if(iw(ii) == 1) cycle
       muii = mu(ii)

       a = sigma(ii,indw(:))
       alpha = sigma(ii,ii)

       ceS = dot_product(a, Ai_eS)
       cmuS = dot_product(a, Ai_muS)

       lambda_ii = (1.0d0 - ceS)/(eS_Ai_muS*(1.0d0 - ceS) - &
            eS_Ai_eS*(muii - cmuS))

       if( lambda_ii >= lambda_c) cycle
       if( lambda_ii > lambda_new) then
          lambda_new = lambda_ii
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
    a = sigma(ii,indw(1:k))
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
         print *,'*** ii0 missing:  ii0=', ii0
       end if
       error_flag = 2
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
       error_flag = 3
       return
       !~ stop
    end if
  end subroutine check_inverse

  subroutine get_muP(iw, k, indw, lambda, muP, sig1)
    ! finds the value of muP at a given lambda in the restricted set

    ! iw(i) = 1, i=1..nn  for i in S, and 0 for i not in S
    ! ind(i), i=1..k is the meaning of the indices of Ai

    implicit none
    integer, dimension(nn),  intent(in) :: iw
    integer,                 intent(in) :: k
    integer, dimension(k),   intent(in) :: indw
    real(8),                 intent(in) :: lambda
    real(8),                intent(out) :: muP, sig1

    real(8), dimension(k) :: muS, eS, Ai_eS, Ai_muS
    real(8) :: eS_Ai_eS, eS_Ai_muS, muS_Ai_muS, aux, gamma
    integer :: i

    eS = 1.0d0
    Ai_eS = matmul( Ai(1:k,1:k), eS)
    eS_Ai_eS = sum(Ai_eS(1:k))

    do i = 1, k
       muS(i) = mu(indw(i))
    end do
    Ai_muS = matmul( Ai(1:k,1:k), muS)
    eS_Ai_muS = sum(Ai_muS(1:k))
    muS_Ai_muS = dot_product(muS(1:k), Ai_muS(1:k))

    muP = lambda*( muS_Ai_muS - eS_Ai_muS**2/eS_Ai_eS ) + eS_Ai_muS/eS_Ai_eS

    gamma = (1.0d0 - lambda*eS_Ai_muS)/eS_Ai_eS
    aux = lambda**2 * muS_Ai_muS  +  gamma**2 * eS_Ai_eS + &
         2.0d0 * lambda * gamma * eS_Ai_muS
    sig1 = sqrt(aux)
  end subroutine get_muP

  subroutine get_weight(iw, k, indw, w0, w1, w0mu, w1mu)
    ! finds the weight: weight = w0 + lambda*w1 
    ! returns w0 and w1 in order to make an interpoolation
    !  muP = w0mu + lambda*w1mu

    ! iw(i) = 1, i=1..nn  for i in S, and 0 for i not in S
    ! ind(i), i=1..k is the meaning of the indices of Ai

    implicit none
    integer, dimension(nn),  intent(in) :: iw
    integer,                 intent(in) :: k
    integer, dimension(k),   intent(in) :: indw
    real(8), dimension(nn), intent(out) :: w0, w1
    real(8),                intent(out) :: w0mu, w1mu

    real(8), dimension(k) :: muS, eS, Ai_eS, Ai_muS
    real(8) :: eS_Ai_eS, eS_Ai_muS, muS_Ai_muS, aux
    integer :: i

    eS = 1.0d0
    Ai_eS = matmul( Ai(1:k,1:k), eS)
    eS_Ai_eS = sum(Ai_eS(1:k))

    do i = 1, k
       muS(i) = mu(indw(i))
    end do
    Ai_muS = matmul( Ai(1:k,1:k), muS)
    eS_Ai_muS = sum(Ai_muS(1:k))
    muS_Ai_muS = dot_product(muS(1:k), Ai_muS(1:k))

    w0 = 0.
    w1 = 0.

    aux = 1.0d0/eS_Ai_eS
    do i = 1, k
       w0(indw(i)) = aux*Ai_eS(i)
    end do

    aux = eS_Ai_muS/eS_Ai_eS
    do i = 1, k
       w1(indw(i)) = Ai_muS(i) - aux*Ai_eS(i)
    end do

    w0mu = eS_Ai_muS/eS_Ai_eS
    w1mu = muS_Ai_muS - eS_Ai_muS**2/eS_Ai_eS

  end subroutine get_weight

  subroutine get_weight1(k, indw, lambda, weight1)
    ! finds the weight: weight = w0 + lambda*w1 

    ! iw(i) = 1, i=1..nn  for i in S, and 0 for i not in S
    ! ind(i), i=1..k is the meaning of the indices of Ai

    implicit none
    integer,                 intent(in) :: k
    integer, dimension(k),   intent(in) :: indw
    real(8),                 intent(in) :: lambda
    real(8), dimension(nn), intent(out) :: weight1

    real(8), dimension(k) :: muS, eS, Ai_eS, Ai_muS
    real(8) :: eS_Ai_eS, eS_Ai_muS, muS_Ai_muS, gamma
    real(8) :: auxv(1), err, aux
    integer :: i

    eS = 1.0d0
    Ai_eS = matmul( Ai(1:k,1:k), eS)
    eS_Ai_eS = sum(Ai_eS(1:k))

    do i = 1, k
       muS(i) = mu(indw(i))
    end do
    Ai_muS = matmul( Ai(1:k,1:k), muS)
    eS_Ai_muS = sum(Ai_muS(1:k))
    muS_Ai_muS = dot_product(muS(1:k), Ai_muS(1:k))

    gamma = (1.0d0 - lambda*eS_Ai_muS)/eS_Ai_eS

    weight1 = 0.
    do i = 1, k 
       weight1(indw(i)) = lambda*Ai_muS(i) + gamma*Ai_eS(i)
    end do

    auxv = minval(weight1)
    err = abs( sum(weight1) - 1.0d0)

    if( auxv(1) < -1.0d-13 .or. err > 1.0d-13) then
      if( print_on /= 0 ) then
        print *,'**** get_weight1: ', auxv, err
      end if
    end if

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

end SUBROUTINE calcturningpoints
