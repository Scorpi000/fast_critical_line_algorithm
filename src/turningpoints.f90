! mex turningpoints.f90 calcturningpointsgen.f90 calcturningpoints.f90

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  use mexf90                    ! API function definitions
  ! calculates turning points of the efficient frontier
  ! 
  ! calling from matlab
  ! [mu_tp, sigma_tp, weights_tp, lambda_tp] = turningpointsgen(mu, sigma[, lower, upper][, max_tps])
  ! dimensions of arguments:
  ! mu_tp, sigma_tp, lambda_tp: 1*T
  ! weights_tp: n*T
  ! mu, lower, upper: n*1
  ! Sigma: n*n
  ! max_tps: scalar
  !
  ! n: number of assets
  ! T: number of turning points
  ! max_tps: maximal space allocated for turning points, default: 2*n, if T exceeds max_tps 
  !          the function is interrupted and an error message returned
  ! mu_tp: mu at each turning point
  ! sigma_tp: sigma at each turning point
  ! weights_tp: weights at each turning point
  ! lambda_tp: lambda at each turning point
  ! mu: returns
  ! sigma: variance-covariance matrix
  ! lower, upper: lower and upper bounds for assets

!  implicit none

  integer, intent(in) :: nlhs, nrhs ! number of left (right) hand side arguments
  integer, intent(in), dimension(*) :: prhs ! pointer to ...
  integer, intent(out), dimension(*) :: plhs
  integer :: error_flag


  integer :: nn,n_tps
  INTEGER :: max_tps
  integer, pointer :: mu,Sigma,mm,ss,ll,xx,lower,upper
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mu_tps, sig_tps, lambda_tps
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: weight_tps
  INteger :: bounds_given, max_tps_given
  double precision :: sum_of_bounds
  
  character(len=150) :: text

  ! Check input arguments
  if(nrhs < 2 .OR. nrhs > 5) then
     call mexErrMsgTxt('Function turningpoints: Wrong number of input arguments, required arguents: &
                        &mu, Sigma[, lower, upper][, max_tps].');
  end if
  
  if(nrhs == 4 .OR. nrhs == 5) then
    bounds_given = 1
  else
    bounds_given = 0
  end if
  if(nrhs == 3 .OR. nrhs == 5) then
    max_tps_given = 1
  else
    max_tps_given = 0
  end if

  mu => mxGetPr(prhs(1))       
  ! Get data and size of the input matrix
  nn = mxGetM(prhs(1))
  Sigma => mxGetPr(prhs(2))
  if( bounds_given /= 0 ) then
    lower => mxGetPr(prhs(3))
    upper => mxGetPr(prhs(4))
  end if
  IF(mxGetN(prhs(2)) /= nn .OR. mxGetM(prhs(2)) /= nn .OR. mxGetN(prhs(1)) /= 1 ) THEN
     CALL mexErrMsgTxt('Function turningpoints: Matrix size mismatch in arguments.')
     !~ CALL mexErrMsgTxt('text2')
  END IF
  if( bounds_given /= 0 ) then
    IF( mxGetM(prhs(3)) /= nn .OR. mxGetM(prhs(4)) /= nn .OR. mxGetN(prhs(3)) /= 1 .OR. mxGetN(prhs(4)) /= 1 ) THEN
      CALL mexErrMsgTxt('Function turningpoints: Lower or upper bound has wrong dimensions.')
    end if
    call specialsum(lower, nn, sum_of_bounds)
    if(  sum_of_bounds > 1) then
      call mexErrMsgTxt('Function turningpoints: Portfolio not feasible: sum of lower bounds > 1')
    end if
    call specialsum(upper, nn, sum_of_bounds)
    if( sum_of_bounds < 1) then
      call mexErrMsgTxt('Function turningpoints: Portfolio not feasible: sum of upper bounds < 1')
    end if
  end if
  IF(max_tps_given == 0) THEN
     max_tps = 2*nn
  ELSE
    if( bounds_given /= 0 ) then
      max_tps = mxGetScalar(prhs(5))
    else
      max_tps = mxGetScalar(prhs(5))
    end if
  END IF
  
  ALLOCATE(mu_tps(max_tps), sig_tps(max_tps), lambda_tps(max_tps), weight_tps(nn,max_tps))
  !~ CALL dummyturning(mu, Sigma, nn, mu_tps, sig_tps, weight_tps, lambda_tps, n_tps, error_flag)
  if( bounds_given /= 0 ) then
    CALL calcturningpointsgen(mu, Sigma, nn, mu_tps, sig_tps, weight_tps, lambda_tps, n_tps, error_flag, max_tps, lower, upper)
  else
    CALL calcturningpoints(mu, Sigma, nn, mu_tps, sig_tps, weight_tps, lambda_tps, n_tps, error_flag, max_tps)
  end if
  
  if( error_flag == 4 ) then
    call mexErrMsgTxt('Function turningpoints: Error during execution: maximal space allocated for turning points exceeded.')
  elseif( error_flag /= 0 ) then
    !~ call mexErrMsgTxt('Unknown error. (Error code ' // achar(48+error_flag) // ').') !achar(i) ichar(c)
    write(text,'("Function turningpoints: Unknown error. (Error code ",i2,")")') error_flag
    call mexErrMsgTxt(text)
  endif
  
  
  ! Create output matrix
  if(nlhs >= 1) then
    plhs(1) = mxCreateDoubleMatrix(1,n_tps,0)
    mm=>mxGetPr(plhs(1))
  endif
  if(nlhs >= 2) then
    plhs(2) = mxCreateDoubleMatrix(1,n_tps,0)
    ss=>mxGetPr(plhs(2))
  endif
  if(nlhs >= 3) then
    plhs(3) = mxCreateDoubleMatrix(nn,n_tps,0)
    xx=>mxGetPr(plhs(3))
  endif
  if(nlhs >= 4) then
    plhs(4) = mxCreateDoubleMatrix(1,n_tps,0)
    ll=>mxGetPr(plhs(4))
  endif
  
    
  CALL assignvalues(mm, ss, xx, ll, mu_tps, sig_tps, weight_tps, lambda_tps, n_tps, nn, nlhs)
  deallocate(mu_tps, sig_tps, lambda_tps, weight_tps)

end subroutine mexFunction

subroutine dummyturning(mu, sigma, nn, mu_tps, sig_tps, weight_tps, lambda_tps, n_tps, error_flag)

  implicit none

  integer, INTENT(IN) :: nn
  integer :: error_flag
  DOUBLE PRECISION, DIMENSION(nn), INTENT(IN) :: mu
  DOUBLE PRECISION, DIMENSION(nn,nn), INTENT(IN) :: sigma
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(2*nn) :: mu_tps, sig_tps, lambda_tps
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(nn,2*nn) :: weight_tps
  INTEGER, INTENT(OUT) :: n_tps
  
  INTEGER :: i, j

  error_flag = 0
  ! Now do something useful...
  n_tps=3
  DO i=1, n_tps
    mu_tps(i)=2*i
    sig_tps(i)=3*i
    lambda_tps(i)=4*i
    DO j=1, nn
      weight_tps(j,i) = i+j
    END DO
  END DO

end subroutine dummyturning

SUBROUTINE assignvalues(mm, ss, xx, ll, mu_tps, sig_tps, weight_tps, lambda_tps, n_tps, nn, nlhs)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n_tps, nn, nlhs
  DOUBLE PRECISION, DIMENSION(n_tps) :: mm, ss, ll, mu_tps, sig_tps, lambda_tps
  DOUBLE PRECISION, DIMENSION(nn,n_tps) :: xx
  DOUBLE PRECISION, DIMENSION(nn,2*nn) :: weight_tps
  
  if( nlhs >= 1 ) then
    mm(1:n_tps)=mu_tps(1:n_tps)
  endif
  if( nlhs >= 2 ) then
    ss(1:n_tps)=sig_tps(1:n_tps)
  endif
  if( nlhs >= 3 ) then
    xx(1:nn,1:n_tps)=weight_tps(1:nn,1:n_tps)
  endif
  if( nlhs >= 4 ) then
    ll(1:n_tps)=lambda_tps(1:n_tps)
  endif

END SUBROUTINE assignvalues

subroutine specialsum(mexpointer, nn, retvalue)
  implicit none
  integer, intent(in) :: nn
  double precision, dimension(nn) :: mexpointer
  double precision, intent(out) :: retvalue
  
  retvalue = sum(mexpointer)
end subroutine specialsum
