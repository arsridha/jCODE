module stat_ibm_1d

  ! External modules
  use stat
  use string
  
  implicit none

  ! Global variables
  integer :: nStatPoints, nStatPoints2
  real(WP) :: xmin, xmax, dx, ymin, ymax, dy
  real(WP), dimension(:), allocatable :: grid1D, TKE_OLD, TKE_NEW, RHO_OLD, RHO_NEW, U_OLD,  &
       U_NEW
 
end module stat_ibm_1d


! =========================!
! Initialize 1D statistics !
! =========================!
subroutine stat_ibm_1d_setup

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use parser
  use parallel
  use simulation_flags
  use geometry
  use grid

  implicit none

  ! Local variables
  integer :: i
  
  ! Return if IBM is not used
  if (.not. useIBM) return

  ! Get stat size
  call parser_read('stat ibm bin size', nStatPoints, globalGridSize(1))

  ! Allocate arrays
  allocate(grid1D(nStatPoints))
  allocate(TKE_OLD(nStatPoints)); TKE_OLD = 0.0_WP
  allocate(TKE_NEW(nStatPoints)); TKE_NEW = 0.0_WP
  allocate(RHO_OLD(nStatPoints)); RHO_OLD = 0.0_WP
  allocate(RHO_NEW(nStatPoints)); RHO_NEW = 0.0_WP
  allocate(U_OLD(nGridPoints)); U_OLD = 0.0_WP
  allocate(U_NEW(nGridPoints)); U_NEW = 0.0_WP

  ! Generate the 1D grid
  xmin = minval(coordinates(:,1))
  call parallel_min(xmin)
  xmax = maxval(coordinates(:,1))
  call parallel_max(xmax)
  dx = (xmax - xmin) / real(nStatPoints-1, WP)
  do i = 1, nStatPoints
     grid1D(i) = xmin + real(i-1,WP) * dx
  end do
  
  return
end subroutine stat_ibm_1d_setup

! ======================!
! Cleanup 1D statistics !
! ======================!
subroutine stat_ibm_1d_cleanup

  ! Internal modules
  use stat_ibm_1d
  
  implicit none

  if (allocated(grid1D)) deallocate(grid1D)
  if (allocated(TKE_OLD)) deallocate(TKE_OLD)
  if (allocated(TKE_NEW)) deallocate(TKE_NEW)
  if (allocated(RHO_OLD)) deallocate(RHO_OLD)
  if (allocated(RHO_NEW)) deallocate(RHO_NEW)
  if (allocated(U_OLD)) deallocate(U_OLD)
  if (allocated(U_NEW)) deallocate(U_NEW)

  return
end subroutine stat_ibm_1d_cleanup

! ===================== !
! Compute 1D statistics !
! ===================== !
subroutine stat_ibm_1d_compute

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use parallel
  use simulation_flags
  use geometry
  use ibm
  use grid
  use state
  use grid_functions
  use grid_levelset
  use math

  implicit none

  ! Local variables
  integer :: i, j
  real(WP), dimension(nStatPoints, nDimensions) :: favreVelocity
  real(WP), dimension(nStatPoints) :: meanRho, volume, totalVolume

  ! Return if IBM is not used
  if (.not. useIBM) return

  ! Update old value
  TKE_OLD = TKE_NEW
  RHO_OLD = RHO_NEW
  U_OLD = U_NEW

  ! Compute mean stats
  meanRho = 0.0_WP
  favreVelocity = 0.0_WP
  totalVolume = 0.0_WP
  volume = 0.0_WP
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     volume(j) = volume(j) + gridNorm(i,1)
     totalVolume(j) = totalVolume(j) + primitiveGridNorm(i,1)
     meanRho(j) = meanRho(j) + gridNorm(i,1) * conservedVariables(i,1)
     favreVelocity(j,1:nDimensions) = favreVelocity(j,1:nDimensions) + gridNorm(i, 1) *      &
          conservedVariables(i,2:nDimensions+1)
  end do
  call parallel_sum(meanRho)
  call parallel_sum(totalVolume)
  
  do i = 1, nDimensions
     call parallel_sum(favreVelocity(:,i))
     favreVelocity(:,i) = favreVelocity(:,i) / meanRho
  end do
  
  ! Update new <alpha> <rho> K
  TKE_NEW = 0.0_WP
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     TKE_NEW(j) = TKE_NEW(j) + 0.5_WP * conservedVariables(i, 1) * sum((                      &
          velocity(i,1:nDimensions) - favreVelocity(j,1:nDimensions))**2) * gridNorm(i, 1)
     U_NEW(i) = conservedVariables(i,2) 
  end do
  call parallel_sum(TKE_NEW)
  TKE_NEW = TKE_NEW / totalVolume
  !U_NEW = U_NEW            
  
  return
end subroutine stat_ibm_1d_compute


! =================== !
! Read ibm statistics !
! =================== !
subroutine stat_ibm_1d_read

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use parallel
  use ibm

  implicit none

  ! Nothing to do
  
  return
end subroutine stat_ibm_1d_read


! ================== !
! Read 1D statistics !
! ================== !
subroutine stat_ibm_1d_write

  ! Internal modules
  use stat_ibm_1d

  ! External modules
  use time_info
  use parallel
  use fileio
  use math
  use geometry
  use ibm
  use grid
  use state
  use grid_functions
  use grid_levelset
  use first_derivative
  use state_jacobian
  use state_functions

  implicit none

  ! Local variables
  character(len=str_medium) :: filename
  integer  :: i, j, n, iunit, ierror
  real(WP) :: D, F(3), pPrime
  real(WP), dimension(nDimensions) :: uDoublePrime, objectVelocity
  real(WP), dimension(nDimensions**2) :: sigmaPrime
  !real(WP), dimension(nStatPoints, nDimensions) :: meanVelocity, favreVelocity
  real(WP), dimension(nStatPoints, 3) :: meanVelocity, favreVelocity
  real(WP), dimension(nStatPoints, nDimensions**2) :: meanSigma, avgFlux
  real(WP), dimension(nStatPoints) :: meanRho, meanUU, meanVV, meanWW, meanUV, meanUW,       &
       meanVW, alpha, volume, totalVolume, meanP, dragP, dragV, meanUprime, pDivU, sigmaDivU,&
       dkdt, mean_gradP, mean_gradV, uPrimeP, uPrimeSigma, UdRhoUdt, VdrhoVdt, avgdpUdx,     &
       udPIdx, uprime_dpdx, avgdsigmaUdx, convectiveterm, avgrhoUUU, dragP_II, x_pdivU,      &
       dragP_III, rhs_term_avg, numericalDissipation, avginviscFlux, avgviscFlux, mean_sound,&
       mean_Viscosity, meanT, forceP, forceV, pprimeUprime
  logical :: procHasObject
  real(WP), dimension(nGridPoints, nUnknowns) :: rhs_terms, dissipationSource
  real(WP), dimension(nGridPoints,nDimensions) :: gradI, gradP, ptimesI, Utilde, dpUdx, pU,  &
       uSigma, dpdx, gradUprime, Udoubleprime_store, dsigmaUdx, convective_term2, gradrhoUUU
  real(WP), dimension(nGridPoints, nUnknowns, nDimensions) :: fluxes1, fluxes2, nfluxes1,    &
       nfluxes2, nfluxes1temp, nfluxes2temp
  real(WP), dimension(nGridPoints, nDimensions, nDimensions) :: pressureStress, temp1,       &
       temp2, viscousStress
  real(WP), dimension(nGridPoints) :: drhoUdt, drhovdt, dkdx, rhoUUU
  real(WP), dimension(nGridPoints, nDimensions**2) :: gradUtilde

  ! Return if IBM is not used
  if (.not. useIBM) return
  
  !Compute gradients of pressure
  ptimesI(:,1) = pressure(:,1) * indicatorFunction(:,1)
  call gradient(ptimesI(:,1), gradP)
  call gradient(pressure(:,1), dpdx)

  call state_rhs(FORWARD, rhs_terms)
 
 ! Compute convective term
  select case (nDimensions)
  case (1)
     fluxes1(:,1,1) = conservedVariables(:,2) * velocity(:,1)
  case (2)
     fluxes1(:,1,1) = conservedVariables(:,2) * velocity(:,1)
     fluxes1(:,2,1) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,1,2) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,2,2) = conservedVariables(:,3) * velocity(:,2)
     pressureStress(:,1,1) = -pressure(:,1) 
     pressureStress(:,2,2) = -pressure(:,1) 
     viscousStress(:,1,1) = stressTensor(:,1)
     viscousStress(:,1,2) = stressTensor(:,2)
     viscousStress(:,2,1) = stressTensor(:,3)
     viscousStress(:,2,2) = stressTensor(:,4)
   case (3)
     fluxes1(:,1,1) = conservedVariables(:,2) * velocity(:,1)
     fluxes1(:,1,2) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,1,3) = conservedVariables(:,2) * velocity(:,3)
     fluxes1(:,2,1) = conservedVariables(:,2) * velocity(:,2)
     fluxes1(:,2,2) = conservedVariables(:,3) * velocity(:,2)
     fluxes1(:,2,3) = conservedVariables(:,3) * velocity(:,3)
     fluxes1(:,3,1) = conservedVariables(:,2) * velocity(:,3)
     fluxes1(:,3,2) = conservedVariables(:,3) * velocity(:,3)
     fluxes1(:,3,3) = conservedVariables(:,4) * velocity(:,3)
     pressureStress(:,1,1) = -pressure(:,1)
     pressureStress(:,2,2) = -pressure(:,1)
     pressureStress(:,3,3) = -pressure(:,1)
     viscousStress(:,1,1) = stressTensor(:,1)
     viscousStress(:,1,2) = stressTensor(:,2)
     viscousStress(:,1,3) = stressTensor(:,3)
     viscousStress(:,2,1) = stressTensor(:,4)
     viscousStress(:,2,2) = stressTensor(:,5)
     viscousStress(:,2,3) = stressTensor(:,6)
     viscousStress(:,3,1) = stressTensor(:,7)
     viscousStress(:,3,2) = stressTensor(:,8)
     viscousStress(:,3,3) = stressTensor(:,9)
  end select
  call transform_fluxes(fluxes1, fluxes2)
  do i = 1, nDimensions
     call first_derivative_apply(i, fluxes2(:,:,i))
  end do
  
  ! Transform stress tensor from Cartesian to contravariant form
  call transform_fluxes(pressureStress, temp1)
  call transform_fluxes(viscousStress, temp2)

  ! Take derivatives of stress tensor
  do i = 1, nDimensions
     call first_derivative_apply(i, temp1(:,:,i))
     call first_derivative_apply(i, temp2(:,:,i))
  end do
  ! Compute gradient of indicator Function
  gradI=0.0_WP
  !call gradient(indicatorFunction(:,1),gradI)

  dissipationSource = 0.0_WP
  ! Compute dissipation
  call add_dissipation(FORWARD,dissipationSource)
  ! Compute boundary source term
  call boundary_sources(FORWARD, dissipationSource)

  ! Check
  ! Compute inviscid flux
  call compute_cartesian_inviscid_fluxes(conservedVariables, velocity, pressure(:,1),        &
       nfluxes1)
  ! Compute viscous fluc
  call compute_cartesian_viscous_fluxes(velocity, stressTensor, heatFlux,                    &
       nfluxes2, speciesFlux, enthalpyFlux)
  call transform_fluxes(nfluxes1, nfluxes1temp)
  call transform_fluxes(nfluxes2, nfluxes2temp)
  ! Take derivatives of fluxes
  do i = 1, nDimensions
     call first_derivative_apply(i, nfluxes1temp(:,:,i))
     call first_derivative_apply(i, nfluxes2temp(:,:,i))
  end do
    
  ! Compute mean stats
  meanRho = 0.0_WP
  meanVelocity = 0.0_WP
  favreVelocity = 0.0_WP
  totalVolume = 0.0_WP
  volume = 0.0_WP
  meanSigma = 0.0_WP
  meanP = 0.0_WP
  mean_gradP = 0.0_WP
  mean_gradV = 0.0_WP
  mean_sound = 0.0_WP
  mean_Viscosity = 0.0_WP
  meanT = 0.0_WP
  do i = 1, nGridPoints
    j = 1 + nint((coordinates(i,1) - xmin) / dx)
     ! Sum up grid volumes
     volume(j) = volume(j) + gridNorm(i,1)
     totalVolume(j) = totalVolume(j) + primitiveGridNorm(i,1)
     ! Mean stats
     meanRho(j) = meanRho(j) + gridNorm(i,1) * conservedVariables(i,1)                     
     meanVelocity(j,1:nDimensions) = meanVelocity(j,1:nDimensions) +                         &
          gridNorm(i,1) * velocity(i,1:nDimensions)
     favreVelocity(j,1:nDimensions) = favreVelocity(j,1:nDimensions) +                       &
          gridNorm(i,1) * conservedVariables(i,2:nDimensions+1)
     meanP(j) = meanP(j) + gridNorm(i,1) * pressure(i,1)
     meanSigma(j,:) = meanSigma(j,:) + gridNorm(i,1) * stressTensor(i,:)
     ! Average of gradient of P and V
     mean_gradP(j) = mean_gradP(j) + temp1(i,1,1) * gridNorm(i,1) * jacobian(i,1)
     mean_gradV(j) = mean_gradV(j) + temp2(i,1,1) * gridNorm(i,1) * jacobian(i,1)

     ! Mean sound speed
     mean_sound(j) = mean_sound(j) + gridNorm(i,1) * sqrt(ratioOfSpecificHeats *             &
          specificVolume(i,1) * pressure(i,1))
     mean_Viscosity(j)  = mean_Viscosity(j) + gridNorm(i,1) * dynamicViscosity(i,1)
  end do
  
  ! Sum them over procs
  call parallel_sum(volume)
  call parallel_sum(totalVolume)
  call parallel_sum(meanRho)
  call parallel_sum(meanP)
  call parallel_sum(mean_gradP)
  call parallel_sum(mean_gradV)
  call parallel_sum(mean_sound)
  call parallel_sum(mean_Viscosity)
  call parallel_sum(meanT)
  do i = 1, nDimensions
     call parallel_sum(meanVelocity(:,i))
     call parallel_sum(favreVelocity(:,i))
  end do
  do i = 1, nDimensions**2
     call parallel_sum(meanSigma(:,i))
  end do

  ! Normalize
  do i = 1, nStatPoints
     alpha(i) = volume(i) / totalVolume(i)
     meanRho(i) = meanRho(i) / volume(i)
     meanP(i) = meanP(i) / totalVolume(i)                 ! Changing here for momentum
     mean_gradP(i) = mean_gradP(i) / totalVolume(i)
     mean_gradV(i) = mean_gradV(i) / totalVolume(i)
     meanVelocity(i,:) = meanVelocity(i,:) / volume(i)
     favreVelocity(i,:) = favreVelocity(i,:) / volume(i) / meanRho(i)
     meanSigma(i,:) = meanSigma(i,:) / totalvolume(i)
     mean_sound(i) = mean_sound(i) / volume(i)
     mean_Viscosity(i) = mean_Viscosity(i)/ volume(i)
     meanT(i) = meanT(i)/volume(i)
  end do
 
  drhoUdt = (U_NEW - U_OLD ) / timeStepSize
  drhovdt = (U_NEW - U_OLD)  / timeStepSize

  ! Compute variance
  meanUU = 0.0_WP; meanVV = 0.0_WP; meanWW = 0.0_WP
  meanUV = 0.0_WP; meanUW = 0.0_WP; meanVW = 0.0_WP
  meanUprime = 0.0_WP             
  meanUPrime = 0.0_WP; sigmaDivU = 0.0_WP
  pDivU = 0.0_WP
  dragP = 0.0_WP; dragV = 0.0_WP
  UdrhoUdt=0.0_WP; VdRhoVdt=0.0_WP 
  udpIdx = 0.0_WP
  uprime_dpdx = 0.0_WP
  uPrimeP = 0.0_WP
  uPrimeSigma = 0.0_WP
  dragP_II = 0.0_WP
  x_pdivU = 0.0_WP
  dragP_III = 0.0_WP
  avgFlux = 0.0_WP
  rhs_term_avg = 0.0_WP
  numericalDissipation = 0.0_WP
  avginviscFlux = 0.0_WP
  avgviscFlux = 0.0_WP
  forceP = 0.0_WP; forceV = 0.0_WP
  pprimeUprime = 0.0_WP
  ! Get the gradient of favre Velocity
  do i=1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     Utilde(i,1:nDimensions) = favreVelocity(j,1:nDimensions)
  end do
  call gradient(Utilde, gradUtilde)
  ! Compute higher order statistics
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     ! Get Reynolds & Favre fluctuations
     uDoublePrime = velocity(i,1:nDimensions) - favreVelocity(j,1:nDimensions)
     pPrime = pressure(i,1) - meanP(j)
     sigmaPrime = stressTensor(i,:) - meanSigma(j,:)

     Udoubleprime_store(i,1) = uDoublePrime(1)
     
     ! Reynolds stresses
     meanUU(j) = meanUU(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1)**2
     if (nDimensions .gt. 1) then
        meanVV(j) = meanVV(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(2)**2
        meanUV(j) = meanUV(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1) *  &
             uDoublePrime(2)
     end if
     if (nDimensions .gt. 2) then
        meanWW(j) = meanWW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(3)**2
        meanUW(j) = meanUW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(1) *  &
             uDoublePrime(3)
        meanVW(j) = meanVW(j) + gridNorm(i,1) * conservedVariables(i,1) * uDoublePrime(2) *  &
             uDoublePrime(3)
     end if
     rhs_term_avg(j) = rhs_term_avg(j) + gridNorm(i,1) * rhs_terms(i,2) * uDoublePrime(1)
     ! Dissipation
     numericalDissipation(j) = numericalDissipation(j) + gridNorm(i,1) * (sum(uDoublePrime * &
          dissipationSource(i,2:nDimensions+1))  + sum(uDoublePrime**2 * 0.5_WP *           &
          dissipationSource(i,1)))  

     ! Other terms
     meanUprime(j) = meanUprime(j) + uDoublePrime(1) * gridNorm(i,1)
     pU(i,1) = pPrime * uDoublePrime(1) * indicatorFunction(i,1)                                   
     uSigma(i,1) = sum(uDoublePrime(1) * sigmaPrime(1:nDimensions)) * indicatorFunction(i,1)
     rhoUUU(i) = conservedVariables(i,1) * sum(uDoublePrime ** 2) * uDoublePrime(1)  *       &
          indicatorFunction(i,1)

     pDivU(j) = pDivU(j) + pPrime * (velocityGradient(i,1) - gradUtilde(i,1)) * gridNorm(i,1)
     sigmaDivU(j) = sigmaDivU(j) + ( sum(sigmaPrime * velocityGradient(i,:)) -               &
          sum(sigmaPrime * gradUtilde(i,:)) ) * gridNorm(i,1)

     dragV(j) = dragV(j) - primitiveGridNorm(i,1) * (uDoublePrime(1) * gradI(i,1) *          &
          sum(sigmaPrime(1:nDimensions)))
     forceV(j) = forceV(j) - primitiveGridNorm(i,1) * (gradI(i,1) *                          &
          sum(sigmaPrime(1:nDimensions)))
     !sigmaDivU(l) = sigmaDivU(l) + ((sigmaPrime(1) *                                        &
     !     velocityGradient(i,1)) - sigmaPrime(1) * meanDuDx(l)) *                           &
     !     gridNorm(i,1) 
     select case (nDimensions)
     case (1)
        avgFlux(j,1) = avgFlux(j,1) + uDoublePrime(1) * fluxes2(i,1,1) * gridNorm(i,1) * jacobian(i,1) 
     case (2)
        pDivU(j) = pDivU(j) + pPrime * (velocityGradient(i,4) - gradUtilde(i,4)) * gridNorm(i,1)
        dragV(j) = dragV(j) - primitiveGridNorm(i,1) * (uDoublePrime(2) * gradI(i,2) *       &
             sum(sigmaPrime(3:4)))
        forceV(j) =  forceV(j) - primitiveGridNorm(i,1) * (gradI(i,2) *                      &
          sum(sigmaPrime(3:4)))
        !avgFlux(j,1) = avgFlux(j,1) + uDoublePrime(1) * fluxes2(i,1,1) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,2) = avgFlux(j,2) + uDoublePrime(1) * fluxes2(i,1,2) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,3) = avgFlux(j,3) + uDoublePrime(2) * fluxes2(i,2,1) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,4) = avgFlux(j,4) + uDoublePrime(2) * fluxes2(i,2,2) * gridNorm(i,1) * jacobian(i,1)
        !sigmaDivU(l) = sigmaDivU(l) + (2.0_WP*(sigmaPrime(2) *                              &
        !  velocityGradient(i,2)) - sigmaPrime(1) * meanDuDx(l)) *                           &
        !  gridNorm(i,1)

     case (3)
        pDivU(j) = pDivU(j) + pPrime * (velocityGradient(i,5) + velocityGradient(i,9) -        &
             gradUtilde(i,5) -  gradUtilde(i,9) ) * gridNorm(i,1)
        dragV(j) = dragV(j) - primitiveGridNorm(i,1) * ( (uDoublePrime(2) * gradI(i,2) *          &
             sum(sigmaPrime(4:6))) + (uDoublePrime(3) * gradI(i,3) * sum(sigmaPrime(7:9))) )
        forceV(j) = forceV(j) - primitiveGridNorm(i,1) * (gradI(i,2) *                          &
          sum(sigmaPrime(4:6))+ gradI(i,3) * sum(sigmaPrime(7:9)))
        avgFlux(j,1) = avgFlux(j,1) + uDoublePrime(1) * fluxes2(i,1,1) * gridNorm(i,1) * jacobian(i,1)
        avgFlux(j,2) = avgFlux(j,2) + uDoublePrime(1) * fluxes2(i,1,2) * gridNorm(i,1) * jacobian(i,1)
        avgFlux(j,3) = avgFlux(j,3) + uDoublePrime(1) * fluxes2(i,1,3) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,4) = avgFlux(j,4) + uDoublePrime(2) * fluxes2(i,2,1) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,5) = avgFlux(j,5) + uDoublePrime(2) * fluxes2(i,2,2) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,6) = avgFlux(j,6) + uDoublePrime(2) * fluxes2(i,2,3) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,7) = avgFlux(j,7) + uDoublePrime(3) * fluxes2(i,3,1) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,8) = avgFlux(j,8) + uDoublePrime(3) * fluxes2(i,3,2) * gridNorm(i,1) * jacobian(i,1)
        !avgFlux(j,9) = avgFlux(j,9) + uDoublePrime(3) * fluxes2(i,3,3) * gridNorm(i,1) * jacobian(i,1)

     end select

     ! Get velocity of associated object
     ! if (ibm_move) then
     !    n = objectIndex(i)
     !    objectVelocity = object(n)%velocity(1:nDimensions)
     ! else
     !    objectVelocity = 0.0_WP
     ! end if
        
     ! Drag production
     dragP(j) = dragP(j) + pPrime * primitiveGridNorm(i,1) * sum(uDoublePrime * gradI(i,:))
     forceP(j) = forceP(j) + pPrime * primitiveGridNorm(i,1) * sum(gradI(i,:))
     !dragV(j) = dragV(j) - primitiveGridNorm(i,1) * (uDoublePrime(1) * gradI(i,1) *          &
     !     sum(sigmaPrime(1:nDimensions)) + gradI(i,2) * uDoublePrime(2) *                    &
     !     sum(sigmaPrime(3:4)))                   

     VdRhoVdt(j) = VdRhoVdt(j) + drhoVdt(i) * primitiveGridNorm(i,1)
     UdrhoUdt(j) = UdrhoUdt(j) + uDoublePrime(1) * drhoUdt(i) * gridNorm(i,1)

     ! Stage I
     uPrimeP(j) = uPrimeP(j) + ( uDoublePrime(1) * dpdx(i,1) )  * gridNorm(i,1)
     uPrimeSigma(j) =  uPrimeSigma(j) + uDoubleprime(1) * sum(temp2(i,1,1:nDimensions) *     &
          jacobian(i,1)) * gridNorm(i,1)
     
     ! Stage II
     udPIdx(j) = udPIdx(j) + uDoubleprime(1) * gradP(i,1) * primitiveGridNorm(i,1)
     dragP_II(j) = dragP_II(j) + pressure(i,1) * primitiveGridNorm(i,1) * (uDoublePrime(1) * &
          gradI(i,1))

     ! Stage III
     uprime_dpdx(j) = uPrime_dpdx(j) + uDoubleprime(1) * mean_gradP(j) * gridNorm(i,1)
     x_pdivU(j) = x_pdivU(j) + pPrime * (velocityGradient(i,1) - gradUtilde(i,1)) *          &
          gridNorm(i,1)
     dragP_III(j) = dragP_III(j) + pPrime * primitiveGridNorm(i,1) * (uDoublePrime(1) *      &
          gradI(i,1))

     dkdx(i) = indicatorFunction(i,1) * conservedVariables(i,1) * favreVelocity(j,1) *       &
          0.5_WP * sum(uDoublePrime**2)

     !pprimeUprime(i) = pprimeUprime(i) + pPrime * sum(uDoublePrime) * gridNorm(i,1) 
     
  end do
  call parallel_sum(meanUU); meanUU = meanUU / meanRho / volume
  call parallel_sum(meanVV); meanVV = meanVV / meanRho / volume
  call parallel_sum(meanWW); meanWW = meanWW / meanRho / volume
  call parallel_sum(meanUV); meanUV = meanUV / meanRho / volume
  call parallel_sum(meanUW); meanUW = meanUW / meanRho / volume
  call parallel_sum(meanVW); meanVW = meanVW / meanRho / volume
  call parallel_sum(meanUprime); meanUprime = meanUprime / volume

  call parallel_sum(pDivU); pDivU = pDivU / totalVolume                                                    
  call parallel_sum(sigmaDivU); sigmaDivU = sigmaDivU / totalVolume
  call parallel_sum(UdRhoUdt); UdRhoUdt = UdRhoUdt / totalVolume
  call parallel_sum(VdRhoVdt); VdRhoVdt = VdRhoVdt / totalVolume

  call parallel_sum(uPrimeP); uPrimeP = uPrimeP / totalVolume
  call parallel_sum(uPrimeSigma); uPrimeSigma = uPrimeSigma / totalVolume
  call parallel_sum(udPIdx); udPIdx = udPIdx / totalVolume
  call parallel_sum(uPrime_dpdx); uprime_dpdx = uprime_dpdx /totalvolume
  call parallel_sum(dragP_II); dragP_II = dragP_II / totalVolume
  call parallel_sum(x_pDivU); x_pDivU =x_pDivU / totalVolume
  call parallel_sum(dragP_III); dragP_III =dragP_III / totalVolume
  call parallel_sum(rhs_term_avg); rhs_term_avg = rhs_term_avg / totalVolume
  call parallel_sum(numericalDissipation); numericalDissipation = numericalDissipation /     &
       totalVolume
  call parallel_sum(pprimeUprime); pprimeUprime = pprimeUprime / volume
  
  ! Compute gradients of transport terms
  call gradient(Udoubleprime_store(:,1), gradUprime)
  call gradient(pU(:,1), dpUdx)
  call gradient(uSigma(:,1), dsigmaUdx)
  call gradient(dkdx, convective_term2)
  call gradient(rhoUUU, gradrhoUUU)

  ! Normalise
  do i =1,nDimensions**2
     call parallel_sum(avgFlux(:,i))
  end do
  do i =1,nStatpoints
     avgFlux(i,:) = avgFlux(i,:) / totalvolume(i)
  end do
  avgdpUdx = 0.0_WP
  avgdsigmaUdx = 0.0_WP
  convectiveterm = 0.0_WP
  avgrhoUUU = 0.0_WP
  do i=1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     avgdpUdx(j) = avgdpUdx(j) + dpUdx(i,1) *primitiveGridNorm(i,1)
     avgdsigmaUdx(j) = avgdsigmaUdx(j) + dsigmaUdx(i,1) *primitiveGridNorm(i,1)
     convectiveterm(j) = convectiveterm(j) + convective_term2(i,1) * primitiveGridNorm(i,1)
     avgrhoUUU(j) = avgrhoUUU(j) + 0.5_WP * gradrhoUUU(i,1) * primitiveGridNorm(i,1)
  end do

  call parallel_sum(avgdpUdx); avgdpUdx = avgdpUdx / totalVolume
  call parallel_sum(convectiveterm); convectiveterm = convectiveterm/totalVolume
  call parallel_sum(avgrhoUUU); avgrhoUUU = avgrhoUUU / totalVolume
  call parallel_sum(dragP); dragP = dragP / totalVolume
  call parallel_sum(dragV); dragV = dragV / totalVolume
  call parallel_sum(avgdsigmaUdx); avgdsigmaUdx = avgdsigmaUdx / totalVolume
  call parallel_sum(forceP); forceP = forceP / totalVolume
  call parallel_sum(forceV); forceV = forceV / totalVolume
  ! Root process writes
  if (iRank .eq. iRoot) then

     ! ------- WRITE THE ASCI FILE ------
     filename = "stat/stat-ibm"
     write(filename, '(2A,I8.8,A)') trim(filename), "-", timestep, ".txt"
     iunit = iopen()
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(A,1ES20.12)') "t=",time
     write(iunit,'(5a20)') 'X','D','Fx','Fy','Fz'
     F = 0.0_WP
     do n = 1, nObjects
        D = (2.0_WP * real(nDimensions, WP) * object(n)%volume / pi)                         &
             ** (1.0_WP / real(nDimensions, WP))
        F(1:nDimensions) =  object(n)%pForce(1:nDimensions) +  object(n)%vForce(1:nDimensions)
        write(iunit,'(10000ES20.12)') object(n)%position(1), D, F(1), F(2), F(3)
     end do
     close(iclose(iunit)) 
     !UdRhoUdt = (RHO_NEW - RHO_OLD) / timeStepSize
     dkdt= (TKE_NEW - TKE_OLD) / (timeStepSize)
     filename= "stat/stat-tke"
     write(filename, '(2A,I8.8,A)') trim(filename), "-", timestep, ".txt"
     iunit = iopen()
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(A,1ES21.12)') "t=",time
     write(iunit,'(50a21)') "X", "alpha", "Volume", "<U>", "<V>", "<W>", "<U'U'>", "<V'V'>", & 
           "<W'W'>", "<U'V'>", "<V'W'>", "<U'W'>", "dragP", "dragV", "<rho>", "<rhoUUU>",    &
           "pU", "uSigma", "pdivU", "sigmaDivU", "<P>", "<Sigma>", "<u''>", "Utilde",        &
           "Vtilde", "Wtilde", "dkdt", "flux_x", "flux_y", "<dP/dx>","dV/dx", "udrhoUdt",    &
           "uPrimeP", "uPrimesigma", "<dpUdx>", "uPdIdx", "uprime_dpdx", "VdrhoVdt",         &
           "avgdsigmaUdx", "convection", "avgrhoUUU", "dragP_II", "x_pDivU", "dragP_III",    &
           "rhs_term_avg", "Dissipation", "Sound", "Viscosity", "Temperature",      & ! forceP, forceV
            "p'u''"
     do i = 1, nStatPoints
        write(iunit,'(10000ES21.12E3)') grid1D(i), alpha(i), totalVolume(i), meanVelocity(i,1),   &
             meanVelocity(i,2), meanVelocity(i,3), meanUU(i), meanVV(i), &
             meanWW(i), meanUV(i), meanVW(i), meanUW(i),  &
             dragP(i), dragV(i), meanRho(i), rhoUUU(i),     &
             avgFlux(i,3), avgFlux(i,4), pDivU(i), sigmaDivU(i), &
             meanP(i), meanSigma(i,1), meanUPrime(i), &
             favreVelocity(i,1), favreVelocity(i,2), favreVelocity(i,3), dkdt(i), &
             avgFlux(i,1), avgFlux(i,2),       &
             mean_gradP(i), mean_gradV(i), &
             UdrhoUdt(i), uPrimeP(i),           &
             uPrimeSigma(i), avgdpUdx(i), udPIdx(i), uprime_dpdx(i),      &
             VdrhoVdt(i), avgdsigmaUdx(i), convectiveterm(i), &
             avgrhoUUU(i), dragP_II(i), x_pDivU(i), dragP_III(i), &
             rhs_term_avg(i), numericalDissipation(i), mean_sound(i), mean_Viscosity(i),&
             meanT(i), pprimeUprime(i) ! forceP(i), forceV(i), 
     end do
     close(iclose(iunit))
     
  end if

  ! Log
  call monitor_log("IBM STATISTICS FILE WRITTEN")
  
  return
end subroutine stat_ibm_1d_write

