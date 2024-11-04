module stat_el_ibm

  ! External modules
  use stat
  use string
  
  implicit none

  ! Global variables
  integer :: nStatPoints
  real(WP) :: xmin, xmax, dx, ymin, ymax, dy
  real(WP), dimension(:), allocatable :: grid1D
 
end module stat_el_ibm


! =========================!
! Initialize 1D statistics !
! =========================!
subroutine stat_el_ibm_setup

  ! Internal modules
  use stat_el_ibm

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
  if (.not. useIBM .and. .not. useParticles) return

  ! Get stat size
  call parser_read('stat ibm bin size', nStatPoints, globalGridSize(1))

  ! Allocate arrays
  allocate(grid1D(nStatPoints))
!   allocate(TKE_OLD(nStatPoints)); TKE_OLD = 0.0_WP
!   allocate(TKE_NEW(nStatPoints)); TKE_NEW = 0.0_WP
!   allocate(RHO_OLD(nStatPoints)); RHO_OLD = 0.0_WP
!   allocate(RHO_NEW(nStatPoints)); RHO_NEW = 0.0_WP
!   allocate(U_OLD(nGridPoints)); U_OLD = 0.0_WP
!   allocate(U_NEW(nGridPoints)); U_NEW = 0.0_WP

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
end subroutine stat_el_ibm_setup

! ======================!
! Cleanup 1D statistics !
! ======================!
subroutine stat_el_ibm_cleanup

  ! Internal modules
  use stat_el_ibm
  
  implicit none

  if (allocated(grid1D)) deallocate(grid1D)
!   if (allocated(TKE_OLD)) deallocate(TKE_OLD)
!   if (allocated(TKE_NEW)) deallocate(TKE_NEW)
!   if (allocated(RHO_OLD)) deallocate(RHO_OLD)
!   if (allocated(RHO_NEW)) deallocate(RHO_NEW)
!   if (allocated(U_OLD)) deallocate(U_OLD)
!   if (allocated(U_NEW)) deallocate(U_NEW)

  return
end subroutine stat_el_ibm_cleanup

! ===================== !
! Compute 1D statistics !
! ===================== !
subroutine stat_el_ibm_compute

  ! Internal modules
  use stat_el_ibm

  ! External modules
  use parallel
  use simulation_flags
  use geometry
  use ibm
  use particle
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
  if (.not. useIBM .and. .not. useParticles) return      
  
  return
end subroutine stat_el_ibm_compute

! =================== !
! Read ibm statistics !
! =================== !
subroutine stat_el_ibm_read

  ! Internal modules
  use stat_el_ibm

  ! External modules
  use parallel
  use ibm

  implicit none

  ! Nothing to do
  
  return
end subroutine stat_el_ibm_read


! ================== !
! Read 1D statistics !
! ================== !
subroutine stat_el_ibm_write

  ! Internal modules
  use stat_el_ibm

  ! External modules
  use time_info
  use parallel
  use fileio
  use math
  use geometry
  use ibm
  use particle
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
  real(WP), dimension(nDimensions) :: objectVelocity, uDoublePrime
  real(WP), dimension(nStatPoints, nDimensions) :: meanVelocity, favreVelocity, meanSigma
  !real(WP), dimension(nStatPoints, 3) :: meanVelocity, favreVelocity
  real(WP), dimension(nStatPoints) :: meanRho, meanP, meanUU, meanVV, meanWW, meanUV, meanUW,&
       meanVW, alpha_ibm, alpha_lpt, volume_ibm, volume_lpt, totalVolume, mean_Viscosity,    &
       numericalDissipation, mean_sound
  logical :: procHasObject
  real(WP), dimension(nGridPoints, nUnknowns) :: dissipationSource

  ! Return if IBM is not used
  if (.not. useIBM .and. .not. useParticles) return
    
  ! Compute mean stats
  meanRho = 0.0_WP
  meanVelocity = 0.0_WP
  favreVelocity = 0.0_WP
  totalVolume = 0.0_WP
  volume_ibm = 0.0_WP
  volume_lpt = 0.0_WP
  mean_sound = 0.0_WP
  mean_Viscosity = 0.0_WP
  meanP = 0.0_WP
  
  do i = 1, nGridPoints
    j = 1 + nint((coordinates(i,1) - xmin) / dx)
     ! Sum up grid volumes
     volume_ibm(j) = volume_ibm(j) + gridNorm(i,1)
     volume_lpt(j) = volume_lpt(j) + primitiveGridNorm(i,1) * (1.0_WP - volumeFraction(i,1))
     totalVolume(j) = totalVolume(j) + primitiveGridNorm(i,1)
     ! Mean stats
     meanRho(j) = meanRho(j) + gridNorm(i,1) * conservedVariables(i,1)                     
     meanVelocity(j,1:nDimensions) = meanVelocity(j,1:nDimensions) +                         &
          gridNorm(i,1) * velocity(i,1:nDimensions)
     favreVelocity(j,1:nDimensions) = favreVelocity(j,1:nDimensions) +                       &
          gridNorm(i,1) * conservedVariables(i,2:nDimensions+1)
     meanP(j) = meanP(j) + gridNorm(i,1) * pressure(i,1)
     meanSigma(j,:) = meanSigma(j,:) + gridNorm(i,1) * stressTensor(i,:)
     ! Mean sound speed
     mean_sound(j) = mean_sound(j) + gridNorm(i,1) * sqrt(ratioOfSpecificHeats *             &
          specificVolume(i,1) * pressure(i,1))
     mean_Viscosity(j)  = mean_Viscosity(j) + gridNorm(i,1) * dynamicViscosity(i,1)
  end do
  
  ! Sum them over procs
  call parallel_sum(volume_ibm)
  call parallel_sum(volume_lpt)
  call parallel_sum(totalVolume)
  call parallel_sum(meanRho)
  call parallel_sum(meanP)
  call parallel_sum(mean_sound)
  call parallel_sum(mean_Viscosity)
  do i = 1, nDimensions
     call parallel_sum(meanVelocity(:,i))
     call parallel_sum(favreVelocity(:,i))
  end do
  do i = 1, nDimensions**2
     call parallel_sum(meanSigma(:,i))
  end do

  ! Normalize
  do i = 1, nStatPoints
     alpha_ibm(i) = 1.0_WP - (volume_ibm(i) / totalVolume(i))            ! Mean volume fraction for IBM
     alpha_lpt(i) = volume_lpt(i) / totalVolume(i)            ! Mean volume fraction for IBM
     meanRho(i) = meanRho(i) / volume_ibm(i)
     meanP(i) = meanP(i) / totalVolume(i)                 ! Changing here for momentum
     meanVelocity(i,:) = meanVelocity(i,:) / volume_ibm(i)
     favreVelocity(i,:) = favreVelocity(i,:) / volume_ibm(i) / meanRho(i)
     meanSigma(i,:) = meanSigma(i,:) / totalvolume(i)
     mean_sound(i) = mean_sound(i) / volume_ibm(i)
     mean_Viscosity(i) = mean_Viscosity(i)/ volume_ibm(i)
  end do

  ! Compute variance
  meanUU = 0.0_WP; meanVV = 0.0_WP; meanWW = 0.0_WP
  meanUV = 0.0_WP; meanUW = 0.0_WP; meanVW = 0.0_WP  
  numericalDissipation = 0.0_WP

  ! Compute higher order statistics
  do i = 1, nGridPoints
     j = 1 + nint((coordinates(i,1) - xmin) / dx)
     uDoublePrime = velocity(i,1:nDimensions) - favreVelocity(j,1:nDimensions)
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
     ! Dissipation
     numericalDissipation(j) = numericalDissipation(j) + gridNorm(i,1) * (sum(uDoublePrime * &
          dissipationSource(i,2:nDimensions+1))  + sum(uDoublePrime**2 * 0.5_WP *           &
          dissipationSource(i,1)))       
  end do
  call parallel_sum(meanUU); meanUU = meanUU / meanRho / volume_ibm
  call parallel_sum(meanVV); meanVV = meanVV / meanRho / volume_ibm
  call parallel_sum(meanWW); meanWW = meanWW / meanRho / volume_ibm
  call parallel_sum(meanUV); meanUV = meanUV / meanRho / volume_ibm
  call parallel_sum(meanUW); meanUW = meanUW / meanRho / volume_ibm
  call parallel_sum(meanVW); meanVW = meanVW / meanRho / volume_ibm
  call parallel_sum(numericalDissipation); numericalDissipation = numericalDissipation /     &
       totalVolume

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
     
     filename= "stat/stat-el_ibm"
     write(filename, '(2A,I8.8,A)') trim(filename), "-", timestep, ".txt"
     iunit = iopen()
     open (iunit, file=adjustl(trim(filename)), form="formatted",iostat=ierror)
     write(iunit,'(A,1ES21.12)') "t=",time
     write(iunit,'(22a21)') "X", "alphaibm", "alphalpt", "Volume", "<U>", "<V>", "<W>",      &
          "<U'U'>", "<V'V'>", "<W'W'>", "<U'V'>", "<V'W'>", "<U'W'>", "<rho>", "<P>",        &
          "<Sigma_1>", "Utilde", "Vtilde", "Wtilde", "Dissipation", "Sound", "Viscosity"
     do i = 1, nStatPoints
        write(iunit,'(10000ES21.12E3)') grid1D(i), alpha_ibm(i), alpha_lpt(i), totalVolume(i), &
             meanVelocity(i,1), meanVelocity(i,2), meanVelocity(i,3), meanUU(i), meanVV(i),  &
             meanWW(i), meanUV(i), meanVW(i), meanUW(i), meanRho(i), meanP(i),               &
             meanSigma(i,1), favreVelocity(i,1), favreVelocity(i,2), favreVelocity(i,3),     &
             numericalDissipation(i), mean_sound(i), mean_Viscosity(i)
     end do
     close(iclose(iunit))
     
  end if

  ! Log
  call monitor_log("IBM-LPT STATISTICS FILE WRITTEN")
  
  return
end subroutine stat_el_ibm_write

