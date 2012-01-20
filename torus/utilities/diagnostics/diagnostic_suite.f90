!Perform a range of nebular diagnostics
!Author: Thomas Haworth
!Created: 20/01/2012
!Last Updated: 20/01/2012

program diagnostic_suite

implicit none

integer :: menuChoice
integer :: submenu
integer :: subsubmenu
integer :: ratioID, tempratioID, neratioID
real :: neMax
real :: neMin
real :: n_e
real :: ne_guess
real :: Te_guess
real :: lineRatio
integer :: numFiles
real :: temperature
real :: TempLineRatio
real :: NeLineRatio
real :: neArray(45)
real :: neRatioArray(45)

menuChoice = 0
submenu = 0
subsubmenu = 0

do
   call display_navigate_menu(menuChoice, submenu, subsubmenu)
   
   if(menuChoice == 1) then
      tempratioID = submenu
      neratioID = subsubmenu

      call setup_density_data(neArray, neRatioArray)

      print *, "Enter the temperature sensitive line ratio value:"
      read(*,*) TemplineRatio
      print *, "Enter the electron density sensitive line ratio value:"
      read(*,*) NelineRatio
      call calculate_temperature_and_ne(tempratioID, neratioID, tempLineRatio, &
           NeLineRatio, neArray, neRatioArray)
      
   else if(menuChoice == 2) then
      ratioID = submenu
      print *, "Enter the electron density"
      read(*,*) n_e
      print *, "Enter the line ratio value:"
      read(*,*) lineRatio
      call temp_calc(n_e, lineRatio, ratioID, temperature) 
   end if

   submenu = 0
   subsubmenu = 0
   menuChoice = 0

end do

contains

recursive subroutine display_navigate_menu(choice, submenu, subsubmenu)
  implicit none
  integer :: choice
  integer :: submenu
  integer :: subsubmenu

  if(choice == 0 .and. submenu == 0 .and. subsubmenu == 0) then

     print *, "****************"
     print *, "Diagnostic Suite"
     print *, "****************"
     print *, " "
     print *, "Main Menu - enter a number from the options below"
     print *, " "
     print *, "1. Get Ne and Te from two line ratios "
     print *, "2. Get temperature from a given electron density and ratio"
     print *, "3. Exit "
     print *, " "
     print *, "choice: "
     
     read(*,*) choice
     
     if(choice > 0 .and. choice < 3) then
        call display_navigate_menu(choice, submenu, subsubmenu)
     else if (choice == 3) then
        call exit_program(0)
    else if (choice /= 2) then
        print *, "Please enter a valid choice: ""1"", ""2"" or ""3"""
     end if

  else if(choice == 1 .or. choice == 2) then
     if(submenu == 0) then
        print *, "************************************************"
        print *, "Diagnostic Suite: Temperature Sensitive Ratios"
        print *, "************************************************"
        print *, " "
        print *, "Here you can generate temperature vs electron number density plots "
        print *, "for temperature sensitive line ratios."
        print *, " "
        print *, "Choose a ratio from the options below"
        print *, " "
        print *, "1. [OIII] (5007+4959)/4363 "
        print *, "2. [NII] (6583+6548)/5755 "
        print *, "3. [NeIII] (3968+3869)/3343 "
        print *, "4. [SIII] (9532+9069)/6312 "
        print *, "5. Exit "
        print *, " "
        print *, "choice: "
        
        read(*,*) submenu
        
        if(submenu > 0 .and. submenu < 5) then
           !all if ok
           if(choice == 1) then
              !also need the ne sensitive line choice
              call display_navigate_menu(choice, submenu, subsubmenu)
           end if
        else if(submenu == 5) then
           call exit_program(0)
        else
           print *, "Please enter a valid choice: ""1"",""2"", ""3"", ""4""or ""5"""
        end if
     else if(submenu > 0 .and. submenu < 5) then
        print *, "***************************************************"
        print *, "Diagnostic Suite: Electron Density Sensitive Ratios"
        print *, "***************************************************"
        print *, " "
        print *, "Here you can generate temperature vs electron number density plots "
        print *, "for temperature sensitive line ratios."
        print *, " "
        print *, "Choose a ratio from the options below"
        print *, " "
        print *, "1. [O II] (3729/3726) "
        print *, "2. Exit "
        print *, " "
        print *, "choice: "
        
        read(*,*) subsubmenu
        
        if(subsubmenu == 1) then

        else if(submenu == 2) then
           call exit_program(0)
        else
           print *, "Please enter a valid choice: ""1"" or ""2"""
        end if

     end if

  else
     print *, "Hard coded error in menu setup"
     call exit_program(0)
  end if

end subroutine

subroutine exit_program(code)
  integer :: code

  print *, "Exiting"
  stop
end subroutine


subroutine setup_density_data(neArray, neRatioArray)
  implicit none

  integer, parameter :: n_oii_lines=45
  character(len=34), parameter :: data_path="electron_density_data/wang_OII.dat"
  integer :: ier
  integer :: counter
  real, intent(out) :: neArray(n_oii_lines), neRatioArray(n_oii_lines)
  real :: ne_err_plus(n_oii_lines), ne_err_minus(n_oii_lines)
  real :: neRatio_err_plus(n_oii_lines), neRatio_err_minus(n_oii_lines)

  print *, "!- Opening electron density data file"
  open(1, file=data_path, status="old", iostat=ier)
  if(ier /= 0) then
     print *, "Error opening "//data_path
     call exit_program(0)
  else
     print *, "!- Opened electron density data file successfully"
  end if
 
  do counter = 1, n_oii_lines
     read (1,*) neArray(counter), ne_err_plus(counter), ne_err_minus(counter), &
          neRatioArray(counter), neRatio_err_plus(counter), neRatio_err_minus(counter)
  end do

  close(1)

end subroutine setup_density_data

subroutine Temp_Calc(n_e, ratio, ratioID, T)
  implicit none
  real :: ratio, n_e, ratio_check, diff, Tmax
  real, intent(inout) :: T
  logical :: converged
  integer :: ratioID

  converged = .false.
  T = 1.
  Tmax = 100000.

  do while(.not. converged)
     if(ratioID == 1) then
        ratio_check = 7.9*exp(3.29e4/T)/(1.+(4.5e-4*n_e/sqrt(T)))
     else if(ratioID == 2) then
        ratio_check = 8.23*exp(2.5e4/T)/(1.+(4.4e-3*n_e/sqrt(T)))
     else if(ratioID == 3) then
        ratio_check = 13.7*exp(4.3e4/T)/(1.+(3.8e-5*n_e/sqrt(T)))
     else if(ratioID == 4) then
        ratio_check = 5.44*exp(2.28e4/T)/(1.+(3.5e-4*n_e/sqrt(T)))
     else
        print *, "invalid ratio ID"
        call exit_program(0)
     end if
     
     diff = ratio / ratio_check
     if(diff .ge. 0.999 .and. diff .lt. 1.001) then
        converged = .true.
        goto 100
     end if
     T = T + 0.1
  enddo
100 continue
  print *, "!- Temperature calculation converged."
  print *, "!- Derived temperature is ", T, "K"
end subroutine Temp_Calc



subroutine Ne_Calc(n_e, ratio, ratioID, neArray, neRatioArray)
  implicit none
  real :: ratio, n_e, ratio_check, diff, Ne_max
  logical :: converged
  integer :: ratioID
  real :: neArray(:)
  real :: neRatioArray(:)
  integer :: i
  real :: dratio, factor, dNe, grad
  
  converged = .false.

  if(ratioID == 1) then
     if(ratio < neRatioArray(45)) then
        print *, "electron density tends to zero"
        print *, "neRatioArray(1) ", neRatioArray(1)
        print *, "ratio ", ratio
        call exit_program(0)
     else if(ratio > neRatioArray(1)) then
        print *, "electron density tends to infinity"
        print *, "neRatioArray(45) ", neRatioArray(45)
        print *, "ratio ", ratio
        call exit_program(0)
     else

        do i = 44, 2, -1
           if(ratio >= neRatioArray(i+1) .and. ratio < neRatioArray(i)) then
              dNe = neArray(i) - neArray(i+1)
              dratio = neRatioArray(i) - neRatioArray(i+1)
              grad = dNe/dratio
!              n_e = 10.**(neArray(i+1) +grad*(ratio-neArray(i+1)))
              n_e = (neArray(i+1) +grad*(ratio-neArray(i+1)))
              goto 100
           end if
        end do
     end if
  else
     print *, "invalid ratio ID"
     call exit_program(0)
  end if
 
100 continue
  print *, "!- Derived electron density is ", n_e, "cm^3", i
end subroutine Ne_Calc


subroutine calculate_temperature_and_ne(tempratioID, neratioID, tempLineRatio, &
     NeLineRatio, neArray, neRatioArray)
  implicit none

  integer :: tempratioID, neratioID
  logical :: globalConverged
  real :: currentTemperature
  real :: currentNe
  real :: tempLineRatio
  real :: NeLineRatio
  real :: oldTemperature
  real :: oldNe
  real :: check1
  real :: check2
  real :: neArray(:)
  real :: neRatioArray(:)

  globalConverged = .false.

!  oldTemperature = 0.0
!  oldNe = 0.0
!  currentTemperature = 0.0
!  currentNe = 0.0

!  do while(.not. globalCOnverged) 
     print *, "pre calc step electron density ", currentNe
     call ne_calc(currentNe, neLineRatio, neratioID, neArray, neRatioArray)
     print *, "post calc step electron density ", currentNe

     print *, "pre calc step temperature ", currentTemperature
     call temp_calc(currentNe, templineRatio, tempratioID, currentTemperature) 
     print *, "post calc step temperature ", currentTemperature
     
  !   check1 = (currentTemperature - oldTemperature)/currentTemperature
  !   check2 = (currentNe - oldNe)/currentNe     

 !    oldTemperature = currentTemperature
 !    oldNe = currentNe

 !    if(abs(check1) < 0.001 .and. check2 < abs(0.001)) then
 !       globalConverged = .true.
 !    else
 !       print *, "check1 ", check1
 !       print *, "check2 ", check2
 !    end if
 ! End do

  print *, " "
  print *, " "
  print *, " "
  print *, "************************"
  print *, "        RESULTS         "
  print *, "************************"
  print *, " "
  print *, "---------------------------------------"
  print *, "Final Temperature: ", currentTemperature
  print *, "Final Electron Density: ", currentNe
  print *, "check1 ", check1
  print *, "check2 ", check2
  print *, "---------------------------------------"
  
end subroutine calculate_temperature_and_ne

end program
