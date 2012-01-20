!Perform a range of line diagnostics
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

menuChoice = 0
submenu = 0
subsubmenu = 0

call display_navigate_menu(menuChoice, submenu, subsubmenu)

if(menuChoice == 1) then
   tempratioID = submenu
   neratioID = submenu
   print *, "n_e guess for temperature sensitive ratio: "
   read(*,*) ne_guess
   print *, "Te guess for temperature sensitive ratio: "
   read(*,*) Te_guess

   call calculate_temperature_and_ne(tempratioID, neratioID, ne_guess, Te_guess)

else if(menuChoice == 2) then
   ratioID = submenu
   print *, "Enter the electron density"
   read(*,*) n_e
   print *, "Enter the line ratio value:"
   read(*,*) lineRatio
   call temp_calc(n_e, lineRatio, ratioID) 
end if

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
     print *, "1. Get Ne and Te from a line ratio "
     print *, "2. I know the electron density and ratio,  want temperature"
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
        print *, "1. ?? "
        print *, "2. ?? "
        print *, "3. ?? "
        print *, "4. ?? "
        print *, "5. Exit "
        print *, " "
        print *, "choice: "
        
        read(*,*) subsubmenu
        
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



subroutine Temp_Calc(n_e, ratio, ratioID)
  implicit none
  real :: ratio, n_e, T, ratio_check, diff, Tmax
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
        print *, "invalid ration ID"
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


subroutine calculate_temperature_and_ne(tempratioID, neratioID, ne_guess, Te_guess)
  implicit none

  integer :: tempratioID, neratioID
  real :: ne_guess, Te_guess
  logical :: globalConverged

  globalConverged = .false.
  
  do while (.not. globalCOnverged) 
     if(tempratioID == 1) then
        ratio_check = 7.9*exp(3.29e4/T)/(1.+(4.5e-4*n_e/sqrt(T)))
     else if(ratioID == 2) then
        ratio_check = 8.23*exp(2.5e4/T)/(1.+(4.4e-3*n_e/sqrt(T)))
     else if(ratioID == 3) then
        ratio_check = 13.7*exp(4.3e4/T)/(1.+(3.8e-5*n_e/sqrt(T)))
     else if(ratioID == 4) then
        ratio_check = 5.44*exp(2.28e4/T)/(1.+(3.5e-4*n_e/sqrt(T)))
     else
        print *, "invalid ration ID"
        call exit_program(0)
     end if

  end do

end subroutine calculate_temperature_and_ne

end program
