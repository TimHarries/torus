module linked_list_class

  !
  ! Class definitions for the small linked list (only in one
  ! direction. 
  ! For now, this stores only one integer value at a node.  Can be
  ! generalaized later if you wish.
  !
  ! Ignores the head value.
  !
  ! WARNING :: Function get (get_value_from_list_seq) should not be used
  !            in a recursive routine!
  !
  !------------------------------------------------------------------------
  !  LOGS:    Created on Jan-20-2003  (R. Kurosawa)
  !
  !
  implicit none

  public :: new,        &  ! Initialize a list
       &    put,        &  ! Put a value in a list
       &    get,        &  ! Get a value from a list
       &    delete,     &  ! Delete a value from a list
       &    delete_all, &  ! Delete all the values from a list
       &    get_n_node     ! Get the number of elements in the list


  private:: int_link_node, int_linked_list, put_value, put_node, put_value_in_list, &
       & get_tail, get_value, get_tail_node, get_value_from_list, get_value_from_list_seq

  
  type linked_list
     private
     integer                  :: n_node  ! number of nodes in the list
     type(link_node), pointer :: head
     type(link_node), pointer :: tail
     type(link_node), pointer :: current ! used for a fast access
  end type linked_list

  
  type link_node
     private
     integer                  :: value      ! data
     type(link_node), pointer :: next       ! next
  end type link_node


  !
  !
  interface new
     module procedure int_link_node
     module procedure int_linked_list
  end interface
  
  interface put
     module procedure put_value_in_list
     module procedure put_value
     module procedure put_node
  end interface

  interface get
     module procedure get_value_from_list
     module procedure get_value_from_list_seq
  end interface
  
contains

  !========================================================
  ! Constructors
  !=======================================================

  !
  !
  subroutine int_link_node(node)
    implicit none
    type(link_node), intent(inout) :: node
    
    node%value = 0
    NULLIFY(node%next)
  end subroutine int_link_node

  !
  !
  subroutine int_linked_list(list)
    implicit none
    type(linked_list), intent(inout) :: list

    list%n_node = 0
    NULLIFY(list%head)
    NULLIFY(list%tail)
    NULLIFY(list%current)
    ALLOCATE(list%head)
    ALLOCATE(list%tail)
    call int_link_node(list%head)
    call int_link_node(list%tail)

    list%tail => list%head ! initially
    
  end subroutine int_linked_list
  


  
  !==============================================================
  ! Basic method for accessing the data
  !===============================================================
  
  !
  ! Assigns a value to a given node
  subroutine put_value(node, value)
    implicit none
    type(link_node), intent(inout) :: node
    integer, intent(in) :: value
    !
    node%value = value
  end subroutine put_value
    

  !
  ! Get data from a node
  function get_value(node) RESULT(out)
    implicit none
    integer :: out
    type(link_node), intent(in) :: node
    !  
    out = node%value
  end function get_value
  
  !
  ! Add one node to the tail of the list
  subroutine put_node(list, node)
    implicit none
    type(linked_list), intent(inout) :: list
    type(link_node), intent(in):: node
    !
    type(link_node), pointer :: pnode ! tmp pointer

    NULLIFY(pnode)
    ALLOCATE(pnode)

    pnode => list%tail

    NULLIFY(pnode%next)
    ALLOCATE(pnode%next)
    pnode%next = node

    list%tail => pnode%next
    
  end subroutine put_node



  !
  ! function to return the tail node in the list

  function get_tail_node(list) RESULT(out)
    implicit none
    type(link_node),pointer  :: out 
    type(linked_list), intent(in) :: list

    !out = get_tail(list%head)
    out => list%tail

  end function get_tail_node

  !
  !
  !
  recursive function get_tail(in_node) RESULT(out)
    implicit none
    type(link_node), pointer :: out
    type(link_node), pointer :: in_node
    
    if (ASSOCIATED(in_node%next)) then  ! goto next node
       !out => get_tail(in_node%next)
       out => get_tail(in_node%next)
    else
       ! this is the tail
       out => in_node
    end if
  end function get_tail



  !
  ! Add a value in the list (at the last position by default)
  !
  subroutine put_value_in_list(list, value)
    implicit none
    type(linked_list), intent(inout) :: list
    integer, intent(in) :: value
    !
    type(link_node)  :: new_node
    
    ! create a new node
    call int_link_node(new_node)

    ! put the value to this node
    call put_value(new_node, value)

    ! put the node in the tail of the list
    call put_node(list, new_node)

    ! counter
    list%n_node = list%n_node + 1
    
  end subroutine put_value_in_list

  
  !
  ! Get the value in the i-th node/position
  !
  function get_value_from_list(list, i) RESULT(out)
    implicit none
    integer :: out
    integer, intent(in) :: i   ! index 
    type(linked_list), intent(in) :: list
    
    type(link_node), pointer :: pnode ! tmp pointer
    
    integer :: j, n
    
    NULLIFY(pnode)
    pnode => list%head
    n = list%n_node
    
    if (i > n) then 
       write(*,*) 'Error:: i > n  in get_value_from_list.'
       write(*,*) '   i = ', i
       write(*,*) '   n = ', n
       stop
    end if


    j = 0
    do while ( i /= j)
       j = j+1
       pnode => pnode%next
    end do
    
    out = pnode%value        
    
  end function get_value_from_list


  !
  ! Get a value in a list sequentially
  ! This is much faster than "get_value_from_list" routine above,
  ! if you are accessing the data sequentially many time.
  function get_value_from_list_seq(list, restart) RESULT(out)
    implicit none
    integer :: out 
    type(linked_list), intent(inout) :: list
    logical, intent(in) :: restart
    !
    !
    !    type(link_node), save, pointer :: pnode

    if (restart) then
       list%current => list%head
    end if
    
    list%current => list%current%next

    out = list%current%value
    
  end function get_value_from_list_seq
  
   
    

  !
  ! print the list
  !
  subroutine print_list(list)
    implicit none
    type(linked_list), intent(in) :: list    

    if (ASSOCIATED(list%head%next)) then
       call print_node_value(list%head%next)
    else
       write(*,*) 'Error:: No element presents in the linked list (print_list)'
       stop
    end if

  contains
    recursive subroutine print_node_value(in_node)
      implicit none
      type(link_node), pointer :: in_node

      if (associated(in_node%next)) then
         write(*,*) in_node%value
         call print_node_value(in_node%next)
      else
         write(*,*) in_node%value
      end if
      
    end subroutine print_node_value
  end subroutine print_list

  
  !
  ! get the number of nodes in the list
  !
  function get_n_node(list) RESULT(out)
    implicit none   
    integer :: out
    type(linked_list), intent(in) :: list
    out = list%n_node
  end function get_n_node



  !
  !
  ! Delete i-th node from a list.
  !
  ! NB: If slow, you should modefy this routine
  ! to access data directly without doing dumb
  ! sequantial index search.
  !
  subroutine delete(list, i)
    implicit none
    type(linked_list), intent(inout) :: list
    integer, intent(in) :: i
    !
    integer :: j, n
    type(link_node), pointer :: pnode      ! tmp pointer
    type(link_node), pointer :: pnode_next ! tmp pointer
    type(link_node), pointer :: pnode_prev ! tmp pointer

    NULLIFY(pnode)
    pnode => list%head
    n = list%n_node

    j =0
    do while (i /= j)
       j = j+1
       if (j> n) then ! something went wrong
          write(*,*) 'Error:: j > n  in delete.'
          write(*,*) 'j=', j
          write(*,*) 'n=', n
          stop
       end if
       pnode_prev => pnode   ! will be used later
       pnode => pnode%next
    end do
    
    pnode_next => pnode%next

    ! Now we have
    ! pnode_prev -->  pnode --> pnode_next
    !
    ! and want to point pre_prev%next to pnode_next, and
    ! delete the memory for pnode.
    ! 

    pnode_prev%next => pnode_next
    
    DEALLOCATE(pnode)
    NULLIFY(pnode)

    ! decrease the total number by one
    list%n_node = list%n_node - 1
    

  end subroutine delete
   


  !
  ! deletes all the nodes in the list
  ! 

  subroutine delete_all(list)
    implicit none
    type(linked_list), intent(inout) :: list
    !

    ! use delete function defined above.
    ! This maybe extremely slow.
    
    integer :: i, n
  
    n = list%n_node

    i =0
    do while (i < n)
       i = i+1
       call delete(list, 1)
    end do

    ! Keep the head for later use!

    ! clean up
    if (ASSOCIATED(list%head)) DEALLOCATE(list%head)
    if (ASSOCIATED(list%tail)) DEALLOCATE(list%tail)

    ! prepare for the next use.
    
    call int_linked_list(list)
    
  end subroutine delete_all
    
    
    
  
end module linked_list_class
  
