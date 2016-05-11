module interface_hierarchy_mod
  use kind_mod
    use vector_mod, only: vector
    use triangle_mesh_mod, only: triangle

    integer, parameter :: backgroundId = 0 ! if no objects are present or tests outside all objects

    type objectIntersection
        integer      :: globalObjectId
        real(double) :: distance
        type(vector) :: normal
    end type objectIntersection


    type localInterfaceHierarchy

        integer                               :: globalObjectId
        type(triangle), pointer, dimension(:) :: triangleList => null() ! if list is empty, test is always true
        type(triangle), pointer, dimension(:) :: normalList => null()

        type(vector) :: reference
        logical      :: referenceIsInside

        type(localInterfaceHierarchy), pointer :: nextNode => null()

    end type localInterfaceHierarchy



    !type localObjectTree
    !
    !end type localObjectDecisionTree



contains

  function returnGlobalIDfromFullMesh(meshList, point) result(id)
    use triangle_mesh_mod, only : triangleMesh, isPointInTriangleMesh
    implicit none
    type(triangleMesh), dimension(:), intent(in) :: meshList
    type(VECTOR) :: point
    integer :: id
    integer :: iMesh
    id = 0
    
    do iMesh = 1, size(meshList)
       if (isPointInTriangleMesh(point, meshList(iMesh))) then
          id = iMesh
          exit
       endif
    enddo
  end function returnGlobalIDfromFullMesh

    function makeLocalInterfaceHierarchy(meshList, boxMin, boxMax, refIn) result(rootNode)

        use vector_mod, only: vector, operator(+), operator(*)
        use triangle_mesh_mod, only: triangleMesh, isAxesAlignedBoxInTriangleMesh, &
                                     isPointInTriangleMesh, intersectAxesAlignedBoxTriangleMesh, &
                                     intersectAxesAlignedBoxTriangle, findTrianglesInAxesAlignedBox

        implicit none
        type(triangleMesh), dimension(:), intent(in) :: meshList
        type(vector), intent(in)                     :: boxMin, boxMax
        type(localInterfaceHierarchy), target        :: rootNode
        type(vector), intent(in), optional           :: refIn

        integer :: iMesh, nMesh
        type(vector) :: reference
        type(localInterfaceHierarchy), pointer :: thisNode

        nMesh = size(meshList)
        !reference = 0.5d0 * (boxMin + boxMax)
        if(present(refIn)) then
            reference = refIn
        else
            reference = randomReference(boxMin,boxMax)
        end if

        thisNode => rootNode

        do iMesh = 1, nMesh

            if(isAxesAlignedBoxInTriangleMesh(boxMin, boxMax, meshList(iMesh))) then ! cell is entirely inside the mesh

                thisNode%globalObjectId = iMesh
                thisNode%reference = reference
                thisNode%referenceIsInside = .true.

                return ! full cell contained in mesh, no need to look for further meshes
            else if(intersectAxesAlignedBoxTriangleMesh(boxMin, boxMax, meshList(iMesh))) then ! cell partially overlaps

                thisNode%globalObjectId = iMesh
                thisNode%reference = reference
                thisNode%referenceIsInside = isPointInTriangleMesh(reference, meshList(iMesh))

                ! find and copy relevant triangles and normals
                call findTrianglesInAxesAlignedBox(boxMin, boxMax, &
                                                   meshList(iMesh)%faces, thisNode%triangleList, &
                                                   meshList(iMesh)%faceNormals, thisNode%normalList)

                allocate(thisNode%nextNode)
                thisNode => thisNode%nextNode
            end if
        end do

        thisNode%globalObjectId = backgroundId
        thisNode%reference = reference
        thisNode%referenceIsInside = .true.

        !thisNode%triangleList => null()
        !thisNode%normalList => null()

        return

    end function makeLocalInterfaceHierarchy


    function refineLocalInterfaceHierarchy(parentRoot, boxMin, boxMax, debug) result(childRoot)

        use vector_mod, only: vector, operator(+), operator(-), operator(*)
        use triangle_mesh_mod, only: findTrianglesInAxesAlignedBox, countIntersectionRayTriangleList

        implicit none
        type(localInterfaceHierarchy), target, intent(in) :: parentRoot
        type(vector), intent(in)                     :: boxMin, boxMax
        logical, optional :: debug
        type(localInterfaceHierarchy), target        :: childRoot
        integer :: countIntersect
        type(vector) :: childReference
        type(localInterfaceHierarchy), pointer :: thisChildNode, thisParentNode

        logical :: writeDebug

        writeDebug = .false.
        if (present(debug)) writeDebug=debug

        thisParentNode => parentRoot
        thisChildNode => childRoot

        !childReference = 0.5d0 * (boxMin + boxMax) + vector(1.d-8, 2.d-8, 2.1d-8)
        childReference = randomReference(boxMin, boxMax)
        if (writedebug) write(*,*) "child reference ",childReference
        do

            if(.not. associated(thisParentNode%nextNode)) then
                thisChildNode%globalObjectId = thisParentNode%globalObjectId
                thisChildNode%reference = childReference
                thisChildNode%referenceIsInside = .true.
                if (writedebug) then
                   write(*,*) "thisParent%nextnode not associated. assigning global object id and returning"
                   write(*,*) "thisChildNode%globalObjectID ",thisChildNode%globalObjectID 
                endif
                return
            endif

            call findTrianglesInAxesAlignedBox(boxMin, boxMax, &
                                               thisParentNode%triangleList, thisChildNode%triangleList, &
                                               thisParentNode%normalList, thisChildNode%normalList)

            countIntersect = countIntersectionRayTriangleList(childReference, thisParentNode%reference - childReference, &
                                                              thisParentNode%triangleList, 1.d0)


            if(size(thisChildNode%triangleList) > 0) then ! childNode shares triangles with parentNode
                thisChildNode%globalObjectId = thisParentNode%globalObjectId

                if (Writedebug)write(*,*) "parent global id is ",thisParentNode%globalObjectId
                if (Writedebug)write(*,*) "thisChildNode%referenceIsInside ",thisChildNode%referenceIsInside
                if (Writedebug)write(*,*) "thisParentNode%referenceIsInside ",thisParentNode%referenceIsInside
!                thisChildNode%globalObjectId = findGlobalObjectID(thisChildNode, childReference) ! TJH

                thisChildNode%reference = childReference
                thisChildNode%referenceIsInside = (modulo(countIntersect,2) == 0 .EQV. thisParentNode%referenceIsInside)

                allocate(thisChildNode%nextNode)
                thisChildNode => thisChildNode%nextNode
                thisParentNode => thisParentNode%nextNode
                if (writedebug) then
                   write(*,*) "child node shares triangles with parent node ",SIZE(thisChildNode%triangleList)
                   write(*,*) "thisChildNode%globalObjectId set to ",thisChildNode%globalObjectId 
                endif
            else ! child does not contain triangles of thisParent
                if(associated(thisChildNode%triangleList)) deallocate(thisChildNode%triangleList)
                if(associated(thisChildNode%normalList)) deallocate(thisChildNode%normalList)

                if(modulo(countIntersect,2) == 0 .EQV. thisParentNode%referenceIsInside) then
                    ! child is entirely inside thisParent, done!
                    thisChildNode%globalObjectId = thisParentNode%globalObjectId
                    thisChildNode%reference = childReference
                    thisChildNode%referenceIsInside = .true.
                    if (writedebug) then
                       write(*,*) "child is entirely inside parent"
                       write(*,*) "thisChildNode%globalObjectId set to ",thisChildNode%globalObjectId 
                       write(*,*) "returning"
                    endif

                    return
                else ! child is entirely outside thisParent, check next parent
                    thisParentNode => thisParentNode%nextNode
                end if

            end if

        end do
        if (writedebug) write(*,*) "returning from end of code"
        return

    end function refineLocalInterfaceHierarchy


    type(vector) function randomReference(boxMin, boxMax)
        use vector_mod, only: vector, operator(+), operator(-), operator(*)
        use random_mod, only: randomNumberGenerator

        type(vector) :: boxMin, boxMax, d
        real(double) :: rx,ry,rz
        real(double), parameter :: radius = 0.6d0

        d = boxMax - boxMin
        randomReference = boxMin + 0.5d0 * d*(1.d0-radius)

        call randomNumberGenerator(getDouble=rx)
        call randomNumberGenerator(getDouble=ry)
        call randomNumberGenerator(getDouble=rz)
        randomReference%x = randomreference%x + radius * d%x * rx
        randomReference%y = randomreference%y + radius * d%y * ry
        randomReference%z = randomreference%z + radius * d%z * rz

        return

    end function randomReference


    subroutine printLocalInterfaceHierarchy(root)
        implicit none
        type(localInterfaceHierarchy), intent(in), target :: root
        type(localInterfaceHierarchy), pointer :: thisNode
        character(len=80) :: text

        thisNode => root
        do while(associated(thisNode))
            if(associated(thisNode%triangleList)) then
                if (thisNode%referenceIsInside) then
                    write (text,'(a,i3,a,i6,a)') 'Object with id ', thisNode%globalObjectId,  ' is tested on ', &
                                size(thisNode%triangleList), ' triangles. (reference is inside)'
                else
                    write (text,'(a,i3,a,i6,a)') 'Object with id ', thisNode%globalObjectId,  ' is tested on ', &
                                size(thisNode%triangleList), ' triangles. (reference is outside)'
                end if
                print *, text
                thisNode => thisNode%nextNode
            else
                if(thisNode%globalObjectId /= backgroundId) then
                    write (text,'(a,i3,a)') 'Object with id ', thisNode%globalObjectId, &
                          ' is tested trivially'
                else
                    text = 'Background'
                end if
                print *, text
                exit
            end if
        end do

    end subroutine printLocalInterfaceHierarchy


    recursive function findGlobalObjectId(thisNode, point) result(globalObjectId)

        use vector_mod, only: operator(-)
        use triangle_mesh_mod, only: countIntersectionRayTriangleList

        implicit none
        type(localInterfaceHierarchy), intent(in) :: thisNode
        type(vector), intent(in) :: point
        integer :: globalObjectId

        integer :: countIntersect

        if(.not. associated(thisNode%triangleList)) then ! empty list; no need to test
            globalObjectId = thisNode%globalObjectId
        else ! triangles are given, test them
            ! Count intersections of line segment point -> reference with given list of triangles
            countIntersect = countIntersectionRayTriangleList(point, thisNode%reference - point, &
                                                              thisNode%triangleList, 1.d0)
            if(modulo(countIntersect,2) == 0 .EQV. thisNode%referenceIsInside) then ! point is inside mesh
                globalObjectId = thisNode%globalObjectId
            else ! point is outside mesh, look at next level in the hierarchy
                if(associated(thisNode%nextNode)) then ! check if this node has a child
                    globalObjectId = findGlobalObjectId(thisNode%nextNode, point)
                else ! 'background' response
                    globalObjectId = backgroundId
                end if
            end if
        end if
        return
    end function findGlobalObjectId


    function combineAllIds(root) result(combined)
        type(localInterfaceHierarchy), intent(in), target :: root
        integer :: combined

        type(localInterfaceHierarchy), pointer :: thisNode

        combined = 0

        thisNode => root

        do while(associated(thisNode))
            combined = combined + 2**thisNode%globalObjectId
            thisNode => thisNode%nextNode
        end do
    end function



    function findFirstIntersection(root, point, direction, currentId) result(intersect)
        use vector_mod, only: operator(+), operator(*)
        use triangle_mesh_mod, only: findFirstIntersectionTriangleList
        implicit none
        type(localInterfaceHierarchy), intent(in), target :: root
        type(vector), intent(in) :: point, direction
        integer, intent(in) :: currentId ! object Id in which the given point lies
        type(objectIntersection) :: intersect

        type(localInterfaceHierarchy), pointer :: thisNode
        type(triangle) :: tempT
        integer :: thisIndex
        real(double) :: thisDist, minDist
        type(vector) :: thisCoords

        intersect%globalObjectId = -1
        intersect%distance = -1.d0
        intersect%normal = vector(0,0,0)

        minDist = -1

        thisNode => root


        !do while(associated(thisNode) .and. associated(thisNode%triangleList))
        do while(associated(thisNode%nextNode))
            call findFirstIntersectionTriangleList(point, direction, thisNode%triangleList, &
                                                   thisIndex, thisDist, thisCoords)
            if(thisIndex > 0) then ! this test hit something
                if(minDist == -1 .or. thisDist < minDist) then ! new hit is closer than previous
                    minDist = thisDist
                    if(thisNode%globalObjectId == currentId) then ! exits current object to child, use that ID
                        !if(associated(thisNode%nextNode)) then ! next object does exist
                        !    intersect%globalObjectId = thisNode%nextNode%globalObjectId
                        !else ! next object does not exist
                        !    intersect%globalObjectId = 0
                        !end if
                        intersect%globalObjectId = thisNode%nextNode%globalObjectId
                    else ! exits current into (grand) parent
                        intersect%globalObjectId = thisNode%globalObjectId
                    end if
                    intersect%distance = thisDist
                    tempT = thisNode%normalList(thisIndex)
                    intersect%normal = thisCoords%x*tempT%a + &
                                       thisCoords%y*tempT%b + &
                                       thisCoords%z*tempT%c
                end if
            end if

            if(thisNode%globalObjectId == currentId) return ! look further than current object

            thisNode => thisNode%nextNode

        end do

    end function findFirstIntersection


    function countInterfaceTriangles(root) result(n)
        ! Returns the total number of triangles in this tree
        implicit none
        type(localInterfaceHierarchy), intent(in), target :: root
        integer :: n

        type(localInterfaceHierarchy), pointer :: thisNode

        n = 0
        thisNode => root

        do while( associated(thisNode) .and. associated(thisNode%triangleList) )
            n = n + size(thisNode%triangleList)
            thisNode => thisNode%nextNode
        end do

    end function countInterfaceTriangles

    subroutine freeTree(root)
        ! Since the tree is a linked-list type structure, its nodes need to be deleted one by one
        ! This function deallocates all children, but not the root itself
        type(localInterfaceHierarchy), intent(inout), target :: root

        type(localInterfaceHierarchy), pointer :: this, next
        this => root%nextNode

        do while(associated(this))
            next => this%nextNode
            deallocate(this)
            this => next
        end do
    end subroutine freeTree


end module interface_hierarchy_mod
