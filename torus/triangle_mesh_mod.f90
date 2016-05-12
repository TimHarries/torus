module triangle_mesh_mod
    use vector_mod, only: vector
    use kind_mod

    type triangle
        ! triangle definition with pointers to its vertices
        type(vector), pointer :: a
        type(vector), pointer :: b
        type(vector), pointer :: c
    end type triangle

    type triangleMesh
        ! 3d object consisting of a collection of points (nodes) in space
        ! and a set of triangles (faces) between them.
        character(len=80) :: label
        character(len=80) :: material

        integer :: nNode
        integer :: nFace
        integer :: nNormal

        logical :: isClosed

        type(vector), pointer   :: nodes(:)
        type(vector), pointer   :: normals(:)
        type(triangle), pointer :: faces(:)
        type(triangle), pointer :: faceNormals(:)

        type(vector) :: boundingBoxMin, boundingBoxMax ! values for axis-aligned bounding box of mesh
    end type triangleMesh


contains



    subroutine printVector(v)
        implicit none
        type(vector) :: v

        print *, v%x, v%y, v%z
    end subroutine printVector


    function readTriangleMesh(fileName) result(meshList)
        use obj_io_mod, only: obj_size, obj_read
        use vector_mod, only: vector, operator(+), operator(-), modulus, &
                              xHat, yHat, zHat
        use messages_mod, only: writeFatal, writeWarning

        implicit none

        character ( len = * )  :: fileName
        type(triangleMesh), pointer, dimension(:) :: meshList
        type(triangleMesh), pointer :: thisMesh

        integer :: iObject, nObject
        character(len=255), dimension(:), pointer :: objectNames
        character(len=255), dimension(:), pointer :: objectMaterials
        integer, dimension(:), pointer :: objStartNode, objStartFace, objStartNormal

        integer :: nNode, nFace, nNormal, maxOrder

        ! arrays for intermediate results
        real(kind = 8), pointer :: rawNodes(:,:)
        real(kind = 8), pointer :: rawNormals(:,:)
        integer(kind=4), pointer :: rawFaces(:,:)
        integer(kind=4), pointer :: rawFaceNormals(:,:)
        integer(kind=4), pointer :: face_order(:)

        type(vector), pointer :: allNodes(:), allNormals(:)
        type(triangle), pointer :: allFaces(:), allFaceNormals(:)

        character (len=255) :: tempBuffer


        type(vector), pointer :: edgeSums(:) ! use for checking if mesh is closed

        integer :: i,j ! looping variable
        integer :: node1, node2
        integer :: nBoundaryNodes

        real(double) :: xMin, xMax, yMin, yMax, zMin, zMax

        call obj_size ( fileName, nNode, nFace, nNormal, maxOrder )

        if(.not. maxOrder == 3) then
            call writeFatal('Given mesh does not consist of triangles exclusively.')
            return
        end if

        ! allocate temporary buffers
        allocate(rawNodes(3, nNode), &
                 face_order(nFace), &
                 rawFaces(3, nFace), &
                 rawNormals(3, nNormal), &
                 rawFaceNormals(3, nFace))

        ! read into temporary buffers
        call obj_read ( fileName, nNode, nFace, nNormal, 3, &
                        rawNodes, face_order, rawFaces, rawNormals, rawFaceNormals, &
                        objectNames, objectMaterials, objStartNode, objStartFace, objStartNormal)

        ! The following 3 loops would be unnecessary if the 3 x N real-arrays could be typecast into type(vector) arrays of size N

        ! allocate final buffers
        allocate(allNodes(nNode), &
                 allNormals(nNormal), &
                 allFaces(nFace), &
                 allFaceNormals(nFace))

        ! copy rawNodes into allNodes
        do i= 1, nNode
            allNodes(i) = vector(rawNodes(1,i), rawNodes(2,i), rawNodes(3,i))
        end do

        ! copy rawNormals
        do i = 1, nNormal
            allNormals(i) = vector(rawNormals(1,i), rawNormals(2,i), rawNormals(3,i))
        end do

        ! construct triangles as pointers to nodes and pointers to normals
        do i = 1, nFace
            allFaces(i)%a => allNodes(rawFaces(1,i))
            allFaces(i)%b => allNodes(rawFaces(2,i))
            allFaces(i)%c => allNodes(rawFaces(3,i))

            allFaceNormals(i)%a => allNormals(rawFaceNormals(1,i))
            allFaceNormals(i)%b => allNormals(rawFaceNormals(2,i))
            allFaceNormals(i)%c => allNormals(rawFaceNormals(3,i))
        end do

        deallocate(rawNodes, rawNormals, rawFaceNormals)

        nObject = size(objectNames)
        allocate(meshList(nObject))

        do iObject = 1, nObject
            thisMesh => meshList(iObject)

            ! assign object label
            thisMesh%label = trim(objectNames(iObject))

            ! assign object material
            thisMesh%material = trim(objectMaterials(iObject))

            ! using starting index of each object to find out the appropriate size
            if(iObject < nObject) then
                thisMesh%nNode = objStartNode(iObject+1) - objStartNode(iObject)
                thisMesh%nFace = objStartFace(iObject+1) - objStartFace(iObject)
                thisMesh%nNormal = objStartNormal(iObject+1) - objStartNormal(iObject)
            else
                thisMesh%nNode = nNode + 1 - objStartNode(iObject)
                thisMesh%nFace = nFace + 1 - objStartFace(iObject)
                thisMesh%nNormal = nNormal + 1 - objStartNormal(iObject)
            end if

            ! set pointers to relevant position in global lists
            thisMesh%nodes => allNodes( objStartNode(iObject) : objStartNode(iObject)+thisMesh%nNode-1 )
            thisMesh%normals => allNormals( objStartNormal(iObject) : objStartNormal(iObject)+thisMesh%nNormal-1 )
            thisMesh%faces => allFaces( objStartFace(iObject) : objStartFace(iObject)+thisMesh%nFace-1 )
            thisMesh%faceNormals => allFaceNormals( objStartFace(iObject) : objStartFace(iObject)+thisMesh%nFace-1 )

            ! Project on axis aligned unit vectors to obtain the bounding box
            call projectVerticesMinMax(thisMesh%nodes, xHat, xMin, xMax)
            call projectVerticesMinMax(thisMesh%nodes, yHat, yMin, yMax)
            call projectVerticesMinMax(thisMesh%nodes, zHat, zMin, zMax)
            thisMesh%boundingBoxMin = vector(xMin, yMin, zMin)
            thisMesh%boundingBoxMax = vector(xMax, yMax, zMax)

            ! Check if mesh is closed, by checking that no vertex lies on a boundary
            ! At internal points, the sum of incoming+outgoing edges vanishes.
            allocate(edgeSums(thisMesh%nNode))
            edgeSums = vector(0,0,0)

            do i = 1,thisMesh%nFace
                do j = 1,3
                    node1 = rawFaces(j,objStartFace(iObject)+i-1) - objStartNode(iObject) + 1
                    node2 = rawFaces(1+modulo(j,3), objStartFace(iObject)+i-1) - objStartNode(iObject) + 1
                    edgeSums(node1) = edgeSums(node1) + thisMesh%nodes(node1) - thisMesh%nodes(node2)
                    edgeSums(node2) = edgeSums(node2) + thisMesh%nodes(node1) - thisMesh%nodes(node2)
                end do
            end do

            nBoundaryNodes = 0
            do i = 1,thisMesh%nNode
                if(modulus(edgeSums(i)) > 1d-12) then
                    nBoundaryNodes = nBoundaryNodes + 1
                end if
            end do

            if(nBoundaryNodes == 0) then
                thisMesh%isClosed = .true.
            else
                thisMesh%isClosed = .false.
            end if

            if(.not. thisMesh%isClosed) then
                write (tempBuffer, '(a,a,a,a,a)') "Object labeled '", thisMesh%label, "' in ", &
                    fileName, " is not closed. This could lead to unexpected results"
                call writeWarning(tempBuffer)
            end if

            deallocate(edgeSums)

        end do

        deallocate(rawFaces)

        return

    end function readTriangleMesh


    subroutine printTriangleMeshSummary(mesh, level)
        ! Prints a summary of the triangle mesh, at given verbosity level
        use messages_mod
        use vector_mod, only: reportVector
        implicit none
        type(triangleMesh) :: mesh
        integer, optional :: level
        integer :: l

        character(len=80) :: temp

        if(present(level)) then
            l = level
        else
            l = 0
        end if

        call writeInfo("--- TRIANGLE MESH SUMMARY ---", l)
        write (temp, '(a,a)') "label: ", trim(mesh%label)
        call writeInfo(temp, l)

        if(mesh%isClosed) then
            call writeInfo("topology: closed", l)
        else
            call writeInfo("topology: open", l)
        end if
        call writeFormatted("(a,1I8)", "no. nodes:", mesh%nNode, l)
        call writeFormatted("(a,1I6)", "no. normals:", mesh%nNormal, l)
        call writeFormatted("(a,1I8)", "no. faces:", mesh%nFace, l)
        call writeInfo("contained in axes aligned bounding box:", l)
        write (temp, '(a,3F10.5)') "  min: ", mesh%boundingBoxMin%x, mesh%boundingBoxMin%y, mesh%boundingBoxMin%z
        call writeInfo(temp, l)
        write (temp, '(a,3F10.5)') "  max: ", mesh%boundingBoxMax%x, mesh%boundingBoxMax%y, mesh%boundingBoxMax%z
        call writeInfo(temp, l)

    end subroutine printTriangleMeshSummary

    subroutine printMesh(mesh)
        ! prints the mesh points and coordinates of all triangles
        ! (triangles are stored as pointers to the vertices, so indices cannot be obtained)
        implicit none
        type(triangleMesh) :: mesh
        integer :: i

        print *, "nodes:"
        print *, "        idx            x                y                z"
        do i = 1, mesh%nNode
            print *, i, mesh%nodes(i)%x, mesh%nodes(i)%y, mesh%nodes(i)%z
        end do

        print *, ""
        print *, "faces:"
        do i = 1, mesh%nFace
            print *, 'face', i
            call printVector(mesh%faces(i)%a)
            call printVector(mesh%faces(i)%b)
            call printVector(mesh%faces(i)%c)
        end do

        print *, ""
        call printTriangleMeshSummary(mesh)

    end subroutine printMesh

    logical function isPointInAxesAlignedBox(point, boxMin, boxMax)
        ! Tests whether 'point' is inside axes-aligned box given by boxMin and boxMax
        use vector_mod, only: vector

        implicit none
        type(vector), intent(in) :: point, boxMin, boxMax

        isPointInAxesAlignedBox = .false.

        if(point%x < boxMin%x) then
            return
        else if(point%x > boxMax%x) then
            return
        else if(point%y < boxMin%y) then
            return
        else if(point%y > boxMax%y) then
            return
        else if(point%z < boxMin%z) then
            return
        else if(point%z > boxMax%z) then
            return
        else
            isPointInAxesAlignedBox = .true.
            return
        end if
    end function isPointInAxesAlignedBox

    real(double) function averageZ(mesh)
        type(triangleMesh) :: mesh
        integer :: i
        averageZ = 0.d0
        do i = 1, SIZE(mesh%faces)
           averageZ = averageZ + mesh%faces(i)%a%z
           averageZ = averageZ + mesh%faces(i)%b%z
           averageZ = averageZ + mesh%faces(i)%c%z
        enddo
        averageZ = averageZ / dble(3*size(mesh%faces))
      end function averageZ


    logical function intersectAxesAlignedBoxTriangle(boxMin, boxMax, &
                                                     triangleA, triangleB, triangleC) &
                                                     result(intersect)
        ! Intersection between Axis-Aligned Box and triangle
        ! The box is given by its lower and upper corners boxMin and boxMax
        ! The triangle is given by the points A, B, C
        ! Returns true if the two objects intersect, false otherwise
        !
        ! Based on: http://stackoverflow.com/questions/17458562/efficient-aabb-triangle-intersection-in-c-sharp
        !
        use vector_mod, only: vector, dotProd, crossProd, operator(-), &
                              xHat, yHat, zHat ! unit vectors

        implicit none
        type(vector) :: boxMin, boxMax
        type(vector) :: triangleA, triangleB, triangleC

        type(vector), dimension(8) :: vertexBox
        type(vector), dimension(3) :: normalBox
        type(vector), dimension(3) :: vertexTriangle
        type(vector), dimension(3) :: edgesTriangle
        type(vector) :: normalTriangle
        type(vector) :: projectAxis

        real(double) :: projectMinBox, projectMaxBox
        real(double) :: projectTriangle, projectMinTriangle, projectMaxTriangle
        integer :: i,j


        intersect = .false.

        ! check axis aligned bounding box of triangle with box
        if(max(triangleA%x, triangleB%x, triangleC%x) < boxMin%x) then
            return
        else if(min(triangleA%x, triangleB%x, triangleC%x) > boxMax%x) then
            return
        else if(max(triangleA%y, triangleB%y, triangleC%y) < boxMin%y) then
            return
        else if(min(triangleA%y, triangleB%y, triangleC%y) > boxMax%y) then
            return
        else if(max(triangleA%z, triangleB%z, triangleC%z) < boxMin%z) then
            return
        else if(min(triangleA%z, triangleB%z, triangleC%z) > boxMax%z) then
            return
        end if

        vertexBox = verticesFromAxesAlignedBox(boxMin, boxMax)
        normalTriangle = crossProd(triangleB-triangleA, triangleC-triangleA) ! compute a normal vector of the triangle

        ! project everything on the triangle's normal vector
        projectTriangle = dotProd(normalTriangle, triangleA)
        call projectVerticesMinMax(vertexBox, normalTriangle, projectMinBox, projectMaxBox)
        if(projectTriangle < projectMinBox .or. projectMaxBox < projectTriangle) then
            return
        end if

        normalBox = (/ xHat, yHat, zHat /)
        vertexTriangle = (/ triangleA, triangleB, triangleC /)
        edgesTriangle = (/ triangleB-triangleA, triangleC-triangleB, triangleA-triangleC /)

        ! all normals resulting from combinations between box' edges and triangle's edges
        do i = 1,3
            do j = 1,3
                projectAxis = crossProd(normalBox(i), edgesTriangle(j))
                call projectVerticesMinMax(vertexBox, projectAxis, projectMinBox, projectMaxBox)
                call projectVerticesMinMax(vertexTriangle, projectAxis, projectMinTriangle, projectMaxTriangle)
                if (projectMaxTriangle < projectMinBox .or. projectMaxBox < projectMinTriangle) then
                    return
                end if
            end do
        end do

        ! if all foregoing tests where negative, they must intersect
        intersect = .true.
        return

    end function intersectAxesAlignedBoxTriangle

    function verticesFromAxesAlignedBox(boxMin, boxMax) result(vertices)
        ! returns an array of 8 type(vectors) containing the vertices
        ! of the axis aligned bouding box defined by boxMin and boxMax
        use vector_mod, only: vector
        implicit none
        type(vector), intent(in) :: boxMin, boxMax
        type(vector), dimension(8) :: vertices

        vertices = (/ vector(boxMin%x, boxMin%y, boxMin%z), & ! make array with vertices of box
                      vector(boxMin%x, boxMin%y, boxMax%z), &
                      vector(boxMin%x, boxMax%y, boxMin%z), &
                      vector(boxMin%x, boxMax%y, boxMax%z), &
                      vector(boxMax%x, boxMin%y, boxMin%z), &
                      vector(boxMax%x, boxMin%y, boxMax%z), &
                      vector(boxMax%x, boxMax%y, boxMin%z), &
                      vector(boxMax%x, boxMax%y, boxMax%z) /)
        return
    end function verticesFromAxesAlignedBox


    subroutine projectVerticesMinMax(vertices, axis, minProject, maxProject)
        ! computes the projection (dot product) of each vertex in vertices with given axis
        ! returns the minimal and maximal value
        use vector_mod, only: vector, dotProd
        type(vector), dimension(:), intent(in) :: vertices
        type(vector), intent(in) :: axis
        real(double) :: minProject, maxProject

        integer :: iVertex, nVertex
        real(double) :: p

        nVertex = size(vertices, 1)

        minProject = dotProd(vertices(1), axis)
        maxProject = minProject

        do iVertex = 2, nVertex
            p = dotProd(vertices(iVertex), axis)
            minProject = min(p, minProject)
            maxProject = max(P, maxProject)
        end do
    end subroutine projectVerticesMinMax


    logical function intersectAxesAlignedBoxTriangleMesh(boxMin, boxMax, mesh) result(intersect)
        ! Tests whether given AABox intersects with the given triangle mesh

        implicit none
        type(vector), intent(in) :: boxMin, boxMax
        type(triangleMesh), intent(in) :: mesh

        integer :: i

        ! Test all triangles
        intersect = .false.
        do i= 1, mesh%nFace
            !if(isPointInAxesAlignedBox(mesh%nodes(iNode), boxMin, boxMax)) then
            if(intersectAxesAlignedBoxTriangle(boxMin, boxMax, mesh%faces(i)%a, mesh%faces(i)%b, mesh%faces(i)%c)) then
                intersect = .true.
                return
            end if
        end do

        return
    end function intersectAxesAlignedBoxTriangleMesh

    integer function countIntersectionRayTriangleList(point, direction, list, limit) result(intersectCount)
        ! Shoot a ray in given direction and count the number of triangles it intersects
        ! When enabled, optional argument limit>0, does not count intersections beyond limit*direction
        use vector_mod, only: vector, intersectLineTriangle

        implicit none
        type(vector), intent(in)                 :: point, direction
        type(triangle), intent(in), dimension(:) :: list
        real(kind=8), intent(in), optional       :: limit

        logical hasLimit
        real(kind=8) :: limitVal, distance

        integer :: iTriangle, nTriangle

        nTriangle = size(list)
        if(nTriangle == 0) then
            intersectCount = 0
            return
        end if

        if(present(limit)) then
            hasLimit = .true.
            limitVal = limit
        else
            hasLimit = .false.
            limitVal = 0
        end if


        intersectCount = 0
        do iTriangle = 1, nTriangle
            if( intersectLineTriangle(point, direction, .true., &
                                      list(iTriangle)%a, list(iTriangle)%b, list(iTriangle)%c, &
                                      distance)) then ! check for intersection with triangle
                if(.not. hasLimit .or. distance < limitVal) then ! if limited, check if distance is not beyond limit
                    intersectCount = intersectCount + 1
                end if
            end if
        end do

        return
    end function countIntersectionRayTriangleList

    integer function countIntersectionRayTriangleMesh(point, direction, mesh, limit) result(intersectCount)
        ! Shoot a ray in given direction and, for given mesh, count the number of triangles it intersects
        ! Wrapper for countIntersectionRayTriangleList

        implicit none
        type(vector), intent(in) :: point, direction
        type(triangleMesh), intent(in) :: mesh
        real(kind=8), intent(in), optional :: limit

        if(present(limit)) then
            intersectCount = countIntersectionRayTriangleList(point, direction, mesh%faces, limit)
        else
            intersectCount = countIntersectionRayTriangleList(point, direction, mesh%faces)
        end if
        return
    end function countIntersectionRayTriangleMesh



    logical function isPointInTriangleMesh(point, mesh, probe)
        ! Tests whether given point is 'inside' given triangle mesh
        ! If mesh is closed, returns true is point is in the interior of mesh
        ! If mesh is open, returns .... FIXME

        use vector_mod, only: vector, operator(-), zHat, dotProd

        implicit none
        type(vector), intent(in) :: point
        type(triangleMesh), intent(in) :: mesh
        type(vector), intent(in), optional :: probe

        type(vector) :: axis
        real(double) :: pointProject, boxProjectMin, boxProjectMax

        integer :: intersectCount

        if(present(probe)) then
            axis = probe - point
        else
            axis = zHat
        end if

        ! Use bouding box to obtain a quick answer
        if(mesh%isClosed) then ! mesh is closed
            if(.not. isPointInAxesAlignedBox(point, mesh%boundingBoxMin, mesh%boundingBoxMax) ) then ! check AA-bounding box
                isPointInTriangleMesh = .false.
                return
            end if
        else ! mesh is open
            ! see if the axis separates the point and the bounding box
            pointProject = dotProd(point, axis)
            call projectVerticesMinMax( verticesFromAxesAlignedBox(mesh%boundingBoxMin, mesh%boundingBoxMax), &
                                        axis, boxProjectMin, boxProjectMax )

            if(pointProject < boxProjectMin) then
                isPointInTriangleMesh = .true.
                return
            else if(boxProjectMax < pointProject) then
                isPointInTriangleMesh = .false.
                return
            end if
        end if
        ! Bounding box is not conclusive.

        ! Shoot ray and count intersection with triangle
        intersectCount = countIntersectionRayTriangleMesh(point, axis, mesh)

        ! Even number implies point is outside mesh, odd implies inside
        if(modulo(intersectCount, 2) == 0) then
            isPointInTriangleMesh = .false.
        else
            isPointInTriangleMesh = .true.
        end if

        return
    end function isPointInTriangleMesh

    logical function isAxesAlignedBoxInTriangleMesh(boxMin, boxMax, mesh)
        ! Tests whether given AABox is entirely inside given triangle mesh
        use vector_mod, only: operator(+), operator(*)

        implicit none
        type(vector),       intent(in) :: boxMin, boxMax
        type(triangleMesh), intent(in) :: mesh

        type(vector) :: centre
        centre = 0.5d0*(boxMin + boxMax)

        ! box is entirely in mesh if one point is inside the mesh and
        ! none of them mesh's faces intersect the box
        isAxesAlignedBoxInTriangleMesh = isPointInTriangleMesh(centre, mesh) .and. &
                                         .not. intersectAxesAlignedBoxTriangleMesh(boxMin, boxMax, mesh)
        return
    end function isAxesAlignedBoxInTriangleMesh


    subroutine findFirstIntersectionTriangleList(point, direction, list, &         ! input
                                                 minIndex, minDistance, minCoords) ! output
        ! Given ray and list of triangles, finds the index of the first triangle hit
        ! by the ray and the distance to this point.
        ! If no intersection exists, minIndex = 0 and minDistance = 0
        ! Optionally sets minCoords: the homogeneous coordinates of the relevant intersection
        use vector_mod, only: vector, intersectLineTriangle

        implicit none
        type(vector),   intent(in)               :: point, direction
        type(triangle), intent(in), dimension(:) :: list
        integer,        intent(out)              :: minIndex
        real(kind=8),   intent(out)              :: minDistance
        type(vector),   intent(out), optional    :: minCoords

        real(kind=8) :: distance
        type(vector) :: minC, coords

        integer :: iTriangle, nTriangle


        minIndex = 0
        minDistance = -1 ! silly placeholder


        nTriangle = size(list, 1)
        do iTriangle = 1, nTriangle
            if( intersectLineTriangle(point, direction, .true., &
                                      list(iTriangle)%a, list(iTriangle)%b, list(iTriangle)%c, &
                                      distance, coords)) then ! check for intersection with triangle
                if(minDistance < 0 .or. distance < minDistance) then ! if first hit or shorter distance
                    minIndex = iTriangle
                    minDistance = distance
                    minC = coords
                end if
            end if
        end do

        if(minIndex == 0) minDistance = 0
        if(present(minCoords)) minCoords = minC

    end subroutine findFirstIntersectionTriangleList



    subroutine findTrianglesInAxesAlignedBox(boxMin, boxMax, trianglesIn, trianglesOut, &
                                             normalsIn, normalsOut) ! optional
        ! Given an axes-aligned bounding box and a list of triangles, returns an array of
        ! all triangles that intersect the given box.
        ! Optionally returns a list of triangle normals correpsonding to the found triangles
        use vector_mod, only: vector

        implicit none
        type(vector), intent(in) :: boxMin, boxMax
        type(triangle), pointer, dimension(:), intent(in) :: trianglesIn
        type(triangle), pointer, dimension(:), intent(out) :: trianglesOut
        type(triangle), pointer, dimension(:), intent(in), optional :: normalsIn
        type(triangle), pointer, dimension(:), intent(out), optional :: normalsOut

        integer :: i, triangleCount, aux
        type(triangle) :: thisT
        type(triangle), pointer, dimension(:) :: tempList

        ! allocate memory
        allocate(trianglesOut(0))
        if(present(normalsOut)) allocate(normalsOut(0))

        triangleCount = 0

        do i = 1, size(trianglesIn)
            thisT = trianglesIn(i)
            if(intersectAxesAlignedBoxTriangle(boxMin, boxMax, thisT%a, thisT%b, thisT%c)) then
                aux = appendTriangleArray(triangleCount, trianglesOut, thisT)
                if(present(normalsOut)) then
                    aux = appendTriangleArray(triangleCount, normalsOut, normalsIn(i))
                end if
                triangleCount = triangleCount + 1
            end if
        end do

        ! Copy to temporary array and back again to prevent memory leakage
        allocate(tempList(triangleCount))

        tempList = trianglesOut(1:triangleCount)
        deallocate(trianglesOut)
        allocate(trianglesOut(triangleCount))
        trianglesOut = tempList

        if(present(normalsOut)) then
            tempList = normalsOut(1:triangleCount)
            deallocate(normalsOut)
            allocate(normalsOut(triangleCount))
            normalsOut = tempList
        end if

        deallocate(tempList)

    contains

        integer function appendTriangleArray(curSize, array, val) result(newSize)
            ! Appends an element to given array, increasing the array size necessary
            ! If the array is too small, an array of twice the original size is allocated
            ! and existing elements are copied into the new array
            type(triangle), pointer, dimension(:) :: array, tmpArray
            type(triangle), intent(in) :: val
            integer, intent(in) :: curSize

            if(size(array) == curSize) then
                allocate( tmpArray(max(1,2*curSize)) )
                tmpArray(1:curSize) = array
                deallocate(array)
                array => tmpArray
            end if
            newSize = curSize + 1
            array(newSize) = val
        end function appendTriangleArray

    end subroutine findTrianglesInAxesAlignedBox


!
!    type(vector) function sumNormals(mesh) result(total)
!        ! Use Gauss' divergence theorem to see if the mesh is closed
!        use vector_mod, only: vector, crossProd, &
!                              operator(+), operator(-), operator(*)
!
!        implicit none
!        type(triangleMesh) :: mesh
!
!        integer :: iFace
!
!        total = vector(0,0,0)
!
!        do iFace = 1, mesh%nFace
!            total = total + crossProd(mesh%faces(iFace)%b - mesh%faces(iFace)%a, &
!                                      mesh%faces(iFace)%c - mesh%faces(iFace)%a)
!        end do
!
!        total = total * 0.5d0
!
!        return
!    end function sumNormals
!
!    type(vector) function findProbeForOpenTriangleMesh(mesh) result(probe)
!        use vector_mod, only: vector
!
!        implicit none
!        type(triangleMesh) :: mesh
!
!        type(vector) :: axis
!        real :: minProject, maxProject
!
!        axis = sumNormals(mesh)
!        call projectVerticesMinMax(mesh%nodes, axis, minProject, maxProject)
!
!    end function findProbeForOpenTriangleMesh
!

end module triangle_mesh_mod


