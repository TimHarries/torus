

program main
    use messages_mod
    use triangle_mesh_mod, only: triangleMesh, readTriangleMesh, printTriangleMeshSummary
    use vector_mod, only: vector, operator(+), operator(-), operator(*)

    use interface_hierarchy_mod

    use object_media_mod



    !implicit none
    !type(triangleMesh), pointer :: meshList(:)
    type(triangleMesh), allocatable :: meshList(:)
    !type(localInterfaceHierarchy), pointer :: pTree
    type(localInterfaceHierarchy) :: root, sub
    !type(vector), dimension(2) :: box = (/ vector(0,0,-1), vector(1,1,0) /)
    type(vector), dimension(2) :: box = (/ vector(-1,-1,-1), vector(1,1,0) /)
    type(vector), dimension(2) :: ray = (/ vector(0,0,0), vector(0.2,0.5,-0.5) /)

    type(objectIntersection) :: isect

    type(medium) :: vcm

    integer :: i, curType

    writeoutput = .true.
    verbositylevel = 5
    outputInfo = .true.
    outputWarnings = .true.

    vcm = createNewMedium('vcm', (/ 0.d0, 1.d0 /) )
    print *, vcm%label

    !meshList => readTriangleMesh("MESH02.obj")
    meshList = readTriangleMesh("MESH02.obj")

    meshList = meshList(size(meshList):1:-1)

    do i = 1, size(meshList)
        call printTriangleMeshSummary(meshList(i),0)
    end do


    print *, 'Local hierarchy:'
    !pTree => makeLocalInterfaceHierarchy( meshList, box(1), box(2))
    root = makeLocalInterfaceHierarchy( meshList, box(1), box(2))
    call printLocalInterfaceHierarchy(root)
    print *, 'combined ids', combineAllIds(root)
    print *, 'triangles in tree:', countInterfaceTriangles(root)

    curType = findGlobalObjectId(root, ray(1))
    print *, 'type at ', ray(1), ':', curType
    isect = findFirstIntersection(root, ray(1), ray(2), curType)

    print *, 'entering object: ', isect%globalObjectId
    print *, 'at posistion: ', ray(1) + ray(2)*isect%distance
    print *, 'and  normal: ', isect%normal

    sub = refineLocalInterfaceHierarchy( root, vector(-1,-1,-0.2), vector(-0.5,-0.5,0))
    call printLocalInterfaceHierarchy(sub)

    print *, 'combined ids', combineAllIds(sub)

    call freeTree(root)

end program main
