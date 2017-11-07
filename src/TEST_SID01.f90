

program main
    use messages_mod
    use triangle_mesh_mod
    use vector_mod, only: vector, intersectLineTriangle

    use object_tree_mod



    !implicit none
    type(triangleMesh), pointer :: ml(:)
    type(triangleMesh) :: m
    type(vector) :: v


    type(vector) :: p,d, a,b,c
    logical :: intersect
    real(kind=8) :: distance
    type(vector) :: coords

    type(localTriangleMeshHierarchy), pointer :: pTree


    type(vector), dimension(3) :: t = (/ vector(12,9,9), vector(9,12,9), vector(19,19,20) /)
    type(vector), dimension(2) :: box = (/ vector(-10,-10,-10), vector(10,10,10) /)


    writeoutput = .true.
    verbositylevel = 5
    outputInfo = .true.
    outputWarnings = .true.

    ml => readTriangleMesh("MESH01.obj")
    m = ml(1)

    call printTriangleMeshSummary(m,0)

    p = vector(0,0,-2.d-16)
    d = vector(-1,-1,-3.d-16)
    a = vector(1,0,0)
    b = vector(1,1,0)
    c = vector(0,1,0)

    intersect = intersectLineTriangle(p,d,.false. ,a,b,c, distance,coords)
    print *, 'line-triangle test:', intersect, distance

    print *, 'triangle-box test:'
    print *, intersectAxesAlignedBoxTriangle(box(1), box(2), t(1), t(2), t(3))
    print *, intersectAxesAlignedBoxTriangle(box(1), box(2), m%faces(1)%a, m%faces(1)%b, m%faces(1)%c)

    intersect = isAxesAlignedBoxInTriangleMesh(box(1), box(2), m)
    print *, 'box-in-mesh test:', intersect

    intersect = isPointInTriangleMesh(vector(3,1,1), m)
    print *, 'point in mesh:', intersect


    intersect = intersectAxesAlignedBoxTriangleMesh(vector(0,0,0), vector(2,1,1), m)
    print *, 'box-mesh intersect:', intersect


    print *, 'local hierarchy'
    pTree => makeLocalTriangleMeshHierarchy( (/m/), vector(0,0,0), vector(2,1,1))
    print *, 'objectId', pTree%globalObjectId
    print *, 'reference inside: ', pTree%referenceIsInside
    if(associated(pTree%triangleList)) then
        print *, 'num triangles inside: ', size(pTree%triangleList)
    else
        print *, 'no triangles listed'
    end if

    !print *, findGlobalObjectId(pTree, vector(5,5,5))


end program main
