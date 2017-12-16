import PyNormaliz

ineq = [ [ 0,0,1 ], [ 0,1,0 ],[1,0,0] ];

cone = PyNormaliz.NmzCone(["cone", ineq ])

PyNormaliz.NmzResult( cone, "HilbertBasis" )

cone2 = PyNormaliz.NmzResult( cone, "IntegerHull" )

PyNormaliz.NmzResult( cone2, "HilbertBasis" )
