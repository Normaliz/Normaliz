import PyQNormaliz_cpp
C = PyQNormaliz_cpp.NmzCone(number_field="min_poly (a2-2) embedding 1.4+/-0.1",cone=[[[1],[0,1]],[[1],[-1]]])
PyQNormaliz_cpp.NmzCompute( C, [ "SupportHyperplanes" ] )
PyQNormaliz_cpp.NmzResult( C, "ExtremeRays" )

def rat_handler(list):
    return list[0]/list[1]

PyQNormaliz_cpp.NmzResult( C, "ExtremeRays", RationalHandler=rat_handler, NumberfieldElementHandler=tuple, VectorHandler=tuple, MatrixHandler=tuple )
