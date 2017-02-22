;2017/02/16
;Gonna need this one for transforming geo to mag. It's got your solution:
;http://jsoc.stanford.edu/~jsoc/keywords/Chris_Russel/Geophysical%20Coordinate%20Transformations.htm#s3
;Gonna need this one to brush up on Euler rotations: https://en.wikipedia.org/wiki/Euler_angles
;
;More enlightenment can be had in JOURNAL__20170221__GEO_TO_MAG_AND_VICE_VERSA__JUST_CHECK_BASIS_VECTOR_AND_COORDINATE_TRANSFORMS,
;where I think I show that the solution on the Stanford website above is wrong
PRO JOURNAL__20170216__GEO_TO_MAG_VECTOR_TRANSFORM__CHECK_INVERSION_OF_MATRIX_IN_NOTEBOOK,Z1Y2Z3=z1y2z3,BRO=bro

  COMPILE_OPT IDL2

  ;;east longitude of magnetic dipole, in degrees
  alpha   = (-1.)*69.761 * !DTOR

  ;;colatitude of geomagnetic pole
  beta    = 11.435 * !DTOR

  ;;third angle(?)
  gamma   = 0. * !DTOR

  ;;Angle things
  c1      = COS(alpha)
  s1      = SIN(alpha)
  
  c2      = COS(beta)
  s2      = SIN(beta)

  c3      = COS(gamma)
  s3      = SIN(gamma)

  ;;Assume rotation about Z first, then rotation about Y, then rotation about Z a second time
  z1y2z3  = [[c1*c2*c3 - s1*s3, (-1.)*c3*s1 - c1*c2*s3, c1*s2], $
             [c1*s3 + c2*c3*s1,       c1*c3 - c2*s1*s3, s1*s2], $
             [     (-1.)*c3*s2,                  s2*s3,    c2]]

  bro    = INVERT(z1y2z3)

  PRINT,bro

  ;; matrix = TRANSPOSE([[0.338,0.9192,-0.198], $
  ;;                     [-0.938,0.345,0.], $
  ;;                     [0.0683,0.1857,0.98]])

  ;; bro    = INVERT(matrix)

END
