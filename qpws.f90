!	This module defines quadratic potential with a step and its derivatives
!	with respect to phi
	MODULE potential

	IMPLICIT NONE

! 	Values of parameters involved

	DOUBLE PRECISION, PARAMETER :: m = 7.147378d-6, alpha=1.606d-3, &
	phi0=14.67d0, Deltaphi=311.0d-4
		
!	Initial conditions
		
	DOUBLE PRECISION, PARAMETER :: phii = 16.5d0, e1i = 2.0d-8, &
	phiNi = -(2.0d0*e1i)**0.5d0
  
	CONTAINS

	DOUBLE PRECISION FUNCTION V(phi)

	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(IN) :: phi

	V = (m*m*phi*phi/2.0d0)*(1.0d0+alpha*tanh((phi-phi0)/Deltaphi))

	END FUNCTION V

	DOUBLE PRECISION FUNCTION DV(phi)

	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(IN) :: phi

	DV = m*m*phi*(1.0d0 + alpha*tanh((phi -phi0)/Deltaphi)) &
	+ (m*m*phi*phi/2.0d0)*(alpha/Deltaphi)*(1.0d0/((cosh((phi -phi0)/Deltaphi))**2.0d0))

	END FUNCTION DV

	DOUBLE PRECISION FUNCTION DDV(phi) 

	IMPLICIT NONE

	DOUBLE PRECISION, INTENT(IN) :: phi

	DOUBLE PRECISION :: DDV1, DDV2  

	DDV1 = m*m*(1.0d0 + alpha*tanh((phi -phi0)/Deltaphi)) &
	+ (2.0d0*m*m*phi)*(alpha/Deltaphi)*(1.0d0/((cosh((phi -phi0)/Deltaphi))**2.0d0)) 

	DDV2 = - (m*m*phi*phi)*(alpha/(Deltaphi*Deltaphi)) &
	*(sinh((phi -phi0)/Deltaphi)/((cosh((phi -phi0)/Deltaphi))**3.0d0))

	DDV = DDV1 + DDV2	

	END FUNCTION DDV

	END MODULE potential
