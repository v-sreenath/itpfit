!	This code evaluates scalar-tensor cross-correlations and tensor bi-spectrum 
!	for any single field canonical scalar field inflationary model, for an 
!	arbitrary configuration of k. While evaluating three point functions, k1  
!	is fixed to be kp and k3/k1 and k2/k1 are varied from 0 to 1 and from 0.5 
!	to 1 respectively.
        PROGRAM powerspectra

	Use potential

	IMPLICIT NONE

!	Variable decalarations

	INTEGER :: I, J, Iheps, INEic, INEshs, JNEic, JNEshs, L, L2, L3

	INTEGER, PARAMETER :: NOSB = 1*10**8,NOSP = NOSB/4, Nk = 100

	DOUBLE PRECISION, PARAMETER :: NEi = 0.0d0, NEf = 69.0d0, Pi = 3.14159265359d0  

	DOUBLE PRECISION, PARAMETER :: kappa= (1.0d0/50.0d0)
	DOUBLE PRECISION :: NE, hb, NEend, NEic, NEshs
        DOUBLE PRECISION :: ai, Hi

	DOUBLE PRECISION, DIMENSION(1:NOSB) :: phi, phiN, phiNN, zNz, H, e1, e2, e2N

        DOUBLE PRECISION :: kphi1, phitemp1, kphi2, phitemp2, kphi3, &
	phitemp3, kphi4, phitemp4, kphi5, phitemp5, kphi6 
 	DOUBLE PRECISION ::  kphiN1, phiNtemp1, kphiN2, phiNtemp2, &
	kphiN3, phiNtemp3, kphiN4, phiNtemp4, kphiN5, phiNtemp5, kphiN6

	DOUBLE PRECISION, PARAMETER :: kp = 2.0d-3, kmin = 2.0d-5, kmax = kp
	DOUBLE PRECISION, PARAMETER :: ymin = log10(kmin), ymax = log10(kmax)
	DOUBLE PRECISION :: k,  k1, k2, k3, hp, y, y2, y2min, y2max, y3, y3min, y3max, hy, hy2,hy3
	DOUBLE PRECISION :: PSk1, PTk1, PSk2, PTk2, PSk3, PTk3, r

	COMPLEX*16, DIMENSION(1:NOSP) :: Rk1, Rk1N, Rk2, Rk2N, Rk3, Rk3N 
	COMPLEX*16 :: Rk1Nic, Rk1NNic, Rk2Nic, Rk2NNic, Rk3Nic, Rk3NNic	

	COMPLEX*16 :: kRk11, Rk1temp1, kRk12, Rk1temp2, kRk13, &
	Rk1temp3, kRk14, Rk1temp4, kRk15, Rk1temp5, kRk16
	COMPLEX*16 :: kRk1N1, Rk1Ntemp1, kRk1N2, Rk1Ntemp2, kRk1N3, &
	Rk1Ntemp3, kRk1N4, Rk1Ntemp4, kRk1N5, Rk1Ntemp5, kRk1N6
	COMPLEX*16 :: kRk21, Rk2temp1, kRk22, Rk2temp2, kRk23, &
	Rk2temp3, kRk24, Rk2temp4, kRk25, Rk2temp5, kRk26
	COMPLEX*16 :: kRk2N1, Rk2Ntemp1, kRk2N2, Rk2Ntemp2, kRk2N3, &
	Rk2Ntemp3, kRk2N4, Rk2Ntemp4, kRk2N5, Rk2Ntemp5, kRk2N6
	COMPLEX*16 :: kRk31, Rk3temp1, kRk32, Rk3temp2, kRk33, &
	Rk3temp3, kRk34, Rk3temp4, kRk35, Rk3temp5, kRk36
	COMPLEX*16 :: kRk3N1, Rk3Ntemp1, kRk3N2, Rk3Ntemp2, kRk3N3, &
	Rk3Ntemp3, kRk3N4, Rk3Ntemp4, kRk3N5, Rk3Ntemp5, kRk3N6

	COMPLEX*16, DIMENSION(1:NOSP) :: hk1, hk1N, hk2, hk2N, hk3, hk3N 
	COMPLEX*16 :: hk1Nic, hk1NNic, hk2Nic, hk2NNic, hk3Nic, hk3NNic

	COMPLEX*16 :: khk11, hk1temp1, khk12, hk1temp2, khk13, &
	hk1temp3, khk14, hk1temp4, khk15, hk1temp5, khk16
	COMPLEX*16 :: khk1N1, hk1Ntemp1, khk1N2, hk1Ntemp2, khk1N3, &
	hk1Ntemp3, khk1N4, hk1Ntemp4, khk1N5, hk1Ntemp5, khk1N6
	COMPLEX*16 :: khk21, hk2temp1, khk22, hk2temp2, khk23, &
	hk2temp3, khk24, hk2temp4, khk25, hk2temp5, khk26
	COMPLEX*16 :: khk2N1, hk2Ntemp1, khk2N2, hk2Ntemp2, khk2N3, &
	hk2Ntemp3, khk2N4, hk2Ntemp4, khk2N5, hk2Ntemp5, khk2N6
	COMPLEX*16 :: khk31, hk3temp1, khk32, hk3temp2, khk33, &
	hk3temp3, khk34, hk3temp4, khk35, hk3temp5, khk36
	COMPLEX*16 :: khk3N1, hk3Ntemp1, khk3N2, hk3Ntemp2, khk3N3, &
	hk3Ntemp3, khk3N4, hk3Ntemp4, khk3N5, hk3Ntemp5, khk3N6

	COMPLEX*16 :: CI
	COMPLEX*16 :: cGsss13I, cGsss13Ic, cGsss2I, cGsss2Ic, cGsss4I, &
	cGsss4Ic, cGsss56I, cGsss56Ic, &
        Gsss13, Gsss2, Gsss4, Gsss56, Gsss7, Gsss47, k6Gsss13, k6Gsss2, &
	k6Gsss4, k6Gsss56, k6Gsss7, k6Gsss47, &
	fnl13, fnl2, fnl4, fnl56, fnl7, fnl47
	COMPLEX*16 :: cGgggI, cGgggIc, Gggg, k6Gggg, hnl
	COMPLEX*16 :: cGssg1I, cGssg1Ic, Gssg1, k6Gssg1, cGssg2I, cGssg2Ic, Gssg2, k6Gssg2, &
	cGssg3I, cGssg3Ic, Gssg3, k6Gssg3, CnlR1, CnlR2, CnlR3, cnls, Gssg
	COMPLEX*16 :: cGsgg1I, cGsgg1Ic, Gsgg1, k6Gsgg1, cGsgg2I, cGsgg2Ic, Gsgg2, k6Gsgg2, &
	cGsgg3I, cGsgg3Ic, Gsgg3, k6Gsgg3, Cnlg1, Cnlg2, Cnlg3,cnlg, Gsgg

	hb = (NEf - NEi)/NOSB

!	Print *, "This code evalautes scalar-temsor crosscorrelations and tensor &
!	bispectrum for a quadratic potential with a step"  
!	Print *, "The initial value of the field is:", phii
!	Print *, "The initial 'velocity' of the field is:", phiNi

! 	Solving for the background using Butcher's fifth order Runge-Kutta method

        phi(1) = phii
        phiN(1) = phiNi
	Hi =(2.0d0*V(phii)/(6.0d0 - phiNi*phiNi))**0.5d0 
	NE = NEi
	
	Do I=1,NOSB-1,1 	

	kphi1= phiN(I)
	phitemp1 = phi(I) + (1.0d0/4.0d0)*hb*kphi1 
	kphiN1 =  -(6.0d0 - phiN(I)*phiN(I))*(phiN(I)/2.0d0) &
	- ((6.0d0 - phiN(I)*phiN(I))/2.0d0)*(DV(phi(I))/V(phi(I))) 
	phiNtemp1 = phiN(I) + (1.0d0/4.0d0)*hb*kphiN1 

	kphi2 = phiNtemp1
	phitemp2 =  phi(I) + (1.0d0/8.0d0)*hb*kphi1 + (1.0d0/8.0d0)*hb*kphi2
	kphiN2 =  -(6.0d0 - phiNtemp1*phiNtemp1)*(phiNtemp1/2.0d0) &
	- ((6.0d0 - phiNtemp1*phiNtemp1)/2.0d0)*(DV(phitemp1)/V(phitemp1))  
	phiNtemp2 = phiN(I) + (1.0d0/8.0d0)*hb*kphiN1 + (1.0d0/8.0d0)*hb*kphiN2

	kphi3 = phiNtemp2
	phitemp3 = phi(I) - (1.0d0/2.0d0)*hb*kphi2 + hb*kphi3 
	kphiN3 = -(6.0d0 - phiNtemp2*phiNtemp2)*(phiNtemp2/2.0d0) &
	- ((6.0d0 - phiNtemp2*phiNtemp2)/2.0d0)*(DV(phitemp2)/V(phitemp2)) 
	phiNtemp3 = phiN(I) - (1.0d0/2.0d0)*hb*kphiN2 + hb*kphiN3

	kphi4 = phiNtemp3 
	phitemp4 =  phi(I) + (3.0d0/16.0d0)*hb*kphi1 + (9.0d0/16.0d0)*hb*kphi4 
	kphiN4 = -(6.0d0 - phiNtemp3*phiNtemp3)*(phiNtemp3/2.0d0) &
	- ((6.0d0 - phiNtemp3*phiNtemp3)/2.0d0)*(DV(phitemp3)/V(phitemp3))
 	phiNtemp4 = phiN(I) + (3.0d0/16.0d0)*hb*kphiN1 + (9.0d0/16.0d0)*hb*kphiN4 

	kphi5 = phiNtemp4
	phitemp5 =  phi(I) - (3.0d0/7.0d0)*hb*kphi1 + (2.0d0/7.0d0)*hb*kphi2 &
	+ (12.0d0/7.0d0)*hb*kphi3 - (12.0d0/7.0d0)*hb*kphi4 + (8.0d0/7.0d0)*hb*kphi5
	kphiN5 = -(6.0d0 - phiNtemp4*phiNtemp4)*(phiNtemp4/2.0d0) &
	- ((6.0d0 - phiNtemp4*phiNtemp4)/2.0d0)*(DV(phitemp4)/V(phitemp4))
	phiNtemp5 =  phiN(I) - (3.0d0/7.0d0)*hb*kphiN1 + (2.0d0/7.0d0)*hb*kphiN2 &
	+ (12.0d0/7.0d0)*hb*kphiN3 - (12.0d0/7.0d0)*hb*kphiN4 + (8.0d0/7.0d0)*hb*kphiN5

	kphi6 = phiNtemp5
	kphiN6 = -(6.0d0 - phiNtemp5*phiNtemp5)*(phiNtemp5/2.0d0) &
	- ((6.0d0 - phiNtemp5*phiNtemp5)/2.0d0)*(DV(phitemp5)/V(phitemp5))

	phi(I + 1) = phi(I) + (hb/90.0d0)*(7.0d0*kphi1 + 32.0d0*kphi3 &
	+ 12.0d0*kphi4 + 32.0d0*kphi5 + 7.0d0*kphi6)
	phiN(I + 1) = phiN(I) + (hb/90.0d0)*(7.0d0*kphiN1 + 32.0d0*kphiN3 &
	+ 12.0d0*kphiN4 + 32.0d0*kphiN5 + 7.0d0*kphiN6)

	phiNN(I) = -(6.0d0 - phiN(I)*phiN(I))*(phiN(I)/2.0d0) &
	- ((6.0d0 - phiN(I)*phiN(I))/2.0d0)*(DV(phi(I))/V(phi(I)))  
	zNz(I) = 1.0d0 + (phiNN(I)/phiN(I))
	H(I) = (2.0d0*V(phi(I))/(6.0d0 - phiN(I)*phiN(I)))**0.5d0 
	e1(I) = 0.5d0*phiN(I)*phiN(I)     
	e2(I)= -6.0d0 + 2.0d0*e1(I) - (2.0d0*DV(phi(I))/(H(I)*H(I)*phiN(I)))
	e2N(I)=-((2.0d0*DDV(phi(I)))/(H(I)*H(I))) + 12.0d0*e1(I) - 3.0d0*e2(I) &
	- 4.0d0*e1(I)*e1(I) + 5.0d0*e1(I)*e2(I) - 0.5d0*e2(I)*e2(I)

!	Write(10,*) NE, phi(I) 
!	Write(20,*) NE, phiN(I) 
!	Write(30,*) NE, H(I)/Hi
!	Write(45,*) NE, e1(I)
!	Write(50,*) NE, phi(I), phiN(I), phiNN(I), H(I), e1(I) 

!	Write(40,*) NE, e2(I)

!       Evolving the scalar and tensor perturbations  

	If (Abs(e1(I)-1.0d0).le.1.0d-3) then 
        NEend = NE 
	exit 
	end if

        NE = NEi+I*hb 

	end do

	Print *,"Inflation ends at the e-fold:", NEend

  	Print *, "The value of the field when inflation ends is:", phi(I + 1) 

	Iheps=Int(((NEend-50.0d0)-NEi)/hb)

!       The pivot scale and the value of the scale factor at NEi

        ai = kp/(Exp(NEend-50.0d0)*H(Iheps))

	Print *, "The initial value of the scale factor is:", ai

! 	Evaluating the power spectrum

	hp = (NEf - NEi)/NOSP
	hy=(ymax-ymin)/Nk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	Finding NEic and NEshs

	JNEic = NOSP 

	Do J=1,NOSP-1,1 
	NE = NEi+(J-1)*hp 
	If (J .gt. JNEic) go to 12
	If (Abs((kmin/(ai*Exp(NE)*H(4*J-3)))-1000.0d0).le.1.0d-1) then
	NEic = NE 
	JNEic = J
	INEic=4*JNEic-3
	end if

12  	If (Abs((kmax/(ai*Exp(NE)*H(4*J-3)))-1.0d-5).le.1.0d-3) then
        NEshs = NE
 	JNEshs = J
	INEshs=4*JNEshs-3
	exit 
	end if

	end do

!	Print *, NEic
!	Print *, NEshs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!Evolving fixed mode k1 using Butcher's fifth order Runge-Kutta method

	CI=DCMPLX(0,1.0d0)
	k1 = kp
	
	Do J = JNEic,JNEshs-1,1 
	NE = NEi+(J-1)*hp 
	
	If (J .eq. JNEic) then
	Rk1Nic = DCMPLX(((2.0d0*k1)**(-0.5d0))*(1.0d0/(ai*Exp(NEic)*phiN(INEic))),0.0d0)
	Rk1NNic = DCMPLX(-((2.0d0*k1)**(-0.5d0))*(zNz(INEic)/(ai*Exp(NEic)*phiN(INEic))), &
        -(((k1/2.0d0)**(0.5d0))/(ai*ai*Exp(NEic)*Exp(NEic)*H(INEic)*phiN(INEic))))

	Rk1(JNEic) = Rk1Nic
	Rk1N(JNEic) = Rk1NNic

	hk1Nic = DCMPLX(((2.0d0)**(0.5d0))*((2.0d0*k1)**(-0.5d0))*(1.0d0/(ai*Exp(NEic))),0.0d0)
	hk1NNic = DCMPLX(-((2.0d0)**(0.5d0))*((2.0d0*k1)**(-0.5d0))*(1.0d0/(ai*Exp(NEic))), &
        -((2.0d0)**(0.5d0))*(((k1/2.0d0)**(0.5d0))/(ai*ai*Exp(NEic)*Exp(NEic)*H(INEic))))

	hk1(JNEic) = hk1Nic
	hk1N(JNEic) = hk1NNic
	end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	kRk11 = Rk1N(J) 
	kRk1N1 = -(1.0d0 - e1(4*J-3) + 2.0d0*zNz(4*J-3))*Rk1N(J)-((k1/(ai*Exp(NE)*H(4*J-3)))**2.0d0)*Rk1(J)
	Rk1temp1 = Rk1(J) + (1.0d0/4.0d0)*hp*kRk11
	Rk1Ntemp1 = Rk1N(J) + (1.0d0/4.0d0)*hp*kRk1N1 

	kRk12 = Rk1Ntemp1
	kRk1N2 = -(1.0d0 - e1(4*J-2) + 2.0d0*zNz(4*J-2))*Rk1Ntemp1 -((k1/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*Rk1temp1
	Rk1temp2 =  Rk1(J) + (1.0d0/8.0d0)*hp*kRk11 + (1.0d0/8.0d0)*hp*kRk12 
	Rk1Ntemp2 = Rk1N(J) + (1.0d0/8.0d0)*hp*kRk1N1 + (1.0d0/8.0d0)*hp*kRk1N2 

	kRk13 = Rk1Ntemp2
	kRk1N3 = -(1.0d0 - e1(4*J-2) + 2.0d0*zNz(4*J-2))*Rk1Ntemp2 - ((k1/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*Rk1temp2
	Rk1temp3 = Rk1(J) - (1.0d0/2.0d0)*hp*kRk12  + hp*kRk13  
	Rk1Ntemp3 = Rk1N(J) - (1.0d0/2.0d0)*hp*kRk1N2  + hp*kRk1N3  

	kRk14 = Rk1Ntemp3 
	kRk1N4 = -(1.0d0 - e1(4*J-1) + 2.0d0*zNz(4*J-1))*Rk1Ntemp3 - ((k1/(ai*Exp(NE+0.5d0*hp)*H(4*J-1)))**2.0d0)*Rk1temp3
	Rk1temp4 = Rk1(J) + (3.0d0/16.0d0)*hp*kRk11  +  (9.0d0/16.0d0)*hp*kRk14  
	Rk1Ntemp4 = Rk1N(J) + (3.0d0/16.0d0)*hp*kRk1N1  +  (9.0d0/16.0d0)*hp*kRk1N4

	kRk15 = Rk1Ntemp4 
	kRk1N5 = -(1.0d0 - e1(4*J) + 2.0d0*zNz(4*J))*Rk1Ntemp4 - ((k1/(ai*Exp(NE+0.75d0*hp)*H(4*J)))**2.0d0)*Rk1temp4
	Rk1temp5 = Rk1(J) - (3.0d0/7.0d0)*hp*kRk11  +  (2.0d0/7.0d0)*hp*kRk12 &
	+  (12.0d0/7.0d0)*hp*kRk13  -  (12.0d0/7.0d0)*hp*kRk14 + (8.0d0/7.0d0)*hp*kRk15
	Rk1Ntemp5 = Rk1N(J) - (3.0d0/7.0d0)*hp*kRk1N1  +  (2.0d0/7.0d0)*hp*kRk1N2 &
	+  (12.0d0/7.0d0)*hp*kRk1N3 - (12.0d0/7.0d0)*hp*kRk1N4 + (8.0d0/7.0d0)*hp*kRk1N5

	kRk16 = Rk1Ntemp5
	kRk1N6 = -(1.0d0 - e1(4*J+1) + 2.0d0*zNz(4*J+1))*Rk1Ntemp5 - ((k1/(ai*Exp(NE+1.0d0*hp)*H(4*J+1)))**2.0d0)*Rk1temp5

	Rk1(J + 1) = Rk1(J) + (hp/90.0d0)*(7.0d0*kRk11 + 32.0d0*kRk13 &
	+ 12.0d0*kRk14 + 32.0d0*kRk15 + 7.0d0*kRk16)
	Rk1N(J + 1) = Rk1N(J) + (hp/90.0d0)*(7.0d0*kRk1N1 + 32.0d0*kRk1N3 &
	+ 12.0d0*kRk1N4 + 32.0d0*kRk1N5 + 7.0d0*kRk1N6)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	khk11 = hk1N(J)
	khk1N1 = -(3.0d0 - e1(4*J-3))*hk1N(J)-((k1/(ai*Exp(NE)*H(4*J-3)))**2.0d0)*hk1(J)
	hk1temp1 = hk1(J) + (1.0d0/4.0d0)*hp*khk11 
	hk1Ntemp1 = hk1N(J) + (1.0d0/4.0d0)*hp*khk1N1 

	khk12 = hk1Ntemp1
	khk1N2 = -(3.0d0 - e1(4*J-2))*hk1Ntemp1 -((k1/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*hk1temp1
	hk1temp2 =  hk1(J) + (1.0d0/8.0d0)*hp*khk11 + (1.0d0/8.0d0)*hp*khk12
	hk1Ntemp2 = hk1N(J) + (1.0d0/8.0d0)*hp*khk1N1 + (1.0d0/8.0d0)*hp*khk1N2

	khk13 = hk1Ntemp2 
	khk1N3 = -(3.0d0 - e1(4*J-2))*hk1Ntemp2 - ((k1/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*hk1temp2
	hk1temp3 = hk1(J) - (1.0d0/2.0d0)*hp*khk12 + hp*khk13
	hk1Ntemp3 = hk1N(J) - (1.0d0/2.0d0)*hp*khk1N2 + hp*khk1N3

	khk14 = hk1Ntemp3 
	khk1N4 = -(3.0d0 - e1(4*J-1))*hk1Ntemp3 - ((k1/(ai*Exp(NE+0.5d0*hp)*H(4*J-1)))**2.0d0)*hk1temp3
	hk1temp4 = hk1(J) + (3.0d0/16.0d0)*hp*khk11  +  (9.0d0/16.0d0)*hp*khk14  
	hk1Ntemp4 = hk1N(J) + (3.0d0/16.0d0)*hp*khk1N1  +  (9.0d0/16.0d0)*hp*khk1N4

	khk15 = hk1Ntemp4
	khk1N5 = -(3.0d0 - e1(4*J))*hk1Ntemp4 - ((k1/(ai*Exp(NE+0.75d0*hp)*H(4*J)))**2.0d0)*hk1temp4
	hk1temp5 = hk1(J) - (3.0d0/7.0d0)*hp*khk11  +  (2.0d0/7.0d0)*hp*khk12 &
	+  (12.0d0/7.0d0)*hp*khk13  -  (12.0d0/7.0d0)*hp*khk14 + (8.0d0/7.0d0)*hp*khk15
	hk1Ntemp5 = hk1N(J) - (3.0d0/7.0d0)*hp*khk1N1  +  (2.0d0/7.0d0)*hp*khk1N2 &
	+  (12.0d0/7.0d0)*hp*khk1N3 - (12.0d0/7.0d0)*hp*khk1N4 + (8.0d0/7.0d0)*hp*khk1N5

	khk16 = hk1Ntemp5
	khk1N6 = -(3.0d0 - e1(4*J+1))*hk1Ntemp5 - ((k1/(ai*Exp(NE+1.0d0*hp)*H(4*J+1)))**2.0d0)*hk1temp5

	hk1(J + 1) = hk1(J) + (hp/90.0d0)*(7.0d0*khk11 + 32.0d0*khk13 &
	+ 12.0d0*khk14 + 32.0d0*khk15 + 7.0d0*khk16)
	hk1N(J + 1) = hk1N(J) + (hp/90.0d0)*(7.0d0*khk1N1 + 32.0d0*khk1N3 &
	+ 12.0d0*khk1N4 + 32.0d0*khk1N5 + 7.0d0*khk1N6)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!	write(201,*)NE,hk1(J),hk1N(J)

	end do

	PSk1=((k1**3.0d0)/(2.0d0*(Pi**2.0d0)))*((Rk1(JNEshs))*DCONJG(Rk1(JNEshs)))
	
	PTk1=4.0d0*((k1**3.0d0)/(2.0d0*(Pi**2.0d0)))*((hk1(JNEshs))*DCONJG(hk1(JNEshs)))

	write(51,*) k1, PSk1, PTk1
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	Here, we evolve the other two modes (k2 and k3) and also integrate over 
!	N using Boole's rule to calculate three point cross-correlations.
	y3min = log10(kmin/k1)
	y3max = log10(kmax/k1)
	hy3 = (y3max-y3min)/Nk

        Do L3 = 1, Nk, 1
        y3 = y3min + hy3*L3 !-1)
 	k3=k1*10.0d0**y3

! k3
	Rk3Nic = DCMPLX(((2.0d0*k3)**(-0.5d0))*(1.0d0/(ai*Exp(NEic)*phiN(INEic))),0.0d0)
	Rk3NNic = DCMPLX(-((2.0d0*k3)**(-0.5d0))*(zNz(INEic)/(ai*Exp(NEic)*phiN(INEic))), &
        -(((k3/2.0d0)**(0.5d0))/(ai*ai*Exp(NEic)*Exp(NEic)*H(INEic)*phiN(INEic))))

	Rk3(JNEic) = Rk3Nic
	Rk3N(JNEic) = Rk3NNic

	hk3Nic = DCMPLX(((2.0d0)**(0.5d0))*((2.0d0*k3)**(-0.5d0))*(1.0d0/(ai*Exp(NEic))),0.0d0)
	hk3NNic = DCMPLX(-((2.0d0)**(0.5d0))*((2.0d0*k3)**(-0.5d0))*(1.0d0/(ai*Exp(NEic))), &
        -((2.0d0)**(0.5d0))*(((k3/2.0d0)**(0.5d0))/(ai*ai*Exp(NEic)*Exp(NEic)*H(INEic))))

	hk3(JNEic) = hk3Nic
	hk3N(JNEic) = hk3NNic

	Do J = JNEic,JNEshs-1,1 
	NE = NEi+(J-1)*hp 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	kRk31 = Rk3N(J) 
	kRk3N1 = -(1.0d0 - e1(4*J-3) + 2.0d0*zNz(4*J-3))*Rk3N(J)-((k3/(ai*Exp(NE)*H(4*J-3)))**2.0d0)*Rk3(J)
	Rk3temp1 = Rk3(J) + (1.0d0/4.0d0)*hp*kRk31
	Rk3Ntemp1 = Rk3N(J) + (1.0d0/4.0d0)*hp*kRk3N1 

	kRk32 = Rk3Ntemp1
	kRk3N2 = -(1.0d0 - e1(4*J-2) + 2.0d0*zNz(4*J-2))*Rk3Ntemp1 -((k3/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*Rk3temp1
	Rk3temp2 =  Rk3(J) + (1.0d0/8.0d0)*hp*kRk31 + (1.0d0/8.0d0)*hp*kRk32 
	Rk3Ntemp2 = Rk3N(J) + (1.0d0/8.0d0)*hp*kRk3N1 + (1.0d0/8.0d0)*hp*kRk3N2 

	kRk33 = Rk3Ntemp2
	kRk3N3 = -(1.0d0 - e1(4*J-2) + 2.0d0*zNz(4*J-2))*Rk3Ntemp2 - ((k3/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*Rk3temp2
	Rk3temp3 = Rk3(J) - (1.0d0/2.0d0)*hp*kRk32  + hp*kRk33  
	Rk3Ntemp3 = Rk3N(J) - (1.0d0/2.0d0)*hp*kRk3N2  + hp*kRk3N3  

	kRk34 = Rk3Ntemp3 
	kRk3N4 = -(1.0d0 - e1(4*J-1) + 2.0d0*zNz(4*J-1))*Rk3Ntemp3 - ((k3/(ai*Exp(NE+0.5d0*hp)*H(4*J-1)))**2.0d0)*Rk3temp3
	Rk3temp4 = Rk3(J) + (3.0d0/16.0d0)*hp*kRk31  +  (9.0d0/16.0d0)*hp*kRk34  
	Rk3Ntemp4 = Rk3N(J) + (3.0d0/16.0d0)*hp*kRk3N1  +  (9.0d0/16.0d0)*hp*kRk3N4

	kRk35 = Rk3Ntemp4 
	kRk3N5 = -(1.0d0 - e1(4*J) + 2.0d0*zNz(4*J))*Rk3Ntemp4 - ((k3/(ai*Exp(NE+0.75d0*hp)*H(4*J)))**2.0d0)*Rk3temp4
	Rk3temp5 = Rk3(J) - (3.0d0/7.0d0)*hp*kRk31  +  (2.0d0/7.0d0)*hp*kRk32 &
	+  (12.0d0/7.0d0)*hp*kRk33  -  (12.0d0/7.0d0)*hp*kRk34 + (8.0d0/7.0d0)*hp*kRk35
	Rk3Ntemp5 = Rk3N(J) - (3.0d0/7.0d0)*hp*kRk3N1  +  (2.0d0/7.0d0)*hp*kRk3N2 &
	+  (12.0d0/7.0d0)*hp*kRk3N3 - (12.0d0/7.0d0)*hp*kRk3N4 + (8.0d0/7.0d0)*hp*kRk3N5

	kRk36 = Rk3Ntemp5
	kRk3N6 = -(1.0d0 - e1(4*J+1) + 2.0d0*zNz(4*J+1))*Rk3Ntemp5 &
	- ((k3/(ai*Exp(NE+1.0d0*hp)*H(4*J+1)))**2.0d0)*Rk3temp5

	Rk3(J + 1) = Rk3(J) + (hp/90.0d0)*(7.0d0*kRk31 + 32.0d0*kRk33 &
	+ 12.0d0*kRk34 + 32.0d0*kRk35 + 7.0d0*kRk36)
	Rk3N(J + 1) = Rk3N(J) + (hp/90.0d0)*(7.0d0*kRk3N1 + 32.0d0*kRk3N3 &
	+ 12.0d0*kRk3N4 + 32.0d0*kRk3N5 + 7.0d0*kRk3N6)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	khk31 = hk3N(J)
	khk3N1 = -(3.0d0 - e1(4*J-3))*hk3N(J)-((k3/(ai*Exp(NE)*H(4*J-3)))**2.0d0)*hk3(J)
	hk3temp1 = hk3(J) + (1.0d0/4.0d0)*hp*khk31 
	hk3Ntemp1 = hk3N(J) + (1.0d0/4.0d0)*hp*khk3N1 

	khk32 = hk3Ntemp1
	khk3N2 = -(3.0d0 - e1(4*J-2))*hk3Ntemp1 -((k3/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*hk3temp1
	hk3temp2 =  hk3(J) + (1.0d0/8.0d0)*hp*khk31 + (1.0d0/8.0d0)*hp*khk32
	hk3Ntemp2 = hk3N(J) + (1.0d0/8.0d0)*hp*khk3N1 + (1.0d0/8.0d0)*hp*khk3N2

	khk33 = hk3Ntemp2 
	khk3N3 = -(3.0d0 - e1(4*J-2))*hk3Ntemp2 - ((k3/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*hk3temp2
	hk3temp3 = hk3(J) - (1.0d0/2.0d0)*hp*khk32 + hp*khk33
	hk3Ntemp3 = hk3N(J) - (1.0d0/2.0d0)*hp*khk3N2 + hp*khk3N3

	khk34 = hk3Ntemp3 
	khk3N4 = -(3.0d0 - e1(4*J-1))*hk3Ntemp3 - ((k3/(ai*Exp(NE+0.5d0*hp)*H(4*J-1)))**2.0d0)*hk3temp3
	hk3temp4 = hk3(J) + (3.0d0/16.0d0)*hp*khk31  +  (9.0d0/16.0d0)*hp*khk34  
	hk3Ntemp4 = hk3N(J) + (3.0d0/16.0d0)*hp*khk3N1  +  (9.0d0/16.0d0)*hp*khk3N4

	khk35 = hk3Ntemp4
	khk3N5 = -(3.0d0 - e1(4*J))*hk3Ntemp4 - ((k3/(ai*Exp(NE+0.75d0*hp)*H(4*J)))**2.0d0)*hk3temp4
	hk3temp5 = hk3(J) - (3.0d0/7.0d0)*hp*khk31  +  (2.0d0/7.0d0)*hp*khk32 &
	+  (12.0d0/7.0d0)*hp*khk33  -  (12.0d0/7.0d0)*hp*khk34 + (8.0d0/7.0d0)*hp*khk35
	hk3Ntemp5 = hk3N(J) - (3.0d0/7.0d0)*hp*khk3N1  +  (2.0d0/7.0d0)*hp*khk3N2 &
	+  (12.0d0/7.0d0)*hp*khk3N3 - (12.0d0/7.0d0)*hp*khk3N4 + (8.0d0/7.0d0)*hp*khk3N5

	khk36 = hk3Ntemp5
	khk3N6 = -(3.0d0 - e1(4*J+1))*hk3Ntemp5 - ((k3/(ai*Exp(NE+1.0d0*hp)*H(4*J+1)))**2.0d0)*hk3temp5

	hk3(J + 1) = hk3(J) + (hp/90.0d0)*(7.0d0*khk31 + 32.0d0*khk33 &
	+ 12.0d0*khk34 + 32.0d0*khk35 + 7.0d0*khk36)
	hk3N(J + 1) = hk3N(J) + (hp/90.0d0)*(7.0d0*khk3N1 + 32.0d0*khk3N3 &
	+ 12.0d0*khk3N4 + 32.0d0*khk3N5 + 7.0d0*khk3N6)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end do

	PSk3=((k3**3.0d0)/(2.0d0*(Pi**2.0d0)))*((Rk3(JNEshs))*DCONJG(Rk3(JNEshs)))

	PTk3=4.0d0*((k3**3.0d0)/(2.0d0*(Pi**2.0d0)))*((hk3(JNEshs))*DCONJG(hk3(JNEshs)))

!	write(53,*) k3, PSk3, PTk3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(k3 .le. kp/2.0d0)then
	y2min = log10((k1-k3)/k1)
	y2max = log10(k1/k1)
	else
	y2min = log10(k3/k1)
	y2max = log10(k1/k1)
	end if
	hy2 = (y2max-y2min)/Nk

	Do L2 = 1, Nk, 1

	y2 = y2min + hy2*(L2-1)
 	k2=k1*10.0d0**y2

!	write(202,*)y2min, y2max, y2, k2

	Rk2Nic = DCMPLX(((2.0d0*k2)**(-0.5d0))*(1.0d0/(ai*Exp(NEic)*phiN(INEic))),0.0d0)
	Rk2NNic = DCMPLX(-((2.0d0*k2)**(-0.5d0))*(zNz(INEic)/(ai*Exp(NEic)*phiN(INEic))), &
        -(((k2/2.0d0)**(0.5d0))/(ai*ai*Exp(NEic)*Exp(NEic)*H(INEic)*phiN(INEic))))

	Rk2(JNEic) = Rk2Nic
	Rk2N(JNEic) = Rk2NNic

	hk2Nic = DCMPLX(((2.0d0)**(0.5d0))*((2.0d0*k2)**(-0.5d0))*(1.0d0/(ai*Exp(NEic))),0.0d0)
	hk2NNic = DCMPLX(-((2.0d0)**(0.5d0))*((2.0d0*k2)**(-0.5d0))*(1.0d0/(ai*Exp(NEic))), &
        -((2.0d0)**(0.5d0))*(((k2/2.0d0)**(0.5d0))/(ai*ai*Exp(NEic)*Exp(NEic)*H(INEic))))

	hk2(JNEic) = hk2Nic
	hk2N(JNEic) = hk2NNic

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGgggI = -7.0d0*(1.0d0/4.0d0)*CI*(k1**2.0d0+k2**2.0d0+k3**2.0d0) &
	*(ai*Exp(NEic)/H(INEic))*DCONJG(hk1(JNEic)) &
	*DCONJG(hk2(JNEic))*DCONJG(hk3(JNEic)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

	cGgggIc = 7.0d0*(1.0d0/4.0d0)*CI*(k1**2.0d0+k2**2.0d0+k3**2.0d0) &
	*(ai*Exp(NEic)/H(INEic))*hk1(JNEic)*hk2(JNEic) &
	*hk3(JNEic)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ssg
	cGssg1I = -7.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NEic)*e1(INEic)/H(INEic))*DCONJG(hk3(JNEic)) &
	*DCONJG(Rk1(JNEic))*DCONJG(Rk2(JNEic))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))
	cGssg1Ic = 7.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NEic)*e1(INEic)/H(INEic))*hk3(JNEic) &
	*Rk1(JNEic)*Rk2(JNEic)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg2I = 7.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NEic)) &
	*H(INEic)*e1(INEic)*e1(INEic)*DCONJG(hk3(JNEic))*DCONJG(Rk1N(JNEic)) &
	*DCONJG(Rk2N(JNEic))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))
	cGssg2Ic = -7.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NEic)) &
	*H(INEic)*e1(INEic)*e1(INEic)*hk3(JNEic)*Rk1N(JNEic)*Rk2N(JNEic) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg3I = 7.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NEic))*H(INEic)*e1(INEic)*e1(INEic) &
	*((k2/k1)*DCONJG(hk3N(JNEic))*DCONJG(Rk1N(JNEic))*DCONJG(Rk2(JNEic)) &
	+(k1/k2)*DCONJG(hk3N(JNEic))*DCONJG(Rk2N(JNEic))*DCONJG(Rk1(JNEic))) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))
	cGssg3Ic = -7.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NEic))*H(INEic)*e1(INEic)*e1(INEic) &
	*((k2/k1)*hk3N(JNEic)*Rk1N(JNEic)*Rk2(JNEic) &
	+(k1/k2)*hk3N(JNEic)*Rk2N(JNEic)*Rk1(JNEic)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sgg
	cGsgg1I = 7.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NEic))*H(INEic)*e1(INEic) &
	*DCONJG(hk2N(JNEic))*DCONJG(hk3N(JNEic))*DCONJG(Rk1(JNEic)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))
	cGsgg1Ic = -7.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NEic))*H(INEic)*e1(INEic)&
	*hk2N(JNEic)*hk3N(JNEic)*Rk1(JNEic)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg2I = -7.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NEic)*e1(INEic)/H(INEic)) &
	*DCONJG(hk2(JNEic))*DCONJG(hk3(JNEic))*DCONJG(Rk1(JNEic))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))
	cGsgg2Ic = 7.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NEic)*e1(INEic)/H(INEic)) &
	*hk2(JNEic)*hk3(JNEic)*Rk1(JNEic)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg3I = -7.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NEic))*H(INEic)*e1(INEic) &
	*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*DCONJG(hk3(JNEic))*DCONJG(hk2N(JNEic)) &
	*DCONJG(Rk1N(JNEic))+((k3*k3-k2*k2-k1*k1)/(k1*k1))*DCONJG(hk2(JNEic)) &
	*DCONJG(hk3N(JNEic))*DCONJG(Rk1N(JNEic)))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))
	cGsgg3Ic = 7.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NEic))*H(INEic)*e1(INEic) &
	*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*hk3(JNEic)*hk2N(JNEic)*Rk1N(JNEic) &
	+((k3*k3-k2*k2-k1*k1)/(k1*k1))*hk2(JNEic)*hk3N(JNEic) &
	*Rk1N(JNEic))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEic)*H(INEic)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Do J = JNEic,JNEshs-1,1 
	NE = NEi+(J-1)*hp 
!k2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	kRk21 = Rk2N(J) 
	kRk2N1 = -(1.0d0 - e1(4*J-3) + 2.0d0*zNz(4*J-3))*Rk2N(J) &
	- ((k2/(ai*Exp(NE)*H(4*J-3)))**2.0d0)*Rk2(J)
	Rk2temp1 = Rk2(J) + (1.0d0/4.0d0)*hp*kRk21
	Rk2Ntemp1 = Rk2N(J) + (1.0d0/4.0d0)*hp*kRk2N1 

	kRk22 = Rk2Ntemp1
	kRk2N2 = -(1.0d0 - e1(4*J-2) + 2.0d0*zNz(4*J-2))*Rk2Ntemp1 &
	-((k2/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*Rk2temp1
	Rk2temp2 =  Rk2(J) + (1.0d0/8.0d0)*hp*kRk21 + (1.0d0/8.0d0)*hp*kRk22 
	Rk2Ntemp2 = Rk2N(J) + (1.0d0/8.0d0)*hp*kRk2N1 + (1.0d0/8.0d0)*hp*kRk2N2 

	kRk23 = Rk2Ntemp2
	kRk2N3 = -(1.0d0 - e1(4*J-2) + 2.0d0*zNz(4*J-2))*Rk2Ntemp2 &
	- ((k2/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*Rk2temp2
	Rk2temp3 = Rk2(J) - (1.0d0/2.0d0)*hp*kRk22  + hp*kRk23  
	Rk2Ntemp3 = Rk2N(J) - (1.0d0/2.0d0)*hp*kRk2N2  + hp*kRk2N3  

	kRk24 = Rk2Ntemp3 
	kRk2N4 = -(1.0d0 - e1(4*J-1) + 2.0d0*zNz(4*J-1))*Rk2Ntemp3 &
	- ((k2/(ai*Exp(NE+0.5d0*hp)*H(4*J-1)))**2.0d0)*Rk2temp3
	Rk2temp4 = Rk2(J) + (3.0d0/16.0d0)*hp*kRk21  +  (9.0d0/16.0d0)*hp*kRk24  
	Rk2Ntemp4 = Rk2N(J) + (3.0d0/16.0d0)*hp*kRk2N1  +  (9.0d0/16.0d0)*hp*kRk2N4

	kRk25 = Rk2Ntemp4 
	kRk2N5 = -(1.0d0 - e1(4*J) + 2.0d0*zNz(4*J))*Rk2Ntemp4 &
	- ((k2/(ai*Exp(NE+0.75d0*hp)*H(4*J)))**2.0d0)*Rk2temp4
	Rk2temp5 = Rk2(J) - (3.0d0/7.0d0)*hp*kRk21  +  (2.0d0/7.0d0)*hp*kRk22 &
	+  (12.0d0/7.0d0)*hp*kRk23  -  (12.0d0/7.0d0)*hp*kRk24 + (8.0d0/7.0d0)*hp*kRk25
	Rk2Ntemp5 = Rk2N(J) - (3.0d0/7.0d0)*hp*kRk2N1  +  (2.0d0/7.0d0)*hp*kRk2N2 &
	+  (12.0d0/7.0d0)*hp*kRk2N3 - (12.0d0/7.0d0)*hp*kRk2N4 + (8.0d0/7.0d0)*hp*kRk2N5

	kRk26 = Rk2Ntemp5
	kRk2N6 = -(1.0d0 - e1(4*J+1) + 2.0d0*zNz(4*J+1))*Rk2Ntemp5 &
	- ((k2/(ai*Exp(NE+1.0d0*hp)*H(4*J+1)))**2.0d0)*Rk2temp5

	Rk2(J + 1) = Rk2(J) + (hp/90.0d0)*(7.0d0*kRk21 + 32.0d0*kRk23 &
	+ 12.0d0*kRk24 + 32.0d0*kRk25 + 7.0d0*kRk26)
	Rk2N(J + 1) = Rk2N(J) + (hp/90.0d0)*(7.0d0*kRk2N1 + 32.0d0*kRk2N3 &
	+ 12.0d0*kRk2N4 + 32.0d0*kRk2N5 + 7.0d0*kRk2N6)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	khk21 = hk2N(J)
	khk2N1 = -(3.0d0 - e1(4*J-3))*hk2N(J)-((k2/(ai*Exp(NE)*H(4*J-3)))**2.0d0)*hk2(J)
	hk2temp1 = hk2(J) + (1.0d0/4.0d0)*hp*khk21 
	hk2Ntemp1 = hk2N(J) + (1.0d0/4.0d0)*hp*khk2N1 

	khk22 = hk2Ntemp1
	khk2N2 = -(3.0d0 - e1(4*J-2))*hk2Ntemp1 &
	- ((k2/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*hk2temp1
	hk2temp2 =  hk2(J) + (1.0d0/8.0d0)*hp*khk21 + (1.0d0/8.0d0)*hp*khk22
	hk2Ntemp2 = hk2N(J) + (1.0d0/8.0d0)*hp*khk2N1 + (1.0d0/8.0d0)*hp*khk2N2

	khk23 = hk2Ntemp2 
	khk2N3 = -(3.0d0 - e1(4*J-2))*hk2Ntemp2 &
	- ((k2/(ai*Exp(NE+0.25d0*hp)*H(4*J-2)))**2.0d0)*hk2temp2
	hk2temp3 = hk2(J) - (1.0d0/2.0d0)*hp*khk22 + hp*khk23
	hk2Ntemp3 = hk2N(J) - (1.0d0/2.0d0)*hp*khk2N2 + hp*khk2N3

	khk24 = hk2Ntemp3 
	khk2N4 = -(3.0d0 - e1(4*J-1))*hk2Ntemp3 &
	- ((k2/(ai*Exp(NE+0.5d0*hp)*H(4*J-1)))**2.0d0)*hk2temp3
	hk2temp4 = hk2(J) + (3.0d0/16.0d0)*hp*khk21  +  (9.0d0/16.0d0)*hp*khk24  
	hk2Ntemp4 = hk2N(J) + (3.0d0/16.0d0)*hp*khk2N1  +  (9.0d0/16.0d0)*hp*khk2N4

	khk25 = hk2Ntemp4
	khk2N5 = -(3.0d0 - e1(4*J))*hk2Ntemp4 &
	- ((k2/(ai*Exp(NE+0.75d0*hp)*H(4*J)))**2.0d0)*hk2temp4
	hk2temp5 = hk2(J) - (3.0d0/7.0d0)*hp*khk21 +  (2.0d0/7.0d0)*hp*khk22 &
	+  (12.0d0/7.0d0)*hp*khk23  -  (12.0d0/7.0d0)*hp*khk24 + (8.0d0/7.0d0)*hp*khk25
	hk2Ntemp5 = hk2N(J) - (3.0d0/7.0d0)*hp*khk2N1  +  (2.0d0/7.0d0)*hp*khk2N2 &
	+  (12.0d0/7.0d0)*hp*khk2N3 - (12.0d0/7.0d0)*hp*khk2N4 + (8.0d0/7.0d0)*hp*khk2N5

	khk26 = hk2Ntemp5
	khk2N6 = -(3.0d0 - e1(4*J+1))*hk2Ntemp5 - ((k2/(ai*Exp(NE+1.0d0*hp)*H(4*J+1)))**2.0d0)*hk2temp5

	hk2(J + 1) = hk2(J) + (hp/90.0d0)*(7.0d0*khk21 + 32.0d0*khk23 &
	+ 12.0d0*khk24 + 32.0d0*khk25 + 7.0d0*khk26)
	hk2N(J + 1) = hk2N(J) + (hp/90.0d0)*(7.0d0*khk2N1 + 32.0d0*khk2N3 &
	+ 12.0d0*khk2N4 + 32.0d0*khk2N5 + 7.0d0*khk2N6)	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	If (((REAL(J-JNEic))/4.0d0)-AINT((REAL(J-JNEic))/4.0d0).eq. 0.0d0) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ggg
	cGgggI = cGgggI + 14.0d0*(-1.0d0/4.0d0)*CI*(ai*Exp(NE)/H(4*J-3)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*DCONJG(hk1(J)) &
	*DCONJG(hk2(J))*DCONJG(hk3(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGgggIc = cGgggIc + 14.0d0*(1.0d0/4.0d0)*CI*(ai*Exp(NE)/H(4*J-3)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*hk1(J)*hk2(J) &
	*hk3(J)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ssg
	cGssg1I = cGssg1I - 14.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NE)*e1(4*J-3)/H(4*J-3))*DCONJG(hk3(J)) &
	*DCONJG(Rk1(J))*DCONJG(Rk2(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg1Ic = cGssg1Ic +14.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NE)*e1(4*J-3)/H(4*J-3))*hk3(J) &
	*Rk1(J)*Rk2(J)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg2I = cGssg2I + 14.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*e1(4*J-3)*DCONJG(hk3(J))*DCONJG(Rk1N(J))*DCONJG(Rk2N(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg2Ic = cGssg2Ic - 14.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*e1(4*J-3)*hk3(J)*Rk1N(J)*Rk2N(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg3I = cGssg3I + 14.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NE))*H(4*J-3) &
	*e1(4*J-3)*e1(4*J-3)*((k2/k1)*DCONJG(hk3N(J))*DCONJG(Rk1N(J))*DCONJG(Rk2(J)) &
	+(k1/k2)*DCONJG(hk3N(J))*DCONJG(Rk2N(J))*DCONJG(Rk1(J))) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg3Ic = cGssg3Ic - 14.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NE))*H(4*J-3) &
	*e1(4*J-3)*e1(4*J-3)*((k2/k1)*hk3N(J)*Rk1N(J)*Rk2(J) &
	+(k1/k2)*hk3N(J)*Rk2N(J)*Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sgg
	cGsgg1I = cGsgg1I + 14.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*DCONJG(hk2N(J))*DCONJG(hk3N(J))*DCONJG(Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg1Ic = cGsgg1Ic - 14.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*hk2N(J)*hk3N(J)*Rk1(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg2I = cGsgg2I - 14.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NE) &
	*e1(4*J-3)/H(4*J-3))*DCONJG(hk2(J))*DCONJG(hk3(J))*DCONJG(Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg2Ic = cGsgg2Ic + 14.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NE) &
	*e1(4*J-3)/H(4*J-3))*hk2(J)*hk3(J)*Rk1(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg3I = cGsgg3I - 14.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*DCONJG(hk3(J)) &
	*DCONJG(hk2N(J))*DCONJG(Rk1N(J))+((k3*k3-k2*k2-k1*k1)/(k1*k1))*DCONJG(hk2(J)) &
	*DCONJG(hk3N(J))*DCONJG(Rk1N(J)))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg3Ic = cGsgg3Ic + 14.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*hk3(J)*hk2N(J)*Rk1N(J) &
	+((k3*k3-k2*k2-k1*k1)/(k1*k1))*hk2(J)*hk3N(J) &
	*Rk1N(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	elseif(((REAL(J-JNEic))/2.0d0)-AINT((REAL(J-JNEic))/2.0d0).eq. 0.0d0 &
        .and. ((REAL(J-JNEic))/4.0d0)-AINT((REAL(J-JNEic))/4.0d0).ne. 0.0d0) then

	cGgggI = cGgggI + 12.0d0*(-1.0d0/4.0d0)*CI*(ai*Exp(NE)/H(4*J-3)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*DCONJG(hk1(J)) &
	*DCONJG(hk2(J))*DCONJG(hk3(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGgggIc = cGgggIc + 12.0d0*(1.0d0/4.0d0)*CI*(ai*Exp(NE)/H(4*J-3)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*hk1(J)*hk2(J) &
	*hk3(J)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ssg
	cGssg1I = cGssg1I - 12.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NE)*e1(4*J-3)/H(4*J-3))*DCONJG(hk3(J)) &
	*DCONJG(Rk1(J))*DCONJG(Rk2(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg1Ic = cGssg1Ic + 12.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NE)*e1(4*J-3)/H(4*J-3))*hk3(J) &
	*Rk1(J)*Rk2(J)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg2I = cGssg2I + 12.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*e1(4*J-3)*DCONJG(hk3(J))*DCONJG(Rk1N(J))*DCONJG(Rk2N(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg2Ic = cGssg2Ic - 12.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*e1(4*J-3)*hk3(J)*Rk1N(J)*Rk2N(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg3I = cGssg3I + 12.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NE))*H(4*J-3) &
	*e1(4*J-3)*e1(4*J-3)*((k2/k1)*DCONJG(hk3N(J))*DCONJG(Rk1N(J))*DCONJG(Rk2(J)) &
	+(k1/k2)*DCONJG(hk3N(J))*DCONJG(Rk2N(J))*DCONJG(Rk1(J))) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg3Ic = cGssg3Ic - 12.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NE))*H(4*J-3) &
	*e1(4*J-3)*e1(4*J-3)*((k2/k1)*hk3N(J)*Rk1N(J)*Rk2(J) &
	+(k1/k2)*hk3N(J)*Rk2N(J)*Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sgg
	cGsgg1I = cGsgg1I + 12.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*DCONJG(hk2N(J))*DCONJG(hk3N(J))*DCONJG(Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg1Ic = cGsgg1Ic - 12.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*hk2N(J)*hk3N(J)*Rk1(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg2I = cGsgg2I - 12.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NE) &
	*e1(4*J-3)/H(4*J-3))*DCONJG(hk2(J))*DCONJG(hk3(J))*DCONJG(Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg2Ic = cGsgg2Ic + 12.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NE) &
	*e1(4*J-3)/H(4*J-3))*hk2(J)*hk3(J)*Rk1(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg3I = cGsgg3I - 12.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*DCONJG(hk3(J)) &
	*DCONJG(hk2N(J))*DCONJG(Rk1N(J))+((k3*k3-k2*k2-k1*k1)/(k1*k1))*DCONJG(hk2(J)) &
	*DCONJG(hk3N(J))*DCONJG(Rk1N(J)))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg3Ic = cGsgg3Ic + 12.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*hk3(J)*hk2N(J)*Rk1N(J) &
	+((k3*k3-k2*k2-k1*k1)/(k1*k1))*hk2(J)*hk3N(J) &
	*Rk1N(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


	else

	cGgggI = cGgggI + 32.0d0*(-1.0d0/4.0d0)*CI*(ai*Exp(NE)/H(4*J-3)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*DCONJG(hk1(J)) &
	*DCONJG(hk2(J))*DCONJG(hk3(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGgggIc = cGgggIc + 32.0d0*(1.0d0/4.0d0)*CI*(ai*Exp(NE)/H(4*J-3)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*hk1(J)*hk2(J) &
	*hk3(J)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ssg
	cGssg1I = cGssg1I - 32.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NE)*e1(4*J-3)/H(4*J-3))*DCONJG(hk3(J)) &
	*DCONJG(Rk1(J))*DCONJG(Rk2(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg1Ic = cGssg1Ic + 32.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NE)*e1(4*J-3)/H(4*J-3))*hk3(J) &
	*Rk1(J)*Rk2(J)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg2I = cGssg2I + 32.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*e1(4*J-3)*DCONJG(hk3(J))*DCONJG(Rk1N(J))*DCONJG(Rk2N(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg2Ic = cGssg2Ic - 32.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*e1(4*J-3)*hk3(J)*Rk1N(J)*Rk2N(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg3I = cGssg3I + 32.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NE))*H(4*J-3) &
	*e1(4*J-3)*e1(4*J-3)*((k2/k1)*DCONJG(hk3N(J))*DCONJG(Rk1N(J))*DCONJG(Rk2(J)) &
	+(k1/k2)*DCONJG(hk3N(J))*DCONJG(Rk2N(J))*DCONJG(Rk1(J))) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGssg3Ic = cGssg3Ic - 32.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NE))*H(4*J-3) &
	*e1(4*J-3)*e1(4*J-3)*((k2/k1)*hk3N(J)*Rk1N(J)*Rk2(J) &
	+(k1/k2)*hk3N(J)*Rk2N(J)*Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sgg
	cGsgg1I = cGsgg1I + 32.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*DCONJG(hk2N(J))*DCONJG(hk3N(J))*DCONJG(Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg1Ic = cGsgg1Ic - 32.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*hk2N(J)*hk3N(J)*Rk1(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg2I = cGsgg2I - 32.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NE) &
	*e1(4*J-3)/H(4*J-3))*DCONJG(hk2(J))*DCONJG(hk3(J))*DCONJG(Rk1(J)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg2Ic = cGsgg2Ic + 32.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NE) &
	*e1(4*J-3)/H(4*J-3))*hk2(J)*hk3(J)*Rk1(J) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg3I = cGsgg3I - 32.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*DCONJG(hk3(J)) &
	*DCONJG(hk2N(J))*DCONJG(Rk1N(J))+((k3*k3-k2*k2-k1*k1)/(k1*k1))*DCONJG(hk2(J)) &
	*DCONJG(hk3N(J))*DCONJG(Rk1N(J)))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))
	cGsgg3Ic = cGsgg3Ic + 32.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NE)) &
	*H(4*J-3)*e1(4*J-3)*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*hk3(J)*hk2N(J)*Rk1N(J) &
	+((k3*k3-k2*k2-k1*k1)/(k1*k1))*hk2(J)*hk3N(J) &
	*Rk1N(J))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NE)*H(4*J-3)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end if  

	end do    

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGgggI = cGgggI + 7.0d0*(-1.0d0/4.0d0)*CI*(ai*Exp(NEshs)/H(INEshs)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*DCONJG(hk1(JNEshs)) &
	*DCONJG(hk2(JNEshs))*DCONJG(hk3(JNEshs))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))
	cGgggIc = cGgggIc + 7.0d0*(1.0d0/4.0d0)*CI*(ai*Exp(NEshs)/H(INEshs)) &
	*(k1**2.0d0+k2**2.0d0+k3**2.0d0)*hk1(JNEshs)*hk2(JNEshs) &
	*hk3(JNEshs)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!ssg
	cGssg1I = cGssg1I - 7.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NEshs)*e1(INEshs)/H(INEshs))*DCONJG(hk3(JNEshs)) &
	*DCONJG(Rk1(JNEshs))*DCONJG(Rk2(JNEshs))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))
	cGssg1Ic = cGssg1Ic + 7.0d0*2.0d0*CI*(k1*k2)*(ai*Exp(NEshs)*e1(INEshs)/H(INEshs))*hk3(JNEshs) &
	*Rk1(JNEshs)*Rk2(JNEshs)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg2I = cGssg2I + 7.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NEshs)) &
	*H(INEshs)*e1(INEshs)*e1(INEshs) &
	*DCONJG(hk3(JNEshs))*DCONJG(Rk1N(JNEshs))*DCONJG(Rk2N(JNEshs)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))
	cGssg2Ic = cGssg2Ic -7.0d0*(CI/2.0d0)*(k3*k3/(k1*k2))*(ai*ai*ai*Exp(3.0d0*NEshs)) &
	*H(INEshs)*e1(INEshs)*e1(INEshs) &
	*hk3(JNEshs)*Rk1N(JNEshs)*Rk2N(JNEshs) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGssg3I = cGssg3I + 7.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NEshs))*H(INEshs)*e1(INEshs)*e1(INEshs) &
	*((k2/k1)*DCONJG(hk3N(JNEshs))*DCONJG(Rk1N(JNEshs))*DCONJG(Rk2(JNEshs)) &
	+(k1/k2)*DCONJG(hk3N(JNEshs))*DCONJG(Rk2N(JNEshs))*DCONJG(Rk1(JNEshs))) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))
	cGssg3Ic = cGssg3Ic - 7.0d0*(CI/2.0d0)*(ai*ai*ai*Exp(3.0d0*NEshs))*H(INEshs)*e1(INEshs)*e1(INEshs) &
	*((k2/k1)*hk3N(JNEshs)*Rk1N(JNEshs)*Rk2(JNEshs) &
	+(k1/k2)*hk3N(JNEshs)*Rk2N(JNEshs)*Rk1(JNEshs)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sgg
	cGsgg1I = cGsgg1I + 7.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NEshs))*H(INEshs)*e1(INEshs) &
	*DCONJG(hk2N(JNEshs))*DCONJG(hk3N(JNEshs))*DCONJG(Rk1(JNEshs)) &
	*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))
	cGsgg1Ic = cGsgg1Ic - 7.0d0*(CI/4.0d0)*(ai*ai*ai*Exp(3.0d0*NEshs))*H(INEshs)*e1(INEshs)&
	*hk2N(JNEshs)*hk3N(JNEshs)*Rk1(JNEshs)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg2I = cGsgg2I - 7.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NEshs)*e1(INEshs)/H(INEshs)) &
	*DCONJG(hk2(JNEshs))*DCONJG(hk3(JNEshs))*DCONJG(Rk1(JNEshs))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))
	cGsgg2Ic = cGsgg2Ic + 7.0d0*(CI/8.0d0)*(k1*k1-k2*k2-k3*k3)*(ai*Exp(NEshs)*e1(INEshs)/H(INEshs)) &
	*hk2(JNEshs)*hk3(JNEshs)*Rk1(JNEshs)*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	cGsgg3I = cGsgg3I -7.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NEshs))*H(INEshs)*e1(INEshs) &
	*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*DCONJG(hk3(JNEshs))*DCONJG(hk2N(JNEshs)) &
	*DCONJG(Rk1N(JNEshs))+((k3*k3-k2*k2-k1*k1)/(k1*k1))*DCONJG(hk2(JNEshs)) &
	*DCONJG(hk3N(JNEshs))*DCONJG(Rk1N(JNEshs)))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))
	cGsgg3Ic = cGsgg3Ic + 7.0d0*(CI/8.0d0)*(ai*ai*ai*Exp(3.0d0*NEshs))*H(INEshs)*e1(INEshs) &
	*(((k2*k2-k3*k3-k1*k1)/(k1*k1))*hk3(JNEshs)*hk2N(JNEshs)*Rk1N(JNEshs) &
	+((k3*k3-k2*k2-k1*k1)/(k1*k1))*hk2(JNEshs)*hk3N(JNEshs) &
	*Rk1N(JNEshs))*Exp(-kappa*(k1+k2+k3)/(3.0d0*ai*Exp(NEshs)*H(INEshs)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gggg = (2.0d0*hp/45.0d0)*(hk1(JNEshs)*hk2(JNEshs)*hk3(JNEshs)*cGgggI &
	+ DCONJG(hk1(JNEshs))*DCONJG(hk2(JNEshs))*DCONJG(hk3(JNEshs))*cGgggIc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gssg1 = (2.0d0*hp/45.0d0)*(hk3(JNEshs)*Rk1(JNEshs)*Rk2(JNEshs)*cGssg1I &
	+ DCONJG(hk3(JNEshs))*DCONJG(Rk1(JNEshs))*DCONJG(Rk2(JNEshs))*cGssg1Ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gssg2 = (2.0d0*hp/45.0d0)*(hk3(JNEshs)*Rk1(JNEshs)*Rk2(JNEshs)*cGssg2I &
	+ DCONJG(hk3(JNEshs))*DCONJG(Rk1(JNEshs))*DCONJG(Rk2(JNEshs))*cGssg2Ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gssg3 = (2.0d0*hp/45.0d0)*(hk3(JNEshs)*Rk1(JNEshs)*Rk2(JNEshs)*cGssg3I &
	+ DCONJG(hk3(JNEshs))*DCONJG(Rk2(JNEshs))*DCONJG(Rk1(JNEshs))*cGssg3Ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gssg = Gssg1 + Gssg2 + Gssg3 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gsgg1 = (2.0d0*hp/45.0d0)*(hk2(JNEshs)*hk3(JNEshs)*Rk1(JNEshs)*cGsgg1I &
	+ DCONJG(hk2(JNEshs))*DCONJG(hk3(JNEshs))*DCONJG(Rk1(JNEshs))*cGsgg1Ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gsgg2 = (2.0d0*hp/45.0d0)*(hk2(JNEshs)*hk3(JNEshs)*Rk1(JNEshs)*cGsgg2I &
	+ DCONJG(hk2(JNEshs))*DCONJG(hk3(JNEshs))*DCONJG(Rk1(JNEshs))*cGsgg2Ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gsgg3 = (2.0d0*hp/45.0d0)*(hk2(JNEshs)*hk3(JNEshs)*Rk1(JNEshs)*cGsgg3I &
	+ DCONJG(hk2(JNEshs))*DCONJG(hk3(JNEshs))*DCONJG(Rk1(JNEshs))*cGsgg3Ic)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	Gsgg = Gsgg1 + Gsgg2 + Gsgg3 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	PSk2=((k2**3.0d0)/(2.0d0*(Pi**2.0d0)))*((Rk2(JNEshs))*DCONJG(Rk2(JNEshs)))
	
	PTk2=4.0d0*((k2**3.0d0)/(2.0d0*(Pi**2.0d0)))*((hk2(JNEshs))*DCONJG(hk2(JNEshs)))

	write(52,*) k2, PSk2, PTk2

	hnl= -((4.0d0/(2.0d0*Pi*Pi))**2.0d0)*(((k1*k2*k3)**3.0d0)*Gggg/(2.0d0*((k1**3.0d0) &
	*PTk2*PTk3+(k2**3.0d0)*PTk1*PTk3+(k3**3.0d0)*PTk1*PTk2)))

	cnls= -(4.0d0/((2.0d0*Pi*Pi)**2.0d0))*(((k1*k2*k3)**3.0d0)*Gssg/((k1**3.0d0) &
	*PSk2*PTk3+(k2**3.0d0)*PSk1*PTk3))

	cnlg= -(4.0d0/((2.0d0*Pi*Pi)**2.0d0))*(((k1*k2*k3)**3.0d0)*Gsgg/((k2**3.0d0) &
	*PSk1*PTk3+(k3**3.0d0)*PSk1*PTk2))


	write(351,*) k3/k1, k2/k1,DBLE(hnl)
	write(352,*) k3/k1, k2/k1, DBLE(cnls)
	write(353,*) k3/k1, k2/k1, DBLE(cnlg)

	end do !k2 loop

	end do	!k3 loop
!	print *,"precision of dp =",precision(Abs(2.0d0))
!	print*,"range of dp =", range(Abs(2.0d0))

	END PROGRAM powerspectra
