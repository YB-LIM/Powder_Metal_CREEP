C#################################################################
C                 
C     Duva and Crow, 1994 creep model implemented with CREEP subroutine
C
C
C     Author: Youngbin LIM
C     Contact: lyb0684@naver.com
C
C################################################################
C
C     STATEV(1) : True (logarithmic) strain in 1-direction
C     STATEV(2) : True (logarithmic) strain in 2-direction
C     STATEV(3) : True (logarithmic) strain in 3-direction
C     STATEV(4) : Volumetric strain (Sum of STATEV 1~3)
C     STATEV(5) : Relative density     
C     DECRA(1)  : Creep strain increment 
C     DESWA(1)  : Swelling strain increment
C
      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(1),STATEV(5),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C
C    Get strain in 1,2,3 direction to calculated volumetric strain
C
      CALL GETVRM('LE',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,MATLAYO,
     1 LACCFLA)
C
      STATEV(1) = ARRAY(1)
      STATEV(2) = ARRAY(2)
      STATEV(3) = ARRAY(3)	  
C
      RETURN    
      END
C
C       
      SUBROUTINE CREEP(DECRA,DESWA,STATEV,SERD,EC,ESW,P,QTILD,
     1 TEMP,DTEMP,PREDEF,DPRED,TIME,DTIME,CMNAME,LEXIMP,LEND,
     2 COORDS,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION DECRA(5),DESWA(5),STATEV(5),PREDEF(1),DPRED(1),
     1 TIME(3),COORDS(3),EC(2),ESW(2)
C
C
      LEXIMP = 0        ! Flag for explicit time integration
      LEND = 0          ! Calculation is performed at the start of the increment
      D_ZERO = 0.78     ! Initial relative density
      AN = 1.53         ! Constant n
      A0 = 95.67        ! Constant A0
      RR = 8.31         ! Universal Gas constant in J/(K*mol)
      ACTE = 150000.0   ! Activation energy in J/mol
C        
      SMEAN = -1.0*P    ! Hydrostatic (mean) stress
C                    
      STATEV(4) = STATEV(1) + STATEV(2) + STATEV(3)   ! Volumetric strain
C
      VOL = STATEV(4)           
      STATEV(5) = D_ZERO*EXP(-VOL)  ! Calculate relative density from volumetric strain
      DENS = STATEV(5)  ! DENS = Relative density
C
      SMALLA = (1.0+(2.0/3.0)*(1.0-DENS))/(DENS**(2.0*AN/(AN+1.0)))
C
      B1 = AN*(1.0-DENS)
      B2 = (1.0-(1.0-DENS)**(1.0/AN))**AN
      B3 = 2.0/(AN+1.0)
      B4 = 3.0/(2.0*AN)
C 
      SMALLB = ((B1/B2)**(B3))*(B4**2.0)
C
      SS = SMALLA*QTILD**(2.0) + SMALLB*SMEAN**(2.0) ! Square of S, i.e., S^2
C
      AA = A0*EXP(-ACTE/(RR*TEMP))
C
      CEEQ_ABS = ABS(AA*SMALLA*QTILD*(SS**((AN-1.0)/2))*DTIME) ! Calculate Absolute creep strain increment
C
C     Handle division by zero for creep strain parts
C
      IF (CEEQ_ABS > 1.0E-10) THEN
         DECRA(1) = AA*SMALLA*QTILD*(SS**((AN-1.0)/2))*DTIME
         DECRA(4) = -DECRA(1)*((AN-1.0)*SMALLB*SMEAN/SS)
         DECRA(5) = DECRA(1)*(1.0/QTILD+(AN-1.0)*SMALLA*QTILD/SS)		 
      ELSE
         DECRA(1) = 0.0
         DECRA(4) = 0.0
         DECRA(5) = 0.0
      END IF
C
      CESW_ABS = ABS(AA*SMALLB*SMEAN*(SS**((AN-1.0)/2))*DTIME) ! Calculate Absolute swelling strain increment	     
C
C     Handle division by zero for swelling strain parts
C
      IF (CESW_ABS > 1.0E-10) THEN
         DESWA(1) = AA*SMALLB*SMEAN*(SS**((AN-1.0)/2))*DTIME
         DESWA(4) = -DESWA(1)*(1.0/SMEAN+(AN-1.0)*SMALLB*SMEAN/SS)
         DESWA(5) = DESWA(1)*((AN-1.0)*SMALLA*QTILD/SS)		 
      ELSE
         DESWA(1) = 0.0
         DESWA(4) = 0.0
         DESWA(5) = 0.0
      END IF
C
      write (7,*) 'DECRA(1) = ',DECRA(1)
      write (7,*) 'DECRA(4) = ',DECRA(4)
      write (7,*) 'DECRA(5) = ',DECRA(5)
      write (7,*) 'DESWA(1) = ',DESWA(1)
      write (7,*) 'DESWA(4) = ',DESWA(4)
      write (7,*) 'DESWA(5) = ',DESWA(5)	  
C	  
      RETURN
      END