!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TVDLB( BCX,CX,CRLB,DZX,FCL,FCLP,FLB,SLX,WCLZ,YLX,
     &  INDX,K,KLX,N,NBT,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Boundary condition subroutine for TVD, Superbee and Upwind
!     solute transport for lower surfaces (i.e., bottom, south, west ).
!
!     INDX =  0:  Dirichlet
!     INDX =  1:  Outflow
!     INDX = -1:  Inflow
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 17, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE SOLTN
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CX(LFDC,LSOLU),DZX(*),SLX(LSV,LFD),YLX(LFDC,LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TVDLB'
      WCLZ = 0.D+0
      IF( INDEX.LT.0 ) THEN
        IF( FLB.GE.ZERO ) WCLZ = BCX*FCL*FLB
      ELSEIF( INDEX.EQ.0 ) THEN
        IF( FLB.GE.ZERO ) THEN
          WCLZ = BCX*FCL*FLB
        ELSEIF( FLB.LT.ZERO .AND. K.LT.KLX ) THEN
          FCLT = YLX(NBT,NSL)/(SLX(2,NBT)*PORD(2,NBT)+SMALL)
          R = ((CX(N,NSL)*FCLP-CX(NBT,NSL)*FCLT)
     &      /(BCX*FCL-CX(N,NSL)*FCLP+SMALL))
     &      *(DZX(K)/(DZX(K)+DZX(K+1)))
          THETA = FLIMIT( R,CRLB,ISLC(1) )
          WCLZ = BCX*FLB*THETA*FCL
     &      + CX(N,NSL)*FLB*(1.D+0-THETA)*FCLP
         ELSEIF( FLB.LT.ZERO .AND. K.EQ.KLX ) THEN
          WCLZ = CX(N,NSL)*FLB*FCLP
         ENDIF
      ELSEIF( INDEX.GT.0 ) THEN
        IF( FLB.LT.ZERO .AND. K.LT.KLX ) THEN
          FCLT = YLX(NBT,NSL)/(SLX(2,NBT)*PORD(2,NBT)+SMALL)
          R = ((CX(N,NSL)*FCLP-CX(NBT,NSL)*FCLT)
     &      /(BCX*FCL-CX(N,NSL)*FCLP+SMALL))
     &      *(DZX(K)/(DZX(K)+DZX(K+1)))
          THETA = FLIMIT( R,CRLB,ISLC(1) )
          WCLZ = BCX*FLB*THETA*FCL
     &      + CX(N,NSL)*FLB*(1.D+0-THETA)*FCLP
        ELSEIF( FLB.LT.ZERO .AND. K.EQ.KLX ) THEN
          WCLZ = CX(N,NSL)*FLB*FCLP
        ENDIF
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TVDLB group  ---
!
      RETURN
      END

!----------------------Subroutine--------------------------------------!
!
      SUBROUTINE TVDUB( BCX,CX,CRLT,DZX,FCL,FCLP,FLT,SLX,WCLZ,YLX,
     &  INDX,K,KLX,N,NBB,NSL )
!
!-------------------------Disclaimer-----------------------------------!
!
!     This material was prepared as an account of work sponsored by
!     an agency of the United States Government. Neither the
!     United States Government nor the United States Department of
!     Energy, nor Battelle, nor any of their employees, makes any
!     warranty, express or implied, or assumes any legal liability or
!     responsibility for the accuracy, completeness, or usefulness
!     of any information, apparatus, product, software or process
!     disclosed, or represents that its use would not infringe
!     privately owned rights.
!
!----------------------Acknowledgement---------------------------------!
!
!     This software and its documentation were produced with Government
!     support under Contract Number DE-AC06-76RLO-1830 awarded by the
!     United Department of Energy. The Government retains a paid-up
!     non-exclusive, irrevocable worldwide license to reproduce,
!     prepare derivative works, perform publicly and display publicly
!     by or for the Government, including the right to distribute to
!     other Government contractors.
!
!---------------------Copyright Notices--------------------------------!
!
!            Copyright Battelle Memorial Institute, 1996
!                    All Rights Reserved.
!
!----------------------Description-------------------------------------!
!
!     Boundary condition subroutine for TVD, Superbee and Upwind
!     solute transport for upper surfaces (i.e., top, north, east ).
!
!     INDX =  0:  Dirichlet
!     INDX =  1:  Outflow
!     INDX = -1:  Inflow
!
!----------------------Authors-----------------------------------------!
!
!     Written by MD White, Battelle, December, 1999.
!     Last Modified by Mark White, Battelle, December 17, 1999.
!
!----------------------Fortran 90 Modules------------------------------!
!
      USE GLB_PAR
      USE TRNSPT
      USE SOLTN
      USE FDVP
      USE CONST
!
!----------------------Implicit Double Precision-----------------------!
!
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!
!----------------------Type Declarations-------------------------------!
!
      REAL*8 CX(LFDC,LSOLU),DZX(*),SLX(LSV,LFD),YLX(LFDC,LSOLU)
!
!----------------------Executable Lines--------------------------------!
!
      ISUB_LOG = ISUB_LOG+1
      SUB_LOG(ISUB_LOG) = '/TVDUB'
      WCLZ = 0.D+0
      IF( INDX.LT.0 ) THEN
        IF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
          FCLB = YLX(NBB,NSL)/(SLX(2,NBB)*PORD(2,NBB)+SMALL)
          R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &      /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &      *(DZX(K)/(DZX(K)+DZX(K-1)))
          THETA = FLIMIT( R,CRLT,ISLC(1) )
          WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &      + BCX*FLT*THETA*FCL
        ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
          WCLZ =  C(N,NSL)*FLT*FCLP
        ENDIF
      ELSEIF( INDX.EQ.0 ) THEN
        IF( FLT.LT.ZERO ) THEN
          WCLZ = BCX*FCL*FLT
        ELSEIF( FLT.GE.ZERO .AND. K.GT.1 ) THEN
          FCLB = YLX(NBB,NSL)/(SLX(2,NBB)*PORD(2,NBB)+SMALL)
          R = ((C(N,NSL)*FCLP-C(NBB,NSL)*FCLB)
     &      /(BCX*FCL-C(N,NSL)*FCLP+SMALL))
     &      *(DZX(K)/(DZX(K)+DZX(K-1)))
          THETA = FLIMIT( R,CRLT,ISLC(1) )
          WCLZ =  C(N,NSL)*FLT*(1.D+0-THETA)*FCLP
     &      + BCX*FLT*THETA*FCL
        ELSEIF( FLT.GE.ZERO .AND. K.EQ.1 ) THEN
          WCLZ =  C(N,NSL)*FLT*FCLP
        ENDIF
      ELSEIF( INDX.GT.0 ) THEN
        IF( FLT.LT.ZERO ) WCLZ = BCX*FCL*FLT
      ENDIF
      ISUB_LOG = ISUB_LOG-1
!
!---  End of TVDUB group  ---
!
      RETURN
      END

