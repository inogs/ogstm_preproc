      SUBROUTINE comp_ssh16_relax()
!CC---------------------------------------------------------------------
!CC
!CC                       ROUTINE div
!CC                     ***************
!CC
!CC  Purpose :
!CC  --------
!CC	compute the now horizontal divergence of the velocity field.
!CC
!C   Method :
!C   -------
!C	The now divergence is given by :
!C       * s-coordinate ('key_s_coord' defined)
!C         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
!C       * z-coordinate (default key)
!C         hdivn = 1/(e1t*e2t) [ di(e2u  un) + dj(e1v  vn) ]
!C
!C      Apply lateral boundary conditions on hdivn through a call to
!C      routine mpplnk ('key_mpp' defined) or lbc.
!C
!C      macro-tasked on horizontal slab (jk-loop :  1, jpk-1)
!C
!C   Input :
!C   ------
!C      argument
!C              ktask           : task identificator
!C              kt              : time step
!C      common
!C            /comcoh/   	: scale factors
!C            /comtsk/   	: multitasking
!C            /comnow/		: present fields (now)
!C
!C   Output :
!C   -------
!C      common
!C	      /comnow/ hdivn	: now horizontal divergence
!C
!C   External : mpplnk or lbc
!C   ---------
!C
!C   Modifications :
!C   --------------
!C      original :  87-06 (P. Andrich, D. L Hostis)
!C      additions : 91-11 (G. Madec)
!C                : 93-03 (M. Guyon) symetrical conditions
!C                : 96-01 (G. Madec) s-coordinates
!C                : 97-06 (G. Madec) lateral boundary cond., lbc
!C----------------------------------------------------------------------
!C parameters and commons
!C ======================
 USE modulo16
      IMPLICIT NONE

!C----------------------------------------------------------------------
!C local declarations
      INTEGER ji16, jj16, jk,i
      INTEGER jpi16m1,jpj16m1,jpk16m1
      INTEGER iter
      REAL(8) Tau,Tau_v
!
      REAL(8) zwu16, zwv16, zww16, zww16up, inv_fact
!
      jpi16m1=jpi16-1
      jpj16m1=jpj16-1
!
      un16(jpi16,:,1)=0; un16(:,jpj16,1)=0;
      vn16(jpi16,:,1)=0; vn16(:,jpj16,1)=0;
      wn16(jpi16,:,1)=0; wn16(:,jpj16,1)=0;
!

            e3t16bb=e3t16
            e3u16bb=e3u16
            e3v16bb=e3v16

            e3t16b=e3t16
            e3u16b=e3u16
            e3v16b=e3v16

            un16bb =un16
            vn16bb =vn16
            
            un16b =un16
            vn16b =vn16

        DO jj16 = 1,jpj16
          DO ji16 = 1,jpi16
               ssh16t_R8(ji16,jj16) = REAL(e3t16(ji16,jj16,1),8)
               ssh16u_R8(ji16,jj16) = REAL(e3u16(ji16,jj16,1),8)
               ssh16v_R8(ji16,jj16) = REAL(e3v16(ji16,jj16,1),8)
          END DO  
        END DO  

        Tau = 1
        Tau_v = 0.1

        DO iter =1,1000


            e3t16bb=e3t16b
            e3u16bb=e3u16b
            e3v16bb=e3v16b

            e3t16b=e3t16
            e3u16b=e3u16
            e3v16b=e3v16

            un16bb =un16b
            vn16bb =vn16b

            un16b =un16
            vn16b =vn16

            DO jj16 = 2,jpj16m1
               DO ji16 = 2,jpi16m1

                  wn16(ji16,jj16,1) = wn16(ji16,jj16,2)  &
                        - e3t16(ji16,jj16,1)*hdivn16(ji16,jj16,1)

                  e3u16(ji16,jj16,1)   = 0.75*e3u16b(ji16,jj16,1)   + 0.25*e3u16bb(ji16,jj16,1)   + umask16(ji16,jj16,1)  *(+1.)*Tau*sign(1.,un16b(ji16,jj16,1)  * wn16(ji16,jj16,1))
                  e3v16(ji16,jj16,1)   = 0.75*e3v16b(ji16,jj16,1)   + 0.25*e3v16bb(ji16,jj16,1)   + vmask16(ji16,jj16,1)  *(+1.)*Tau*sign(1.,vn16b(ji16,jj16,1)  * wn16(ji16,jj16,1))
 !                e3u16(ji16-1,jj16,1) = 0.75*e3u16b(ji16-1,jj16,1) + 0.25*e3u16bb(ji16-1,jj16,1) + umask16(ji16-1,jj16,1)*(-1.)*Tau*sign(1.,un16b(ji16-1,jj16,1)* wn16(ji16,jj16,1))
 !                e3v16(ji16,jj16-1,1) = 0.75*e3v16b(ji16,jj16-1,1) + 0.25*e3v16bb(ji16,jj16-1,1) + vmask16(ji16,jj16-1,1)*(-1.)*Tau*sign(1.,vn16b(ji16,jj16-1,1)* wn16(ji16,jj16,1))

                  e3t16(ji16,jj16,1)  = 0.25*(e3u16(ji16,jj16,1)+e3u16(ji16+1,jj16,1) + e3v16(ji16,jj16,1)+e3v16(ji16,jj16+1,1))

                  un16(ji16,jj16,1)   = 0.75*un16b(ji16,jj16,1)     + 0.25*un16bb(ji16,jj16,1)   + umask16(ji16,jj16,1)  *(+1.)*Tau_v*sign(1.,un16b(ji16,jj16,1)  * wn16(ji16,jj16,1))
                  vn16(ji16,jj16,1)   = 0.75*vn16b(ji16,jj16,1)     + 0.25*vn16bb(ji16,jj16,1)   + vmask16(ji16,jj16,1)  *(+1.)*Tau_v*sign(1.,vn16b(ji16,jj16,1)  * wn16(ji16,jj16,1))
 !                un16(ji16-1,jj16,1) = 0.75*un16b(ji16-1,jj16,1)   + 0.25*un16bb(ji16-1,jj16,1) + umask16(ji16-1,jj16,1)*(-1.)*Tau_v*sign(1.,un16b(ji16-1,jj16,1)* wn16(ji16,jj16,1))
 !                vn16(ji16,jj16-1,1) = 0.75*vn16b(ji16,jj16-1,1)   + 0.25*vn16bb(ji16,jj16-1,1) + vmask16(ji16,jj16-1,1)*(-1.)*Tau_v*sign(1.,vn16b(ji16,jj16-1,1)* wn16(ji16,jj16,1))

               END DO  
            END DO  
            write(*,*) '-----------------'
            write(*,*) 'iter ', iter
            write(*,*) 'wn16 ', SUM(wn16)
            write(*,*) '-----------------'


            call div16(1)

        END DO  

        DO jj16 = 1,jpj16
            DO ji16 = 1,jpi16
               ssh16u_R8(ji16,jj16) = REAL(e3u16(ji16,jj16,1),8)
               ssh16v_R8(ji16,jj16) = REAL(e3v16(ji16,jj16,1),8)
               e3t16(ji16,jj16,1)   = 0.25*(e3u16(ji16,jj16,1)+e3u16(ji16+1,jj16,1) + e3v16(ji16,jj16,1)+e3v16(ji16,jj16+1,1))
               ssht16(ji16,jj16)    = e3t16(ji16,jj16,1)
               sshu16(ji16,jj16)    = e3u16(ji16,jj16,1)
               sshv16(ji16,jj16)    = e3v16(ji16,jj16,1)
            END DO  
        END DO  
!
!
!
      RETURN
      END
