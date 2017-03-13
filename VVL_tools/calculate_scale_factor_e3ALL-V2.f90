! Fortran version of calculate_scale_factor_e3ALL-V2.py
! before to introduce it in ogstm/src/IO/forcing_phys.f90
program SCALEFACTORS

USE netcdf
IMPLICIT NONE

character*1024 maskfile, forcingT, outfilename
integer jpi,jpj,jpk,jpim1, jpjm1
integer ji,jj,jk
logical IS_INGV_E3T
double precision somma_e3t_0, s0,s1,s2

real(8), allocatable, dimension(:,:,:) :: tmask, umask, vmask
real(8), allocatable, dimension(:,:,:) :: e3t_0,e3t, diff_e3t, e3u, e3v, e3w, e3u_0, e3v_0, e3w_0
real(8), allocatable, dimension(:,:)   :: ssh, e1u,e2u,e1v,e2v,e1t,e2t
real(8), allocatable, dimension(:,:)   :: e1u_x_e2u, e1v_x_e2v, e1t_x_e2t

maskfile    = 'mesh_mask_cut.nc'
forcingT    = 'T20140501-12:00:00.nc'
outfilename = 'E20140501.nc'

call getDIMENSION(maskfile,'x',jpi)
call getDIMENSION(maskfile,'y',jpj)
call getDIMENSION(maskfile,'z',jpk)

jpim1=jpi-1
jpjm1=jpj-1

allocate(tmask   (jpk,jpj,jpi))
allocate(umask   (jpk,jpj,jpi))
allocate(vmask   (jpk,jpj,jpi))
allocate(e3t_0   (jpk,jpj,jpi))
allocate(e3u_0   (jpk,jpj,jpi))
allocate(e3v_0   (jpk,jpj,jpi))
allocate(e3w_0   (jpk,jpj,jpi))
allocate(e3t     (jpk,jpj,jpi))
allocate(e3u     (jpk,jpj,jpi))
allocate(e3v     (jpk,jpj,jpi))
allocate(e3w     (jpk,jpj,jpi))
allocate(diff_e3t(jpk,jpj,jpi))



allocate(e1u(jpj,jpi))
allocate(e1v(jpj,jpi))
allocate(e1t(jpj,jpi))
allocate(e2u(jpj,jpi))
allocate(e2v(jpj,jpi))
allocate(e2t(jpj,jpi))
allocate(ssh(jpj,jpi))
allocate(e1u_x_e2u(jpj,jpi))
allocate(e1v_x_e2v(jpj,jpi))
allocate(e1t_x_e2t(jpj,jpi))

write(*,*) 'Start reading'


call readNetCDF_3dvar(maskfile,'tmask',jpi,jpj,jpk, tmask)
call readNetCDF_3dvar(maskfile,'umask',jpi,jpj,jpk, umask)
call readNetCDF_3dvar(maskfile,'vmask',jpi,jpj,jpk, vmask)
call readNetCDF_3dvar(maskfile,'e3t_0',jpi,jpj,jpk, e3t_0)
call readNetCDF_3dvar(maskfile,'e3u_0',jpi,jpj,jpk, e3u_0)
call readNetCDF_3dvar(maskfile,'e3v_0',jpi,jpj,jpk, e3v_0)
call readNetCDF_3dvar(maskfile,'e3w_0',jpi,jpj,jpk, e3w_0)

write(*,*) 'End reading 3D'

call readNetCDF_2dvar(maskfile,'e1u',jpi,jpj, e1u)
call readNetCDF_2dvar(maskfile,'e2u',jpi,jpj, e2u)

call readNetCDF_2dvar(maskfile,'e1v',jpi,jpj, e1v)
call readNetCDF_2dvar(maskfile,'e2v',jpi,jpj, e2v)

call readNetCDF_2dvar(maskfile,'e1t',jpi,jpj, e1t)
call readNetCDF_2dvar(maskfile,'e2t',jpi,jpj, e2t)

call readNetCDF_2dvar(forcingT,'sossheig',jpi,jpj, ssh)


write(*,*) 'Start calculation'

ssh = ssh*tmask(1,:,:)

IS_INGV_E3T = .false.
     if (.not.IS_INGV_E3T) then



          e3t = e3t_0

          DO ji= 1,jpi
          DO jj= 1,jpj
          somma_e3t_0 =0.0
          if (tmask(1,jj,ji).eq.1) then
              DO jk=1,jpk
                   somma_e3t_0 = somma_e3t_0 + e3t_0(jk,jj,ji)
                   e3t(jk,jj,ji)  = e3t_0(jk,jj,ji) * (1.0 + ssh(jj,ji)/somma_e3t_0)
              ENDDO
          endif
          ENDDO
          ENDDO

      write(*,*) 'step 1'
         do ji=1,jpi
         do jj=1,jpj
            e1u_x_e2u(jj,ji) = e1u(jj,ji)*e2u(jj,ji)
            e1v_x_e2v(jj,ji) = e1v(jj,ji)*e2v(jj,ji)
            e1t_x_e2t(jj,ji) = e1t(jj,ji)*e2t(jj,ji)
         enddo
         enddo


         diff_e3t = e3t - e3t_0

         e3u = 0.0
         e3v = 0.0


         DO ji = 1,jpim1
         DO jj = 1,jpjm1
         DO jk = 1,jpk
             s0= e1t_x_e2t(jj,ji )  * diff_e3t(jk,jj,ji)
             s1= e1t_x_e2t(jj,ji+1) * diff_e3t(jk,jj,ji+1)
             s2= e1t_x_e2t(jj+1,ji) * diff_e3t(jk,jj+1,ji)
             e3u(jk,jj,ji) = 0.5*(umask(jk,jj,ji)/(e1u_x_e2u(jj,ji)) * (s0 + s1))
             e3v(jk,jj,ji) = 0.5*(vmask(jk,jj,ji)/(e1v_x_e2v(jj,ji)) * (s0 + s2))
         ENDDO
         ENDDO
         ENDDO


         DO ji = 1,jpi
         DO jj = 1,jpj
         DO jk = 1,jpk
             e3u(jk,jj,ji) = e3u_0(jk,jj,ji) + e3u(jk,jj,ji)
             e3v(jk,jj,ji) = e3v_0(jk,jj,ji) + e3v(jk,jj,ji)
         ENDDO
         ENDDO
         ENDDO

    write(*,*) 'step 4'

         DO ji = 1,jpi
         DO jj = 1,jpj
             e3w(1,jj,ji) = e3w_0(1,jj,ji) + diff_e3t(1,jj,ji)
         ENDDO
         ENDDO

         DO ji = 1,jpi
         DO jj = 1,jpj
         DO jk = 2,jpk
             !e3w(jk,jj,ji) = e3w_0(jk,jj,ji) + (1.0 - 0.5*tmask(jk,jj,ji)) * (diff_e3t(jk-1,jj,ji) + 0.5*tmask(jk,jj,ji)*diff_e3t(jk,jj,ji)
             if (tmask(jk,jj,ji).eq.1) then
                  e3w(jk,jj,ji) = e3w_0(jk,jj,ji) + 0.5*( diff_e3t(jk-1,jj,ji) + diff_e3t(jk,jj,ji))
             else
                 e3w(jk,jj,ji) = e3w_0(jk,jj,ji) + diff_e3t(jk-1,jj,ji)
             endif

         ENDDO

         ENDDO
         ENDDO

     endif ! IS_INGV_E3T



call PREPARE_BOX_FILE_T(outfilename)

call MODIFY_NC_3D(outfilename,'e3u',jpi,jpj,jpk,e3u)
call MODIFY_NC_3D(outfilename,'e3v',jpi,jpj,jpk,e3v)
call MODIFY_NC_3D(outfilename,'e3w',jpi,jpj,jpk,e3w)
call MODIFY_NC_3D(outfilename,'e3t',jpi,jpj,jpk,e3t)


CONTAINS

!***************************************************************************
SUBROUTINE getDIMENSION(fileNetCDF,stringname,n)
use netcdf
implicit none
character  fileNetCDF*(*)
integer n
integer DIMid,ncid,stat
character  dim_name*30, stringname*(*)
integer counter

counter = 0
stat = nf90_open(fileNetCDF, nf90_nowrite, ncid)        ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_inq_dimid (ncid, stringname, DIMid)         ; call handle_err1(stat,counter,fileNetCDF)
stat = nf90_Inquire_Dimension (ncid, DIMid, dim_name, n); call handle_err1(stat,counter,fileNetCDF)
stat = nf90_close(ncid)                                 ; call handle_err1(stat,counter,fileNetCDF)
END SUBROUTINE getDIMENSION

!****************************************************************************
SUBROUTINE readNetCDF_3dvar(fileNetCDF,varname,jpi,jpj,jpk,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*) ,varname*(*)
integer ncid, stat, VARid
integer jpi,jpj,jpk
real(8) MATRIX(jpk,jpj,jpi)
integer counter, i,j,k
real copy_in(jpi,jpj,jpk)


counter = 0

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid) ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_inq_varid (ncid, varname, VARid)     ;
call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_get_var (ncid,VARid,copy_in)
call handle_err2(stat, fileNetCDF,varname)       ; call handle_err1(stat,counter,fileNetCDF )
stat = nf90_close(ncid)                          ; call handle_err1(stat,counter,fileNetCDF )


      DO i=1,jpi
        DO j=1,jpj
          DO k=1,jpk
            MATRIX(k,j,i) = real(copy_in(i,j,k),8)
          ENDDO
        ENDDO
      ENDDO

end SUBROUTINE readNetCDF_3dvar

!****************************************************************************

SUBROUTINE readNetCDF_2dvar(fileNetCDF,varname,jpi,jpj,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer mycount
integer ncid, stat, VARid
integer jpi,jpj, i,j
real(8) MATRIX(jpj,jpi)

real copy_in(jpi,jpj)

stat = nf90_open(fileNetCDF, nf90_nowrite, ncid); call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, varname, VARid)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_get_var (ncid,VARid,copy_in)
call handle_err2(stat, fileNetCDF,varname)      ; call handle_err1(stat,mycount,fileNetCDF )

stat = nf90_close(ncid)                         ; call handle_err1(stat,mycount,fileNetCDF )


      DO i=1,jpi
        DO j=1,jpj
            MATRIX(j,i) = real(copy_in(i,j),8)
        ENDDO
      ENDDO


END SUBROUTINE readNetCDF_2dvar


!****************************************************************************
subroutine handle_err1(status,mycount, FileNetCDF)
USE netcdf
integer status,mycount
character fileNetCDF*(*)

mycount =mycount+1
if(status .ne. nf90_NoErr)  then
   write(*,*) 'netcdf call',mycount,'with status = ',status
   write(*,*)  'file :', fileNetCDF
   write(*,*) nf90_strerror(status)
   write(*,*) 'Stopped'

endif
end


        subroutine handle_err2(status,fileNetCDF,varname)
        USE netcdf
        integer status
        character fileNetCDF*(*) ,varname*(*)
        if(status .ne. nf90_NoErr)  then
           write(*,*) 'ERROR in Var = ', varname, ' file :', fileNetCDF
        endif

        end subroutine handle_err2


    LOGICAL FUNCTION PREPARE_BOX_FILE_T(fileNetCDF)
        use netcdf
        implicit none

        character fileNetCDF*(*)
        integer s, nc,G, counter
        integer timid, depid, yid, xid
        integer idvartime, idgdept, idphit, idlamt, idT, idS, idR, idW,idH,idE,ide3t
        integer idHt, idHu, idHv

        G     = nf90_global

        s = nf90_create(fileNetCDF, NF90_CLOBBER, nc)
        ! *********** GLOBAL ********************
        s = nf90_put_att(nc, G, 'Convenctions'     , 'OPA')

        ! *********** DIMENSIONS ****************
        s= nf90_def_dim(nc,'x'           , im,  xid)
        s= nf90_def_dim(nc,'y'           , jm,  yid)
        s= nf90_def_dim(nc,'deptht'      , km,depid)
        s= nf90_def_dim(nc,'time_counter', NF90_UNLIMITED,timid)


        ! ********** VARIABLES *****************
        s = nf90_def_var(nc,'time_counter',  nf90_double,(/timid/),           idvartime)
        s = nf90_def_var(nc,'deptht',        nf90_float, (/depid/),             idgdept)
        s = nf90_def_var(nc,'nav_lat',       nf90_float, (/xid,yid/),            idphit)
        s = nf90_def_var(nc,'nav_lon',       nf90_float, (/xid,yid/),            idlamt)
        s = nf90_def_var(nc,'e3u' ,          nf90_float, (/xid,yid,depid,timid/),   ide3u)
        s = nf90_def_var(nc,'e3v' ,          nf90_float, (/xid,yid,depid,timid/),   ide3v)
        s = nf90_def_var(nc,'e3w' ,          nf90_float, (/xid,yid,depid,timid/),   ide3w)
        s = nf90_def_var(nc,'e3t' ,          nf90_float, (/xid,yid,depid,timid/),   ide3t)


        s = nf90_put_att(nc,idvartime,'units', UnitsTime )

        s = nf90_put_att(nc,idgdept,'units'        ,'m')
        s = nf90_put_att(nc,idgdept,'positive'     ,'down')
        s = nf90_put_att(nc,idphit, 'units'        ,'degrees_north')
        s = nf90_put_att(nc,idphit, 'long_name'    ,'Latitude')
        s = nf90_put_att(nc,idlamt, 'units'        ,'degrees_east')
        s = nf90_put_att(nc,idlamt, 'long_name'    ,'Longitude')



        s =nf90_enddef(nc)

        counter=0




        s= nf90_close(nc)
        PREPARE_BOX_FILE_T = .true.
    END FUNCTION PREPARE_BOX_FILE_T

SUBROUTINE MODIFY_NC_3D(fileNetCDF,Varname,im,jm,km,MATRIX)
use netcdf
implicit none
character fileNetCDF*(*), Varname*(*)
integer ncid, stat, VARid
integer jpi, jpj, jpk, i,j,k
real(8) MATRIX(jpk,jpj,jpi)
real copy_in(jpi,jpj,jpk)
integer mycount


      DO i=1,jpi
        DO j=1,jpj
          DO k=1,jpk
            copy_in(i,j,k) = real(MATRIX(k,i,j),4)
          ENDDO
        ENDDO
      ENDDO

mycount = 0;
stat = nf90_open( fileNetCDF, NF90_WRITE, ncid) ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_inq_varid (ncid, Varname, VARid)
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )
stat = nf90_put_var(ncid, VARid,  copy_in  )
call handle_err2(stat, fileNetCDF, Varname)     ; call handle_err1(stat,mycount,fileNetCDF )

stat=nf90_close(ncid)                           ; call handle_err1(stat,mycount,fileNetCDF )

END SUBROUTINE MODIFY_NC_3D



end PROGRAM SCALEFACTORS



