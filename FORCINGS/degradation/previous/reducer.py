import numpy as np
import netCDF4


class mesh():
    def __init__(self,filename=None, isFreeSurface=True):
        if filename is not None:
            NCin = netCDF4.Dataset(filename,'r')
            self.jpi=NCin.dimensions['x'].size
            self.jpj=NCin.dimensions['y'].size
            self.jpk   =NCin.dimensions['z'].size

            self.glamt  = np.array(NCin.variables['glamt'][0,0,:,:])
            self.glamu  = np.array(NCin.variables['glamu'][0,0,:,:])
            self.glamv  = np.array(NCin.variables['glamv'][0,0,:,:])
            self.gphit  = np.array(NCin.variables['gphit'][0,0,:,:])
            self.gphiu  = np.array(NCin.variables['gphiu'][0,0,:,:])
            self.gphiv  = np.array(NCin.variables['gphiv'][0,0,:,:])


            self.tmask  = np.array(NCin.variables['tmask'][0,:,:,:]).astype(np.bool)
            self.umask  = np.array(NCin.variables['umask'][0,:,:,:]).astype(np.bool)
            self.vmask  = np.array(NCin.variables['vmask'][0,:,:,:]).astype(np.bool)

            if isFreeSurface:
                self.e3t_0  = np.array(NCin.variables['e3t_0'][0,:,:,:])#_0 a 1/24
                self.e3u_0  = np.array(NCin.variables['e3u_0'][0,:,:,:])
                self.e3v_0  = np.array(NCin.variables['e3v_0'][0,:,:,:])
                self.e3w_0  = np.array(NCin.variables['e3w_0'][0,:,:,:])

            self.e1t  = np.array(NCin.variables['e1t'][0,0,:,:])
            self.e1u  = np.array(NCin.variables['e1u'][0,0,:,:])
            self.e1v  = np.array(NCin.variables['e1v'][0,0,:,:])
            self.e2t  = np.array(NCin.variables['e2t'][0,0,:,:])
            self.e2u  = np.array(NCin.variables['e2u'][0,0,:,:])
            self.e2v  = np.array(NCin.variables['e2v'][0,0,:,:])
            self.gdept  = np.array(NCin.variables['gdept'][0,:,0,0])
            self.gdepw  = np.array(NCin.variables['gdepw'][0,:,0,0])
            self.nav_lev  = np.array(NCin.variables['nav_lev'])
            NCin.close()
        else:
            self.jpi=None
            self.jpj=None
            self.jpk=None

            self.glamt  = None
            self.glamu  = None
            self.glamv  = None
            self.gphit  = None
            self.gphiu  = None
            self.gphiv  = None


            self.gmask  = None
            self.gmask  = None
            self.gmask  = None
            self.e3t_0  = None
            if isFreeSurface:
                self.e3u_0  = None
                self.e3v_0  = None
                self.e3w_0  = None
            self.e1t  = None
            self.e1u  = None
            self.e1v  = None
            self.e2t  = None
            self.e2u  = None
            self.e2v  = None
            self.gdept  = None
            self.gdepw  = None
            self.nav_lev  = None
            self.finecells = None


class reducer():
    def __init__(self, isFreeSurface=True):
        self.time = 1
        self.z_a = 1
        self.x_a = 1
        self.y_a = 1
        self.isFreeSurface = isFreeSurface
        self.jpi = None
        self.jpj = None
        self.jpk = None

    def read_fine(self,filename):
        self.Mfine = mesh(filename=filename,isFreeSurface=self.isFreeSurface)




    def reducing_by_boxes_2d_sum(self, M2d, I, J, direction=1):
        ''' Valid for e2u, e2v etc'''
        jpi = len(I)
        jpj = len(J)
        step = self.step
        m2d = np.zeros((jpj,jpi),np.double)
        for ji in range(jpi):
            for jj in range(jpj):
                j=J[jj]
                i=I[ji]
                if direction ==1:
                    m2d[jj,ji] = np.sum(M2d[j+step-1,i:i+step])
                if direction ==2:
                    m2d[jj,ji] = np.sum(M2d[j:j+step,i+step-1])
        return m2d

    def reducing_by_boxes_3d_mean(self, M3d, W3d, I, J, direction='u'):
        ''' Valid for e3u, e3v '''
        jpi = len(I)
        jpj = len(J)
        jpk = self.jpk
        step = self.step
        m3d = np.zeros((jpk,jpj,jpi),np.double)

        for jk in range(jpk):
            for ji in range(jpi):
                for jj in range(jpj):
                    j=J[jj]
                    i=I[ji]
                    if direction =="u":
                        e = M3d[jk,j:j+step,i+step-1]
                        w = W3d[j:j+step,i+step-1]
                        ew = (e*w).sum()
                        m3d[jk,jj,ji] = ew/w.sum()
                    if direction =="v":
                        e  = M3d[jk,j+step-1,i:i+step]
                        w  = W3d[j+step-1,i:i+step]
                        ew = (e*w).sum()
                        m3d[jk,jj,ji] = ew/w.sum()
        return m3d


    def reducing_by_boxes_double(self, M3d, I, J):
        '''
        valid for e3u_0, e3v_0 
        mean over box
        meglio fare solo i punti d'acqua?
        '''
        jpi = self.jpi
        jpj = self.jpj
        jpk = self.jpk
        step = self.step
        m3d = np.zeros((jpk,jpj,jpi),np.double)
        for ji in range(jpi):
            for jj in range(jpj):
                j=J[jj]
                i=I[ji]
                for jk in range(jpk):
                    m3d[jk,jj,ji] = np.mean(M3d[jk,j:j+step,i:i+step])
        return m3d



    def reducing_by_boxes_bool(self,Bool3d, I, J, direction='u'):
        '''
        Reduction of  mask arrays
        '''
        jpi = self.jpi
        jpj = self.jpj
        jpk = self.jpk
        step = self.step
        bool3d = np.zeros((jpk,jpj,jpi),np.bool)
        count3d= np.zeros((jpk,jpj,jpi),np.int32)
        for ji in range(jpi):
            for jj in range(jpj):
                j=J[jj]
                i=I[ji]
                if direction == 'u':
                    for jk in range(jpk):
                        bool3d[jk,jj,ji] = np.any(Bool3d[jk,j:j+step,i+step-1])
                if direction == 'v':
                    for jk in range(jpk):
                        bool3d[jk,jj,ji] = np.any(Bool3d[jk,j+step-1,i:i+step])
                if direction == 't':
                    for jk in range(jpk):
                        bool3d[jk,jj,ji] = np.any(Bool3d[jk,j:j+step,i:i+step])
                        count3d[jk,jj,ji]= np.sum(Bool3d[jk,j:j+step,i:i+step])
        return bool3d, count3d

#     def reducing_by_boxes_2d_sum(self, M2d, I, J, direction="u"):
#         ''' Valid for e2u, e2v etc'''
#     
#         m2d = np.zeros((jpj,jpi),np.double)
#         for ji in range(jpi):
#             for jj in range(jpj):
#                 j=J[jj]
#                 i=I[ji]
#                 if direction =="u":
#                     m2d[jj,ji] = np.sum(M2d[j:j+step,i+step-1])
#                 if direction =="v":
#                     m2d[jj,ji] = np.sum(M2d[j+step-1,i:i+step])
#         return m2d
    def reducing_by_boxes_2d_mean(self, M2d, I, J, direction="lam"):
        '''
        Valid for glamu, glamt
        lam = longitude index
        phy = latitude  index
         U box  + + *     V box    * * *       T box:   * * *  
                + + *              + + +                * * *
                + + *              + + +                * * *
        '''
        step = self.step
        m2d = np.zeros((self.jpj,self.jpi),np.double)
        for ji in range(self.jpi):
            for jj in range(self.jpj):
                j=J[jj]
                i=I[ji]
                if direction =="lam":
                    m2d[jj,ji] = np.mean(M2d[j:j+step,  i+step-1])
                if direction =="phi":
                    m2d[jj,ji] = np.mean(M2d[j+step-1   ,i:i+step])
#               if direction =="t":
#                    m2d[jj,ji] = np.mean(M2d[j:j+step, i:i+step])
    
        return m2d

    def reduce(self, step):
        '''
        Generates a mask object reduced by step
        Arguments:
        * step * integer
        
        The coarse mesh matches exactly the fine one at East and North.
        Returns:
        * m *  mymesh object
        * I *  numpy array of integers, such as coarse_lon[k] = fine_lon[I[k]]
               I[k] is the left index of the box.
        * J *  idem, for latitude.
        '''
        
        print "Start Mesh reduction"
        m = mesh(filename=None,isFreeSurface=self.isFreeSurface)
        self.step = step

        m.jpi, i_start = divmod(self.Mfine.jpi,step)
        m.jpj, j_start = divmod(self.Mfine.jpj,step)
        m.jpk = self.Mfine.jpk

        self.jpi = m.jpi
        self.jpj = m.jpj
        self.jpk = m.jpk

        I_glo = np.arange(self.Mfine.jpi)
        J_glo = np.arange(self.Mfine.jpj)
        I = I_glo[i_start::step]
        J = J_glo[j_start::step]


        m.tmask, m.finecells = self.reducing_by_boxes_bool(self.Mfine.tmask, I, J, 't')
        m.umask,  _          = self.reducing_by_boxes_bool(self.Mfine.umask, I, J, 'u')
        m.vmask,  _          = self.reducing_by_boxes_bool(self.Mfine.vmask, I, J, 'v')
        
        
        if self.isFreeSurface:
            print "Start e3 section"
            m.e3t_0 = self.reducing_by_boxes_double(self.Mfine.e3t_0, I, J)
            m.e3u_0 = self.reducing_by_boxes_3d_mean(self.Mfine.e3u_0, self.Mfine.e2u, I, J, 'u')
#           m.e3u_0 = self.reducing_by_boxes_double(self.Mfine.e3u_0, I, J)
            m.e3v_0 = self.reducing_by_boxes_3d_mean(self.Mfine.e3v_0, self.Mfine.e1v, I, J, 'v')
#           m.e3v_0 = self.reducing_by_boxes_double(self.Mfine.e3v_0, I, J)
            m.e3w_0 = self.reducing_by_boxes_double(self.Mfine.e3w_0, I, J)
        
        print "Start e1 e2 section"
        m.e1t = self.reducing_by_boxes_2d_sum(self.Mfine.e1t, I, J, 1)
        m.e1u = self.reducing_by_boxes_2d_sum(self.Mfine.e1u, I, J, 1)
        m.e1v = self.reducing_by_boxes_2d_sum(self.Mfine.e1v, I, J, 1)
        m.e2t = self.reducing_by_boxes_2d_sum(self.Mfine.e2t, I, J, 2)
        m.e2u = self.reducing_by_boxes_2d_sum(self.Mfine.e2u, I, J, 2)
        m.e2v = self.reducing_by_boxes_2d_sum(self.Mfine.e2v, I, J, 2)

        print "Start phi lam section"
        m.glamt = self.reducing_by_boxes_2d_mean(self.Mfine.glamt, I, J,'lam')
        m.glamu = self.reducing_by_boxes_2d_mean(self.Mfine.glamu, I, J,'lam')
        m.glamv = self.reducing_by_boxes_2d_mean(self.Mfine.glamv, I, J,'lam')
        m.gphit = self.reducing_by_boxes_2d_mean(self.Mfine.gphit, I, J,'phi')
        m.gphiu = self.reducing_by_boxes_2d_mean(self.Mfine.gphiu, I, J,'phi')
        m.gphiv = self.reducing_by_boxes_2d_mean(self.Mfine.gphiv, I, J,'phi')

        m.gdept   = self.Mfine.gdept
        m.gdepw   = self.Mfine.gdepw
        m.nav_lev = self.Mfine.nav_lev
        return m,I,J

    def nearest_extend3D(self,M3d):
        '''
        Argument
        * M3d * a numpy array of dimension (jpk, jpj,jpi)

        Returns:
        * out * a numpy array of dimension (jpk,jpj+2,jpi+2),
        with frame of one cell all around, in longitude and latitude
        '''
        jpk, jpj, jpi = M3d.shape
        out = np.ones((jpk,jpj+2,jpi+2), dtype=M3d.dtype)
        out[:,1:-1, 1:-1] = M3d
        out[:, 1:-1 ,0] = M3d[:,:, 0]
        out[:, 1:-1,-1] = M3d[:,:,-1]
        out[:, 0 ,1:-1] = M3d[:, 0,:]
        out[:,-1 ,1:-1] = M3d[:,-1,:]
        return out

    def nearest_extend1D(self,array):
        nP=len(array)+2
        out = np.zeros((nP),dtype=array.dtype)
        out[0]    = array[ 0]
        out[1:-1] = array
        out[-1]   = array[-1]
        return out

    def nearest_extend2D(self,M2d):
        '''
        Argument:
        * M2d * a numpy array having dimension (jpj,jpi)

        Returns:
        * out * a numpy array having dimension (jpj+2, jpi+2)
                with a frame of one cell all around, obtainded by
                nearest extrapolation
        '''

        jpj, jpi = M2d.shape
        out = np.zeros((jpj+2,jpi+2), dtype=M2d.dtype)
        out[1:-1, 1:-1] = M2d
        out[: ,0] = self.nearest_extend1D(M2d[:, 0])
        out[:,-1] = self.nearest_extend1D(M2d[:,-1])
        out[0, :] = self.nearest_extend1D(M2d[ 0,:])
        out[-1,:] = self.nearest_extend1D(M2d[-1,:])
        return out
    def extrap_2d(self,M2d):
        '''
        Argument:
        * M2d * a numpy array having dimension (jpj,jpi)

        Returns:
        * out * a numpy array having dimension (jpj+2, jpi+2)
                with a frame of one cell all around, obtainded by
                linear extrapolation
        The idea is:

        F(n+1)  = f(n) + ( f(n) - (f(n-1)) )  = 2*f(n) - f(n-1)
        F(-1)   = f(0) - ( f(1) - f(0)     )  = 2*f(0) - f(1)
        '''
        jpj, jpi = M2d.shape
        out = np.zeros((jpj+2,jpi+2), dtype=M2d.dtype)
        out[1:-1, 1:-1] = M2d
        out[:, 0] = self.nearest_extend1D( 2*M2d[:, 0] - M2d[:, 1] )
        out[:,-1] = self.nearest_extend1D( 2*M2d[:,-1] - M2d[:,-2] )
        out[0, :] = self.nearest_extend1D( 2*M2d[0, :] - M2d[1, :] )
        out[-1,:] = self.nearest_extend1D( 2*M2d[-1,:] - M2d[-2,:] )
        return out

    def expand(self,m,I,J):
        '''
        Argument:
        * m * a mymesh object, having dimensions jpk,jpj,jpi
        * I * numpy array[jpi] of indexes
        * J * numpy array[jpi] of indexes
        Returns:

        * M * a mymesh object, having dimensions jpk, jpj+2, jpi+2
              obtained by generating a frame of one cell all around,
              in longitude and latitude.
        * I * numpy array of indexes
        '''
        M = mesh(filename=None, isFreeSurface=self.isFreeSurface)
        M.jpi = m.jpi+2
        M.jpj = m.jpj+2
        M.jpk = m.jpk

        I_exp = np.ones((M.jpi),np.int32)*(-999)
        J_exp = np.ones((M.jpj),np.int32)*(-999)
        I_exp[1:-1]=I
        J_exp[1:-1]=J

        # here pure frame of zeros
        M.tmask = np.zeros((M.jpk,M.jpj,M.jpi),np.bool)
        M.umask = np.zeros((M.jpk,M.jpj,M.jpi),np.bool)
        M.vmask = np.zeros((M.jpk,M.jpj,M.jpi),np.bool)
        M.tmask[:,1:-1, 1:-1] = m.tmask
        M.umask[:,1:-1, 1:-1] = m.umask
        M.vmask[:,1:-1, 1:-1] = m.vmask
        M.finecells = np.zeros((M.jpk,M.jpj,M.jpi),np.int32)
        M.finecells[:,1:-1, 1:-1] = m.finecells
        #------------------------------------
        if self.isFreeSurface:
            M.e3t_0 = self.nearest_extend3D(m.e3t_0)
            M.e3u_0 = self.nearest_extend3D(m.e3u_0)
            M.e3v_0 = self.nearest_extend3D(m.e3v_0)
            M.e3w_0 = self.nearest_extend3D(m.e3w_0)
        M.e1t = self.nearest_extend2D(m.e1t)
        M.e1u = self.nearest_extend2D(m.e1u)
        M.e1v = self.nearest_extend2D(m.e1v)
        M.e2t = self.nearest_extend2D(m.e2t)
        M.e2u = self.nearest_extend2D(m.e2u)
        M.e2v = self.nearest_extend2D(m.e2v)


        M.glamu = self.extrap_2d(m.glamu)
        M.glamv = self.extrap_2d(m.glamv)
        M.glamt = self.extrap_2d(m.glamt)
        M.gphiu = self.extrap_2d(m.gphiu)
        M.gphiv = self.extrap_2d(m.gphiv)
        M.gphit = self.extrap_2d(m.gphit)

        M.gdept   = m.gdept
        M.gdepw   = m.gdepw
        M.nav_lev = m.nav_lev

        return M,I_exp,J_exp


    def cut_med(self,m,lon_cut=-8.875,depth_cut=0,biscay_land=True):
        '''
        Cuts the atlantic buffer and unuseful bottom of the Med Sea

        Arguments :
        *  m          * a mesh object
        * lon_cut     * float indicating the first longitude of returned mesh
        * depth_cut   * integer, number of layers of bottom which are land points
        * biscay_land * logical flag; if True, Biscay Gulf will be set as land
        
        Returns :
        * M * mesh object, ready for ogstm
        * lon_ind_cut  * integer
        '''
        lon_ind_cut =np.argmin(np.abs(m.glamt[0,:]-lon_cut))

        M = mesh(filename=None, isFreeSurface=self.isFreeSurface)

        M.jpi = m.jpi - lon_ind_cut
        M.jpj = m.jpj
        M.jpk = m.jpk - depth_cut



        M.tmask = m.tmask[:M.jpk,:,lon_ind_cut:]
        M.umask = m.umask[:M.jpk,:,lon_ind_cut:]
        M.vmask = m.vmask[:M.jpk,:,lon_ind_cut:]

        M.tmask[:,:,0] = False
        M.vmask[:,:,0] = False
        M.tmask[:,:,0] = False
        M.finecells = m.finecells[:M.jpk,:,lon_ind_cut:]
        if self.isFreeSurface:
            M.e3t_0 = m.e3t_0[:M.jpk,:,lon_ind_cut:]
            M.e3u_0 = m.e3u_0[:M.jpk,:,lon_ind_cut:]
            M.e3v_0 = m.e3v_0[:M.jpk,:,lon_ind_cut:]
            M.e3w_0 = m.e3w_0[:M.jpk,:,lon_ind_cut:]
        M.e1t = m.e1t[:,lon_ind_cut:]
        M.e1u = m.e1u[:,lon_ind_cut:]
        M.e1v = m.e1v[:,lon_ind_cut:]
        M.e2t = m.e2t[:,lon_ind_cut:]
        M.e2u = m.e2u[:,lon_ind_cut:]
        M.e2v = m.e2v[:,lon_ind_cut:]


        M.glamu = m.glamu[:,lon_ind_cut:]
        M.glamv = m.glamv[:,lon_ind_cut:]
        M.glamt = m.glamt[:,lon_ind_cut:]
        M.gphiu = m.gphiu[:,lon_ind_cut:]
        M.gphiv = m.gphiv[:,lon_ind_cut:]
        M.gphit = m.gphit[:,lon_ind_cut:]
        M.gdept   = m.gdept[:M.jpk]
        M.gdepw   = m.gdepw[:M.jpk]
        M.nav_lev = m.nav_lev[:M.jpk]

        if biscay_land:
            bool1 = (M.glamt < 0)  &   (M.gphit > 42)
            bool2 = (M.glamt < -6.)  & (M.gphit > 37.25)
            land = bool1 | bool2

            for k in range(M.jpk):
                tmask = M.tmask[k,:,:]
                umask = M.umask[k,:,:]
                vmask = M.vmask[k,:,:]
                tmask[land] = False
                umask[land] = False
                vmask[land] = False
        return M, lon_ind_cut


    def dumpfile(self, m, outfile):
        ncOUT=netCDF4.Dataset(outfile,"w");
        
        ncOUT.createDimension('x',m.jpi);
        ncOUT.createDimension('y',m.jpj);
        ncOUT.createDimension('z',m.jpk);
        ncOUT.createDimension('time',self.time)
        
        ncOUT.createDimension('x_a',self.x_a);
        ncOUT.createDimension('y_a',self.y_a);
        ncOUT.createDimension('z_a',self.z_a);
        

        ncvar    = ncOUT.createVariable('e1t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[0,0,:,:] = m.e1t   ;
        ncvar    = ncOUT.createVariable('e1u'   ,'d',('time','z_a', 'y', 'x')   ); ncvar[0,0,:,:] = m.e1u   ;
        ncvar    = ncOUT.createVariable('e1v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[0,0,:,:] = m.e1v   ;

        ncvar    = ncOUT.createVariable('e2t'   ,'d',('time','z_a', 'y', 'x')  ) ; ncvar[0,0,:,:] = m.e2t   ;
        ncvar    = ncOUT.createVariable('e2u'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[0,0,:,:] = m.e2u   ;
        ncvar    = ncOUT.createVariable('e2v'   ,'d',('time','z_a', 'y', 'x'))   ; ncvar[0,0,:,:] = m.e2v   ;     
 
        if self.isFreeSurface:
            ncvar    = ncOUT.createVariable('e3t_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[0,:,:,:] = m.e3t_0
            ncvar    = ncOUT.createVariable('e3u_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[0,:,:,:] = m.e3u_0
            ncvar    = ncOUT.createVariable('e3v_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[0,:,:,:] = m.e3v_0
            ncvar    = ncOUT.createVariable('e3w_0'   ,'d',('time','z', 'y', 'x'))     ; ncvar[0,:,:,:] = m.e3w_0   
        
                 
 
        ncvar    = ncOUT.createVariable('gdept' ,'d',('time','z', 'y_a', 'x_a')) ; ncvar[0,:,0,0] = m.gdept
        ncvar    = ncOUT.createVariable('gdepw' ,'d',('time','z', 'y_a', 'x_a')) ; ncvar[0,:,0,0] = m.gdepw
    
        ncvar    = ncOUT.createVariable('glamt'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[0,0,:,:] = m.glamt 
        ncvar    = ncOUT.createVariable('glamu'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[0,0,:,:] = m.glamu      
        ncvar    = ncOUT.createVariable('glamv'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[0,0,:,:] = m.glamv 
    
        ncvar    = ncOUT.createVariable('gphit'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[0,0,:,:] = m.gphit 
        ncvar    = ncOUT.createVariable('gphiu'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[0,0,:,:] = m.gphiu      
        ncvar    = ncOUT.createVariable('gphiv'   ,'d',('time','z_a', 'y', 'x')) ; ncvar[0,0,:,:] = m.gphiv 
        ncvar    = ncOUT.createVariable('nav_lat','f',('y','x'))                 ; ncvar[:] = m.gphit
        ncvar    = ncOUT.createVariable('nav_lev' ,'f',('z',))                   ; ncvar[:] = m.nav_lev
        ncvar    = ncOUT.createVariable('nav_lon','f',('y','x'))                 ; ncvar[:] = m.glamt
        ncvar    = ncOUT.createVariable('tmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = m.tmask 
        ncvar    = ncOUT.createVariable('umask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = m.umask 
        ncvar    = ncOUT.createVariable('vmask' ,'d',('time','z', 'y', 'x') )    ; ncvar[:] = m.vmask
        ncvar    = ncOUT.createVariable('finecells' ,'i',('time','z', 'y', 'x') ); ncvar[:] = m.finecells
        ncOUT.close()

    def nearest_2d_cutting_along_slices(self, M2d,I,J):
        '''
        Performs a nearest 2d interp over two regular grids, based on precalculated indexes.
        Returns a reduced 2d array, having dimension (len(J),len(I))
        
        Arguments:
        * M2d *  2d numpy array
        * I   *  integer numpy array
                 M2d longitude[I[k]] is the nearest of M2d_reduced longitude[k]
        * J   *  integer numpy array, same meaning.
        
        
        Returns: 
        * M2d_reduced * 2d numpy array 
        
        '''
        jpi = len(I)
        jpj = len(J)
        M2d_reduced = np.zeros((jpj,jpi), dtype=M2d.dtype)
        for ji in range(jpi):
            for jj in range(jpj):
                M2d_reduced[jj,ji] = M2d[J[jj],I[ji]]
        return M2d_reduced

if __name__=="__main__":
    import numpy as np
    isFreeSurface=True
    inputmesh="meshmask.nc"
    import pickle

    # just to save time for tests
    #fid=open("temp.dat")
    #LIST=pickle.load(fid)
    #[R, Mred, I, J ] = LIST
    #fid.close()

    R =reducer(isFreeSurface)
    R.read_fine(inputmesh)
    passo=4
    Mred,I,J = R.reduce(passo)
    Mexp,I,J = R.expand(Mred,I,J)
    outputmesh="meshmask_6125.nc"
    R.dumpfile(Mexp, outputmesh)
    np.savetxt('South_West_I_indexes.txt',I,'%d') # of mesh not expanded
    np.savetxt('South_West_J_indexes.txt',J,'%d')




    
