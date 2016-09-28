# -*- coding: utf-8 -*-
import numpy as np
import netCDF4 as nc
import logging
from bclib.obj_lib.lateral_obj import lateral_bc
from bclib.obj_lib.river_obj import river_data
from bclib.obj_lib.bmask_obj import bmesh
from bclib.obj_lib.submesh_obj import sub_mesh


class mesh:
    """
        Class mesh

    """

    def __init__(self,elaboration_obj):# ncf_mesh,ncf_submesh,ncf_bounmesh):
        self.input_data = elaboration_obj
        self.path = self.input_data.file_mask
        self._extract_information()


        if self.input_data.use_as_libray == False:

            logging.info("BC standalone activate")

            if self.input_data.active_atm:
                logging.info("Atmosphere enabled")
                self.submesh = sub_mesh(self,self.input_data.file_submask)
                self.submesh.atmosphere()
            else:
                logging.info("Atmosphere disabled")

            if self.input_data.active_bmask:
                logging.info("Bmask enabled")
                self.bounmesh = bmesh(self,self.input_data.file_bmask)
                self.generate_bounmask()

            else:
                logging.info("Bmask disabled")

            if self.input_data.active_gib:
                logging.info("Gib enabled")
                self.gibilterra = lateral_bc(self,self.input_data.file_nutrients)

            else:
                logging.info("Gib disabled")

            if self.input_data.active_river:
                logging.info("River enabled")
                self.river = river_data(self,self.input_data.file_river,self.input_data.file_runoff)
                self.river.map_contribute_on_sea()

            else:
                logging.info("River disabled")

            if self.input_data.active_river or self.input_data.active_gib :
                self.bc()

            logging.info("end")

        else:
            logging.info("BC lib as library mode activated")
            logging.info("Mesh partially builded")
            logging.info("You must kwon what will do!")


    def _extract_information(self):
        try:
            self.ncfile = nc.Dataset(self.path, 'r')
        except:
            print("MASKFILE NOT FOUND")
            exit()
        for i in self.ncfile.dimensions:
            #print(self.ncfile.dimensions[i])
            setattr(self, self.ncfile.dimensions[i].name, self.ncfile.dimensions[i].size)
        for i in self.ncfile.variables:
            b = self.ncfile.variables[i][:].copy()
            setattr(self, i, b)
        self.ncfile.close()

        self.tmask_dimension = self.tmask.shape




    def generate_bounmask(self):

        """ This fuction generate bounmask """

        if self.input_data.active_bmask == True :
            #input var
            bm = self.bounmesh
            bm.vnudg = self.input_data.variables
            rdpmin = self.input_data.rdpmin
            rdpmax = self.input_data.rdpmax
            #print(type(self.y),type(self.x))
            self.glamt = self.glamt.reshape(int(self.y),int(self.x))
            bm.nudg = len(bm.vnudg)
            bm.jpk = self.tmask_dimension[1]
            bm.jpjglo = self.tmask_dimension[2]
            bm.jpiglo = self.tmask_dimension[3]
            #print(bm.nudg,bm.jpk,bm.jpjglo,bm.jpiglo)
            bm.resto = np.zeros((bm.nudg,bm.jpk,bm.jpjglo,bm.jpiglo));

            for jk in range(0,bm.jpk):

                for jn in range(0,bm.nudg):
                    for jj in range(0,bm.jpjglo):
                        for ji in range(0,bm.jpiglo):

                            if (self.glamt[jj][ji] < bm.vnudg[jn][1]):

                                bm.resto[jn,jk,jj,ji]=1./(rdpmin*86400.);

                for jn in range(0,bm.nudg):
                    for jj in range(0,bm.jpjglo):
                        for ji in range(0,bm.jpiglo):
                            if (self.glamt[jj][ji] > bm.vnudg[jn][1]) and (self.glamt[jj][ji] <= self.input_data.end_nudging):
                                reltim = rdpmin + (rdpmax-rdpmin)*(self.glamt[jj][ji]-bm.vnudg[jn][1])/(self.input_data.end_nudging-bm.vnudg[jn][1]);
                                bm.resto[jn,jk,jj,ji] = 1./(reltim*86400.);


            bm.resto[:,self.tmask[0] == 0] = 1.e+20;
            count = 0
            bm.idx = np.zeros((bm.jpk,bm.jpjglo,bm.jpiglo),dtype=np.int)
            bm.water_points = np.sum(self.tmask)
            bm.idx_inv = np.zeros((bm.water_points,3),dtype=np.int);

            for jk in range(bm.jpk):
                for jj in range(bm.jpjglo):
                    for ji in range(bm.jpiglo):
                        if self.tmask[0,jk,jj,ji] == 1.0:
                            bm.idx[jk,jj,ji] = count+1;
                            bm.idx_inv[count,0]=jk+1;
                            bm.idx_inv[count,1]=jj+1;
                            bm.idx_inv[count,2]=ji+1;
                            count=count+1;


            bm.idx[self.tmask[0] == 0] = 0;
            self.bounmesh.write_netcdf()
            logging.info("bounmesh generation ended")




        else :
            logging.info("bounmesh generation disabled")



    def bc(self):

        logging.info("BC generation file start")

        if not ("bounmesh" in vars()):
            self.bounmesh = bmesh(self,self.input_data.file_bmask)
            self.bounmesh.load_bounmask()
            self.bounmesh.idx = self.bounmesh.index[0]
            #print(type(self.bounmesh.idx),self.bounmesh.idx.shape,type(self.bounmesh.index),self.bounmesh.index.shape)
            self.bounmesh.vnudg = self.input_data.variables
            self.bounmesh.rdpmin = self.input_data.rdpmin
            self.bounmesh.rdpmax = self.input_data.rdpmax
            self.glamt = self.glamt.reshape(int(self.y),int(self.x))
            self.bounmesh.nudg = len(self.bounmesh.vnudg)
            self.bounmesh.jpk = self.tmask_dimension[1]
            self.bounmesh.jpjglo = self.tmask_dimension[2]
            self.bounmesh.jpiglo = self.tmask_dimension[3]
            self.bounmesh.resto = np.zeros((self.bounmesh.nudg,self.bounmesh.jpk,self.bounmesh.jpjglo,self.bounmesh.jpiglo));
            self.bounmesh.resto[0][:] = self.bounmesh.reN1p[:]
            self.bounmesh.resto[1][:] = self.bounmesh.reN3n[:]
            self.bounmesh.resto[2][:] = self.bounmesh.reO2o[:]
            self.bounmesh.resto[3][:] = self.bounmesh.reN5s[:]
            self.bounmesh.resto[4][:] = self.bounmesh.reO3c[:]
            self.bounmesh.resto[5][:] = self.bounmesh.reO3h[:]
            self.bounmesh.resto[6][:] = self.bounmesh.reN6r[:]

        jpk = self.tmask_dimension[1]
        jpj = self.tmask_dimension[2]
        jpi = self.tmask_dimension[3]
        n_coast_cell = self.river.river_georef.shape[0]
        area = self.e1t * self.e2t
        area = area[0,0][:]

        index = self.bounmesh.idx
        print(self.river.river_years)
        if self.input_data.active_gib :
            isNudg = []#np.zeros(aux.shape,dtype=int)
            isNudg.append([])
            isNudg.append([])
            isNudg.append([])
            isNudg.append([])
            isNudg.append([])
            isNudg.append([])
            isNudg.append([])

            for jn in range(self.bounmesh.nudg):
                aux = self.bounmesh.resto[jn][:]

                print(jn)
                for k in range(jpk):
                    for j in range(jpj):
                        for i in range(jpi):
                            if (aux[k,j,i]!=0 and aux[k,j,i] < 1e+19):   # da controllare
                                isNudg[jn].append([k,j,i])

        for yr in (range(self.input_data.simulation_start_time,
                            self.input_data.simulation_end_time)):
            logging.info(str(yr))

            if self.input_data.active_gib :
                for time in range(4):
                    name_file = self.input_data.dir_out+"/GIB_"+str(yr)+self.gibilterra.season[time]+".nc"
                    ncfile = nc.Dataset(name_file, 'w')
                    for jn in range(self.bounmesh.nudg):

                        npi=len(isNudg[jn]);
                        idx=np.zeros((npi),dtype="int32");
                        data=np.zeros(npi);
                        count = 0
                        if jn == 0:
                            for cord in isNudg[jn]:
                                idx[count] = index[cord[0],cord[1],cord[2]];
                                data[count] = self.gibilterra.phos[time,cord[0],cord[1],cord[2]];
                                count = count +1;


                            ncfile.createDimension('gib_idxt_N1p',count)
                            idx_n1p = ncfile.createVariable('gib_idxt_N1p', 'i4', ('gib_idxt_N1p',))
                            idx_n1p[:] = idx[:]
                            data_n1p = ncfile.createVariable('gib_N1p', 'f',('gib_idxt_N1p',))
                            data_n1p[:] = data[:]

                        if jn == 1:
                            for cord in isNudg[jn]:
                                idx[count] = index[cord[0],cord[1],cord[2]];
                                data[count] = self.gibilterra.ntra[time,cord[0],cord[1],cord[2]];
                                count = count +1;

                            ncfile.createDimension('gib_idxt_N3n',count)
                            idx_n3n = ncfile.createVariable('gib_idxt_N3n', 'i', ('gib_idxt_N3n',))
                            idx_n3n[:] = idx[:]
                            data_n3n = ncfile.createVariable('gib_N3n', 'f',('gib_idxt_N3n',))
                            data_n3n[:] = data[:]

                        if jn == 2:
                            for cord in isNudg[jn]:
                                idx[count] = index[cord[0],cord[1],cord[2]];
                                data[count] = self.gibilterra.dox[time,cord[0],cord[1],cord[2]];
                                count = count +1;

                            ncfile.createDimension('gib_idxt_O2o',count)
                            idx_o2o = ncfile.createVariable('gib_idxt_O2o', 'i', ('gib_idxt_O2o',))
                            idx_o2o[:] = idx[:]
                            data_o2o = ncfile.createVariable('gib_O2o', 'f',('gib_idxt_O2o',))
                            data_o2o[:] = data[:]

                        if jn == 3:
                            for cord in isNudg[jn]:
                                idx[count] = index[cord[0],cord[1],cord[2]];
                                data[count] = self.gibilterra.sica[time,cord[0],cord[1],cord[2]];
                                count = count +1;


                            ncfile.createDimension('gib_idxt_N5s',count)
                            idx_n3n = ncfile.createVariable('gib_idxt_N5s', 'i', ('gib_idxt_N5s',))
                            idx_n3n[:] = idx[:]
                            data_n5s = ncfile.createVariable('gib_N5s', 'f',('gib_idxt_N5s',))
                            data_n5s[:] = data[:]

                        if jn == 4:
                            for cord in isNudg[jn]:
                                idx[count] = index[cord[0],cord[1],cord[2]];
                                data[count] = self.gibilterra.dic[time,cord[0],cord[1],cord[2]];
                                count = count +1;

                            ncfile.createDimension('gib_idxt_O3c',count)
                            idx_o3c = ncfile.createVariable('gib_idxt_O3c', 'i', ('gib_idxt_O3c',))
                            idx_o3c[:] = idx[:]
                            data_o3c = ncfile.createVariable('gib_O3c', 'f',('gib_idxt_O3c',))
                            data_o3c[:] = data[:]

                        if jn == 5:
                            for cord in isNudg[jn]:
                                idx[count] = index[cord[0],cord[1],cord[2]];
                                data[count] = self.gibilterra.alk[time,cord[0],cord[1],cord[2]];
                                count = count +1;

                            ncfile.createDimension('gib_idxt_O3h',count)
                            idx_o3h = ncfile.createVariable('gib_idxt_O3h', 'i', ('gib_idxt_O3h',))
                            idx_o3h[:] = idx[:]
                            data_o3h = ncfile.createVariable('gib_O3h', 'f',('gib_idxt_O3h',))
                            data_o3h[:] = data[:]

                        if jn == 6:
                            for cord in isNudg[jn]:
                                idx[count] = index[cord[0],cord[1],cord[2]];
                                data[count] = 0.0025;
                                count = count +1;

                            ncfile.createDimension('gib_idxt_N6r',count)
                            idx_n6r = ncfile.createVariable('gib_idxt_N6r', 'i', ('gib_idxt_N6r',))
                            idx_n6r[:] = idx[:]
                            data_n6r = ncfile.createVariable('gib_N6r', 'f',('gib_idxt_N6r',))
                            data_n6r[:] = data[:]


                    ncfile.close()

                    logging.info("--finish nutrients netcdf write")

            jpt_riv = 12
            #print("n_coast_cell",n_coast_cell)
            index_riv_a = np.zeros((  n_coast_cell,), dtype = np.int);
            position = np.zeros((  n_coast_cell,3), dtype = np.int);
            phos_riv_a  = np.zeros((  n_coast_cell,  jpt_riv));
            ntra_riv_a  = np.zeros((  n_coast_cell,  jpt_riv));
            sili_riv_a  = np.zeros((  n_coast_cell,  jpt_riv));
            alka_riv_a  = np.zeros((  n_coast_cell,  jpt_riv));
            dicc_riv_a  = np.zeros((  n_coast_cell,  jpt_riv));
            w= 1.0e+12;
            t = 1/(365 * 86400)
            n = 1/14;
            totN = 0;
            p = 1/31;
            totP = 0;
            s = 1/28;
            totS = 0;
            totA = 0;
            totD = 0;
            #print(index.shape)
            for jc in range(n_coast_cell):

                jj  = int(self.river.river_georef[jc,1]);
                ji  = int(self.river.river_georef[jc,2]);
                #print(ji,jj)
                Vol2cells = area[jj,ji]*(self.e3t[0,0,jj,ji]+self.e3t[0,0,jj,ji]);
                cn = w*t*n/Vol2cells;
                cp = w*t*p/Vol2cells;
                cs = w*t*s/Vol2cells;
                ca = w*t  /Vol2cells;
                cc = w*t  /Vol2cells;

                totN = totN + sum(self.river.river_data["DIN_KTperYR_NOBLS"][str(yr)][jc,:],2)/12;
                totP = totP + sum(self.river.river_data["DIP_KTperYR_NOBLS"][str(yr)][jc,:],2)/12;
                totS = totS + sum(self.river.river_data["DIS_KTperYR_NOBLS"][str(yr)][jc,:],2)/12;
                totA = totA + sum(self.river.river_data["ALK_GmolperYR_NOBLS"][str(yr)][jc,:],2)/12;
                totD = totD + sum(self.river.river_data["DIC_KTperYR_NOBLS"][str(yr)][jc,:],2)/12;
                index_riv_a[jc]  = index[0,jj,ji];
                position[jc] = [0,jj,ji]
                ntra_riv_a[jc,:] = self.river.river_data["DIN_KTperYR_NOBLS"][str(yr)][jc,:]*cn;
                phos_riv_a[jc,:] = self.river.river_data["DIP_KTperYR_NOBLS"][str(yr)][jc,:]*cp;
                sili_riv_a[jc,:] = self.river.river_data["DIS_KTperYR_NOBLS"][str(yr)][jc,:]*cs;
                alka_riv_a[jc,:] = self.river.river_data["ALK_GmolperYR_NOBLS"][str(yr)][jc,:]*ca;
                dicc_riv_a[jc,:] = self.river.river_data["DIC_KTperYR_NOBLS"][str(yr)][jc,:]*cc;
                # index_riv_a[jc2]  = index[1,jj,ji];
                # ntra_riv_a[jc2,:] = self.river.river_runoff_data["no3_kt_yr"][str(yr)][jc,:]*cn;
                # phos_riv_a[jc2,:] = self.river.river_runoff_data["po4_kt_yr"][str(yr)][jc,:]*cp;
                # sili_riv_a[jc2,:] = self.river.river_runoff_data["sic_kt_yr"][str(yr)][jc,:]*cs;
                # alka_riv_a[jc2,:] = self.river.river_runoff_data["alk_Gmol_yr"][str(yr)][jc,:]*ca;
                # dicc_riv_a[jc2,:] = self.river.river_runoff_data["dic_kt_yr"][str(yr)][jc,:]*cc;

            # print(index_riv_a)
            idxt_riv = np.sort(index_riv_a);
            # print(idxt_riv)
            ix = index_riv_a.argsort();
            #print("-----")
            # print(ntra_riv_a)
            n3n_riv = ntra_riv_a[ix,:];
            # print("!!!")
            # print(n3n_riv)
            # print(n3n_riv.shape)
            # print(n_coast_cell)
            pos_riv = position[ix,:]
            n1p_riv = phos_riv_a[ix,:];
            o3h_riv = alka_riv_a[ix,:];
            o3c_riv = dicc_riv_a[ix,:];
            n5s_riv = sili_riv_a[ix,:];
            count_riv = n_coast_cell;

            for mth in range(12):
                name_file = self.input_data.dir_out+"/TIN_%d%02d15-00:00:00.nc" %(yr, mth+1)

                ncfile = nc.Dataset(name_file, 'w')
                ncfile.createDimension("riv_idxt",count_riv)
                ncfile.createDimension("cords",3)
                riv_idxt_riv = ncfile.createVariable('riv_idxt', 'i4', ('riv_idxt',))
                # print("------------")
                # print(n3n_riv[:,mth])
                # print(n1p_riv[:,mth])
                # print(n5s_riv[:,mth])
                # print(o3c_riv[:,mth])
                # print(o3h_riv[:,mth])

                riv_idxt_riv[:] = idxt_riv[:]
                riv_pos_riv = ncfile.createVariable('position', 'i4', ('riv_idxt','cords'))
                riv_pos_riv[:,:] = pos_riv[:,:]
                riv_a_n3n = ncfile.createVariable('riv_N3n', 'f4', ('riv_idxt',))
                riv_a_n3n[:] = n3n_riv[:,mth]
                riv_a_n1p = ncfile.createVariable('riv_N1p', 'f4', ('riv_idxt',))
                riv_a_n1p[:] = n1p_riv[:,mth]
                riv_a_n5s = ncfile.createVariable('riv_N5s', 'f4', ('riv_idxt',))
                riv_a_n5s[:] = n5s_riv[:,mth]
                riv_a_o3c = ncfile.createVariable('riv_O3c', 'f4', ('riv_idxt',))
                riv_a_o3c[:] = o3c_riv[:,mth]
                riv_a_o3h = ncfile.createVariable('riv_O3h', 'f4', ('riv_idxt',))
                riv_a_o3h[:] = o3h_riv[:,mth]
                ncfile.close()
            logging.info("--finish river netcdf write")
