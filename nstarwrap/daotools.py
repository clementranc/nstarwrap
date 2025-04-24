import ctypes as ct
import numpy as np
import pandas as pd
import pickle
import shutil
import sys

class ListStar:
    """
    """
    def __init__(self,
            file: str,
            path: str = './',
            **kwargs):

        self.path = path
        self.name = file
        self.extension = file.split('.')[-1]
        self._load_list_stars()
        
    def _load_list_stars(self):
        """
        """
        filename = f'{self.path}/{self.name}'
        if self.extension == 'lst':
            colnames = ['id', 'x', 'y', 'mag', 'err1', 'err2']
        elif self.extension == 'coo':
            colnames = ['id', 'x', 'y', 'mag', 'err1', 'err2', 'last']
        else:
            colnames = ['id', 'x', 'y', 'mag', 'err1', 'err2']
        col = range(len(colnames))
        fmt = dict()
        [fmt.update({a: np.float64}) for a in colnames]
        fmt.update({'id': np.int64})
        self.table = pd.read_table(filename, sep=r"\s+", names=colnames, 
                                       usecols=col, dtype=fmt, skiprows=3)
        
        # Load 3 first lines of input file
        three_first_lines = ''
        file = open(filename, 'r')
        for i in range(3):
            three_first_lines = f'{three_first_lines}{file.readline()}'
        file.close()
        self.header = three_first_lines
        
    def to_ds9(self,
               filename: str,
               path: str = None,
               markersize: int = 10,
               markercolor: int = 'green',
               linewidth: float = 1,
               marker: str = 'circle',
               text_offset = [0, 10],
               **kwargs):
        """
        """
        
        if path == None: path = self.path

        opts = dict({
            'marker': marker,
            'markersize': markersize,
            'markercolor': markercolor,
            'linewidth': linewidth,
            'offset': text_offset,
        })
        label_short = ['ms', 'mc', 'lw']
        label_long = ['markersize', 'markercolor', 'linewidth']
        defaults_rank = [1, 2, 3]
        for i in range(len(label_short)):
            if label_short[i] in kwargs:
                if not kwargs[label_short[i]] == self.to_ds9.__defaults__[defaults_rank[i]]:
                    opts.update({label_long[i] : kwargs[label_short[i]]})

        # Create DS9 regions to check the position of the lens and source
        txt = ('# Region file format: DS9 version 4.1\n'
               'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
               'image\n')
        
        for i in range(len(self.table)):
            star = self.table.iloc[[i], [self.table.columns.get_loc('x'), 
                                         self.table.columns.get_loc('y')]].values[0]
            starid = self.table.iloc[[i], [self.table.columns.get_loc('id')]].values[0]
            offset = 25
            textxy = star + [opts['offset'][0], opts['markersize'] + opts['offset'][1]]
            if opts['marker'] == 'point':
                txt = f"{txt}{opts['marker']}({star[0]:.3f},{star[1]:.3f}) # point=x {opts['markersize']} color={opts['markercolor']}\n"
            else:
                txt = f"{txt}{opts['marker']}({star[0]:.3f},{star[1]:.3f},{opts['markersize']}) # width={opts['linewidth']} color={opts['markercolor']}\n"
            txt = f"{txt}text {textxy[0]:.3f} {textxy[1]:.3f} # text={{{starid[0]:d}}} width={opts['linewidth']} color={opts['markercolor']}\n"

        fname = f'{path}/{filename}'
        if len(fname) > 4:
            if not fname[-4:] == '.reg': fname = f'{fname}.reg'
        else: fname = f'{fname}.reg'
        file = open(fname, 'w')
        file.write(txt)
        file.close()


    def to_coo(self,
               filename: str,
               path: str = '.',
               **kwargs):
        """Write a coo file from a list of stars.
        """
        input1 = self.table
        txt = f'{self.header}'
        for i in range(len(input1)):
            txt = (f"{txt}{int(input1['id'].values[i]):7d} {input1['x'].values[i]:8.2f} {input1['y'].values[i]:8.2f} "
                   f"{input1['mag'].values[i]:8.3f} {input1['err1'].values[i]:8.3f} {input1['err2'].values[i]:8.3f} "
                   f"{input1['last'].values[i]:8.3f}\n")

        fn = f'{path}/{filename}'
        print("Create file:")
        file = open(fn, 'w')
        file.write(txt)
        file.close()
        print(f"{fn}")


    def find_my_stars(self,
               my_stars: list,
               inplace: bool = False,
               **kwargs):
        """Cross-match with a table of star coordinates.


            my_stars = [[1093.82, 1171.43, 'Unrelated South West'],
                        [1108.85, 1166.67, 'Unrelated South East'],
                        [1112.51, 1183.30, 'Target']]


        """
        if not inplace:
            input1 = self.table.copy()
        else:
            input1 = self.table

        input1['distance'] = input1['x']
        for i in range(len(my_stars)):
            input1['distance'] = np.sqrt(np.power(input1['x'] - my_stars[i][0], 2) + np.power(input1['y'] - my_stars[i][1], 2))
            mask = input1['distance'] == np.min(input1['distance'])

            if my_stars[i][2] == "Target":
                tmp = input1.sort_values('distance')
                print('TARGET (Lens + Source): distances are in pixels')
                print(tmp[['id', 'x', 'y', 'distance']].head(4).to_string(index=False))
            else:
                print(f"{my_stars[i][2]:20} --> {input1.loc[mask, 'id'].values[0]:7d}  [distance: {np.min(input1['distance']):.1e} pixels]")

        if not inplace:
            return input1

def load_star_list_file(file_name = 'image.lst', path = ".", show=1, header=3):
    """Load a list of file, convention of DAOPHOT-II *.lst files."""

    # --- PARAMETERS ---
    file_input1 = f'{path}/{file_name}'
    colnames = ['id', 'x', 'y', 'mag', 'err1', 'err2']
    col = range(len(colnames))
    fmt = dict()
    [fmt.update({a: np.float64}) for a in colnames]
    fmt.update({'id': np.int64})
    input1 = pd.read_table(file_input1, sep=r"\s+", names=colnames, 
                               usecols=col, dtype=fmt, skiprows=header,
                               comment='#')

    if show > 0:
        print(f'Input list of stars in file {file_input1}:')
        print(input1.head(show))

    # Load 3 first lines of input file
    three_first_lines = ''
    file = open(file_name, 'r')
    for i in range(header):
        three_first_lines = f'{three_first_lines}{file.readline()}'
    file.close()

    return input1, three_first_lines

def cross_match_files(input1, input2):
    """Cross-match two tables of stars based on pixel coordinates."""

    # Cross-matching stars
    print('Cross-matching the two tables of stars...')
    print(' Old ID      New ID')
    input1['match'] = 0
    input2['distance'] = input2['x']
    for i in range(len(input1)):
        input2['distance'] = np.power(input1.loc[i, 'x'] - input2['x'], 2) + np.power(input1.loc[i, 'y'] - input2['y'], 2)
        mask = input2['distance'] == np.min(input2['distance'])
        input1.loc[i, 'match'] = input2.loc[mask, 'id'].values[0]
        print(f"{input1.loc[i, 'id']:7d} --> {input2.loc[mask, 'id'].values[0]:7d}")


class NstarPythonWrapper:
    """
    """
    def __init__(self,
            image : str = u'image.fits',
            psf : str = u'image.psf',
            group : str = u'image.grp',
            maxcol_daophot : int = 16384,
            watch_daophot : float = 0.0,
            fitrad_daophot : float = 9.0,
            e1_daophot : float = 0.75,
            e2_daophot : float = 5.0,
            fitting_box : list = [-1, -1, -1, -1],
            err_factor : float = 1.0,
            pickle : bool = False,
            **kwargs):

        self.image_name = u"{0}".format(image)
        self.maxcol_daophot = maxcol_daophot
        self.psf_name = u"{0}".format(psf)
        self.group_name = u"{0}".format(group)
        self.pickle = pickle

        self.opt = dict({
            'watch' : ct.c_float(watch_daophot),
            'fitrad' : ct.c_float(fitrad_daophot),
            'e1' : ct.c_float(e1_daophot),
            'e2' : ct.c_float(e2_daophot),
            'box_xmin' : ct.c_int(int(fitting_box[0])),
            'box_xmax' : ct.c_int(int(fitting_box[1])),
            'box_ymin' : ct.c_int(int(fitting_box[2])),
            'box_ymax' : ct.c_int(int(fitting_box[3])),
            'err_factor' : ct.c_double(err_factor)})

        self.daophot_path = self.__find_daophot_path__()
        self.daolib = self.__load_fortran_lib__(self.daophot_path)
        self.__init_types_before_nstar__()

        nmax_fit_box = 3000
        self.residuals = np.zeros(nmax_fit_box, order="F")
        self.residuals_ptr = self.residuals.ctypes.data_as(ct.POINTER(ct.c_double))
        self.sigmas = np.zeros(nmax_fit_box, order="F")
        self.sigmas_ptr = self.sigmas.ctypes.data_as(ct.POINTER(ct.c_double))
        self.model = np.zeros(nmax_fit_box, order="F")
        self.model_ptr = self.model.ctypes.data_as(ct.POINTER(ct.c_double))

        # Table of additional optional parameters
        nmax_opt_params = 4
        self.opt_params = np.zeros(nmax_opt_params, order="F")
        self.opt_params_ptr = self.opt_params.ctypes.data_as(ct.POINTER(ct.c_double))

    def __find_daophot_path__(self):
        """Find the path to DAOPHOT."""

        daophot_path = shutil.which("daophot").split('/')
        daophot_path = [a for a in daophot_path if a != '']
        daophot_path = '/{0}'.format('/'.join(daophot_path[:-1]))

        return daophot_path

    def __assign_chararray__(self, string):
        """Create a character array from a string."""
        arr = np.chararray(len(string), unicode=False, order='F')
        for i in range(len(string)):
            arr[i] = string[i]
        return arr

    def __load_fortran_lib__(self, daophot_path):
        """Load the shared fortran library giving access to DAOPHOT."""
        print(f'{daophot_path}/pywrapper.so')
        try:
            fortlib = ct.CDLL(f'{daophot_path}/pywrapper.so')
        except:
            sys.exit('DAOPHOT shared fortran library "pywrapper" could not be loaded. Abort!')
        return fortlib
        
    def attach(self, image='', pickle=False):
        """Attach the image."""
        pickle = self.pickle | pickle
        ext = self.image_name.split('.')[-1]
        if ext == 'fits':
            fattach = self.daolib.imagesize
            fattach.argtypes=[
                ct.POINTER(ct.c_char), ct.POINTER(ct.c_int),
                ct.POINTER(ct.c_int), ct.c_int]
            nrow = ct.c_int(0)
            ncol = ct.c_int(0)
            if image=='': 
                image = self.image_name
                imagel = ct.c_int(len(image))
            string  = u'{0}'.format(image)
            fname_image = self.__assign_chararray__(string)
            fname_image_ptr = fname_image.ctypes.data_as(ct.POINTER(ct.c_char))
            try:
                fattach(fname_image_ptr,  ct.byref(ncol), ct.byref(nrow), imagel)
            except:
                sys.exit('DAOPHOT ATTACH command failed. Check the image file name and location.')

            self.ncol = ncol
            self.nrow = nrow
            self.load_picture(pickle=pickle)

        elif ext == 'pkl':
            ncol, nrow = self.__load_pickleize_picture__()

        else:
            sys.exit("Image extension {ext} not recognized.")

        return ncol, nrow
        
    def load_picture(self, pickle=False):
        """Load a picture and keep it in memory."""

        floadim = self.daolib.loadimage
        floadim.argtypes=[ct.c_int, ct.POINTER(ct.c_float), ct.c_int, 
            ct.POINTER(ct.c_int)]
        maxcol = ct.c_int(self.maxcol_daophot)
        ier = ct.c_int(-1)
        data = np.zeros((maxcol.value, maxcol.value), order="F")
        data_ptr = data.ctypes.data_as(ct.POINTER(ct.c_float))
        try:
            floadim(self.nrow, data_ptr, maxcol, ct.byref(ier))
        except:
            sys.exit('Could not load the picture in memory. Check file names and paths.')
        
        if (ier.value > 0) | (ier.value <0):
            sys.exit('Could not load the picture in memory. Check file names and paths.')

        self.picture = data
        self.picture_ptr = data_ptr

        if pickle:
            self.__pickleize_picture__()

    def __pickleize_picture__(self):
        
        objects_to_save = dict({
            'picture' : self.picture,
            'ncol' : self.ncol,
            'nrow' : self.nrow,
            })
        fname = f"{self.image_name}.pkl"
        save_file = open(fname, "wb")
        pickle.dump(objects_to_save, save_file)
        save_file.close()

    def __load_pickleize_picture__(self):
        
        file = open(self.image_name, "rb")
        objects = pickle.load(file)
        file.close()

        #print(objects)

        self.picture = objects['picture']
        self.picture_ptr = self.picture.ctypes.data_as(ct.POINTER(ct.c_float))
        self.ncol = objects['ncol']
        self.nrow = objects['nrow']

        return self.ncol, self.nrow

    def __init_types_before_nstar__(self):
        """Defines argument types before running nstar"""
        self.pynstartable = self.daolib.pynstartable
        self.pynstartable.argtypes=[
            ct.c_int, ct.c_int, ct.POINTER(ct.c_float), ct.c_int,  # nrow, ncol, data_ptr, maxcol
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x1, y1
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x2, y2
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x3, y3
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x4, y4
            ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),      # fratio12, fratio13
            ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),      # fratio14, z0
            ct.POINTER(ct.c_double), ct.c_float, ct.c_float,       # chi2, watch, fitrad
            ct.c_float, ct.c_float,                                # e1, e2
            ct.POINTER(ct.c_char), ct.c_int,                       # psf, psfl
            ct.POINTER(ct.c_char), ct.c_int,                       # grp, grpl,
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),            # xmin, xmax
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),            # ymin, ymax
            ct.c_double, ct.POINTER(ct.c_double),                  # rfac, residuals,
            ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),      # sigmas, model
            ct.c_int, ct.POINTER(ct.c_float),                      # number of stars, zpmag Zero-point magnitude
            ct.POINTER(ct.c_float)]                                # optional params

        self.pynstar = self.daolib.pynstar
        self.pynstar.argtypes=[
            ct.c_int, ct.c_int, ct.POINTER(ct.c_float), ct.c_int,  # nrow, ncol, data_ptr, maxcol
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x1, y1
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x2, y2
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x3, y3
            ct.POINTER(ct.c_float), ct.POINTER(ct.c_float),        # x4, y4
            ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),      # fratio12, fratio13
            ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),      # fratio14, z0
            ct.POINTER(ct.c_double), ct.c_float, ct.c_float,       # chi2, watch, fitrad
            ct.c_float, ct.c_float,                                # e1, e2
            ct.POINTER(ct.c_char), ct.c_int,                       # psf, psfl
            ct.POINTER(ct.c_char), ct.c_int,                       # grp, grpl,
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),            # xmin, xmax
            ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),            # ymin, ymax
            ct.c_double,                                           # rfac
            ct.c_int, ct.POINTER(ct.c_float)]                      # number of stars, zpmag Zero-point magnitude

    def nstar_goodness_of_fit(self, theta):
        """Run NSTAR once at given position and flux ratio and return the chi2 and total flux."""
        
        nrow, ncol = self.nrow, self.ncol
        data_ptr = self.picture_ptr
        maxcol = self.maxcol_daophot
        opt_params_ptr = self.opt_params_ptr 

        tof_psf = self.__assign_chararray__(u'{0}'.format(self.psf_name))
        tof_psf_ptr = tof_psf.ctypes.data_as(ct.POINTER(ct.c_char))
        tof_grp = self.__assign_chararray__(u'{0}'.format(self.group_name))
        tof_grp_ptr = tof_grp.ctypes.data_as(ct.POINTER(ct.c_char))
        tof_psfl = ct.c_int(len(self.psf_name))
        tof_grpl = ct.c_int(len(self.group_name))

        xfort3 = ct.c_float(0)
        yfort3 = ct.c_float(0)
        fratio13 = ct.c_double(0)
        xfort4 = ct.c_float(0)
        yfort4 = ct.c_float(0)
        fratio14 = ct.c_double(0)

        self.zpmag = ct.c_float(0)

        if isinstance(theta, list):
            nb_params = len(theta)
        elif isinstance(theta, np.ndarray):
            nb_params = theta.shape[0]
        else:
            sys.exit("Error in nstarwrap: do not recognize input parameter type.")

        if nb_params == 5:  # 2-star fit, 5 parameters
            x1, y1, x2, y2, fratio = theta
            xfort1 = ct.c_float(x1)
            yfort1 = ct.c_float(y1)
            xfort2 = ct.c_float(x2)
            yfort2 = ct.c_float(y2)
            fratio12 = ct.c_double(fratio)
            fit_stars = ct.c_int(2)  # number of stars
        elif nb_params == 8:  # 3-star fit, 8 parameters
            x1, y1, x2, y2, x3, y3, fratio12, fratio13 = theta
            xfort1 = ct.c_float(x1)
            yfort1 = ct.c_float(y1)
            xfort2 = ct.c_float(x2)
            yfort2 = ct.c_float(y2)
            fratio12 = ct.c_double(fratio12)
            xfort3 = ct.c_float(x3)
            yfort3 = ct.c_float(y3)
            fratio13 = ct.c_double(fratio13)
            fit_stars = ct.c_int(3)  # number of stars
        elif nb_params-1 == 11:  # 4-star fit, 11 parameters
            x1,y1,x2,y2,x3,y3,x4,y4,fratio12,fratio13,fratio14, k = theta
            xfort1 = ct.c_float(x1)
            yfort1 = ct.c_float(y1)
            xfort2 = ct.c_float(x2)
            yfort2 = ct.c_float(y2)
            fratio12 = ct.c_double(fratio12)
            xfort3 = ct.c_float(x3)
            yfort3 = ct.c_float(y3)
            fratio13 = ct.c_double(fratio13)
            xfort4 = ct.c_float(x4)
            yfort4 = ct.c_float(y4)
            fratio14 = ct.c_double(fratio14)
            fit_stars = ct.c_int(4)  # number of stars

        chi2 = ct.c_double(-1.0)
        z0 = ct.c_double(-1.0)

        tof_watch = self.opt['watch']
        tof_fitrad = self.opt['fitrad']
        tof_e1 = self.opt['e1']
        tof_e2 = self.opt['e2']
        tof_rfac = self.opt['err_factor']

        self.pynstar(nrow, ncol, data_ptr, maxcol, 
            ct.byref(xfort1), ct.byref(yfort1),
            ct.byref(xfort2), ct.byref(yfort2),
            ct.byref(xfort3), ct.byref(yfort3),
            ct.byref(xfort4), ct.byref(yfort4),
            ct.byref(fratio12), ct.byref(fratio13),
            ct.byref(fratio14), ct.byref(z0), ct.byref(chi2),
            tof_watch, tof_fitrad, tof_e1, tof_e2,
            tof_psf_ptr, tof_psfl,
            tof_grp_ptr, tof_grpl,
            ct.byref(self.opt['box_xmin']), ct.byref(self.opt['box_xmax']),
            ct.byref(self.opt['box_ymin']), ct.byref(self.opt['box_ymax']),
            tof_rfac, fit_stars, ct.byref(self.zpmag), opt_params_ptr)

        self.dof = (self.opt['box_xmax'].value - self.opt['box_xmin'].value + 1)\
            * (self.opt['box_ymax'].value - self.opt['box_ymin'].value + 1)\
            - 3 * fit_stars.value

        return chi2.value, z0.value

    def nstartable_goodness_of_fit(self, theta):
        """Run NSTAR once at given position and flux ratio and return the chi2 and total flux."""
        
        nrow, ncol = self.nrow, self.ncol
        data_ptr = self.picture_ptr
        maxcol = self.maxcol_daophot
        residuals_ptr = self.residuals_ptr
        sigmas_ptr = self.sigmas_ptr
        model_ptr = self.model_ptr

        tof_psf = self.__assign_chararray__(u'{0}'.format(self.psf_name))
        tof_psf_ptr = tof_psf.ctypes.data_as(ct.POINTER(ct.c_char))
        tof_grp = self.__assign_chararray__(u'{0}'.format(self.group_name))
        tof_grp_ptr = tof_grp.ctypes.data_as(ct.POINTER(ct.c_char))
        tof_psfl = ct.c_int(len(self.psf_name))
        tof_grpl = ct.c_int(len(self.group_name))

        xfort3 = ct.c_float(0)
        yfort3 = ct.c_float(0)
        fratio13 = ct.c_double(0)
        xfort4 = ct.c_float(0)
        yfort4 = ct.c_float(0)
        fratio14 = ct.c_double(0)

        self.zpmag = ct.c_float(0)

        if isinstance(theta, list):
            nb_params = len(theta)
        elif isinstance(theta, np.ndarray):
            nb_params = theta.shape[0]
        else:
            sys.exit("Error in nstarwrap: do not recognize input parameter type.")

        if nb_params == 5:  # 2-star fit, 5 parameters
            x1, y1, x2, y2, fratio = theta
            xfort1 = ct.c_float(x1)
            yfort1 = ct.c_float(y1)
            xfort2 = ct.c_float(x2)
            yfort2 = ct.c_float(y2)
            fratio12 = ct.c_double(fratio)
            fit_stars = ct.c_int(2)  # number of stars
        elif nb_params == 8:  # 3-star fit, 8 parameters
            x1, y1, x2, y2, x3, y3, fratio12, fratio13 = theta
            xfort1 = ct.c_float(x1)
            yfort1 = ct.c_float(y1)
            xfort2 = ct.c_float(x2)
            yfort2 = ct.c_float(y2)
            fratio12 = ct.c_double(fratio12)
            xfort3 = ct.c_float(x3)
            yfort3 = ct.c_float(y3)
            fratio13 = ct.c_double(fratio13)
            fit_stars = ct.c_int(3)  # number of stars
        elif nb_params-1 == 11:  # 4-star fit, 11 parameters
            x1,y1,x2,y2,x3,y3,x4,y4,fratio12,fratio13,fratio14, k = theta
            xfort1 = ct.c_float(x1)
            yfort1 = ct.c_float(y1)
            xfort2 = ct.c_float(x2)
            yfort2 = ct.c_float(y2)
            fratio12 = ct.c_double(fratio12)
            xfort3 = ct.c_float(x3)
            yfort3 = ct.c_float(y3)
            fratio13 = ct.c_double(fratio13)
            xfort4 = ct.c_float(x4)
            yfort4 = ct.c_float(y4)
            fratio14 = ct.c_double(fratio14)
            fit_stars = ct.c_int(4)  # number of stars

        chi2 = ct.c_double(-1.0)
        z0 = ct.c_double(-1.0)

        tof_watch = self.opt['watch']
        tof_fitrad = self.opt['fitrad']
        tof_e1 = self.opt['e1']
        tof_e2 = self.opt['e2']
        tof_rfac = self.opt['err_factor']

        self.pynstar(nrow, ncol, data_ptr, maxcol, 
            ct.byref(xfort1), ct.byref(yfort1),
            ct.byref(xfort2), ct.byref(yfort2),
            ct.byref(xfort3), ct.byref(yfort3),
            ct.byref(xfort4), ct.byref(yfort4),
            ct.byref(fratio12), ct.byref(fratio13),
            ct.byref(fratio14), ct.byref(z0), ct.byref(chi2),
            tof_watch, tof_fitrad, tof_e1, tof_e2,
            tof_psf_ptr, tof_psfl,
            tof_grp_ptr, tof_grpl,
            ct.byref(self.opt['box_xmin']), ct.byref(self.opt['box_xmax']),
            ct.byref(self.opt['box_ymin']), ct.byref(self.opt['box_ymax']),
            tof_rfac, residuals_ptr, sigmas_ptr, model_ptr, fit_stars, 
            ct.byref(self.zpmag))

        self.dof = (self.opt['box_xmax'].value - self.opt['box_xmin'].value + 1)\
            * (self.opt['box_ymax'].value - self.opt['box_ymin'].value + 1)\
            - 3 * fit_stars.value

        return chi2.value, z0.value
