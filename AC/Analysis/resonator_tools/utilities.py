import warnings
import numpy as np
import matplotlib.pyplot as plt


def Watt2dBm(x):
    '''
    converts from units of watts to dBm
    '''
    return 10.*np.log10(x*1000.)
    
def dBm2Watt(x):
    '''
    converts from units of watts to dBm
    '''
    return 10**(x/10.) /1000.   
    
class plotting(object):
    '''
    some helper functions for plotting
    Normalisation can be achieved by: 
        S21_vect_norm = np.vectorize(complex)(S21_real, S21_imag)/resonator.fitresults["a"]
    '''
##
## z_data_raw denotes the raw data
## z_data denotes the normalized data
## 
#   def plotall(self):
#       plotallfig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = [10,2])
#       real = self.z_data_raw.real
#       imag = self.z_data_raw.imag
#       real2 = self.z_data_sim.real
#       imag2 = self.z_data_sim.imag
        

#       ax1.plot(real,imag,label='rawdata')
#       ax1.plot(real2,imag2,label='fit')
#       ax1.set(xlabel = 'Re(S21)', ylabel = 'Im(S21)')
#       ax1.legend()
        

#       ax2.plot(self.f_data*1e-9,np.absolute(self.z_data_raw),label='rawdata')
#       ax2.plot(self.f_data*1e-9,np.absolute(self.z_data_sim),label='fit')
#       ax2.set (xlabel = 'f (GHz)', ylabel = '|S21|')
#       ax2.legend()

#       ax3.plot(self.f_data*1e-9,np.angle(self.z_data_raw),label='rawdata')
#       ax3.plot(self.f_data*1e-9,np.angle(self.z_data_sim),label='fit')
#       ax3.set(xlabel = 'f (GHz)', ylabel = 'arg(|S21|)')
#       ax3.legend()
#       plotallfig.suptitle("Raw data/ fit plot")
#       plt.show()

    def plotall(self, title= '', S21_or_S11 = 'S21', fontsize_user = 12):
        #print(self.z_data_raw)
        real = self.z_data_raw.real
        imag = self.z_data_raw.imag
        real2 = self.z_data_sim.real
        imag2 = self.z_data_sim.imag
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3)
        plt.subplot(221)
        plt.suptitle(title)
        plt.plot(real,imag,label='raw data')
        plt.plot(real2,imag2,label='fit')
        if S21_or_S11 == 'S21':
            plt.xlabel('Re(S$_{21}$)', fontsize = fontsize_user)
            plt.ylabel('Im(S$_{21}$)', fontsize = fontsize_user)
        else:
            plt.xlabel('Re(S$_{11}$)', fontsize = fontsize_user)
            plt.ylabel('Im(S$_{11}$)', fontsize = fontsize_user)
        plt.legend()
        plt.subplot(222)
        plt.plot(self.f_data*1e-9,np.absolute(self.z_data_raw),label='raw data')
        plt.plot(self.f_data*1e-9,np.absolute(self.z_data_sim),label='fit')
        plt.xlabel('f (GHz)')
        if S21_or_S11 == 'S21':
            plt.ylabel(r'|S$_{21}$|', fontsize = fontsize_user)
        else:
            plt.ylabel(r'|S$_{11}$|', fontsize = fontsize_user)
        plt.legend()
        plt.subplot(223)
        plt.plot(self.f_data*1e-9,np.angle(self.z_data_raw),label='raw data')
        plt.plot(self.f_data*1e-9,np.angle(self.z_data_sim),label='fit')
        plt.xlabel('f (GHz)')
        if S21_or_S11 == 'S21':
            plt.ylabel(r'arg(|S$_{21}$|)', fontsize = fontsize_user)
        else:
            plt.ylabel(r'arg(|S$_{11}$|)', fontsize = fontsize_user)
        plt.legend()
        plt.subplot(224)
        plt.plot(self.f_data*1e-9,np.unwrap(np.angle(self.z_data_raw)),label='raw data')
        plt.plot(self.f_data*1e-9,np.unwrap(np.angle(self.z_data_sim)),label='fit')
        plt.xlabel('f (GHz)')
        if S21_or_S11 == 'S21':
            plt.ylabel(r'Unwrapped arg(|S$_{21}$|)', fontsize = fontsize_user)
        else:
            plt.ylabel(r'Unwrapped arg(|S$_{11}$|)', fontsize = fontsize_user)
        plt.legend()
        plt.show()
        plt.tight_layout()
        
        
        
    def plotfit(self, title= '', location = 'best', fontsize_user = 12, color_fit = 'g', linestyle_fit ='-', S21_or_S11 = 'S21', max_xticks = 0, max_yticks = 0):
        real = self.z_data_raw.real
        imag = self.z_data_raw.imag
        real2 = self.z_data_sim.real
        imag2 = self.z_data_sim.imag
        plt.plot(self.f_data*1e-9,np.absolute(self.z_data_raw),label='Experiment', color = 'tab:blue')
        plt.plot(self.f_data*1e-9,np.absolute(self.z_data_sim),label='Fit', color = color_fit, linestyle = linestyle_fit)
        plt.xlabel('f (GHz)', fontsize = fontsize_user)
        if S21_or_S11 == 'S21':
            plt.ylabel(r'|S$_{21}$|', fontsize = fontsize_user)
        else:
            plt.ylabel(r'|S$_{11}$|', fontsize = fontsize_user)
        plt.xticks(fontsize= fontsize_user)
        plt.yticks(fontsize= fontsize_user)
        plt.legend(loc = location, fontsize = fontsize_user)
        #plt.title("Magnitude fit")
        plt.title(title)
        if max_xticks:
            plt.locator_params(axis='x', nbins=max_xticks)
        if max_yticks:
            plt.locator_params(axis='y', nbins=max_yticks)
        plt.tight_layout()
        plt.show()

    def plotfitcircle(self, title= '', location = 'best', fontsize_user = 12, color_fit = 'g', linestyle_fit ='-', S21_or_S11 = 'S21', max_xticks = 0, max_yticks = 0):
        real = self.z_data_raw.real
        imag = self.z_data_raw.imag
        real2 = self.z_data_sim.real
        imag2 = self.z_data_sim.imag
        plt.plot(real,imag,label='Experiment', color = 'tab:blue')
        plt.plot(real2,imag2,label='Fit', color = color_fit, linestyle = linestyle_fit)
        if S21_or_S11 == 'S21':
            plt.xlabel('Re(S$_{21}$)', fontsize = fontsize_user)
            plt.ylabel('Im(S$_{21}$)', fontsize = fontsize_user)
        else:
            plt.xlabel('Re(S$_{11}$)', fontsize = fontsize_user)
            plt.ylabel('Im(S$_{11}$)', fontsize = fontsize_user)
        plt.xticks(fontsize= fontsize_user)
        plt.yticks(fontsize= fontsize_user)
        plt.legend(loc = location, fontsize = fontsize_user)
        plt.title(title)
        if max_xticks:
            plt.locator_params(axis='x', nbins=max_xticks)
        if max_yticks:
            plt.locator_params(axis='y', nbins=max_yticks)
        plt.plot(real2,imag2,label='Fit', color = 'tab:green')
        plt.xlabel('Re(S21)')
        plt.ylabel('Im(S21)')
        plt.legend(loc = location)
        plt.title(title)
        plt.tight_layout()
        plt.show()


    def plotcalibrateddata(self):
        real = self.z_data.real
        imag = self.z_data.imag
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3)
        plt.subplot(221)
        plt.plot(real,imag,label='calib')
        plt.xlabel('Re(S21)')
        plt.ylabel('Im(S21)')
        plt.legend()
        plt.subplot(222)
        plt.plot(self.f_data*1e-9,np.absolute(self.z_data),label='calib')
        plt.xlabel('f (GHz)')
        plt.ylabel('|S21|')
        plt.legend()
        plt.subplot(223)
        plt.plot(self.f_data*1e-9,np.angle(self.z_data),label='calib')
        plt.xlabel('f (GHz)')
        plt.ylabel('arg(|S21|)')
        plt.legend()
        plt.show()
        
    def plotrawdata(self):
        real = self.z_data_raw.real
        imag = self.z_data_raw.imag
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3)
        plt.subplot(221)
        plt.plot(real,imag,label='rawdata')
        plt.xlabel('Re(S21)')
        plt.ylabel('Im(S21)')
        plt.legend()
        plt.subplot(222)
        plt.plot(self.f_data*1e-9,np.absolute(self.z_data_raw),label='raw data')
        plt.xlabel('f (GHz)')
        plt.ylabel('|S21|')
        plt.legend()
        plt.subplot(223)
        plt.plot(self.f_data*1e-9,np.angle(self.z_data_raw),label='raw data')
        plt.xlabel('f (GHz)')
        plt.ylabel('arg(|S21|)')
        plt.legend()
        plt.show()
        
    def plot_fine(self, value_low, value_high):
        if value_low > value_high:
            value_high = valu
            value_low = vals
            value_low = np.copy(valu)
            value_high = np.copy(vals)
        real = self.z_data_raw.real
        imag = self.z_data_raw.imag
        real2 = self.z_data_sim.real
        imag2 = self.z_data_sim.imag
        
        plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.3)
        plt.subplot(221)
        cutoff_lw = np.argwhere(np.array(self.f_data)>value_low)[0][0]
        cutoff_hgh = np.argwhere(np.array(self.f_data)<value_high)[:]
        
        print(self.f_data[cutoff_lw], self.f_data[cutoff_hgh[len(cutoff_hgh) - 1][0]])
        plt.plot(real[cutoff_lw:cutoff_hgh[len(cutoff_hgh) - 1][0]],imag[cutoff_lw:cutoff_hgh[len(cutoff_hgh) - 1][0]],label='rawdata')
        plt.plot(real2[cutoff_lw:cutoff_hgh[len(cutoff_hgh) - 1][0]],imag2[cutoff_lw:cutoff_hgh[len(cutoff_hgh) - 1][0]],label='fit')
        plt.xlabel('Re(S21)')
        plt.ylabel('Im(S21)')
        plt.legend()
        plt.subplot(222)
  
        plt.plot(self.f_data*1e-3,np.absolute(self.z_data_raw),label='raw data')
        plt.plot(self.f_data*1e-3,np.absolute(self.z_data_sim),label='fit')
        plt.xlabel('f (kHz)')
        plt.xlim(value_low*1e-3, value_high*1e-3)
        plt.ylabel('|S21|')
        plt.legend()
        plt.subplot(223)
    
        plt.plot(self.f_data*1e-3,np.angle(self.z_data_raw),label='raw data')
        plt.plot(self.f_data*1e-3,np.angle(self.z_data_sim),label='fit')
        plt.xlabel('f (kHz)')
        plt.xlim(value_low*1e-3, value_high*1e-3)
        plt.ylabel('arg(|S21|)')
        plt.legend()
        plt.show()

class save_load(object):
    '''
    procedures for loading and saving data used by other classes
    '''
    def _ConvToCompl(self,x,y,dtype):
        '''
        dtype = 'realimag', 'dBmagphaserad', 'linmagphaserad', 'dBmagphasedeg', 'linmagphasedeg'
        '''
        if dtype=='realimag':
            return x+1j*y
        elif dtype=='linmagphaserad':
            return x*np.exp(1j*y)
        elif dtype=='dBmagphaserad':
            return 10**(x/20.)*np.exp(1j*y)
        elif dtype=='linmagphasedeg':
            return x*np.exp(1j*y/180.*np.pi)
        elif dtype=='dBmagphasedeg':
            return 10**(x/20.)*np.exp(1j*y/180.*np.pi)   
        else: warnings.warn("Undefined input type! Use 'realimag', 'dBmagphaserad', 'linmagphaserad', 'dBmagphasedeg' or 'linmagphasedeg'.", SyntaxWarning)
    
    def add_data(self,f_data,z_data):
        self.f_data = np.array(f_data)
        self.z_data_raw = np.array(z_data)
        
    def cut_data(self,f1,f2):
        def findpos(f_data,val):
            pos = 0
            for i in range(len(f_data)):
                if f_data[i]<val: pos=i
            return pos
        pos1 = findpos(self.f_data,f1)
        pos2 = findpos(self.f_data,f2)
        self.f_data = self.f_data[pos1:pos2]
        self.z_data_raw = self.z_data_raw[pos1:pos2]
        
    def add_fromtxt(self,fname,dtype,header_rows,usecols=(0,1,2),fdata_unit=1.,delimiter=None):
        '''
        dtype = 'realimag', 'dBmagphaserad', 'linmagphaserad', 'dBmagphasedeg', 'linmagphasedeg'
        '''
        data = np.loadtxt(fname,usecols=usecols,skiprows=header_rows,delimiter=delimiter)
        self.f_data = data[:,0]*fdata_unit
        self.z_data_raw = self._ConvToCompl(data[:,1],data[:,2],dtype=dtype)
        
    def add_fromhdf():
        pass
    
    def add_froms2p(self,fname,y1_col,y2_col,dtype,fdata_unit=1.,delimiter=None):
        '''
        dtype = 'realimag', 'dBmagphaserad', 'linmagphaserad', 'dBmagphasedeg', 'linmagphasedeg'
        '''
        if dtype == 'dBmagphasedeg' or dtype == 'linmagphasedeg':
            phase_conversion = 1./180.*np.pi
        else: 
            phase_conversion = 1.
        f = open(fname)
        lines = f.readlines()
        f.close()
        z_data_raw = []
        f_data = []
        if dtype=='realimag':
            for line in lines:
                if ((line!="\n") and (line[0]!="#") and (line[0]!="!")) :
                    lineinfo = line.split(delimiter)
                    f_data.append(float(lineinfo[0])*fdata_unit)
                    z_data_raw.append(np.complex(float(lineinfo[y1_col]),float(lineinfo[y2_col])))
        elif dtype=='linmagphaserad' or dtype=='linmagphasedeg':
            for line in lines:
                if ((line!="\n") and (line[0]!="#") and (line[0]!="!") and (line[0]!="M") and (line[0]!="P")):
                    lineinfo = line.split(delimiter)
                    f_data.append(float(lineinfo[0])*fdata_unit)
                    z_data_raw.append(float(lineinfo[y1_col])*np.exp( np.complex(0.,phase_conversion*float(lineinfo[y2_col]))))
        elif dtype=='dBmagphaserad' or dtype=='dBmagphasedeg':
            for line in lines:
                if ((line!="\n") and (line[0]!="#") and (line[0]!="!") and (line[0]!="M") and (line[0]!="P")):
                    lineinfo = line.split(delimiter)
                    f_data.append(float(lineinfo[0])*fdata_unit)
                    linamp = 10**(float(lineinfo[y1_col])/20.)
                    z_data_raw.append(linamp*np.exp( np.complex(0.,phase_conversion*float(lineinfo[y2_col]))))
        else:
            warnings.warn("Undefined input type! Use 'realimag', 'dBmagphaserad', 'linmagphaserad', 'dBmagphasedeg' or 'linmagphasedeg'.", SyntaxWarning)
        self.f_data = np.array(f_data)
        self.z_data_raw = np.array(z_data_raw)
        
    def save_fitresults(self,fname):
        pass
    


