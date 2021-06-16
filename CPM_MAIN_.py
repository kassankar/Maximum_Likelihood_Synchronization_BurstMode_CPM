class MAIN:
    import numpy as np
    from scipy import integrate, special, signal
    import matplotlib.pyplot as plt
    def qfunc(self,x):
        special = self.special
        np      = self.np 
        return 0.5*special.erfc(x/np.sqrt(2))
    def CREATECPMPULSE(self,pulse,pulse_length,pulse_width,os,img,plot_line_style):
        np        = self.np 
        integrate = self.integrate
        plt       = self.plt
        qfunc     = self.qfunc
        Ts =  1/os
        w  = pulse_width
        if pulse==1 :      
            t          = np.arange(-pulse_length/2,pulse_length/2+Ts,Ts) 
            t0         = 0
            g_t        = (2*w)/((t-t0)**2+w**2)            
            g_t        = g_t*(1/(4*np.pi))
            Cst        = np.sum([g_t])*Ts
            nug_t      = (0.5) / Cst
            g_t        = nug_t*g_t
            q_t        = integrate.cumtrapz(g_t)*Ts
        
        elif pulse==2 :
            t          = np.arange(-pulse_length/2,pulse_length/2+Ts,Ts) 
            Bt         = 0.3
            alpha      = 2*np.pi*Bt/(np.sqrt(np.log(2)))
            gauss      = qfunc(alpha*(t-0.5)) - qfunc(alpha*(t+0.5))
            Cst        = 0.5/(np.sum(gauss)*Ts)
            g_t        = Cst*gauss
            q_t        = integrate.cumtrapz(g_t)*Ts
                            
        elif pulse==3 :
            t          = np.arange(0,pulse_length+Ts,Ts)
            g_t        = (1/(2*pulse_length)*(1- np.cos(2*np.pi*t/(pulse_length))))
            K          = 0.5/(np.sum(g_t)*Ts)
            g_t        = K*g_t
            q_t        = integrate.cumtrapz(g_t)*Ts
        elif pulse==4 :
            t          = np.arange(0,pulse_length+Ts,Ts)
            g_t        = 1/(2*pulse_length)*np.ones((1,np.size(t,0)))
            g_t[0,0]   = 0
            g_t[0,-1]  = 0
            g_t        = g_t.flatten()
            K          = 0.5/(np.sum(g_t)*Ts)
            g_t        = K*g_t
            q_t        = integrate.cumtrapz(g_t)*Ts


        if img > 0:
            fig = plt.figure(img)
            if pulse==1:
                plt.plot(t,4*np.pi*g_t, linestyle = plot_line_style)
            else:
                plt.plot(t,g_t, linestyle = plot_line_style)
                
            plt.xlabel('time')
            plt.ylabel('Freqeuncy pulse g(t)')
            plt.grid()
        else:
            fig=0
        return g_t, q_t, fig
    def Data_Modulation(self,pulse,g_t,modulation_index,pulse_length,os,bits,img, plot_line_style):
        signal    = self.signal
        np        = self.np
        integrate = self.integrate
        plt       = self.plt
        h         = modulation_index
        L         = pulse_length
        Ts        = 1/os

        bits_s = np.hstack((bits[:,None], np.zeros((np.size(bits,0),os-1)))).flatten() # same as matlab upsample (insert zeros (nb:os-1) after each element of the array)
        t_seq  = np.arange(0,np.size(bits,0)+Ts,Ts)
        s      = np.convolve(bits_s,g_t)
        s      = s[0:np.size(t_seq,0)]
        Ph     = integrate.cumtrapz(s,initial=0)*Ts
        if pulse==1 :
            Mod_s  = np.exp(1j*2*2*np.pi*h*Ph)
        else:
            Mod_s  = np.exp(1j*2*np.pi*h*Ph)
        if img > 0:
            fig = plt.figure(img,figsize=(15,5))
            plt.subplot(1,2,1)
            plt.plot(t_seq,np.unwrap(np.angle(Mod_s)), linestyle = plot_line_style)
            plt.xlabel('time t/(T_s)')
            plt.ylabel('Phase \phi(t;\alpha)')
            plt.grid(True)
            plt.subplot(1,2,2)
            # plt.plot(t_seq[0:-1-(np.size(bits,0)-L)*os],Mod_s)
            plt.plot(t_seq,Mod_s, linestyle = plot_line_style)
            plt.xlabel('time t/(T_s)')
            plt.ylabel('Modulated Signal s(t,\alpha)')
            plt.grid(True)
            fig.tight_layout(pad=3.0)
        else:
            fig=0
        return Ph, Mod_s, t_seq, fig
    def awgn(self,s,SNRdB,os=1):
        """
        AWGN channel
        Add AWGN noise to input signal. The function adds AWGN noise vector to signal 's' to generate a resulting signal vector 'r' of specified SNR in dB. It also
        returns the noise vector 'n' that is added to the signal 's' and the power spectral density N0 of noise added
        Parameters:
            s : input/transmitted signal vector
            SNRdB : desired signal to noise ratio (expressed in dB) for the received signal
            os : oversampling factor (applicable for waveform simulation) default L = 1.
        Returns:
            r : received signal vector (r=s+n)
        """
        np = self.np
        gamma = 10**(SNRdB/10) #SNR to linear scale
        if s.ndim==1:# if s is single dimensional vector
            P=os*np.sum(np.abs(s)**2)/len(s) #Actual power in the vector
        else: # multi-dimensional signals like MFSK
            P=os*np.sum(np.sum(abs(s)**2))/len(s) # if s is a matrix [MxN]
        N0=P/gamma # Find the noise spectral density
        if np.isrealobj(s):# check if input is real/complex object type
            n = np.sqrt(N0/2)*np.random.standard_normal(s.shape) # computed noise
        else:
            n = np.sqrt(N0/2)*(np.random.standard_normal(s.shape)+1j*np.random.standard_normal(s.shape))
        r = s + n # received signal
        return r
   
  