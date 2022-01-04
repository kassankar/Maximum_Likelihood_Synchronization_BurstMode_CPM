#%%
import numpy as np
from scipy import integrate, special
import matplotlib.pyplot as plt
from CPM_MAIN_ import MAIN
import random
# %matplotlib widget
# %%EE
# Variables ini
Main          = MAIN()  # Call class MAIN
pulse         = 2;      # 1 -> lorentzian pulse
                        # 2 -> GMSK pulse BT = 0.3
                        # 3 -> LRC pulse
                        # 4 -> LREC pulse
L             = 4;      # Pulse length
                        # 1  -> Full response
                        # >1 -> Partial response
OS            = 2**1;   # Over sampling frequency
Ts            = 1/OS;   # Sampling Time (not symbole time, symbole time is always one)
M             = 2**1;   # M_ary symbols used (2 -> Binary)
h             = 0.5;    # Modulation index
width         = 0.8;    # This variable is used for Lorentzian Pulse only. (Not be used for pulse > 1)
# %%
# frequency pulse
g_t, q_t, fig = Main.CREATECPMPULSE(pulse,L,width,OS,1,'-') # Function return the CPM pulse and phase.
                                                            # g_t = g(t) is the CPM pulse shape.
                                                            # q_t -> is the phase, integral of g_t.
                                                            # fig = figure of g_t
fig
# %%
# Training sequence Phase and Signal
antipodal        = 1  # True->(-1,1), Flase->(0,1)
L0               = 64 # sequence length
# Only when antipodal = 0.
R0               = 2  # R0 number of starting zeroes 
R1               = 6  # R1 number of ones
R2               = 8  # R2 number of zeros at the end
if antipodal !=True:
    Tr           = np.hstack((np.zeros((1,R0)),np.ones((1,R1)),np.zeros((1,R2)))).flatten()
else:
    Tr           = np.hstack((-(M-1)*np.ones((1,int(L0/4))),(M-1)*np.ones((1,int(L0/2))),-(M-1)*np.ones((1,int(L0/4))))).flatten()
Ph, Mod_s, t_seq, fig   = Main.Data_Modulation(pulse,g_t,h,L,OS,Tr,2,'-')

# %%
# Add Tl
if antipodal !=True:
    l             = (L-1)/2
    Tr1           = np.hstack((np.zeros((1,R0)),np.ones((1,R1)),np.zeros((1,R2)))).flatten()
else:
    if L>1:
        Tl            = int((L-1)/2)
    else:
        Tl            = 0
    Tr1                  = np.hstack((-(M-1)*np.ones((1,int(L0/4))),(M-1)*np.ones((1,int(L0/2))),-(M-1)*np.ones((1,int(L0/4))), -(M-1)*np.ones((1,Tl)))).flatten()
Ph, Mod_s, t_seq1, fig   = Main.Data_Modulation(pulse,g_t,h,L,OS,Tr1,0,'--')

# Remove lag Tl
if L>1 :
    Mod_s_lag = Mod_s[Tl*OS+1:]
    Ph_lag    = Ph[Tl*OS+1:]
    t_seq     = t_seq[0:np.size(Mod_s_lag,0)]
else:
    Mod_s_lag = Mod_s
    Ph_lag    = Ph
# %%
# Add channel Noise and offset
# Add time offset
Snr      = np.arange(0,11,1).flatten()
itt      = 10**5
fe       = np.zeros((1,itt)).flatten()
fe_int   = np.zeros((1,itt)).flatten()
te       = np.zeros((1,itt)).flatten()
pe       = np.zeros((1,itt)).flatten()


nfev_int = np.zeros((1,np.size(Snr))).flatten()
ntev     = np.zeros((1,np.size(Snr))).flatten()
npev     = np.zeros((1,np.size(Snr))).flatten()

for I in range(0,np.size(Snr)):
    for J in range(0,itt): 
        tau_min  = -OS/2
        tau_max  = OS/2-1
        k        = random.randint(tau_min, tau_max) # Delay tau (time offset)
        if k>0:
            r    =  np.hstack((np.zeros((1,k)).flatten(),Mod_s_lag)).flatten()
            r    = r[0:np.size(r,0)-k]
        else:
            r    =  np.hstack((Mod_s_lag, np.zeros((1,abs(k))).flatten())).flatten()
            r    = r[abs(k):]
        # Add freq offset
        fd_min = -OS/2
        fd_max = OS/2-0.1
        fd     = random.uniform(fd_min, fd_max) #  fd (frequency offset)
        r      = r*np.exp(1j*2*np.pi*fd*t_seq);

        # Add phase offset theta (\pn 1 to avoid 2*pi rotation)
        theta_min = 0+1
        theta_max = 2*np.pi-1
        theta     = random.uniform(theta_min, theta_max) # theta (phase offset)
        r         = r*np.exp(1j*theta)

        # Add Noise
        snr = Snr[I]
        r   = Main.awgn(r,snr,OS)

        # Compute Lambda1 and Lambda2
        if L==1:
            r1 = np.zeros((1,L0*OS+1),dtype=complex).flatten()
            r2 = np.zeros((1,L0*OS+1),dtype=complex).flatten()
        else:
            r1 = np.zeros((1,L0*OS),dtype=complex).flatten()
            r2 = np.zeros((1,L0*OS),dtype=complex).flatten()

        r1[0:int(L0*OS/4)]                  = r[0:int(L0*OS/4)]
        r1[int(3*L0*OS/4):int(L0*OS)]       = r[int(3*L0*OS/4):int(L0*OS)]*np.exp(-1j*(M-1)*np.pi*h*L0)
        r2[int(L0*OS/4):int(3*L0*OS/4)]     = r[int(L0*OS/4):int(3*L0*OS/4)]*np.exp(1j*(M-1)*np.pi*h*L0/2)

        r1p = r1*np.exp(1j*(M-1)*np.pi*h*t_seq)
        r2p = r2*np.exp(-1j*(M-1)*np.pi*h*t_seq)

        kf         = 2**1
        Nf         = kf*OS*L0
        r1p_fft    = np.fft.fftshift(np.fft.fft(r1p,int(Nf)))
        r2p_fft    = np.fft.fftshift(np.fft.fft(r2p,int(Nf)))
        f          = np.fft.fftshift(np.fft.fftfreq(Nf,Ts))
        X          = np.abs(r1p_fft) + np.abs(r2p_fft)
        idx_max    = np.argmax(X)
        fd_est     = f[idx_max]
        v_1        = idx_max-1
        v0         = idx_max
        v1         = idx_max+1
        if v_1< 0 or v1 >= np.size(f,0):
            fe[J]         =  0
            fe_int[J]     =  0
            te[J]         =  0
            pe[J]         =  0
        else:
            fd_est_int    =  fd_est + 1/(2*L0*kf)*(np.log(X[v_1])-np.log(X[v1])) / \
                        (np.log(X[v_1])+ np.log(X[v1]) - 2*np.log(X[v0]))
            fe[J]         =  (fd_est - fd)**2
            fe_int[J]     =  (fd_est_int-fd)**2
            lambda1       = np.sum(r1p*np.exp(-1j*2*np.pi*(fd_est_int)*t_seq))*Ts
            lambda2       = np.sum(r2p*np.exp(-1j*2*np.pi*(fd_est_int)*t_seq))*Ts
            eps_est       = np.angle(lambda1*np.conj(lambda2))/(2*(M-1)*np.pi*h)
            te[J]         = (eps_est - k/OS)**2
            theta_est     = np.angle(np.exp(-1j*(M-1)*np.pi*h*eps_est)*lambda1 +
                                    np.exp(1j*(M-1)*np.pi*h*eps_est)*lambda2)
            theta_est     = (theta_est) % (2 * np.pi)
            pe[J]         = (theta_est-theta)**2
    nfev_int[I]   = (fe_int[fe_int!=0]).mean()
    ntev[I]       = (te[te != 0]).mean()
    npev[I]       = (pe[pe != 0]).mean()

# Plot
fig = plt.figure(3)
plt.plot(Snr,nfev_int)
plt.yscale("log")
plt.xlabel("Eb/N0")
plt.ylabel("CRB $f_d$")
fig = plt.figure(4)
plt.plot(Snr,ntev)
plt.yscale("log")
plt.xlabel("Eb/N0")
plt.ylabel("CRB ($ \\tau $)")
fig = plt.figure(5)
plt.plot(Snr,npev)
plt.yscale("log")
plt.xlabel("Eb/N0")
plt.ylabel("CRB ($ \\theta $)")
