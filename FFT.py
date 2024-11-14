import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets 
from scipy.fft import fft, ifft, fftfreq
from scipy import signal as sg


# Defining pi for ease and the inital values for the variables
pi = np.pi
points = 100
phase = 0
time_tot = 1
amp = 1
form = 'sin'

# Creating the figures and axes
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(8, 3))
fig.subplots_adjust(hspace=0.4)
fig.subplots_adjust(bottom=0.45)

# Initial definition of time values, t
t = np.linspace(0, 1, points)



# This function creates the input signal at the desired variables
# The conditional statements in this function allow to produce sine, square, or sawtooth
def signal(t, f, phase, amp, noisescale, form):
    noise = noisescale * np.random.normal(0, 1, len(t)) # this creates the noise signal based on the 'noisescale' amplitude
    if form == 'sin':
        return amp * np.sin(2 * pi * f * t + np.deg2rad(phase)) + noise
    elif form == 'square':
        return amp * sg.square(2 * pi * f * t + np.deg2rad(phase)) + noise
    elif form == 'sawtooth':
        return amp * sg.sawtooth(2 * pi * f * t + np.deg2rad(phase)) + noise

# Plotting the signal with the initial values specified in the question
funcPlot, = ax1.plot(t, signal(t, 10, 0, 1, 0, form), color='b')
ax1.set_title('Signal')

# This function outputs the power spectrum of the signal for the positive frequency range
def power(t,f, phase, points, time_tot, amp, noisescale, form):
    freqs = fftfreq(len(t), (time_tot)/points) # this calculates the frequency range produced
    posfreqs = freqs[freqs >= 0] # this line only takes the positive frequency range
    func = signal(t, f, phase, amp, noisescale, form) # calling the signal function
    Ft = fft(func) # calculate the Fourier transform
    spectrum = (np.abs(Ft)**2)/ points**2 # calculate the power spectrum and normalise
    return spectrum[freqs >= 0], posfreqs # only returns the spectrum for the positive frequencies so the sizes are the same with posfreqs

# Claculate the power spectrum and frequencies and plot them
powers, posfreqs = power(t, 10, 0, points, time_tot, amp, 0, form)
powerPlot, = ax2.plot(posfreqs, powers, color='r')
ax2.set_title('Power Spectrum')


# Below are many lines of slider and button initialisation so the comments refer to the attribute each slider alters

# Frequency slider
freqAx = plt.axes([0.12, 0.3, 0.25, 0.03]) 
freqSlider = widgets.Slider(freqAx, 'Frequency (Hz)', 0, 100, valinit=10, valstep=0.1)

# Number of points slider
pointAx = plt.axes([0.12,0.25, 0.25, 0.03])
pointSlider = widgets.Slider(pointAx, r'N$^{o}$ of points', 1, 200, valinit=100, valstep=1)

# Time slider
timeAx = plt.axes([0.12, 0.2, 0.25, 0.03])
timeSlider = widgets.Slider(timeAx, 'Time (s)', 0, 1, valinit=1, valstep=0.01)

# Time windowing sliders
# tll = time limit low 
tllAx = plt.axes([0.12, 0.15, 0.25, 0.03])
tllSlider = widgets.Slider(tllAx, 'T low', 0, 1, valinit=0, valstep=0.01) 
#tlh = time limit high
tlhAx = plt.axes([0.12, 0.1, 0.25, 0.03])
tlhSlider = widgets.Slider(tlhAx, 'T high', 0, 1, valinit=1, valstep=0.01)

# Phase slider
phaseAx = plt.axes([0.12, 0.05, 0.25, 0.03])
phaseSlider = widgets.Slider(phaseAx, r'Phase ($^{o})$', 0, 360, valinit=0, valstep=1)

# Amplitude slider
amplitudeAx = plt.axes([0.5, 0.3, 0.25, 0.03])
amplitudeSlider = widgets.Slider((amplitudeAx), 'Amplitude', 0, 1, valinit=1, valstep=0.05)

# Noise slider
noiseAx = plt.axes([0.5, 0.25, 0.25, 0.03])
noiseSlider = widgets.Slider(noiseAx, 'Noise', 0, 1, valinit=0, valstep=0.05)

# Frequency windowing sliders
# fll = freq limit low
fllAx = plt.axes([0.5,0.2,0.25,0.03])
fllSlider = widgets.Slider(fllAx, 'f low', 0, 100, valinit=0, valstep=1)
# flh = freq limit high
flhAx = plt.axes([0.5,0.15,0.25,0.03])
flhSlider = widgets.Slider(flhAx, 'f high', 0, 100, valinit=50, valstep=1)

# Form button (ie what kind of wave is being produced)
formAx = plt.axes([0.8, 0.25, 0.1, 0.1])
formButton = widgets.RadioButtons(formAx, ('sin','square','sawtooth'))

# Close button (self explanatory)
closeAx = plt.axes([0.8, 0.1, 0.1, 0.1])
closeButton = widgets.Button(closeAx, 'close')

# This function calculates the inverse fourier transform within the specified frequency windowing to reconstruct the signal
def ift(t,f, phase, amp, fll, flh, points, time_tot, noisescale, form):
    func = signal(t, f, phase, amp, noisescale, form)
    Ft = fft(func) # find the fft again but seperately from the power function
    freqs = fftfreq(len(t), (time_tot)/points) # similarly find the frequency range
    # The following two lines allow for windowing in frequency space
    # The absolute value is taken to keep sizes of arrays the same as only using positive side would make these arrays half as large
    # note: the frequencies are set to zero and not completely thrown out once again to keep array sizes the same
    Ft[np.abs(freqs) < fll] = 0 # set frequencies below fll to zero
    Ft[np.abs(freqs) > flh] = 0 # set frequencies above flh to zero
    iFt = ifft(Ft) # evaluate the ifft over our desired range
    return iFt.real # only take the real part of the ifft

# Plot the inverse ft with the initial values
iftPlot, = ax3.plot(t, ift(t, 10, 0, 1 , 0 ,50, 100, 1, 0 ,form), color='g')
ax3.set_title('Inverse FT')

# Plots the vertical lines for the time windowing
tll_line = ax1.axvline(tllSlider.val, color='gray', linestyle='--')
tlh_line = ax1.axvline(tlhSlider.val, color='gray', linestyle='--')

# Plots the vertical lines for the frequency windowing
fll_line = ax2.axvline(fllSlider.val, color='gray', linestyle='--')
flh_line = ax2.axvline(flhSlider.val, color='gray', linestyle='--')

# This is the main function of the code that updates all the plots whenever anything is changed
def update(val):
    # Defining the variables as the values of the sliders
    points = int(pointSlider.val)
    freq = freqSlider.val
    phase = phaseSlider.val
    tll = tllSlider.val
    tlh = tlhSlider.val
    fll = fllSlider.val
    flh = flhSlider.val
    time_tot = timeSlider.val
    amp = amplitudeSlider.val
    noisescale = noiseSlider.val
    form = formButton.value_selected
    
    # This logic here just makes sure the time windowing and the big time truncation agree
    if tlh - tll > time_tot:
        tlh = tll + time_tot
        
    # Definition of time range and full time length to make the big time slider useful
    t = np.linspace(tll, tlh, points)
    tfull = np.linspace(0, time_tot, points)
    
    # Frequency range definition and positive freq selection as before
    freqs = fftfreq(len(t), time_tot/points)
    posfreqs = freqs[freqs >= 0]
    
    # Defining the signal that will be used next
    sgnl = signal(t, freq, phase, amp, noisescale, form)
    
    # Signal plot
    funcPlot.set_ydata(sgnl) 
    funcPlot.set_xdata(t)
    
    # Power spectrum calculation and Power spectrum plot
    powers, posfreqs = power(tfull, freq, phase, points, time_tot, amp, noisescale, form)
    powerPlot.set_ydata(powers)
    powerPlot.set_xdata(posfreqs)
    
    # Inverse ft plot
    iftPlot.set_data(t, ift(t, freq, phase, amp, fll, flh, points, time_tot, noisescale, form))
    
    # Time windowing lines
    tll_line.set_xdata(tll)
    tlh_line.set_xdata(tlh)
    
    # Frequency space windowing lines
    fll_line.set_xdata(fll)
    flh_line.set_xdata(flh)
    
    # Re draw everything
    plt.draw()  
    
# Makes the close button actually close the gui
def close(event):
    plt.close('all')    

# All these lines call the update function whenever any slider or button is pressed
freqSlider.on_changed(update) 
pointSlider.on_changed(update)
timeSlider.on_changed(update)   
tllSlider.on_changed(update)
tlhSlider.on_changed(update)
fllSlider.on_changed(update)
flhSlider.on_changed(update)
phaseSlider.on_changed(update)
amplitudeSlider.on_changed(update)
noiseSlider.on_changed(update)
formButton.on_clicked(update)
closeButton.on_clicked(close) # except for this one which closes the gui