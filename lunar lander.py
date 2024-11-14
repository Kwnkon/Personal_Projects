import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle, Polygon
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider

# Create the figures
fig, (ax1,ax2) = plt.subplots(1,2, figsize = (6,8))
plt.subplots_adjust(0.1,0.25)

# Manually set the position of the subplots to make it visible
ax1.set_position([0.1, 0.2, 0.55, 0.7])
ax2.set_position([0.7, 0.2, 0.25, 0.7])

t = np.linspace(0,10,100) # Array of times for 10s
ax_slider = plt.axes([0.25, 0.1, 0.65, 0.03])
slider = Slider(ax_slider, 'thruster', 0, 3, valinit=0, valstep=0.1)

 
#Starting values
h0 = 100 # height
v0 = 0 # velocity
g = 1 # gravity
f0 = 100 # fuel

# The following line took me too long to figure out considering how obvious it is
dt = t[1] - t[0] # time step

# Create the axes for this game
ax1.set_xlim(0,10)
ax1.set_ylim(0,h0)
ax1.set_xlabel('ground')
ax1.set_ylabel('height / m')


ax2.set_xlim(0, 1) 
ax2.set_ylim(0, 100.1)  # fuel ranges from 0 to the initial amount
ax2.set_title('Fuel Gauge')

# Make a cute rocket shape

rocket_body = Rectangle((4.5, 0), 1, 4, color='gray', zorder=5)  # Rocket body
rocket_tip = Polygon([[4.5, 4], [5.5, 4], [5, 5]], color='red', zorder=6)  # Rocket tip
flame = Polygon([[4.5, 0], [5.5, 0], [5, -2]], color='orange', zorder=4, visible=False)  # Flame shape
ax1.add_patch(rocket_body)
ax1.add_patch(rocket_tip)
ax1.add_patch(flame)


# Little tools that will help us later
#           rocket, = ax1.plot([], [], 'ro', label='rocket') OLD WAY of making my 

f_level, = ax2.plot([], [])
txt1 = ax1.text(0.75, 0.9, '')
txt3 = ax1.text(6, 0.9, '')
#txt4 = ax1.text(0.75, 5, 'fuel: 100') this text box is uneccesary now

txt5 = plt.text(-1, 0.4, 'Instructions: Press Spacebar to use maximum thruster')
txt6 = plt.text(-0.3, 0.35, 'Land with Velocity < 2 to land safely')
txt7 = plt.text(-0.3, 0.3, 'Good Luck!')
# don't worry about txt 2 missing above, it is dynamically changed in the 
# update function (bc otherwise it will obscure the actual game)

# initial values are initial values (shocker!)
v = v0
h = h0
f = f0

# these two functions will allow us to use the spacebar to control the thruster
def key_on(event):
    global thruster
    
    if event.key == ' ':
        thruster = 3
    else:
        thruster = slider.val
        
def key_off(event):
    global thruster
    
    if event.key == ' ':
        thruster = 0
    else:
        thruster = slider.val
        
  #this actually connects the above function to the figures      
fig.canvas.mpl_connect('key_press_event', key_on)
fig.canvas.mpl_connect('key_release_event', key_off)



# this update function actually makes the FuncAnimation work
def update(frame):
    global h,v,f
    #thruster = slider.val
    
    
    # this if loop checks whether we have any fuel left so we can alter our acceleration
    if f > 0:
        a = g - thruster
        f = f -  5 * thruster * dt
        
        # this shows cool flames when thruster is on
        if thruster > 0:
            flame.set_visible(True)
        else:
            flame.set_visible(False)
            
            
        # txt4.set_text('fuel: '+str(round(f,1)))   no longer necessary as fuel gauge is here
        
    elif f < 1:
       # txt4.set_text('Out of Fuel')   same as above
        a = g
        flame.set_visible(False)

        
    # dynamically changing our v and h values
    v = v + a * dt
    h = h - v * dt
    # Make sure that we are above the ground or tell the user whether they crashed or landed
    if h > 0:
        txt1.set_text('Height: '+str(round(h, 1)))
        txt3.set_text('Velocity: '+str(round(v, 1)))
    elif h < 1 and v > 2:
        h = 0
        txt1.set_text('')
        txt2 = ax1.text(0.4, 50, '', fontsize = 50)
        txt2.set_text('CRASH!!!')
        ani.event_source.stop()
    elif h < 1 and v < 2:
        h = h0
        txt1.set_text('')
        txt2 = ax1.text(2.5, 50, '', fontsize = 50)
        txt2.set_text('Landed!')
        ani.event_source.stop()
        
        
    
    # Set the rocket position
    rocket_body.set_y(h)  
    rocket_tip.set_xy([[4.5, h + 4], [5.5, h + 4], [5, h + 5]])
    
    # Update flame position
    flame.set_xy([[4.5, h], [5.5, h], [5, h - 2]])  # Move flames with the rocket
    
    #Actually set the position of the rocket and fuel gauge
    #           rocket.set_data(5, h)  OLD WAY of moving rocket
    f_level.set_data(np.array([0,1]), f * np.ones(2))
    
    
    return rocket_body, rocket_tip, flame, f_level


ani = FuncAnimation(fig, update, frames=len(t), interval=100)

plt.show()


