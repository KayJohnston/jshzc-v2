# JSHZC v2
# Jackie Silver / Kay Johnston

# Things to remember to add:
# Error Checking
# Remove spurious implied accuracy from the sig fig rounding

# Some way of greying out stars which are not present.
# Add a scale.
# Add the full mendax spectrum guesser as an option.
# Better looking sliders.
# Graphic for icy, better graphic for WW.
# Maybe increase font size and rejig graphics canvas size to match.

from tkinter import *
import PIL
from PIL import ImageTk, Image, ImageDraw, ImageFont
import math

class App():

    def __init__(self,master):
        # Lists to hold the primary and secondary stars.
        self.pstars = []
        self.sstars = []
        
        # Create an frame which will hold the star frames with entry boxes and calculated values.
        self.entry_box_frame = Frame(master)
        self.entry_box_frame.pack()

        # Create a frame which will hold the distance entry box and calculated temperatures for that distance.
        self.distance = Frame(master)
        self.distance.pack()

        # Create a frame which will hold the gas giant distance and mass entry boxes and guess at its temperature and type.
        self.gasgiant = Frame(master)
        self.gasgiant.pack()

        # Create a frame which will hold general controls.
        self.controls = Frame(master)
        self.controls.pack()

        # Create a frame which will show the system image.
        self.system = Frame(master)
        self.system.pack()

        # Create a set of star frames.
        self.sf_list = []
        for new_star_frame in range(0,6):
            self.sf_list.append(Star_Frame(self,self.entry_box_frame))
        # Add default values for Sol to the first star frame.
        self.sf_list[0].R_entrybox.entry.delete(0,'end')
        self.sf_list[0].R_entrybox.entry.insert(0,'1')
        self.sf_list[0].T_entrybox.entry.delete(0,'end')
        self.sf_list[0].T_entrybox.entry.insert(0,'5778')

        self.sf_list[0].update_calcs('dummy')

        # Disable distance entry box on the primary star.
        self.sf_list[0].D_entrybox.entry['state']='disabled'

        # Create an entry box for the temperature at given distance calculator.
        self.d_temp = Entry_Box(self.distance,'Distance to check (au)','1',20,8)

        # Create a label to show the calculated temperatures at a given distance.
        self.t_at_d = StringVar()
        self.t_at_d.set('Black Body: ??? General: ??? Icy: ???')
        self.t_at_d_label = Label(self.distance, textvariable = self.t_at_d, width = 68)
        self.t_at_d_label.pack()

        # Bind the temperature at distance entry box to automatic update.
        self.d_temp.entry.bind('<Return>', self.auto_calculate)

        # Create entry boxes for the mass and distance of a gas giant.
        # Default values are for approximately Jupiter's position and mass.
        self.gg_distance = Entry_Box(self.gasgiant,'Gas Giant distance (au)','5.2',20,8)
        self.gg_mass = Entry_Box(self.gasgiant,'Gas Giant Mass (ME)','318',20,8)

        # Create a label to show the calculated temperature and type of the gas giant.
        self.gg_tandt = StringVar()
        self.gg_tandt.set('Temperature: ??? Sudarsky Class: ???')
        self.gg_tandt_label = Label(self.gasgiant, textvariable = self.gg_tandt, width = 40)
        self.gg_tandt_label.pack()

        # Bind the gas giant entry boxes to the automatic update.
        self.gg_distance.entry.bind('<Return>', self.auto_calculate)
        self.gg_mass.entry.bind('<Return>', self.auto_calculate)

        # Create entry boxes for the start and end of the distance range.
        self.d_start = Entry_Box(self.controls,'Scale Distance start (au)','0',20,8)
        self.d_end = Entry_Box(self.controls,'Scale Distance end (au)','10',20,8)

        # Bind the distance range entry boxes to the automatic update.
        self.d_start.entry.bind('<Return>', self.auto_calculate)
        self.d_end.entry.bind('<Return>', self.auto_calculate)

        # Load in a font.
        self.fnt = ImageFont.truetype('Quicksand-Regular.otf', FONTSIZE)

        # Load in planet icons.
        self.pic_mrw = Image.open('mrw.png')
        self.pic_mhz = Image.open('mhz.png')
        self.pic_mhz2 = Image.open('mhz2.png')
        self.pic_elw = Image.open('elw.png')
        self.pic_ww = Image.open('ww.png')
        self.pic_amw = Image.open('amw.png')
        self.pic_icy = Image.open('mhz.png')

        # Create a canvas to show the system image.
        self.sys_canvas = Canvas(self.system, width = XDIM, height = YDIM)
        self.sys_canvas.pack()

        # Updates the calculated results on screen by running the auto calculate function.
        self.auto_calculate('dummy_value')

    def auto_calculate(self,A):
        # In this case A is unused.

        # Clear the star lists before we start reading them in.
        self.pstars = []
        self.sstars = []

        # Go through each star frame, reading in the values and creating a star object.
        for sf in self.sf_list:
            R = sf.R_entrybox.entry.get()
            T = sf.T_entrybox.entry.get()
            D = sf.D_entrybox.entry.get()
            newstar = Star(R,T,D)
            if newstar.D == '0':
                if newstar.R != '0':
                    self.pstars.append(newstar)
            else:
                self.sstars.append(newstar)

        # Read the distance for temperature at distance calculations.
        D = self.d_temp.entry.get()
        D = eval(D)
        D = float(D)
        t_bb = snowliner5(self.pstars,self.sstars,D,0)
        t_gen = snowliner5(self.pstars,self.sstars,D,0.35)
        t_icy = snowliner5(self.pstars,self.sstars,D,0.6)
        
        # ED won't let the temperature drop below 20 K in practice, so catch that.
        if t_gen < 20:
            t_gen = 20
        if t_icy < 20:
            t_icy = 20

        # Add a radius note for the primary.
        primary_radius_au = self.pstars[0].Rau
        primary_radius_au = round_to_n(primary_radius_au,3)
        primary_radius_ls = self.pstars[0].Rls
        primary_radius_ls = round_to_n(primary_radius_ls,3)
        p_rad_text = '     Primary radius: ' + str(primary_radius_au) + ' au / ' + str(primary_radius_ls) + ' ls'

        # Text for the temperature at a distance results.
        text = 'Black Body: ' + str(int(t_bb)) + 'K  General: ' + str(int(t_gen)) + 'K  Icy: ' + str(int(t_icy)) + ' K'

        final = text + p_rad_text
        
        self.t_at_d.set(final)

        # Read the distance and mass for the gas giant calculations.
        D = self.gg_distance.entry.get()
        M = self.gg_mass.entry.get()
        D = eval(D)
        D = float(D)
        M = eval(M)
        M = float(M)
        Tbase = snowliner5(self.pstars,self.sstars,D,0.35)
        Tinc = gginc(M)
        T = Tbase + Tinc
        C = 'I'
        if T > cIbound:
            C = 'II'
        if T > cIIbound:
            C = 'III'
        if T > cIIIbound:
            C = 'IV'
        if T > cIVbound:
            C = 'V'
        text = 'Temperature: ' + str(int(T)) + ' K  Sudarsky Class: ' + C
        self.gg_tandt.set(text)

        # Calculate MRW upper boundary.
        self.mrw_inner = 0 # Change this to account for solar radius?
        self.mrw_innerls = 0
        self.mrw_outer = dfortm2('1103',self.pstars,self.sstars)
        self.mrw_outer = round_to_n(self.mrw_outer,3)
        self.mrw_outerls = int(self.mrw_outer * auls)
        
        # Calculate main habitable zone boundaries.
        self.mhz_inner = dfortm2('315',self.pstars,self.sstars)
        self.mhz_inner = round_to_n(self.mhz_inner,3)
        self.mhz_innerls = int(self.mhz_inner * auls)
        self.mhz_outer = dfortm2('223',self.pstars,self.sstars)
        self.mhz_outer = round_to_n(self.mhz_outer,3)
        self.mhz_outerls = int(self.mhz_outer * auls)
        
        # Calculate ELW boundaries.
        self.elw_inner = dfortm2('278',self.pstars,self.sstars)
        self.elw_inner = round_to_n(self.elw_inner,3)
        self.elw_innerls = int(self.elw_inner * auls)
        self.elw_outer = dfortm2('227',self.pstars,self.sstars)
        self.elw_outer = round_to_n(self.elw_outer,3)
        self.elw_outerls = int(self.elw_outer * auls)
        
        # Calculate WW boundaries.
        self.ww_inner = dfortm2('307',self.pstars,self.sstars)
        self.ww_inner = round_to_n(self.ww_inner,3)
        self.ww_innerls = int(self.ww_inner * auls)
        self.ww_outer = dfortm2('156',self.pstars,self.sstars)
        self.ww_outer = round_to_n(self.ww_outer,3)
        self.ww_outerls = int(self.ww_outer * auls)
        
        # Calculate AMW boundaries.
        self.amw_inner = dfortm2('193',self.pstars,self.sstars)
        self.amw_inner = round_to_n(self.amw_inner,3)
        self.amw_innerls = int(self.amw_inner * auls)
        self.amw_outer = dfortm2('117',self.pstars,self.sstars)
        self.amw_outer = round_to_n(self.amw_outer,3)
        self.amw_outerls = int(self.amw_outer * auls)

        # Calculate approximate Icy boundary.
        self.icy_inner = dfortm2('150',self.pstars,self.sstars)
        self.icy_inner = round_to_n(self.icy_inner,3)
        self.icy_innerls = int(self.icy_inner * auls)
        self.icy_outer = 100000
        self.icy_outerls = (100000 * auls)
        
        # Call a function which updates the image of the system.
        self.update_image()

    def update_image(self):
        # Create a new image in PIL.
        self.pil_image = Image.new('RGBA',(XDIM,YDIM),'black')
        self.draw = ImageDraw.Draw(self.pil_image)

        # Get scale start and end points.
        scale_start = eval(self.d_start.entry.get())
        scale_end = eval(self.d_end.entry.get())
        
        # Draw groups for each planet type.
        self.drawgroup('MRW',self.pic_mrw,10,10,self.mrw_inner,self.mrw_outer,scale_start,scale_end)
        self.drawgroup('MHZ',self.pic_mhz2,10,60,self.mhz_inner,self.mhz_outer,scale_start,scale_end)
        self.drawgroup('ELW',self.pic_elw,10,110,self.elw_inner,self.elw_outer,scale_start,scale_end)
        self.drawgroup('WW',self.pic_ww,10,160,self.ww_inner,self.ww_outer,scale_start,scale_end)
        self.drawgroup('AMW',self.pic_amw,10,210,self.amw_inner,self.amw_outer,scale_start,scale_end)
        self.drawgroup('ICY',self.pic_icy,10,260,self.icy_inner,self.icy_outer,scale_start,scale_end)

        # Draw the scale.
        self.drawscale(10,310,scale_start,scale_end)

        # Convert the image to one that tkinter can use, and draw it to the canvas.
        self.working_image = ImageTk.PhotoImage(self.pil_image)
        self.image_on_canvas = self.sys_canvas.create_image(0, 0, anchor = NW, image = self.working_image)

    def drawscale(self,sx,sy,scale_start,scale_end):
        # Axis line.
        lx_s = sx + 240
        lx_e = sx + 640
        ly_s = ly_e = sy
        self.draw.line(((lx_s,ly_s),(lx_e,ly_e)), width = 2, fill = (255,255,255,255))

        # Draw ticks.
        tick_length = 100
        scale_tick_length_au = (scale_end - scale_start) / 4
##        scale_tick_length_ls = (scale_end - scale_start) * (499 / 4)
        for tick in range(0,5):
            tx_s = tx_e = sx + 240 + (tick * tick_length)
            ty_s = sy
            ty_e = sy + 8
            self.draw.line(((tx_s,ty_s),(tx_e,ty_e)), width = 2, fill = (255,255,255,255))
            D = scale_start + (tick * scale_tick_length_au)
            D = str(D) + ' au'
            self.draw.text((tx_s - 14,ty_e + 2),D,font = self.fnt,fill = (255,255,255,255))
##            D = scale_start + (tick * scale_tick_length_ls)
##            D = str(D) + ' ls'
##            self.draw.text((tx_s - 18,ty_e + 16),D,font = self.fnt,fill = (255,255,255,255))
        

    def drawgroup(self,name,picname,sx,sy,inner,outer,scale_start,scale_end):
        # Paste in a picture of the planet.
        self.pil_image.paste(picname,(sx + 40,sy))

        # Text for the name of this group.
        nx = sx
        ny = sy + 14
        self.draw.text((nx,ny),name,font = self.fnt,fill = (255,255,255,255))

        # Text for the inner boundary, in au.
        innerau = str(round(inner,3)) + ' au'
        ix = sx + 90
        iy = sy + 8
        self.draw.text((ix,iy),innerau,font = self.fnt,fill = (255,255,255,255))

        # Text for the outer boundary, in au.
        if outer < 100000:
            outerau = str(round(outer,3)) + ' au'
        else:
            outerau = 'unlimited'
        ox = sx + 150
        oy = sy + 8
        self.draw.text((ox,oy),outerau,font = self.fnt,fill = (255,255,255,255))

        # Text for the inner boundary, in ls.
        innerls = str(int(inner * auls)) + ' ls'
        ix = sx + 90
        iy = sy + 22
        self.draw.text((ix,iy),innerls,font = self.fnt,fill = (255,255,255,255))

        # Text for the outer boundary, in ls.
        if outer < 100000:
            outerls = str(int(outer * auls)) + ' ls'
        else:
            outerls = 'unlimited'
        ox = sx + 150
        oy = sy + 22
        self.draw.text((ox,oy),outerls,font = self.fnt,fill = (255,255,255,255))

        # Reference line.
        lx_s = sx + 240
        lx_e = sx + 640
        ly_s = ly_e = sy + 20
        self.draw.line(((lx_s,ly_s),(lx_e,ly_e)), width = 17, fill = (31,31,31,255))

        # Actual line.
        length = lx_e - lx_s
        scale_length = scale_end - scale_start
        ratio = length / scale_length
        s_inner = inner - scale_start
        s_outer = outer - scale_start
        ax_s = lx_s + (s_inner * ratio)
        if ax_s < lx_s:
            ax_s = lx_s
        elif ax_s > lx_e:
            ax_s = lx_e
        ax_e = lx_s + (s_outer * ratio)
        if ax_e < lx_s:
            ax_e = lx_s
        elif ax_e > lx_e:
            ax_e = lx_e
        ay_s = ay_e = ly_s
        self.draw.line(((ax_s,ay_s),(ax_e,ay_e)), width = 17, fill = (255,255,255,255))


# Holds entry boxes and calculated data labels for one star.
class Star_Frame(App):

    def __init__(self,parent,master):

        # Create a frame to hold the entry box and calculated label frames.
        self.frame = Frame(master)
        self.frame.pack(fill = 'x')

        # Entry box frame.
        self.entrybox_frame = Frame(self.frame)
        self.entrybox_frame.pack(side = LEFT)

        # Calculated label frame.
        self.calc_frame = Frame(self.frame)
        self.calc_frame.pack()

        # Make entry boxes.
        self.R_entrybox = Entry_Box(self.entrybox_frame, 'Radius (Rsol):','0',10,7)
        self.T_entrybox = Entry_Box(self.entrybox_frame, 'Temp. (K):','0',8,7)
        self.D_entrybox = Entry_Box(self.entrybox_frame, 'Distance (SMA) (au):','0',16,7)

        # Bind entry box key press of return to call the parent's auto calculate function.
        self.R_entrybox.entry.bind('<Return>', parent.auto_calculate, add = '+') # This allows it to call two functions from one keypress.
        self.T_entrybox.entry.bind('<Return>', parent.auto_calculate, add = '+')
        self.D_entrybox.entry.bind('<Return>', parent.auto_calculate, add = '+')

        # Bind entry box key press of return to call internal update function.
        self.R_entrybox.entry.bind('<Return>', self.update_calcs, add = '+')
        self.T_entrybox.entry.bind('<Return>', self.update_calcs, add = '+')
        self.D_entrybox.entry.bind('<Return>', self.update_calcs, add = '+')

        # Make label boxes to show luminosity and spectral class.
        self.calc_luminosity = StringVar()
        self.calc_luminosity.set('Luminosity (Lsol): ???')
        self.luminosity_label = Label(self.calc_frame, textvariable = self.calc_luminosity, width = 20, anchor = 'w')
        self.luminosity_label.pack(side = LEFT)

        self.calc_spectral = StringVar()
        self.calc_spectral.set('Spectral Class: ???')
        self.spectral_label = Label(self.calc_frame, textvariable = self.calc_spectral,width = 20, anchor = 'w')
        self.spectral_label.pack(side = LEFT)

        # Call our internal update function.
        self.update_calcs('dummy')

    # Read the entry boxes and update the calculated labels.
    def update_calcs(self,dummy):
        R = self.R_entrybox.entry.get()
        T = self.T_entrybox.entry.get()
        L = luminosity(R,T)
        S = spectral(R,T)
        self.calc_luminosity.set('Luminosity (Lsol): ' + str(round_to_n(L,3)))
        self.calc_spectral.set('Spectral Class: ' + str(S))


# Entry boxes with an attached label.  
class Entry_Box():

    def __init__(self,master,nametext,default,w1,w2):
        # Create a frame for this entry box.
        self.frame = Frame(master)
        self.frame.pack(side = LEFT)

        # Create a label.
        self.label = Label(self.frame,text = nametext,width = w1)
        self.label.pack(side = LEFT)

        # Create an entry box.
        self.entry = Entry(self.frame, width = w2)
        self.entry.pack(side = LEFT)
        self.entry.insert(0,default)

# A class for spectral type boundaries.
class Boundary():

    def __init__(self,name,low,high):
        self.name = name
        self.low = low
        self.high = high

# Generic class for stars, could expand this with additional data.
class Star():

    def __init__(self,R,T,D):
        self.R = R
        self.T = T
        self.L = luminosity(self.R,self.T)
        self.S = spectral(self.R,self.T)
        self.D = D # Semi-major axis.  Only for secondary stars.
        self.Rau, self.Rls = radiusinls(self.R)

    def report(self):
        print('R:',self.R,'solar radius (about',self.Rls,'ls) T:',self.T,'L:',self.L,'Spectral:',self.S)

# Function to calculate relative luminosity
# Given ED's usual terms of R in solar radii and T in Kelvins.
def luminosity(R,T):
    R = eval(R)
    T = eval(T)
    R = float(R)
    T = int(T)
    Tnormal = T / Tsol
    L = (R ** 2) * (Tnormal ** 4)
    return L

# Function to calculate basic spectral class.
def spectral(R,T):
    R = eval(R)
    T = eval(T)
    R = float(R)
    T = int(T)
    result = 'No match found.'
    for boundary in boundary_list:
        if T >= boundary.low and T <= boundary.high:
            subdivs = (boundary.high - boundary.low) / 10
            steps = 10 - (T - boundary.low) / subdivs
            steps = int(steps)
            if steps == 10:
                steps = 9
            result = boundary.name + str(steps)
    # Check for giant stars.
    if R > 20:
        result += ' (giant)'
    # Check for neutron stars and black holes.
    if T == 0:
        result = 'None'
    elif T >= 70000:
        result = 'Neutron Star'
    return result

# Returns the predicted temperature increase for a gas giant of a given mass.
# Not perfect but good enough to be useful.
def gginc(mass):
    coeff = 0.02
    power = 1.3
    Tinc = coeff * (mass**power)
    return Tinc

def radiusinls(R):
    R = eval(R)
    R = float(R)
    Rau = (R * Rsol) / au
##    Rau = round_to_n(Rau,2)
    Rls = (R * Rsol * auls) / au
##    Rls = round_to_n(Rls,2)
    return Rau, Rls

# Function for rounding to n sig figs.
# http://stackoverflow.com/questions/3410976/how-to-round-a-number-to-significant-figures-in-python
def round_to_n(x, n):
    if not x: return 0
    power = -int(math.floor(math.log10(abs(x)))) + (n - 1)
    factor = (10 ** power)
    return round(x * factor) / factor
    

### Function that gets the distance required for a particular black body temperature.
### Now superceded by the dfortmulti function.
##def dfort(Tbb,R,T):
##    Tbb = eval(Tbb)
##    R = eval(R)
##    T = eval(T)
##    Tbb = int(Tbb)
##    R = float(R)
##    T = int(T)
##    R *= Rsol
##    top = R * (T ** 2)
##    bottom = 2 * (Tbb ** 2)
##    d = top / bottom
##    d /= au
##    return d

### Function that gets the distance required for a particular black body temperature, for multiple stars.
### Currently using primary stars only (i.e. stars that the body is directly orbiting.)
### Now superceded by the dfortm2 function, although I have my doubts.  Not least as to how to spell superceded.
##def dfortmulti(Tbb,pstars):
##    Tbb = eval(Tbb)
##    top = 0
##    for star in pstars:
##        R = eval(star.R)
##        T = eval(star.T)
##        R = float(R)
##        T = int(T)
##        R *= Rsol
##        top += (R ** 2) * (T ** 4)
##    bottom = 4 * (Tbb ** 4)
##    calc = top / bottom
##    distance = (calc ** 0.5)
##    distance /= au
##    return distance

# Function that gets the distance required for a particular black body temperature.
# For multiple stars, per suspected ED mechanics.  Untested as yet.  Don't think the bloody thing is working.
def dfortm2(Tbb,pstars,sstars):
    Tbb = eval(Tbb)
    # Calculate for secondary stars.
    # Distance to the secondaries is invariant; calculate their temperature contribution as a constant.
    RHS = 0
    for star in sstars:
        R = eval(star.R)
        T = eval(star.T)
        D = eval(star.D) # Semi-major axis distance.
        R = float(R)
        T = int(T)
        D = float(D)
        R *= Rsol
        D *= au
        top = (R ** 2) * (T ** 4)
        bottom = 4 * (D ** 2)
        RHS += (top / bottom)
    Tsec = (RHS ** 0.25)

    Tbb = Tbb - Tsec # Could go negative, catch.  And it's weird anyway.  I don't trust this bit.
    
    # Calculate for primary stars.
    top = 0
    for star in pstars:
        R = eval(star.R)
        T = eval(star.T)
        R = float(R)
        T = int(T)
        R *= Rsol
        top += (R ** 2) * (T ** 4)
    bottom = 4 * (Tbb ** 4)
    calc = top / bottom
    distance = (calc ** 0.5)
    distance /= au
    
    return distance

def snowliner5(pstars,sstars,D,albedo):
    # PRIMARY stars are calculated properly.
    # SECONDARY stars are calculated properly together.
    # OVERALL temperature is calculated summatively.
    
    # PRIMARY calculation
    # -------------------
    RHS = 0
    for star in pstars:
        R = eval(star.R)
        R = float(R)
        T = eval(star.T)
        T = int(T)
        R *= Rsol
        D = D * au
        RHS += ((R ** 2) * (T ** 4)) / (4 * (D ** 2))
    RHS *= (1 - albedo)
    tempP = (RHS ** 0.25)
    
    # SECONDARY calculation
    # ---------------------
    RHS = 0
    for star in sstars:
        R = eval(star.R)
        R = float(R)
        T = eval(star.T)
        T = int(T)
        D = eval(star.D) * au
        RHS += ((R ** 2) * (T ** 4)) / (4 * (D ** 2))
    RHS *= (1 - albedo)
    tempS = (RHS ** 0.25)

    # OVERALL sum and return.
    # -----------
    result = tempS + tempP

    return int(result)

# Constants
Tsol = 5778 # Kelvin
Rsol = 696000000 # Metres
au = 1.496e11 # Metres
auls = 499 # Light seconds
cIVbound = 1400 # Kelvin
cIIIbound = 800 # Kelvin
cIIbound = 250 # Kelvin
cIbound = 150 # Kelvin

# Spectral type temperature boundaries.
boundary_list = []
boundary_list.append(Boundary('O',33000,69999.9))
boundary_list.append(Boundary('B',10000,32999.9))
boundary_list.append(Boundary('A',7500,9999.9))
boundary_list.append(Boundary('F',6000,7499.9))
boundary_list.append(Boundary('G',5200,5999.9))
boundary_list.append(Boundary('K',3700,5199.9))
boundary_list.append(Boundary('M',2000,3699.9))
boundary_list.append(Boundary('L',1300,1999.9))
boundary_list.append(Boundary('T',700,1299.9))
boundary_list.append(Boundary('Y',1,699.9))

# Handy variables
XDIM, YDIM = 690, 340 # Size of the main drawing canvas.
FONTSIZE = 12

# Main loop
root = Tk()
root.title('Jackie\'s Hab-Zone Calculator v2')

mainapp = App(root)

root.mainloop()
