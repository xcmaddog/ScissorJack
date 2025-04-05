# import packages here
using SNOW
using Plots


# Definition of structs for components of jack
"""
Diagonal member
This is the struct object for the diagonal members. It assumes that they are all U-bars

t is the thickness of the diagonal member in inches
w is the width of the diagonal member in inches
S_y is the yield strength of the diagonal member in psi
ppi is the price per inch of the diagonal member in USD/inch
"""
struct DiagonalMember
    name::String
    t::Float32
    w::Float32
    S_y::Float32
    ppi::Float32
end

"""
Crossbar
This is the type for the crossbars.
"""
abstract type Crossbar end

"""
CrossbarSquare
This is the struct for the hollow square crossbars. It assumes that they are all square bars.

t is the thickness of the crossbar in inches
w is the width of the crossbar in inches
S_y is the yield strength of the crossbar in psi
E is the modulus of elasticity of the crossbar in psi
ppi is the price per inch of the crossbar in USD/inch
"""
struct CrossbarSquare <: Crossbar
    name::String
    t::Float32
    w::Float32
    S_y::Float32
    E::Float32
    ppi::Float32
end

"""
Crossbar Round
This is the struct for the solid round crossbars. It assumes that they are all solid round bars.

d is the diameter of the crossbar in inches
S_y is the yield strength of the crossbar in psi
E is the modulus of elasticity of the crossbar in psi
ppi is the price per inch of the crossbar in USD/inch
"""
struct CrossbarRound <: Crossbar
    name::String
    d::Float32
    S_y::Float32
    E::Float32
    ppi::Float32
end

"""
Pin
This is the struct object for the pins.

d is the diameter of the pin in inches
S_y is the yield strength of the pin in psi
p is the price for a singular pin in USD
"""
struct Pin
    name::String
    d::Float32
    len::Float32
    S_y::Float32
    price::Float32
end

# definition of Best Solution So Far (BSSF) struct
"""
Solution
This is the struct object for the Best Solution So Far (BSSF)

p is the price
diag is the struct object for the diagonal member
cross is the struct object for the crossbar
pin is the struct object for the pin
len_d is the length of the diagonal members
len_cb is the length of the crossbar
hole_d is the length from the center of the hole to the edge of the member for the diagonal members
hole_cb is the length from the center of the hole to the edge of the member for the crossbar
"""
struct Solution #updated such that the default Crossbar is round 
    price::Float16
    diag::DiagonalMember
    cross::CrossbarRound
    pin::Pin
    len_d::Float32
    len_cb::Float32
    hole_d::Float32
    hole_cb::Float32
end

#definition of struct for the safety profile
"""
SafetyProfile
this is the struct object that holds the factors of safety for the failure modes

diag_tearout is the factor of safety for the diagonal member in tearout
diag_axial is the factor of safety for the axial forces in the diagonal member
diag_bearing is the factor of safety for bearing in the diagonal member
cross_bearing is the factor of safety for bearing in the crossbar
cross_tearout is the factor of safety for tearout in the crossbar
cross_buckling is the factor of safety for buckling in the crossbar
pin_shear is the factor of safety for direct shear in the pin
pin_bearing is the factor of safety for bearing in the pin
"""
struct SafetyProfile
    diag_tearout::Float32
    diag_axial::Float32
    diag_bearing::Float32
    cross_bearing::Float32
    cross_buckling::Float32
    pin_shear::Float32
    pin_bearing::Float32
end


#############################################################################
# Stress functions:
#############################################################################
"""
scissor_forces(applied_force, cross_bar_length; diagonal_length = 8.5)
This function returns the force in the diagonal and cross-bar members of a scissor jack.

applied_force is the tension force applied to the jack in Newtons
cross_bar_length is the length of the cross bar in inches
diagonal_length is the length of the diagonal bar in inches. Its default length is 8.5 inches.
"""
function scissor_forces(applied_force, cross_bar_length; diagonal_length = 8.5)
    alpha = acos(cross_bar_length / (2 * diagonal_length)) #rad
    F_d = (1/2) * (applied_force / sin(alpha)) #lbf
    F_cb = applied_force * cot(alpha) #lbf
    return F_d, F_cb
end

"""
diagonal_bar_bearing_stress(force, diameter, thickness)
This function returns the bearing stress in the diagonal bar

force is the force in the diagonal member
diameter is the diameter of the hole in the diagonal member
thickness is the thickness of the diagonal member
"""
function diagonal_bar_bearing_stress(force, diameter, thickness)
    return (force) / (2 * diameter * thickness)
end

"""
cross_bar_bearing_stress(force, diameter, thickness), FOR SQUARE CROSSBAR
This function returns the bearing stress in the crossbar for a square shape

force is the force in the crossbar
diameter is the diameter of the hole in the crossbar
thickness is the thickness of the crossbar
"""
function cross_bar_bearing_stress(force, diameter, thickness)
    return (force) / (2 * diameter * thickness)
end

"""
cross_bar_bearing_stress_round(force, diameter, diameter_crossbar), FOR ROUND CROSSBAR
This function returns the bearing stress in the crossbar for a solid round

force is the force in the crossbar
diameter is the diameter of the hole in the crossbar
diameter_crossbar is the diameter of the crossbar
"""
function cross_bar_bearing_stress_round(force, diameter, diameter_crossbar)
    return (force) / (diameter * diameter_crossbar)
end

"""
cross_bar_tearout_stress_round(force, diameter_cb, diameter_pin, length)
This function returns the tearout stress in the diagonal bar

force is the force in the diagonal member
diameter_cb is the diameter of the crossbar (assumess solid round cross-section)
diameter_pin is the diameter of the pin
length is the distance from the center of the hole to the edge of the crossbar
"""
#function cross_bar_tearout_stress_round(force, diameter_cb, diameter_pin, length)
#    r_cb = diameter_cb / 2
#    r_pin = diameter_pin / 2
#    theta = atan(r_pin/r_cb)
#    width = 2 * r_cb * (cos(theta))
#area = width * length
#    return force / (2 * area)
#end

"""
pin_bearing_stress(force, diameter, thickness), FOR SQUARE CROSSBAR (also works for h-bar)
This function returns the bearing stress in the pin for a square shape crossbar

force is the force in the crossbar (or the diagonal member)
diameter is the diameter of the pin
thickness is the thickness of the crossbar (or the diagonal member)
"""
function pin_bearing_stress(force, diameter, thickness)
    return (force) / (2 * diameter * thickness)
end

"""
pin_bearing_stress(force, diameter, diameter_crossbar), FOR ROUND CROSSBAR
This function returns the bearing stress in the pin for a solid round crossbar

force is the force in the crossbar
diameter is the diameter of the pin
diameter_crossbar is the diameter of the crossbar
"""
function pin_bearing_stress_round(force, diameter, diameter_crossbar)
    r_cb = diameter_crossbar / 2
    r_pin = diameter / 2
    theta = atan(r_pin/r_cb)
    width = 2 * r_cb * (cos(theta))
    area = (pi * (r_cb ^2) * 4 * ((theta/(2*pi)))) + (2 * r_pin * (width/2)) #includes circular area (formed by sweeping theta for area and then adding triangles that were skipped by sweeps)
    return (force) / (diameter * diameter_crossbar)
end

"""
diagonal_bar_tearout_stress(force, length, thickness)
This function returns the tearout stress in the diagonal bar

force is the force in the diagonal member
length is the distance from the center of the hole to the edge of the diagonal member
thickness is the thickness of the diagonal member
"""
function diagonal_bar_tearout_stress(force, length, thickness)
    return (force) / (4 * length * thickness)
end

"""
diagonal_bar_tensile_stress(force, width, diameter, thickness)
This function returns the tensile stress in the diagonal bar

force is the force in the diagonal member
width is the width of the diagonal member
diameter is the diameter of the hole in the diagonal member
thickness is the thickness of the diagonal member
"""
function diagonal_bar_tensile_stress(force, width, diameter, thickness)
    return (force) / ((width - diameter) * thickness * 2)
end

"""
pin_direct_shear(force, n, diameter)
This function returns the direct shear force in the pin

force is the force in the cross bar
n is the number of shear planes in the pin
diameter is the diameter of the pin
"""
function pin_direct_shear(force, n, diameter)
    return (force) / (n * (pi / 4) * (diameter ^ 2))
end

"""
von_mises_stress(sigma_x, sigma_y, tau_xy)
This function returns the von-mises stress from plane stress

sigma_x is the stress in the x direction
sigma_y is the stress in the y direction
tau_xy is the shear in the xy direction
"""
function von_mises_stress(sigma_x, sigma_y, tau_xy)
    return sqrt((sigma_x^2) + (sigma_y^2) - (sigma_x*sigma_y) + (3*(tau_xy^2)))
end

#############################################################################
# Buckling functions:
#############################################################################

function k(I, A)
    return sqrt(I / A)
end

function l_over_k_1(c, E, s_y)
    return sqrt((2*(pi^2)*c*E) / (s_y))
end

function s_r(l, k)
    return l / k
end

"""
max_load_buckling(c, E, s_y, I, A, l)
This function returns the critical loading force for an Euler or Johnson column

c is the end-condition coefficient
E is the modulus of Elasticity
s_y is the yield strength
I is the moment of inertia
A is the area
l is the length
"""
function max_load_buckling(c, E, s_y, I, A, l)
    K = k(I, A)
    lk1 = l_over_k_1(c, E, s_y)
    S_r = s_r(l, K)
    if S_r > lk1
        return (A * c * (pi^2) * E) / (S_r ^ 2)
    else
        return A * (s_y - (1 / (c * E)) * ((s_y * S_r) / (2 * pi))^2)
    end
end

#############################################################################
# Geometry-related functions:
#############################################################################
"""
function that calculates the jack lifting range in inches given a start/end angle (in degrees, measured from the ground to the diagonal member)
the diagonal member length in inches, and the pin hole distance from the end of the member
"""
function jack_lift_range(angleStart, angleEnd, memberLength, holeDist) 
    heightInit = 2 * sind(angleStart) * (memberLength - holeDist) #calculates vertical component of starting jack position using angles and length
    heightFinal = 2 * sind(angleEnd) * (memberLength - holeDist)  #same as above with maximum jack solution
    liftRange = heightFinal - heightInit
    return liftRange
end

"""
calc_price(pin::Pin, diagonal_member::DiagonalMember, crossbar::Crossbar, diag_len, cb_len)
This function returns the price for the materials to make a scissor jack

pin is a Pin struct object
diagonal_member is a DiagonalMember struct object
crossbar is a Crossbar object
diag_len is the length of the diagonal member in inches
cb_len is the length of the crossbar in inches
"""
function calc_price(pin::Pin, diagonal_member::DiagonalMember, crossbar::Crossbar, diag_len, cb_len)
    pin_price = 4 * pin.price
    diag_price = 4 * (diag_len * diagonal_member.ppi)
    cb_price = cb_len * crossbar.ppi
    return pin_price + diag_price + cb_price
end

#############################################################################
# Safety Factor functions:
#############################################################################
"""
safety_diagonal_tear(applied_force, S_y_d, l_cb, l_d , l, t_d)
This function calulates the factor of safety for tearout failure in the diagonal members (Currently written for U-bar)

applied_force is the tension force applied to the scissor jack
S_y_d is the yield strength of the diagonal member
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
l is the tearout length of the diagonal member (center of hole to edge)
t_d is the thickness of the diagonal member
"""
function safety_diagonal_tear(applied_force, S_y_d, l_cb, l_d , l, t_d)
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    sigma_nom = diagonal_bar_tearout_stress(F_d, l, t_d)
    sigma_von = von_mises_stress(0.0, 0.0, sigma_nom)
    n = S_y_d / sigma_von
    return n
end

"""
safety_diagonal_axial(applied_force, S_y_d, l_cb, l_d, w_d, d, t_d)
This function calulates the factor of safety for axial failure in the diagonal members (Currently written for U-bar)

applied_force is the tension force applied to the scissor jack
S_y_d is the yield strength of the diagonal member
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
w_d is the width of the diagonal member
d is the diameter of the hole
t_d is the thickness of the diagonal member
"""
function safety_diagonal_axial(applied_force, S_y_d, l_cb, l_d, w_d, d, t_d)
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    sigma_nom = diagonal_bar_tensile_stress(F_d, w_d, d, t_d)
    sigma_von = von_mises_stress(sigma_nom, 0.0, 0.0)
    n = S_y_d / sigma_von
    return n
end

"""
safety_diagonal_bearing(applied_force, S_y_d, l_cb, l_d, d, t_d)
This function calulates the factor of safety for bearing failure in the diagonal members (Currently written for U-bar)

applied_force is the tension force applied to the scissor jack
S_y_d is the yield strength of the diagonal member
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
d is the diameter of the hole
t_d is the thickness of the diagonal member
"""
function safety_diagonal_bearing(applied_force, S_y_d, l_cb, l_d, d, t_d)
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    sigma_nom = diagonal_bar_bearing_stress(F_d, d, t_d)
    sigma_von = von_mises_stress(sigma_nom, 0.0, 0.0)
    n = S_y_d / sigma_von
    return n
end

"""
function safety_crossbar_bearing_round(applied_force, S_y_d, l_cb, l_d, d, d_cb)
This function calulates the factor of safety for bearing failure in the crossbar member (Currently written for solid round crossbar)

applied_force is the tension force applied to the scissor jack
S_y_d is the yield strength of the diagonal member
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
d is the diameter of the hole
d_cb is the diameter of crossbar (assumes crossbar is solid round)
"""
function safety_crossbar_bearing_round(applied_force, S_y_cb, l_cb, l_d, d, d_cb)
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    sigma_nom = cross_bar_bearing_stress_round(F_cb, d, d_cb)
    sigma_von = von_mises_stress(sigma_nom, 0.0, 0.0)
    n = S_y_cb / sigma_von
    return n
end

"""
safety_buckling_crossbar(applied_force, S_y_cb, E_cb, l_cb, w_cb, t_cb)
This function calulates the factor of safety for buckling failure in the crossbar (Written for hollow square bar)

applied_force is the tension force applied to the scissor jack
S_y_cb is the yield strength of the crossbar
E_cb is the modulus of elasticity for the crossbar
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
w_cb is the width of the crossbar 
t_cb is the thickness of the crossbar (assumes crossbar is hollow square)
"""
function safety_buckling_crossbar(applied_force, S_y_cb, E_cb, l_cb, l_d, w_cb, t_cb) #for square cross section crossbar
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    A = (w_cb ^ 2) - (w_cb - (2 * t_cb)) ^ 2
    I = ((w_cb ^ 4 / 12)) - ((w_cb - (2 * t_cb)) ^ 4 / 12)
    max_buck = max_load_buckling(C, E_cb, S_y_cb, I, A, l_cb)
    n = max_buck / F_cb
    return n
end

"""
safety_buckling_crossbar_round(applied_force, S_y_cb, E_cb, l_cb, l_d, d_cb)
This function calulates the factor of safety for buckling failure in the crossbar (Written for solid round bar)

applied_force is the tension force applied to the scissor jack
S_y_cb is the yield strength of the crossbar
E_cb is the modulus of elasticity for the crossbar
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
d_cb is the diameter of the crossbar (assumes crossbar is solid round)
"""
function safety_buckling_crossbar_round(applied_force, S_y_cb, E_cb, l_cb, l_d, d_cb, C) #for solid round cross section crossbar
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    A = ((d_cb / 2) ^ 2) * pi
    I = (pi / 64) * ((d_cb) ^ 4)
    max_buck = max_load_buckling(C, E_cb, S_y_cb, I, A, l_cb)
    n = max_buck / F_cb
    return n
end

"""
safety_crossbar_tear_round(applied_force, S_y_cb, l_cb, l_d, d_cb, d_pin, l)
This function calulates the factor of safety for tearout failure in the crossbar members (written for solid round)

applied_force is the tension force applied to the scissor jack
S_y_cb is the yield strength of the crossbar
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
d_cb is the diameter of the crossbar
d_pin is the diameter of the pin
l is the tearout length of the crossbar (center of hole to edge)
"""
#function safety_crossbar_tear_round(applied_force, S_y_cb, l_cb, l_d, d_cb, d_pin, l) #for solid round cross section member
#    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
#    sigma_nom = cross_bar_tearout_stress_round(F_cb, d_cb, d_pin, l)
#    sigma_von = von_mises_stress(0.0, 0.0, sigma_nom)
#    n = S_y_cb / sigma_von
#    return n
#end

"""
safety_pin_shear(applied_force, S_y_p, l_cb, l_d, N, d)
This function calulates the factor of safety for direct shear failure in the pins

applied_force is the tension force applied to the scissor jack
S_y_p is the yield strength of the pin
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
N is the number of shear planes on the pin (should be 2)
d is the diameter of the hole
"""
function safety_pin_shear(applied_force, S_y_p, l_cb, l_d, N, d)
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    tau_nom = pin_direct_shear(F_cb, N, d)
    sigma_von = von_mises_stress(0.0, 0.0, tau_nom)
    n = S_y_p / sigma_von
    return n
end

"""
safety_pin_bearing(applied_force, S_y_p, l_cb, l_d, d, t_cb)
This function calulates the factor of safety for bearing failure in the pin (Currently written for square crossbar)

applied_force is the tension force applied to the scissor jack
S_y_p is the yield strength of the pin
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
d is the diameter of the hole
t_cb is the thickness of the crossbar (assumes crossbar is hollow square)

Note that while there are three members in contact with the pin, two of the three are the same. 
    And this function only checks the bearing stress from the square crossbar
"""
function safety_pin_bearing(applied_force, S_y_p, l_cb, l_d, d, t_cb)
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    sigma_nom = pin_bearing_stress(F_cb, d, t_cb)
    sigma_von = von_mises_stress(sigma_nom, 0.0, 0.0)
    n = S_y_p / sigma_von
    return n
end

"""
safety_pin_bearing_round(applied_force, S_y_p, l_cb, l_d, d, t_cb)
This function calculates the factor of safety for bearing failure in the pin (Written for round crossbar)

applied_force is the tension force applied to the scissor jack
S_y_p is the yield strength of the pin
l_cb is the length of the cross bar (center of hole to center of hole)
l_d is the length of the diagonal member (center of hole to center of hole)
d_p is the diameter of the hole / diameter of the pin
d_cb is the diameter of the crossbar (assumes crossbar is solid circular)
t_d is the thickness of the diagonal member

Note that this function looks at the bearing stress at the contact between both the diagonal members and the crossbar
"""
function safety_pin_bearing_round(applied_force, S_y_p, l_cb, l_d, d_p, d_cb, t_d)
    F_d, F_cb = scissor_forces(applied_force, l_cb, diagonal_length = l_d)
    sigma_nom_diag = pin_bearing_stress(F_d, d_p, t_d)
    sigma_nom_cb = pin_bearing_stress_round(F_cb, d_p, d_cb)
    sigma_nom = max(sigma_nom_diag, sigma_nom_cb)
    sigma_von = von_mises_stress(sigma_nom, 0.0, 0.0)
    n = S_y_p / sigma_nom
    return n
end


"""
make_objective(diagonal_member, crossbar, pin, applied_force, end_condition_constant)
This function makes an objective function with the values from the components and loading conditions hard-coded into it.

diagonal_member is a DiagonalMember struct
crossbar is a Crossbar struct
pin is a Pin struct
applied_force is the tension force applied to the jack in lbf
end_condition_constant is the end condition constant for buckling in the crossbar
"""

function make_objective(diagonal_member::DiagonalMember, crossbar::Crossbar, pin::Pin, applied_force, end_condition_constant)
    function objective_function!(x, g)
        #x is a vector that will have:
            #the length of the diagonal member 
            #the length of the crossbar
            #the distance from the center of the hole to the edge of the member in the diagonal member
            #the distance from the center of the hole to the edge of the member in the crossbar
        #here we will "unpack" the vector into variables for legibility
        diagonal_length = x[1] #inches
        crossbar_length = x[2] #inches
        diagonal_tearout_length = x[3] #inches
        crossbar_tearout_length =x[4] #inches
        
        #Calculate safety factors for failure modes (write to items in vector g)
        #factor of safety for tearout on diagonal member
        g[1] = safety_diagonal_tear(applied_force, diagonal_member.S_y, crossbar_length, diagonal_length,
            diagonal_tearout_length, diagonal_member.t)

        #factor of safety for axial stress on diagonal member
        g[2] = safety_diagonal_axial(applied_force, diagonal_member.S_y, crossbar_length, diagonal_length,
            diagonal_member.w, pin.d, diagonal_member.t)

        #factor of safety for bearing stress on diagonal member
        g[3] = safety_diagonal_bearing(applied_force, diagonal_member.S_y, crossbar_length, diagonal_length,
            pin.d, diagonal_member.t) 

        #factor of safety for bearing stress on crossbar
        g[4] = safety_crossbar_bearing_round(applied_force, diagonal_member.S_y, crossbar_length,
            diagonal_length, pin.d, crossbar.d) 

        #for solid round cross section crossbar
        g[5] = safety_buckling_crossbar_round(applied_force, crossbar.S_y, crossbar.E, crossbar_length, 
            diagonal_length, crossbar.d, end_condition_constant)

        #factor of safety for shear on pin
        g[6] = safety_pin_shear(applied_force, pin.S_y, crossbar_length, diagonal_length, 2, pin.d) 

        #factor of safety for bearing stress on pin
        g[7] = safety_pin_bearing_round(applied_force, pin.S_y, crossbar_length, diagonal_length, pin.d,
            crossbar.d, diagonal_member.t) 
        
        #Check geometry (write to item in vector g)
        #the magnitude of the range of movement the jack has
        g[8] = jack_lift_range(15, 75, diagonal_length + (2 * diagonal_tearout_length),
            2 * diagonal_tearout_length) 
        
        #Check any other constraints (write to item in vector g) (weight?)

        #Calculate price
        price = calc_price(pin, diagonal_member, crossbar, diagonal_length, crossbar_length)
        return price
    end
    return objective_function!
end


"""
combo_optimizer(diagonal_member::DiagonalMember, crossbar::Crossbar | CrossbarRound, pin::Pin, safety::SafetyProfile)
This is a function that takes a given combination of a diagonal member, crossbar, pin, and safety profile returns a Solution struct
"""
function combo_optimizer(diagonal_member::DiagonalMember, crossbar::Crossbar, pin::Pin, safety::SafetyProfile)
    #check that the pin is compatable
    if false #replace "false" with a function that checks that the pin is compatible
        return nothing
    end
        
    # make the objective function
    objective! = make_objective(diagonal_member, crossbar, pin, 1000, 1) #double check the applied force

    #set the optimizer variables
    x0 = [1.0, 1.0, 1.0, 1.0]  # starting point
    lx = [0.01, 0.01, 0.01, 0.01]  # lower bounds on x
    ux = [24.0, 48.0, 6.0, 6.0]  # upper bounds on x
    ng = 8  # number of constraints
    lg = [safety.diag_tearout, safety.diag_axial, safety.diag_bearing, safety.cross_bearing, 
            safety.cross_buckling, safety.pin_shear, safety.pin_bearing, 6.0]  # lower bounds on g
    ug = Inf*ones(ng)  # upper bounds on g
    options = Options(solver=IPOPT())  # choosing IPOPT solver

    #run the optimizer
    xopt, fopt, info = minimize(objective!, x0, ng, lx, ux, lg, ug, options)

    #build the solution
    solution = Solution(fopt, diagonal_member, crossbar, pin, xopt[1], xopt[2], xopt[3], xopt[4])

    #return the solution
    return solution
end


#add all of the components here

#pin components
#pinXX = Pin(name, diameter, length_of_shaft, yield_strength, price)
pin1 = Pin("92314A240",0.19,0.375,80000,8.39)
pin2 = Pin("93190A794",0.625,1,80000,1.71)
pin3 = Pin("92186A546",0.25,1.5,80000,3.98)
pin4 = Pin("91735A409",0.3125,0.5,80000,11.08)
pin5 = Pin("91735A611",0.21875,1.5,80000,11.40)
pin6 = Pin("91259A626", 0.375, 1.25, 140000, 2.54)
pin7 = Pin("91259A468", 0.188, 1.25, 140000, 7.91)
pin8 = Pin("91259A475", 0.313, 0.438, 140000, 8.46)
pin9 = Pin("91259A796", 0.625, 1.25, 140000, 5.86)
pin10 = Pin("91259A583", 0.313, 1.00, 140000, 2.09)
pin11 = Pin("91259A626", 0.375, 1.25, 140000, 2.54)
pin12 = Pin("91259A468", 0.188, 1.25, 140000, 7.91)
pin13 = Pin("91259A475", 0.313, 0.438, 140000, 8.46)
pin14 = Pin("91259A796", 0.625, 1.25, 140000, 5.86)
pin15 = Pin("91259A583", 0.313, 1.00, 140000, 2.09)

pins = (pin1, pin2, pin3, pin4, pin5,pin6,pin7,pin8,pin9,pin10, pin11, pin12, pin13, pin14, pin15)  


#Diagonal member components
#diagXX = DiagonalMember(name, thickness, width, yield_strength, price_per_inch)

diag1 = DiagonalMember("7779T31", 0.125, 0.375, 36000, 0.7758333)
diag2 = DiagonalMember("7779T33", 0.150, 1.00, 36000, 0.464166)
diag3 = DiagonalMember("7779T35", 0.125, 0.75, 36000, 0.4475)
diag4 = DiagonalMember("7779T37", 0.156, 0.5, 36000, 0.2725)
diag5 = DiagonalMember("7779T43", 0.17, 1.41, 36000, 0.9125)
diag6 = DiagonalMember("7779T45", 0.184, 1.584, 36000, 1.1691)
diag7 = DiagonalMember("7779T47", 0.19, 1.75, 36000, 1.46166)
diag8 = DiagonalMember("7779T49", 0.2, 1.92, 36000, 1.7475)
diag9 = DiagonalMember("7779T39", 0.1875, 1, 36000, 0.6066)
diag10 = DiagonalMember("7779T41", 0.25, 0.625, 36000, 0.6775)
diag11 = DiagonalMember("7779T31", 0.125, 0.375, 36000, 0.7758333)
diag12 = DiagonalMember("7779T33", 0.150, 1.00, 36000, 0.464166)
diag13 = DiagonalMember("7779T35", 0.125, 0.750, 36000, 0.4475)
diag14 = DiagonalMember("7779T37", 0.156, 0.500, 36000, 0.2725)
diag15 = DiagonalMember("7779T43", 0.170, 1.41, 36000, 0.9125)

diags = (diag1,diag2,diag3,diag4,diag5,diag6,diag7,diag8,diag9,diag10, diag11, diag12, diag13, diag14, diag15)

#Crossbar components
#cbXX = CrossbarRound(name, diameter, yield_strength, modulus_of_elasticity, price_per_inch)

cb1 = CrossbarRound("9210K16",1,25000,225000,4.23916)
cb2 = CrossbarRound("9210K12",0.5,25000,225000,1.67)
cb3 = CrossbarRound("6818T54",0.5625,130000,29000000,1.6708)
cb4 = CrossbarRound("6818T53",0.5,130000,29000000,1.365)
cb5 = CrossbarRound("6818T27",1.5,130000,29000000,4.31)
cb6 = CrossbarRound("8935K336", 0.375, 95000, 29000000, 0.298333)
cb7 = CrossbarRound("8935K24", 3.000, 95000, 29000000, 11.248333)
cb8 = CrossbarRound("8935K11", 1.500, 95000, 29000000, 3.44554)
cb9 = CrossbarRound("8920K231", 1.000, 54000, 29000000, 1.015833)
cb10 = CrossbarRound("8920K195", 0.7500, 54000, 29000000, 0.6341666)
cb11 = CrossbarRound("8935K336", 0.375, 95000, 29000000, 0.298333)
cb12 = CrossbarRound("8935K24", 3.000, 95000, 29000000, 11.248333)
cb13 = CrossbarRound("8935K11", 1.500, 95000, 29000000, 3.44554)
cb14 = CrossbarRound("8920K231", 1.000, 54000, 29000000, 1.015833)
cb15 = CrossbarRound("8920K195", 0.7500, 54000, 29000000, 0.6341666)
    
cbs = (cb1, cb2, cb3, cb4, cb5, cb6, cb7, cb8, cb9, cb10, cb11, cb12, cb13, cb14, cb15)

#define the safety profile 
#SafetyProfile(diag_tearout, diag_axial, diag_bearing, cross_bearing, cross_buckling, pin_shear, pin_bearing)
safety_profile = SafetyProfile(3.0, 2.0, 2.0, 2.0, 4.0, 2.0, 2.0)


#initialize the Best Solultion So Far (bssf)
bssf = Solution(1000.0, diag1, cb1, pin1, 1.0, 1.0, 1.0, 1.0)

#go through the design space
for pin in pins 
    for diag in diags 
        for cb in cbs
            #check that components are compatible

            
            #run the optimization 
            solution = combo_optimizer(diag, cb, pin, safety_profile)

            #compare to bssf
            if solution.price < bssf.price
                global bssf = solution
            end

        end
    end
end
println(bssf)