using Plots

#diameter of hole in diagonal member and cross bar
d = 0.25 #in
#thickness of diagonal member
t_d = 0.125 #in
#thickness of crossbar
t_cb = 0.06 #in
#length from center of hole to edge in diagonal member and cross bar
l = 0.4 #in
# width of diagonal member
w = 1.25 #in
#number of shear planes on pin
n = 2
# Applied Force
F_A = 700 #lbf
# length of diagonal member
l_diag = 8.5 #inches
#start and end angles, degrees
theta_start = 0
theta_end = 45
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

function plot_diagonal(l, d)
    plot(l, d)
    plot!(xlabel = "Cross-bar length (in)")
    plot!(ylabel = "Force diagonal / Force applied")
    plot!(legend = false)
    plot!(grid = false)
    savefig("C:\\Users\\xcmad\\Documents\\Homework\\MEEN372\\ScissorJack\\diagonalFig.png")
end

function plot_cross_bar(l, c)
    plot(l, c)
    plot!(xlabel = "Cross-bar length (in)")
    plot!(ylabel = "Force cross-bar / Force applied")
    plot!(legend = false)
    plot!(grid = false)
    savefig("C:\\Users\\xcmad\\Documents\\Homework\\MEEN372\\ScissorJack\\crossBarFig.png")
end

#force = 100.0 #N
#lengths = range(start = 0.1, stop = 16, length = 100) #inches

#d_ratios = zeros(Float32, length(lengths))
#c_ratios = zeros(Float32, length(lengths))

#for i in range(1,length(lengths))
#    d, cb = scissor_forces(force, lengths[i])
#    d_ratios[i] = d/force 
#    c_ratios[i] = cb/force 
#end

#plot_diagonal(lengths, d_ratios)
#plot_cross_bar(lengths, c_ratios)


function diagonal_bar_bearing_stress(force, diameter, thickness)
    return (force) / (2 * diameter * thickness)
end

function cross_bar_bearing_stress(force, diameter, thickness)
    return (force) / (2 * diameter * thickness)
end

function pin_bearing_stress(force, diameter, thickness)
    return (force) / (2 * diameter * thickness)
end

function diagonal_bar_tearout_stress(force, length, thickness)
    return (force) / (4 * length * thickness)
end

function diagonal_bar_tensile_stress(force, width, diameter, thickness)
    return (force) / ((width - diameter) * thickness * 2)
end

function pin_direct_shear(force, n, diameter)
    return (force) / (n * (pi / 4) * (diameter ^ 2))
end

function von_mises_stress(sigma_x, sigma_y, tau_xy)
    return sqrt((sigma_x^2) + (sigma_y^2) - (sigma_x*sigma_y) + (3*(tau_xy^2)))
end

"""
function that calculates the jack lifting range in inches given a start/end angle (in degrees, measured from the ground to the diagonal member)
the diagonal member length in inches, and the pin hole distance from the end of the member
"""
function jack_lift_range(angleStart, angleEnd, memberLength, holeDist) 
    heightInit = 2 * sind(angleStart) * (memberLength - holeDist) #calculates vertical component of starting jack position using angles and length
    heightFinal = 2 * sind(angleEnd) * (memberLength - holeDist)  #same as above with maximum jack sonition
    liftRange = heightFinal - heightInit
    return liftRange
end

#Forces in members ****************************************
F_d, F_cb = scissor_forces(F_A, 13.02)

println("The force in the crossbar is: " * string(F_cb) * " lbf.")
println("The force in the diagonal member is: " * string(F_d) * " lbf.")

#diagonal member*******************************************
#tearout
dt = diagonal_bar_tearout_stress(F_d, l, t_d)
dt = von_mises_stress(0.0, 0.0, dt)
dt = dt/1000
println("diagonal bar tearout stress: " * string(dt) * " kpsi")
#axial/tensile
da = diagonal_bar_tensile_stress(F_d, w, d, t_d)
da = von_mises_stress(da, 0.0, 0.0)
da = da / 1000
println("diagonal bar axial stress: " * string(da) * " kpsi")
#bearing
db = diagonal_bar_bearing_stress(F_d, d, t_d)
db = von_mises_stress(db, 0.0, 0.0)
db = db / 1000
println("diagonal bar bearing stress: " * string(db) * " kpsi")
#crossbar **************************************************
#bearing
cb = cross_bar_bearing_stress(F_cb, d, t_cb)
cb = von_mises_stress(cb, 0.0, 0.0)
cb = cb / 1000
println("crossbar bearing stress: " * string(cb) * " kpsi")
#pin ****************************************************
#shear
ps = pin_direct_shear(F_cb, n, d)
ps = von_mises_stress(0.0, 0.0, ps)
ps = ps / 1000
println("pin shear stress: " * string(ps) * " kpsi")
#bearing
pb = pin_bearing_stress(F_cb, d, t_cb)
pb = von_mises_stress(pb, 0.0, 0.0)
pb = pb / 1000
println("pin bearing stress: " * string(pb) * " kpsi")

#jack lift range
lr = jack_lift_range(theta_start, theta_end, l_diag, l)
println("jack lift range: " * string(lr) * " inches")