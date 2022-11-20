"""
Author: Sergio Lopez  Banos
GitHub: LopezBanos
Exercise Title: Chaotic Pendulum
Date: October 2022
Description: In this example I study the chaotic pendulum restricted to '
2 dimensions.
"""

function Chaotic_Pendulum(g::Float64, l1::Float64, l2::Float64, m1::Float64, m2::Float64, dt::Float64, time::Int64, θ1_ini::Float64, θ2_ini::Float64, ω1_ini::Float64, ω2_ini::Float64) 

    # Initialize angles and angular velocities
    θ1 = θ1_ini
    θ2 = θ2_ini
    ω1 = ω1_ini
    ω2 = ω2_ini
    x1 = l1*sin(θ1)     
    y1 = -l1*cos(θ1)
    x2 = x1 + l2*sin(θ2)
    y2 = y1 - l2*cos(θ2)

    list_θ1 = [θ1]
    list_θ2 = [θ2]
    list_ω1 = [ω1]
    list_ω2 = [ω2]
    list_x1 = [x1]
    list_x2 = [x2]
    list_y1 = [y1]
    list_y2 = [y2]
    t_list = [0.0]

    # Euler Method
    t = 0.0
    while (t<time)
        t = t + dt
        dω1_dt = (-g*(2*m1 + m2)*sin(θ1) - m2*g*sin(θ1- 2*θ2)-2*sin(θ1 - θ2)*m2*((ω2^2)*l2 + (ω1^2)*l1*cos(θ1 - θ2)))/(l1*(2*m1 + m2 - m2*cos(2*θ1 - 2*θ2)))
        dω2_dt = (2*sin(θ1 - θ2)*((ω1^2)*l1*(m1 + m2) + g*(m1 + m2)*cos(θ1) + (ω2^2)*l2*m2*cos(θ1 - θ2)))/(l1*(2*m1 + m2 - m2*cos(2*θ1 - 2*θ2)))
        
        # Angular Velocities
        ω1 = ω1 + dω1_dt*dt
        ω2 = ω2 + dω2_dt*dt
        # Angles
        θ1 = θ1 + ω1*dt
        θ2 = θ2 + ω2*dt

        # Position
        x1 = l1*sin(θ1)     
        y1 = -l1*cos(θ1)
        x2 = x1 + l2*sin(θ2)
        y2 = y1 - l2*cos(θ2)

        append!(t_list, t)
        append!(list_θ1, θ1)
        append!(list_θ2, θ2)
        append!(list_ω1, ω1)
        append!(list_ω2, ω2)
        append!(list_x1, x1)
        append!(list_x2, x2)
        append!(list_y1, y1)
        append!(list_y2, y2)
    end
    return list_θ1,list_θ2,list_ω1,list_ω2,list_x1,list_x2,list_y1,list_y2,t_list
end 



# PARAMETERS
P1_LENGTH = 1.0
P2_LENGTH = 1.5
M1 = 1.3
M2 = 2.2
GRAVITY = 9.81
TIME = 50 #sec
STEP = 0.01
θ1_INI = -3.1415/2
θ2_INI = 1.2
ω1_INI = 0.0
ω2_INI = 0.0

data = Chaotic_Pendulum(GRAVITY, P1_LENGTH, P2_LENGTH, M1, M2, STEP, TIME, θ1_INI, θ2_INI, ω1_INI, ω2_INI)

# Plotting Section
using Plots
# PHASE SPACE ANGLES
plot(data[1], data[2], label="θ2 vs θ1")
plot!(xlab="θ1", ylab="θ2")
title!("Angle Space Phase")

# PHASE SPACE ANGULAR VELOCITIES
plot(data[3], data[4], label="ω2 vs ω1")
plot!(xlab="ω1", ylab="ω2")
title!("Angular Velocity Space Phase")

# Angles with respect time
plot(data[9], data[1], label="θ1(t)")
plot!(data[9], data[2], label="θ2(t)")
plot!(xlab="Time - Seconds", ylab="θ(t)")
title!("θ as Function of Time")

# Angular Velocities with respect time
plot(data[9], data[3], label="ω1(t)")
plot!(data[9], data[4], label="ω2(t)")
plot!(xlab="Time - Seconds", ylab="ω(t)")
title!("ω as Function of Time")

# Phase Space Spatial Coordinates
plot(data[5], data[7], label="XY - Pendulum 1")
plot!(data[6], data[8], label="XY - Pendulum 2")
plot!(xlab="X Coordinate", ylab="Y - Coordinate")
title!("Phase Space of Spatial Coordinates")

# PHASE SPACE
plot(data[1], data[3], label="θ1")
plot!(xlab="θ1", ylab="ω1")
title!("Phase Space Chaotic Pendulum")

plot(data[2], data[4], label="θ2")
plot!(xlab="θ2", ylab="ω2")
title!("Phase Space Chaotic Pendulum")