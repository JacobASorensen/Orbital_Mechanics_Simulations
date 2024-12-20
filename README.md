# Orbital Mechanics Simulations
## Intro

During one of the first few lectures in this class we hardcoded a 3-body problem in python. I decided to generalize the code to N-bodies since I wanted to look at simulations with more than just 3 bodies (Also, ever since my first comp-sci class I dislike manually coding what can be generalized).

Since my personal laptop is rather slow I thought it would be good to try to convert the simulation code to C++. I had ChatGPT convert my N-Body python ODE code, and the RKF4 solver our professor gave us into C++. With Pybind I was amazed at how easy it was to run compiled C++ as a python function.

The 'lectures' and 'hw' directories contain my notes and homework for this class, while the 'cplusplus' directory contains the n-body code along with helpers.

## Project Purpose

### I want to better understand n-body dynamics for n > 2
Specifically I am interested in learning more about Lagrange points, and how satellites will behave in the different types of lagrange points.
I also want to see what types of maneuevers are required to enter and stay within a lagrange point. (How are there asteroids in Earth-Sun L4 and L5?)
### I want to learn and subsequently improve my ability to interface C++ with Python
The simulation code is written in C++, and the main interface is written in python.
I had never before used PyBind. Learning how to use it was a cool learning experience.
The speed of C++ mixed with the ease of use and data processing prowess of python!

### I want to see what sort of information I can glean
Specifically I want to try to simulate a grid of different starting positions within lagrange points, and see where the most stable position is

## Features

N-Body simulator in C++ bound to python with pybind. 

![5-body animation](animations/5_body_anim2.gif)

NBodySim python class which:
- simplifies the interface between the N-body simulator and JPL Horizons API calls (using astroquery)
- provides graphing methods to easily display simulation results (using matplotlib)

![Jupiter and Earth](animations/Jupiter_Earth.gif)

![Inner Solar System](animations/Inner_Solar_System.gif)

- Translation feature for graphs to see what the orbits look like from a target body (Rotating reference frame option planned to be added to this function)

![Earth Mars and Jupiter translated to Earth](animations/Earth_Focus.gif)


## Features to add

### - mass data from JPL Horizons 
Pull mass data from JPL Horizons rather than using hardcoded values


### - Graph Focus function
Add target object functionality in order to graph rotating reference frames

#### -- A function to add satellites
I would like to be able to add small satellites given orbital ephemerides with reference to the target body

#### -- Preventing satellites from affecting large planetary bodies
If I wanted to add 50 satellites with slightly different orbital parameters, this would add unnecessary computation steps to the simulation due to their inconsequential mass. Therefore I want to add a feature which removes their force effects on other bodies

### - Error Correction From JPL Horizon data
It would be nice to correct the simulation for named bodies at key times based off of JPL Horizons data

#### -- Error analysis compared to JPL Horizons baseline
It would be nice to estimate the error over time of this simulator using JPL Horizons as a baseline

### - some sort of patched conics implementation
Given a focus it would be cool to try out some patched conics

### - 3d or other ways of displaying calculated simulation data


## Runge Kutte Solvers

I grabbed these from wikipedia in order to implement my RKF4(5) solver

![alt text](images/RKF4.png)

![alt text](images/RKF45.png)

![alt text](images/functions.png)

![alt text](images/RK45_Coefficients.png)
