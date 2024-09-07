# Define variables
CXX = c++
CXXFLAGS = -O3 -Wall -shared -std=c++11 -fPIC
PYBIND11_INCLUDES = $(shell python3 -m pybind11 --includes)
EXT_SUFFIX = $(shell python3-config --extension-suffix)
TARGET = Orbit_Sim

# The rule to build the target
$(TARGET)$(EXT_SUFFIX): Orbit_Sim.cpp
	$(CXX) $(CXXFLAGS) $(PYBIND11_INCLUDES) Orbit_Sim.cpp -o $(TARGET)$(EXT_SUFFIX)

# Clean rule
clean:
	rm -f $(TARGET)$(EXT_SUFFIX)
