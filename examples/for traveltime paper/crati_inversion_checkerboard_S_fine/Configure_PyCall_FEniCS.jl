## configure PyCall
using PyCall,Conda
ENV["PYTHON"]="/usr/bin/python";
Pkg.build("PyCall");
Pkg.build("FEniCS");
Pkg.build("Conda");
##
using PyCall
using FEniCS
mesh=RectangleMesh([0,0],[200,100],200,100);
V = FunctionSpace(mesh, "Lagrange",1);
##
fi=pyimport("fenics");
py"""
import numpy as np

def sinpi(x):
    return np.sin(np.pi * x)

"""
py"sinpi"(1)
py"""
import dolfin as dl
"""
