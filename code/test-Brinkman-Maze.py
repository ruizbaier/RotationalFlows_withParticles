from fenics import *


'''
Brinkman flow on a maze-shaped domain

imposing INTLET u and wall u

'''


parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True

fileO = XDMFFile("outputs/Brinkman-Maze2D.xdmf")
fileO.parameters['rewrite_function_mesh']=True
fileO.parameters["functions_share_mesh"] = True
fileO.parameters["flush_output"] = True

mesh = Mesh("meshes/2DMaze.xml")
bdry = MeshFunction("size_t", mesh, "meshes/2DMaze_facet_region.xml")
ds = Measure("ds", subdomain_data=bdry)

wall = 73; outlet=72; inlet=71; 

n = FacetNormal(mesh);

t = 0.0; dt = 0.1; T = 1.; frequencySave = 2;


# ********* model parameters forcing terms ********* #
nu = Constant(0.0001)
K  = Constant(1.e-4)

# ********* Finite dimensional spaces MINI ********* #

Hv = FiniteElement("BDM", mesh.ufl_cell(), 1)
Ht = FiniteElement("CG", mesh.ufl_cell(), 1)
Hq = FiniteElement("DG", mesh.ufl_cell(), 0)
Hh = FunctionSpace(mesh,MixedElement([Hv,Ht,Hq]))

v, th, q = TestFunctions(Hh)
sol = Function(Hh); 
u,  w, p = TrialFunctions(Hh) 

print("....... DoF = ", Hh.dim())
    
# ********* Boundaries, boundary conditions, and initial conditions ********* #

uflow = Expression(('t*200*(x[1]-0.45)*(0.55-x[1])','0'), degree = 2, t = t)
#wflow = Expression('-20*pow(nu,0.5)*(1-2*x[1])', nu = nu, t = t, degree =1)

bcUin    = DirichletBC(Hh.sub(0), uflow, bdry, inlet)
bcUwall  = DirichletBC(Hh.sub(0), Constant((0,0)), bdry, wall)

bcs = [bcUin, bcUwall]


uold = project(Constant((0,0)),Hh.sub(0).collapse())

# ******** Weak forms ************* #

AA = 1/K * dot(u,v)*dx \
     + pow(nu,0.5)*dot(v,curl(w)) * dx \
     + pow(nu,0.5)*dot(u,curl(th)) * dx \
     - p*div(v)*dx - q*div(u)*dx - w*th*dx

BB = 1/K * dot(uold,v)*dx

# ******* Solving *************** #
inc = 0
while (t <= T+dt):
    print("t=%.3f" % t)
    uflow.t = t
    solve(AA == BB, sol, bcs)
    u_h,w_h,p_h = sol.split()
    assign(uold,u_h)

    if (inc % frequencySave == 0):
        u_h.rename("u","u"); fileO.write(u_h,t)
        w_h.rename("w","w"); fileO.write(w_h,t)
        p_h.rename("p","p"); fileO.write(p_h,t)
    t += dt; inc += 1
