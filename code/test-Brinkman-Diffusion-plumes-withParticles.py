'''

* Brinkman - diffusion coupling + particle tracking 

* permeability depends on concentration

* strong coupled form 

       d/dt c + u . grad(c) - div(D*grad(c)) = 0
d/dt u + 1/K(c) u + sqrt(nu)*curl w + grad p = c*(0,1)
                                       div u = 0
                         w - sqrt(nu)*curl u = 0

and particles advected by Brinkman velocity 

* boundary conditions: 
         u.n = 0 everywhere 
           c = 0 on top
           c = 1 on bottom
 D*grad(c).n = 0 on walls

'''

from dolfin import *
from LagrangianParticles import LagrangianParticles
from particle_generators import RandomCircle

parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 3

fileO = XDMFFile("outputs/Brinkman-Diffusion-plumes-and-particles.xdmf")
fileO.parameters['rewrite_function_mesh']=True
fileO.parameters["functions_share_mesh"] = True
fileO.parameters["flush_output"] = True

# ************ Time units ********************* #

time = 0.0; Tfinal = 12.0;  dt = 0.1; 
inc = 0;   frequency = 10;

# ***** Model coefficients and parameters ***** #
D    = Constant(1.e-4)
K0   = Constant(1.e-1)
nu   = Constant(1.e-2)

c0   = Expression("1-x[1]*x[1]+0.01*cos(4*pi*x[0])*sin(pi*x[1])", degree = 1)

f = lambda c : c*Constant((0,1))

mesh = RectangleMesh(Point(0,0),Point(2,1),80,40, 'crossed') 

bdry = MeshFunction("size_t", mesh, 1)
bdry.set_all(0)

bot = 31; top =32; wall= 33

GBot = CompiledSubDomain("near(x[1],0.) && on_boundary")
GTop = CompiledSubDomain("near(x[1],1.) && on_boundary")
GWall = CompiledSubDomain("(near(x[0],0.) || near(x[0],2.)) && on_boundary")
GTop.mark(bdry,top); GBot.mark(bdry,bot); GWall.mark(bdry,wall)
ds = Measure("ds", subdomain_data = bdry)

x = SpatialCoordinate(mesh)
K = lambda c: K0/exp(-9.7044*c + 4.1588*(1.-x[1]))

# ********* Finite dimensional spaces ********* #

Hs = FiniteElement('CG', mesh.ufl_cell(),1)
Hv = FiniteElement('RT', mesh.ufl_cell(),1)
Hq = FiniteElement('DG', mesh.ufl_cell(),0)
Ht = FiniteElement('CG', mesh.ufl_cell(),1)

Hh = FunctionSpace(mesh, MixedElement([Hs,Hv,Hq,Ht]))
Vh = FunctionSpace(mesh, Hv)

print (" ****** Total DoF = ", Hh.dim())

trial = TrialFunction(Hh)
sol = Function(Hh)
c, u, p, w  = split(sol)
s, v, q, th = TestFunctions(Hh)

uold = Function(Hh.sub(1).collapse())
cold = interpolate(c0,Hh.sub(0).collapse())

# BCs
bc1 = DirichletBC(Hh.sub(0), Constant(0), GTop)
bc2 = DirichletBC(Hh.sub(0), Constant(1), GBot)
bcu = DirichletBC(Hh.sub(1), Constant((0,0)), 'on_boundary')
bcs = [bc1,bc2,bcu]


# *********** Insert particles

particle_positions = RandomCircle([1., 0.5], 0.25).generate([50, 50])

lp = LagrangianParticles(Vh)
lp.add_particles(particle_positions)

# ******** Weak form Brinkman-diffusion **** #

FF = 1./dt*(c-cold)*s*dx \
    + dot(u,grad(c))*s*dx \
    + D*dot(grad(c),grad(s))*dx \
    + 1./dt*dot(u-uold,v)*dx \
    + 1./K(c)*dot(u,v)*dx \
    + sqrt(nu)*dot(curl(w),v)*dx \
    - p * div(v) * dx \
    - dot(f(c),v) * dx \
    - q * div(u) * dx \
    - w * th * dx \
    + sqrt(nu)*dot(curl(th),u)*dx

Tang = derivative(FF, sol, trial)
problem = NonlinearVariationalProblem(FF, sol, bcs, J=Tang)
solver  = NonlinearVariationalSolver(problem)
solver.parameters['nonlinear_solver']                    = 'newton'
solver.parameters['newton_solver']['linear_solver']      = 'mumps'
solver.parameters['newton_solver']['absolute_tolerance'] = 1e-7
solver.parameters['newton_solver']['relative_tolerance'] = 1e-7
solver.parameters['newton_solver']['maximum_iterations'] = 25

# ********* Time loop ************* #

while (time <= Tfinal):
    
    print("time = %.3f" % time)
    solver.solve()
    ch,uh,ph,wh = sol.split()

    ## ADVECTING particles
    lp.step(uh, dt=dt)
    
    
    
    assign(cold,ch)    
    assign(uold,uh)
     
    if (inc % frequency == 0):
        ch.rename("c","c"); fileO.write(ch,time)
        uh.rename("u","u"); fileO.write(uh,time)
        ph.rename("p","p"); fileO.write(ph,time)
        wh.rename("w","w"); fileO.write(wh,time)

        ### HERE I don't know how to save the positions into the same Xdmf
        points_list = list(Point(*pp) for pp in particle_positions)
        #lp.rename("part","part"); fileO.write(lp,time)
        #fileO.write(points_list,time)
        XDMFFile("outputs/particles"+str(inc)+".xdmf").write(points_list)
        
    inc += 1; time += dt
    
# ******** End of time loop **** #
