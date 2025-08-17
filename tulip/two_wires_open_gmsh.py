import gmsh
import sys

# Initialize Gmsh
gmsh.initialize()

# Set the model name
gmsh.model.add("two_wires_open_gmsh")

# --- Define Geometry ---

# Parameters for the conductors and the open boundary
center_1_x = 0.025
center_1_y = 0
center_2_x = -0.025
center_2_y = 0
conductor_radius = 0.002
open_boundary_radius = 0.5

# Characteristic length (mesh size)
# A smaller value will result in a finer mesh
mesh_size_conductors = 0.0005
mesh_size_boundary = 0.05

# -- Conductor 1 --
# Create center point and points on the circle
p_center1 = gmsh.model.geo.addPoint(center_1_x, center_1_y, 0, mesh_size_conductors)
p1_1 = gmsh.model.geo.addPoint(center_1_x + conductor_radius, center_1_y, 0, mesh_size_conductors)
p1_2 = gmsh.model.geo.addPoint(center_1_x - conductor_radius, center_1_y, 0, mesh_size_conductors)
p1_3 = gmsh.model.geo.addPoint(center_1_x, center_1_y + conductor_radius, 0, mesh_size_conductors)
p1_4 = gmsh.model.geo.addPoint(center_1_x, center_1_y - conductor_radius, 0, mesh_size_conductors)

# Create circular arcs
c1_1 = gmsh.model.geo.addCircleArc(p1_1, p_center1, p1_3)
c1_2 = gmsh.model.geo.addCircleArc(p1_3, p_center1, p1_2)
c1_3 = gmsh.model.geo.addCircleArc(p1_2, p_center1, p1_4)
c1_4 = gmsh.model.geo.addCircleArc(p1_4, p_center1, p1_1)

# Create a curve loop from the arcs
cl1 = gmsh.model.geo.addCurveLoop([c1_1, c1_2, c1_3, c1_4])

# -- Conductor 2 --
# Create center point and points on the circle
p_center2 = gmsh.model.geo.addPoint(center_2_x, center_2_y, 0, mesh_size_conductors)
p2_1 = gmsh.model.geo.addPoint(center_2_x + conductor_radius, center_2_y, 0, mesh_size_conductors)
p2_2 = gmsh.model.geo.addPoint(center_2_x - conductor_radius, center_2_y, 0, mesh_size_conductors)
p2_3 = gmsh.model.geo.addPoint(center_2_x, center_2_y + conductor_radius, 0, mesh_size_conductors)
p2_4 = gmsh.model.geo.addPoint(center_2_x, center_2_y - conductor_radius, 0, mesh_size_conductors)

# Create circular arcs
c2_1 = gmsh.model.geo.addCircleArc(p2_1, p_center2, p2_3)
c2_2 = gmsh.model.geo.addCircleArc(p2_3, p_center2, p2_2)
c2_3 = gmsh.model.geo.addCircleArc(p2_2, p_center2, p2_4)
c2_4 = gmsh.model.geo.addCircleArc(p2_4, p_center2, p2_1)

# Create a curve loop from the arcs
cl2 = gmsh.model.geo.addCurveLoop([c2_1, c2_2, c2_3, c2_4])

# -- Open Boundary --
# Create center point and points on the circle
p_center_b = gmsh.model.geo.addPoint(0, 0, 0, mesh_size_boundary)
pb_1 = gmsh.model.geo.addPoint(open_boundary_radius, 0, 0, mesh_size_boundary)
pb_2 = gmsh.model.geo.addPoint(-open_boundary_radius, 0, 0, mesh_size_boundary)
pb_3 = gmsh.model.geo.addPoint(0, open_boundary_radius, 0, mesh_size_boundary)
pb_4 = gmsh.model.geo.addPoint(0, -open_boundary_radius, 0, mesh_size_boundary)

# Create circular arcs
cb_1 = gmsh.model.geo.addCircleArc(pb_1, p_center_b, pb_3)
cb_2 = gmsh.model.geo.addCircleArc(pb_3, p_center_b, pb_2)
cb_3 = gmsh.model.geo.addCircleArc(pb_2, p_center_b, pb_4)
cb_4 = gmsh.model.geo.addCircleArc(pb_4, p_center_b, pb_1)

# Create a curve loop from the arcs
cl_boundary = gmsh.model.geo.addCurveLoop([cb_1, cb_2, cb_3, cb_4])

# -- Define Surfaces --
# Create the main surface (vacuum) with holes for the conductors
s_vacuum = gmsh.model.geo.addPlaneSurface([cl_boundary, cl1, cl2])

# Synchronize the CAD model
gmsh.model.geo.synchronize()

# --- Define Physical Groups ---

# Physical group for Conductor 1 (Curve)
gmsh.model.addPhysicalGroup(1, [cl1], name="Conductor_1")
# Physical group for Conductor 2 (Curve)
gmsh.model.addPhysicalGroup(1, [cl2], name="Conductor_2")
# Physical group for the Open Boundary (Curve)
gmsh.model.addPhysicalGroup(1, [cl_boundary], name="OpenBoundary_0")
# Physical group for the Vacuum (Surface)
gmsh.model.addPhysicalGroup(2, [s_vacuum], name="Vacuum")


# --- Generate Mesh ---

# Força o Gmsh a salvar o arquivo no formato da versão 2.2, que é compatível com MFEM
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

# Generate the 2D mesh
gmsh.model.mesh.generate(2)

# --- Save Mesh ---

# Save the mesh to a file
gmsh.write("two_wires_open_gmsh.msh")

# --- Finalize ---

# To visualize the model we can use the Gmsh GUI:
if 'close' not in sys.argv:
    gmsh.fltk.run()

# Finalize Gmsh
gmsh.finalize()