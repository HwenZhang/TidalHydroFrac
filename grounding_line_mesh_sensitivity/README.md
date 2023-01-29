# TidalHydroFrac

Codes for the Amery Tidal Hydrofracture project

This is the tidal-response branch for viscoelastic grounding line migration.
Docker container name: fenics-grounding-line-tidal-response

Rheology:
Incompressible Viscoelastic : A=3e-24, n=3, mu=3.8e9
eta_max = 1e14 Pa s

Inflow: U0 = 9m/y in x-direction
Length = 20km
Height = 500m
Tidal amplitude = 1.0 m

Container name:

docker run -ti --name fenics-grounding-line-mesh-test -w /home/fenics -v /Users/hwenzhang/fenicsproject/grounding_line_mesh_sensitivity:/home/fenics/shared -d -p 127.0.0.1:8888:8888  quay.io/fenicsproject/stable:latest

docker exec -ti -u fenics fenics-grounding-line-mesh-test /bin/bash -l

docker run -ti --name fenics-grounding-line-mesh-test-2 -w /home/fenics -v /Users/hwenzhang/fenicsproject/grounding_line_mesh_sensitivity:/home/fenics/shared -d -p 127.0.0.1:8889:8889  quay.io/fenicsproject/stable:latest

docker exec -ti -u fenics fenics-grounding-line-mesh-test-2 /bin/bash -l



scp -r /Users/eart0487/fenicsproject/grounding_line_tidal_response hwenzhang@163.1.22.143:/Users/hwenzhang/fenicsproject

scp -r hwenzhang@163.1.22.143:/Users/hwenzhang/fenicsproject/grounding_line_tidal_response/results /Users/eart0487/fenicsproject/grounding_line_tidal_response/results

scp -r hwenzhang@163.1.22.143:/Users/hwenzhang/fenicsproject/grounding_line_mesh_sensitivity/results/stokes_tidal_response_U09ma_L20000_A3e_24_n3_0_eta1e14_tide1_00_C1e_7_DX12 /Users/eart0487/fenicsproject/grounding_line_mesh_sensitivity/results

scp -r hwenzhang@163.1.22.143:/Users/hwenzhang/fenicsproject/grounding_line_mesh_sensitivity/results/stokes_tidal_response_U09ma_L20000_A3e_24_n3_0_delta1e_15_tide1_00_C1e_7_DX12_smooth /Users/eart0487/fenicsproject/grounding_line_mesh_sensitivity/results