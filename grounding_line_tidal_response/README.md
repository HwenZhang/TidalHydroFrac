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
fenics-grounding-line-tidal-response A=3e-24

fenics-grounding-line-tidal-response-2 test for relaxation parameter

docker run -ti --name fenics-grounding-line-tidal-response -w /home/fenics -v /Users/hwenzhang/fenicsproject/grounding_line_tidal_response:/home/fenics/shared -d -p 127.0.0.1:8893:8893  quay.io/fenicsproject/stable:latest

docker exec -ti -u fenics fenics-grounding-line-tidal-response /bin/bash -l

docker run -ti --name fenics-grounding-line-tidal-response-2 -w /home/fenics -v /Users/hwenzhang/fenicsproject/grounding_line_tidal_response:/home/fenics/shared -d -p 127.0.0.1:8897:8897  quay.io/fenicsproject/stable:latest

docker exec -ti -u fenics fenics-grounding-line-tidal-response-2 /bin/bash -l

docker run -ti --name fenics-grounding-line-tidal-response-3 -w /home/fenics -v /Users/hwenzhang/fenicsproject/grounding_line_tidal_response:/home/fenics/shared -d -p 127.0.0.1:8891:8891  quay.io/fenicsproject/stable:latest

docker exec -ti -u fenics fenics-grounding-line-tidal-response-3 /bin/bash -l

scp -r /Users/eart0487/fenicsproject/grounding_line_tidal_response hwenzhang@163.1.22.143:/Users/hwenzhang/fenicsproject

scp -r hwenzhang@163.1.22.143:/Users/hwenzhang/fenicsproject/grounding_line_tidal_response/results/stokes_tidal_response_U09ma_L20000_A3e_24_delta6e_19_tide1_00_relaxation0.4 /Users/eart0487/fenicsproject/grounding_line_tidal_response/results

scp -r hwenzhang@163.1.22.143:/Users/hwenzhang/fenicsproject/grounding_line_tidal_response /Users/eart0487/fenicsproject