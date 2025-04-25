# biotracker

<a href="https://zenodo.org/doi/10.5281/zenodo.13121947"><img src="https://zenodo.org/badge/421445542.svg" alt="DOI"></a>

Java-based software to perform particle tracking of biological particles using hydrodynamic models. Originally developed to use output fields from FVCOM in the Scottish west coast region to study the movement of sea lice between fish farm sites.

Built with [OpenJDK 23.0.2](https://jdk.java.net/23/). 

Capabilities:
- Temperature-dependent initialization (egg production)
- Salinity-dependent mortality rates
- Temperature-dependent development
- 3D RK4 advection and diffusion
- Active vertical movement (salinity, light, time triggers)
- Movement between nested FVCOM meshes
- Parallel movement & .nc file reading
- Continual connectivity tracking at two depth ranges

Full functionality is only implemented for FVCOM hydrodynamics. 

See [biotrackR](https://github.com/Sz-Tim/biotrackR) for R-based interface.
