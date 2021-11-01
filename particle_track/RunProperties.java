/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.io.*;
import java.util.Properties;

import static particle_track.ISO_datestr.dateIntParse;

/**
 * @author sa01ta
 */
public class RunProperties {
    String basedir, sitedir, // Location of code and site data files
            datadir, datadirPrefix, datadirSuffix, // Location of default hydrodynamic data, with prefix and suffix for finding annual subdirectories
            datadir2, datadir2Prefix, datadir2Suffix, // Location of secondary (larger domain) hydrodynamic data, with prefix and suffix for finding annual subdirectories 
            mesh1, mesh2, // Full path to the mesh files used describing spatial structure of the hydrodynamic data (
            mesh1Type, mesh2Type, // What type of meshes are being read in (FVCOM or ROMS)
            restartParticles, // Full path to file containing locations of particles for a hot restart (matches last hour of locations file)
            location, sitefile, sitefileEnd, habitat, suffix, species, // Descriptive strings
            coordRef, // Coordinate reference system
            seasonalDensityPath; // Path + filename for month-specific particle start densities; defaults to "" = 1 for all particles


    boolean backwards, // run model backwards? Needs some work on loops to make this work correctly
            timeInterpolate, spatialInterpolate, // interpolate between hydro file data values?
            rk4, // use RK4 numerical integration (alternative is Euler; need about 10 times as many steps)
            parallel, // use multiple cores to speed up run?
            diffusion, variableDiff, // include random walk, use diffusion parameter from hydro output?
            salinityMort, // mortality calculated based on local salinity (sea lice - doesn't presently do anything)?
            endOnArrival, // stop at first suitable habitat site, or simply note arrival and move on?
            setStartDepth, fixDepth, // set particle depth at initiation?
            verticalDynamics, // should vertical dynamics be included? Turns on velocities, diffusion, and swimming behaviour; ONLY implemented for FVCOM
            readHydroVelocityOnly, // read only u,v from hydro files (saves RAM, ignores random extra variables)
            recordPsteps, splitPsteps, // record particle element densities? split by source site?
            recordConnectivity, recordLocations, recordArrivals, // record connectivity? particle locations? arrivals at sites?
            duplicateLastDay, // Should hydro file for last day be duplicated for interpolation puproses during last hour of simulation (false except when in operational mode)    
            checkOpenBoundaries, // Should open boundaries be checked? If reading hydro mesh from file directly, the answer is currently NO (open boundaries treated as closed boundaries).
            verboseSetUp;

    int start_ymd, end_ymd, numberOfDays, // Start and end of run. If numberOfDays = 0, it is ignored and end_ymd is used instead
            releaseScenario, // 0 release all at "releaseTime", 1 continuous release ("nparts" per hour per site)
            nparts, // Number of particles released per site (per hour in releaseScenario == 1
            recordsPerFile1, // Number of records per velocity file (allow two velocity files with different sizes)
            stepsPerStep, // Number of increments between each velocity record (also for time interpolations)
            //depthLayers, // Number of depth layers in hydro output - IDEALLY DEPRECATE
            thresh, // Threshold distance for "settlement" (m)
            behaviour, // Particle behaviour - see Particle.java
            parallelThreads, // Number of threads to use in parallel execution
            minchVersion, // Another element of the filename for hydrodynamic files
            pstepsInterval, connectivityInterval; // Interval in hours between recording element density summaries, and connectivity

    double releaseTime, releaseTimeEnd, viabletime, // Time of particle release (if releaseScenario == "0") and end of particle release (if releaseScenario == 2), Time to attain settlement competency
            dt, // Timestep (s) per record
            D_h, // Horizontal diffuision parameter
            D_hVert, // Vertical diffuision parameter
            mortalityRate, // Hourly mortality rate of particles
            maxParticleAge, // Maximum age for particles. Set to <=0 and it will be ignored.
            viableDegreeDays, maxDegreeDays, // Degree x days to use for settlement viability time and mortality
            sinkingRateMean, sinkingRateStd, // Particle sinking distribution parameters
            vertSwimSpeedMean, vertSwimSpeedStd,
            salinityThreshold,
            startDepth, // Particle initiation depth
            restartParticlesCutoffDays; // when reading the specified restart particles file, cutoff in particle start date to apply (days before start date of run)

    /**
     * @param filename
     */
    public RunProperties(String filename) {
        System.out.println("GETTING PROPERTIES FROM " + filename);
        Properties properties = new Properties();
        try {
            properties.load(new FileInputStream(filename));
        } catch (IOException e) {
            System.err.println("!!!Could not find properties file!!!");
        }

        for (String key : properties.stringPropertyNames()) {
            String value = properties.getProperty(key);
            System.out.println(key + " => " + value);
        }

        basedir = properties.getProperty("basedir", "/home/sa01ta/particle_track/");
        sitedir = properties.getProperty("sitedir", "/home/sa01ta/particle_track/minch_sites/");
        datadir = properties.getProperty("datadir", "/home/sa01da/work/minch2/Archive/");
        datadirPrefix = properties.getProperty("datadirPrefix", "netcdf_");
        datadirSuffix = properties.getProperty("datadirSuffix", "");
        datadir2 = properties.getProperty("datadir2", "/home/sa01da/data/oss/");
        datadir2Prefix = properties.getProperty("datadirPrefix", "");
        datadir2Suffix = properties.getProperty("datadirSuffix", "_OLD");

        seasonalDensityPath = properties.getProperty("seasonalDensityPath", "");

        minchVersion = Integer.parseInt(properties.getProperty("minchVersion", "2"));

        mesh1 = properties.getProperty("mesh1", "/home/sa01ta/particle_track/WestCOMS_mesh.nc");
        mesh2 = properties.getProperty("mesh2", "");
        mesh1Type = properties.getProperty("mesh1Type", "");
        mesh2Type = properties.getProperty("mesh2Type", "");

        checkOpenBoundaries = Boolean.parseBoolean(properties.getProperty("checkOpenBoundaries", "false"));

        restartParticles = properties.getProperty("restartParticles", "");
        restartParticlesCutoffDays = Double.parseDouble(properties.getProperty("restartParticlesCutoffDays", "21"));

        sitefile = properties.getProperty("sitefile", "startlocations.dat");
        sitefileEnd = properties.getProperty("sitefileEnd", sitefile);

        location = properties.getProperty("location", "minch");
        habitat = properties.getProperty("habitat", "");
        suffix = properties.getProperty("suffix", "");

        verboseSetUp = Boolean.parseBoolean(properties.getProperty("verboseSetUp", "true"));

        coordRef = properties.getProperty("coordRef", "WGS84");

        start_ymd = Integer.parseInt(properties.getProperty("start_ymd", "20180101"));
        numberOfDays = Integer.parseInt(properties.getProperty("numberOfDays", "0"));
        if (numberOfDays > 0) {
            int[] startDate = ISO_datestr.dateIntParse(start_ymd);
            ISO_datestr tempIsoDate = new ISO_datestr(startDate[0], startDate[1], startDate[2]);
            for (int i = 1; i < numberOfDays; i++) {
                tempIsoDate.addDay();
            }
            end_ymd = Integer.parseInt(tempIsoDate.getDateStr());
        } else {
            end_ymd = Integer.parseInt(properties.getProperty("end_ymd", "20180102"));
        }
        if (end_ymd < start_ymd) {
            System.err.println("***** End date before start date! *****");
            System.exit(1);
        }

        // Set variable values based on contents of the properties object
        backwards = Boolean.parseBoolean(properties.getProperty("backwards", "false"));
        // DEPRECATED
        //timeInterpolate = Boolean.parseBoolean(properties.getProperty("timeInterpolate","true"));
        //spatialInterpolate = Boolean.parseBoolean(properties.getProperty("spatialInterpolate","true"));

        rk4 = Boolean.parseBoolean(properties.getProperty("rk4", "true"));

        diffusion = Boolean.parseBoolean(properties.getProperty("diffusion", "true"));
        variableDiff = Boolean.parseBoolean(properties.getProperty("variableDiff", "false"));
        salinityMort = Boolean.parseBoolean(properties.getProperty("salinityMort", "false"));
        endOnArrival = Boolean.parseBoolean(properties.getProperty("endOnArrival", "false"));
        // DEPRECATED
        //tidalRelease = Boolean.parseBoolean(properties.getProperty("tidalRelease","false"));
        setStartDepth = Boolean.parseBoolean(properties.getProperty("setStartDepth", "false"));
        fixDepth = Boolean.parseBoolean(properties.getProperty("fixDepth", "false"));
        startDepth = Integer.parseInt(properties.getProperty("startDepth", "0"));

        readHydroVelocityOnly = Boolean.parseBoolean(properties.getProperty("readHydroVelocityOnly", "false"));
        duplicateLastDay = Boolean.parseBoolean(properties.getProperty("duplicateLastDay", "false"));

        parallelThreads = Integer.parseInt(properties.getProperty("parallelThreads", "4"));
        parallel = parallelThreads > 1;
        // DEPRECATED
        //oldOutput = Boolean.parseBoolean(properties.getProperty("oldOutput"));

        recordPsteps = Boolean.parseBoolean(properties.getProperty("recordPsteps", "true"));
        splitPsteps = Boolean.parseBoolean(properties.getProperty("splitPsteps", "true"));
        pstepsInterval = Integer.parseInt(properties.getProperty("pstepsInterval", "24"));
        recordConnectivity = Boolean.parseBoolean(properties.getProperty("recordConnectivity", "true"));
        connectivityInterval = Integer.parseInt(properties.getProperty("connectivityInterval", "24"));
        recordLocations = Boolean.parseBoolean(properties.getProperty("recordLocations", "true"));
        recordArrivals = Boolean.parseBoolean(properties.getProperty("recordArrivals", "true"));

        releaseScenario = Integer.parseInt(properties.getProperty("releaseScenario", "0"));
        nparts = Integer.parseInt(properties.getProperty("nparts", "5"));
        releaseTime = Double.parseDouble(properties.getProperty("releaseTime", "0"));
        releaseTimeEnd = Double.parseDouble(properties.getProperty("releaseTimeEnd", "24"));

        dt = Double.parseDouble(properties.getProperty("dt", "3600"));

        // 21/11/2018 --- Ideally want to remove this. If using two models for hydro 
        // would ideally sense from the hydro file
        recordsPerFile1 = Integer.parseInt(properties.getProperty("recordsPerFile1", "25"));
        //recordsPerFile2 = Integer.parseInt(properties.getProperty("recordsPerFile2","4"));

        stepsPerStep = Integer.parseInt(properties.getProperty("stepsPerStep", "25"));

        // 21/11/2018 --- possible to remove? If using two models for hydro 
        // would ideally sense from the hydro file
        //depthLayers = Integer.parseInt(properties.getProperty("depthLayers","10"));       

        verticalDynamics = Boolean.parseBoolean(properties.getProperty("verticalDynamics", "false"));
        thresh = Integer.parseInt(properties.getProperty("thresh", "500"));
        behaviour = Integer.parseInt(properties.getProperty("behaviour", "1"));

        D_h = Double.parseDouble(properties.getProperty("D_h", "0.1"));
        D_hVert = Double.parseDouble(properties.getProperty("D_hVert", "0.005"));
        mortalityRate = Double.parseDouble(properties.getProperty("mortalityRate", "0.01"));

        species = properties.getProperty("species", "none");
        viabletime = Double.parseDouble(properties.getProperty("viabletime", "86"));
        maxParticleAge = Double.parseDouble(properties.getProperty("maxParticleAge", "336"));

        viableDegreeDays = Double.parseDouble(properties.getProperty("viableDegreeDays", "-1"));
        maxDegreeDays = Double.parseDouble(properties.getProperty("maxDegreeDays", "-1"));
        if (viableDegreeDays != -1) {
            //viabletime = viableDegreeDays*20;
            viabletime = -1;
            System.out.println("viableDegreeDays entered; set viabletime=" + viabletime + " so won't be used at 56N!");
        }
        if (maxDegreeDays != -1) {
            //maxParticleAge = maxDegreeDays*20;
            maxParticleAge = -1;
            System.out.println("maxDegreeDays entered; set maxParticleAge=" + maxParticleAge + " so won't be used at 56N!");
        }
        if ((viableDegreeDays != -1 || maxParticleAge != -1) && readHydroVelocityOnly) {
            System.err.println("readHydroVelocityOnly==true AND trying to use degree-days for development => won't develop or die!");
        }

        vertSwimSpeedMean = Double.parseDouble(properties.getProperty("vertSwimSpeedMean", "0"));
        vertSwimSpeedStd = Double.parseDouble(properties.getProperty("vertSwimSpeedStd", "0"));

        sinkingRateMean = Double.parseDouble(properties.getProperty("sinkingRateMean", "0"));
        sinkingRateStd = Double.parseDouble(properties.getProperty("sinkingRateStd", "0"));

        salinityThreshold = Double.parseDouble(properties.getProperty("salinityThreshold", "0"));

        properties.list(System.out);

    }

}



