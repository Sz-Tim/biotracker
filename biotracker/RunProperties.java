/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

/**
 * @author sa01ta
 */
public class RunProperties {
    String datadir, datadirPrefix, datadirSuffix, // Location of default hydrodynamic data, with prefix and suffix for finding annual subdirectories
            datadir2, datadir2Prefix, datadir2Suffix, // Location of secondary (larger domain) hydrodynamic data, with prefix and suffix for finding annual subdirectories 
            mesh1, mesh2, // Full path to the mesh files used describing spatial structure of the hydrodynamic data
            mesh1Type, mesh2Type, // What type of meshes are being read in (FVCOM or ROMS)
            restartParticles, // Full path to file containing locations of particles for a hot restart (matches last hour of locations file)
            location, location2, sitefile, sitefileEnd, habitat, suffix, species, // Descriptive strings
            siteDensityPath, // Path + filename for daily start densities for each site; defaults to "" = 1 for all particles; col1 = siteNames, col2:N = dates
            daylightPath; // Path + filename for sunrise / sunset hours; defaults to "" = ignore

    boolean backwards, // run model backwards? Needs some work on loops to make this work correctly
            rk4, // use RK4 numerical integration (alternative is Euler; need about 10 times as many steps)
            parallel, // use multiple cores to speed up run?
            diffusion, variableDh, variableDhV, // include random walk, use diffusion parameter from hydro output?
            salinityMort, // mortality calculated based on local salinity
            endOnArrival, // stop at first suitable habitat site, or simply note arrival and move on?
            setStartDepth, // set particle depth at initiation?
            fixDepth,
            swimLightLevel,
            readHydroVelocityOnly, // read only u,v from hydro files (saves RAM, ignores random extra variables)
            recordPsteps, splitPsteps, // record particle element densities? split by source site?
            recordConnectivity, recordLocations, recordArrivals, // record connectivity? particle locations? arrivals at sites?
            recordMovement, recordActivity, // record all movements for a sample of particles? Record sink/swim/float counts within each element and hour?
            recordVertDistr, // record vertical distributions?
            duplicateLastDay, // Should hydro file for last day be duplicated for interpolation purposes during last hour of simulation (false except when in operational mode)
            checkOpenBoundaries, // Should open boundaries be checked? If reading hydro mesh from file directly, the answer is currently NO (open boundaries treated as closed boundaries).
            verboseSetUp,
            needS, needT, needZeta, needK, needLight, needVh, // Internal use: which hydrodynamic variables need to be loaded?
            coordOS, // Coordinate reference system
            FVCOM; // Is mesh FVCOM? If not, not all functions are supported

    ISO_datestr start_ymd, end_ymd;

    int numberOfDays, // Start and end of run. If numberOfDays = 0, it is ignored and end_ymd is used instead
            releaseScenario, // 0 release all at "releaseTime", 1 continuous release ("nparts" per releaseInterval hours per site)
            nparts, // Number of particles released per site (per hour in releaseScenario == 1
            recordsPerFile1, // Number of records per velocity file (allow two velocity files with different sizes)
            stepsPerStep, // Number of increments between each velocity record (also for time interpolations)
            thresh, // Threshold distance for "settlement" (m)
            parallelThreads, // Number of threads to use in parallel execution
            minchVersion, minchVersion2, // Another element of the filename for hydrodynamic files
            pstepsInterval, connectivityInterval,  // Interval in hours between recording element density summaries, connectivity
            vertDistrInterval, vertDistrMax; // Interval in hours for recording vertDistr, max depth for bins (1 m bins from 0 to vertDistrMax)

    double releaseTime, releaseTimeEnd, viabletime, // Time of particle release (if releaseScenario == "0") and end of particle release (if releaseScenario == 2), Time to attain settlement competency
            dt, // Time step (s) per record
            openBoundaryThresh, // distance threshold (m) to open boundary nodes for particles to be ejected
            D_h, // Horizontal diffusion parameter
            D_hVert, // Vertical diffusion parameter
            mortalityRate, // Hourly mortality rate of particles
            maxParticleAge, // Maximum age for particles. Set to <=0 to ignore.
            viableDegreeDays, maxDegreeDays, // Degree x days to use for settlement viability time and mortality
            vertSwimSpeedMean, vertSwimSpeedStd,
            vertSwimSpeedCopepodidMean, vertSwimSpeedCopepodidStd,
            vertSwimSpeedNaupliusMean, vertSwimSpeedNaupliusStd,
            sinkingRateMean, sinkingRateStd, // Particle sinking distribution parameters
            sinkingRateCopepodidMean, sinkingRateCopepodidStd,
            sinkingRateNaupliusMean, sinkingRateNaupliusStd,
            salinityThreshold, salinityThreshMin, salinityThreshMax, // 1) sink below threshold; 2-3) Sanvdik 2020 A3: linear increase in prSink from Max (none sink) to Min (all sink)
            startDepth, // Particle initiation depth
            maxDepth, // maximum particle depth
            pstepsMaxDepth, // maximum depth for recording particle density in psteps output
            releaseInterval, // release frequency in hours
            restartParticlesCutoffDays; // when reading the specified restart particles file, cutoff in particle start date to apply (days before start date of run)

    public RunProperties(String filename) {
        System.out.println("Getting properties from " + filename);
        Properties properties = new Properties();
        try {
            properties.load(new FileInputStream(filename));
        } catch (IOException e) {
            System.err.println("--- Could not find properties file - check filename and working directory ---");
        }
        for (String key : properties.stringPropertyNames()) {
            String value = properties.getProperty(key);
            System.out.println(key + " => " + value);
        }

        // Directories
        datadir = properties.getProperty("datadir", "/media/archiver/common/sa01da-work/WeStCOMS2/Archive/");
        datadirPrefix = properties.getProperty("datadirPrefix", "netcdf_");
        datadirSuffix = properties.getProperty("datadirSuffix", "");
        datadir2 = properties.getProperty("datadir2", "");
        datadir2Prefix = properties.getProperty("datadirPrefix", "");
        datadir2Suffix = properties.getProperty("datadirSuffix", "");

        // Geography, hydrodynamic files, & mesh files
        coordOS = Boolean.parseBoolean(properties.getProperty("coordOS", "true"));
        location = properties.getProperty("location", "westcoms");
        location2 = properties.getProperty("location2", "westcoms");
        minchVersion = Integer.parseInt(properties.getProperty("minchVersion", "2"));
        minchVersion2 = Integer.parseInt(properties.getProperty("minchVersion2", "2"));
        recordsPerFile1 = Integer.parseInt(properties.getProperty("recordsPerFile1", "25"));
        mesh1 = properties.getProperty("mesh1", "/home/sa04ts/FVCOM_meshes/WeStCOMS2_mesh.nc");
        mesh2 = properties.getProperty("mesh2", "");
        mesh1Type = properties.getProperty("mesh1Type", "");
        mesh2Type = properties.getProperty("mesh2Type", "");
        FVCOM = properties.getProperty("mesh1Type").equals("FVCOM");
        checkOpenBoundaries = Boolean.parseBoolean(properties.getProperty("checkOpenBoundaries", "false"));
        openBoundaryThresh = Double.parseDouble(properties.getProperty("openBoundaryThresh", "500"));

        // Sites
        sitefile = properties.getProperty("sitefile", "startlocations.dat");
        sitefileEnd = properties.getProperty("sitefileEnd", sitefile);
        habitat = properties.getProperty("habitat", "");
        suffix = properties.getProperty("suffix", "");
        verboseSetUp = Boolean.parseBoolean(properties.getProperty("verboseSetUp", "true"));

        // Dates
        start_ymd = new ISO_datestr(properties.getProperty("start_ymd", "2090401"));
        numberOfDays = Integer.parseInt(properties.getProperty("numberOfDays", "0"));
        if (numberOfDays > 0) {
            ISO_datestr tempIsoDate = new ISO_datestr(properties.getProperty("start_ymd", "20190401"));
            for (int i = 1; i < numberOfDays; i++) {
                tempIsoDate.addDay();
            }
            end_ymd = new ISO_datestr(tempIsoDate);
        } else {
            end_ymd = new ISO_datestr(properties.getProperty("end_ymd", "20190402"));
            numberOfDays = end_ymd.getDateNum() - start_ymd.getDateNum() + 1;
        }

        // Run parameters
        backwards = Boolean.parseBoolean(properties.getProperty("backwards", "false"));
        parallelThreads = Integer.parseInt(properties.getProperty("parallelThreads", "4"));
        readHydroVelocityOnly = Boolean.parseBoolean(properties.getProperty("readHydroVelocityOnly", "false"));
        duplicateLastDay = Boolean.parseBoolean(properties.getProperty("duplicateLastDay", "false"));
        dt = Double.parseDouble(properties.getProperty("dt", "3600"));
        fixDepth = Boolean.parseBoolean(properties.getProperty("fixDepth", "false"));
        maxDepth = Double.parseDouble(properties.getProperty("maxDepth", "10000"));

        // Release
        restartParticles = properties.getProperty("restartParticles", "");
        restartParticlesCutoffDays = Double.parseDouble(properties.getProperty("restartParticlesCutoffDays", "21"));
        siteDensityPath = properties.getProperty("siteDensityPath", "");
        setStartDepth = Boolean.parseBoolean(properties.getProperty("setStartDepth", "false"));
        startDepth = Integer.parseInt(properties.getProperty("startDepth", "0"));
        releaseScenario = Integer.parseInt(properties.getProperty("releaseScenario", "0"));
        releaseInterval = Double.parseDouble(properties.getProperty("releaseInterval", "1"));
        nparts = Integer.parseInt(properties.getProperty("nparts", "5"));
        releaseTime = Double.parseDouble(properties.getProperty("releaseTime", "0"));
        releaseTimeEnd = Double.parseDouble(properties.getProperty("releaseTimeEnd", "24"));

        // Arrival
        endOnArrival = Boolean.parseBoolean(properties.getProperty("endOnArrival", "false"));
        thresh = Integer.parseInt(properties.getProperty("thresh", "500"));

        // Advection & diffusion
        rk4 = Boolean.parseBoolean(properties.getProperty("rk4", "true"));
        stepsPerStep = Integer.parseInt(properties.getProperty("stepsPerStep", "30"));
        diffusion = Boolean.parseBoolean(properties.getProperty("diffusion", "true"));
        variableDh = Boolean.parseBoolean(properties.getProperty("variableDh", "false"));
        variableDhV = Boolean.parseBoolean(properties.getProperty("variableDhV", "true"));
        D_h = Double.parseDouble(properties.getProperty("D_h", "0.1"));
        D_hVert = Double.parseDouble(properties.getProperty("D_hVert", "0.001"));
        if (fixDepth) {
            D_hVert = 0;
            variableDhV = false;
            System.out.println("fixDepth = true; Setting D_hVert = 0 and variableDiffusion = false.");
        }

        // Behaviour
        species = properties.getProperty("species", "none");
        daylightPath = properties.getProperty("daylightPath", "");
        swimLightLevel = Boolean.parseBoolean(properties.getProperty("swimLightLevel", "false"));
        vertSwimSpeedMean = Double.parseDouble(properties.getProperty("vertSwimSpeedMean", "0"));
        vertSwimSpeedStd = Double.parseDouble(properties.getProperty("vertSwimSpeedStd", "0"));
        vertSwimSpeedCopepodidMean = Double.parseDouble(properties.getProperty("vertSwimSpeedCopepodidMean", "" + vertSwimSpeedMean));
        vertSwimSpeedCopepodidStd = Double.parseDouble(properties.getProperty("vertSwimSpeedCopepodidStd", "" + vertSwimSpeedStd));
        vertSwimSpeedNaupliusMean = Double.parseDouble(properties.getProperty("vertSwimSpeedNaupliusMean", "" + vertSwimSpeedMean));
        vertSwimSpeedNaupliusStd = Double.parseDouble(properties.getProperty("vertSwimSpeedNaupliusStd", "" + vertSwimSpeedStd));
        sinkingRateMean = Double.parseDouble(properties.getProperty("sinkingRateMean", "0"));
        sinkingRateStd = Double.parseDouble(properties.getProperty("sinkingRateStd", "0"));
        sinkingRateCopepodidMean = Double.parseDouble(properties.getProperty("sinkingRateCopepodidMean", "" + sinkingRateMean));
        sinkingRateCopepodidStd = Double.parseDouble(properties.getProperty("sinkingRateCopepodidStd", "" + sinkingRateStd));
        sinkingRateNaupliusMean = Double.parseDouble(properties.getProperty("sinkingRateNaupliusMean", "" + sinkingRateMean));
        sinkingRateNaupliusStd = Double.parseDouble(properties.getProperty("sinkingRateNaupliusStd", "" + sinkingRateStd));
        salinityThreshold = Double.parseDouble(properties.getProperty("salinityThreshold", "0"));
        salinityThreshMin = Double.parseDouble(properties.getProperty("salinityThreshMin", "" + salinityThreshold));
        salinityThreshMax = Double.parseDouble(properties.getProperty("salinityThreshMax", "" + salinityThreshold));

        // Demographics
        salinityMort = Boolean.parseBoolean(properties.getProperty("salinityMort", "true"));
        mortalityRate = Double.parseDouble(properties.getProperty("mortalityRate", "0.01"));
        viabletime = Double.parseDouble(properties.getProperty("viabletime", "-1"));
        maxParticleAge = Double.parseDouble(properties.getProperty("maxParticleAge", "-1"));
        viableDegreeDays = Double.parseDouble(properties.getProperty("viableDegreeDays", "40"));
        maxDegreeDays = Double.parseDouble(properties.getProperty("maxDegreeDays", "150"));
        if (viableDegreeDays != -1) {
            viabletime = -1;
            System.out.println("viableDegreeDays entered; set viabletime=" + viabletime + " so won't be used at 56N!");
        }
        if (maxDegreeDays != -1) {
            maxParticleAge = -1;
            System.out.println("maxDegreeDays entered; set maxParticleAge=" + maxParticleAge + " so won't be used at 56N!");
        }
        if ((viableDegreeDays != -1 || maxParticleAge != -1) && readHydroVelocityOnly) {
            System.err.println("readHydroVelocityOnly==true AND trying to use degree-days for development => won't develop or die!");
        }

        // Output
        recordPsteps = Boolean.parseBoolean(properties.getProperty("recordPsteps", "true"));
        splitPsteps = Boolean.parseBoolean(properties.getProperty("splitPsteps", "true"));
        pstepsInterval = Integer.parseInt(properties.getProperty("pstepsInterval", "24"));
        pstepsMaxDepth = Double.parseDouble(properties.getProperty("pstepsMaxDepth", "10000"));
        recordConnectivity = Boolean.parseBoolean(properties.getProperty("recordConnectivity", "true"));
        connectivityInterval = Integer.parseInt(properties.getProperty("connectivityInterval", "24"));
        recordLocations = Boolean.parseBoolean(properties.getProperty("recordLocations", "true"));
        recordArrivals = Boolean.parseBoolean(properties.getProperty("recordArrivals", "true"));
        recordMovement = Boolean.parseBoolean(properties.getProperty("recordMovement", "false"));
        recordActivity = Boolean.parseBoolean(properties.getProperty("recordElemActivity", "false"));
        recordVertDistr = Boolean.parseBoolean(properties.getProperty("recordVertDistr", "false"));
        vertDistrInterval = Integer.parseInt(properties.getProperty("vertDistrInterval", "1"));
        vertDistrMax = Integer.parseInt(properties.getProperty("vertDistrMax", "20"));

        // hydrodynamic file requirements
        needS = !fixDepth || salinityMort;
        needT = viableDegreeDays > -1;
        needZeta = false;
        needK = !fixDepth && variableDhV;
        needLight = !fixDepth;
        needVh = variableDh;

        properties.list(System.out);
    }


    public void checkForConflictingProperties() {
        boolean conflict = false;
        String conflictingProperties = "Error: Conflicting properties\n";
        if (variableDhV && !diffusion) {
            conflict = true;
            conflictingProperties += "  variableDhV && !diffusion\n";
        }
        if (start_ymd.isLaterThan(end_ymd)) {
            conflict = true;
            conflictingProperties += "  start_ymd.isLaterThan(end_ymd)\n";
        }
        if (!fixDepth && !mesh1Type.equals("FVCOM")) {
            conflict = true;
            conflictingProperties += "  verticalDynamics && !mesh1Type.equals(\"FVCOM\")\n";
        }

        if(conflict) {
            System.err.println(conflictingProperties);
            System.exit(1);
        }
    }
}



