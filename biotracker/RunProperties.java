/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Properties;

/**
 * @author sa01ta
 */
public class RunProperties {
    String mesh0, mesh1, // Full path to the mesh files used describing spatial structure of the hydrodynamic data
            meshType0, meshType1, // What type of meshes are being read in (FVCOM or ROMS)
            hfDir0, hfDirPrefix0, hfDirSuffix0, // Location of default hydrodynamic data, with prefix and suffix for finding annual subdirectories
            hfDir1, hfDirPrefix1, hfDirSuffix1, // Location of secondary (larger domain) hydrodynamic data, with prefix and suffix for finding annual subdirectories
            hfDir2, hfDirPrefix2, hfDirSuffix2, // Location of wave data for Stokes drift, with prefix and suffix for finding annual subdirectories
            hfFilePrefix0, hfFilePrefix1, hfFilePrefix2, // Mesh domain location as in the .nc filenames (westcoms2, etive28)
            restartParticles, // Full path to file containing locations of particles for a hot restart (matches last hour of locations file)
            sitefile, sitefileEnd, habitat, suffix, species, // Descriptive strings
            siteDensityPath, // Path + filename for daily start densities for each site; defaults to "" = 1 for all particles; col1 = siteNames, col2:N = dates
            daylightPath, // Path + filename for sunrise / sunset hours; defaults to "" = ignore
            eggTemp_fn, // Function to relate egg production to temperature; "constant" (default), "linear", "quadratic", or "logistic"
            mortSal_fn; // Function to relate mortality to salinity; "constant" (default), "linear", or "quadratic"

    boolean backwards, // run model backwards? Needs some work on loops to make this work correctly
            rk4, // use RK4 numerical integration (alternative is Euler; need about 10 times as many steps)
            diffusion, variableDh, variableDhV, // include random walk, use diffusion parameter from hydro output?
            endOnArrival, // stop at first suitable habitat site, or simply note arrival and move on?
            fixDepth,
            swimLightLevel,
            readHydroVelocityOnly, // read only u,v from hydro files (saves RAM, ignores random extra variables)
            stokesDrift,
            recordImmature,
            recordPsteps, splitPsteps, // record particle element densities? split by source site?
            recordConnectivity, // record connectivity between specific depths?
            connectImmature, // record connectivity for non-infectious stage?
            recordConnectivityDepth2, // record connectivity at a second set of depths?
            recordLocations, recordArrivals, // record particle locations? arrivals at sites?
            recordMovement, recordActivity, // record all movements for a sample of particles? Record sink/swim/float counts within each element and hour?
            recordVertDistr, // record vertical distributions?
            duplicateLastDay, // Should hydro file for last day be duplicated for interpolation purposes during last hour of simulation (false except when in operational mode)
            checkOpenBoundaries, // Should open boundaries be checked? If reading hydro mesh from file directly, the answer is currently NO (open boundaries treated as closed boundaries).
            verboseSetUp,
            needS, needT, needZeta, needK, needLight, needVh, needStokes, // Internal use: which hydrodynamic variables need to be loaded?
            coordOS, // Coordinate reference system
            FVCOM; // Is mesh FVCOM? If not, not all functions are supported

    ISO_datestr start_ymd, end_ymd;

    int numberOfDays, // Start and end of run. If numberOfDays = 0, it is ignored and end_ymd is used instead
            releaseScenario, // 0 release all at "releaseTime", 1 continuous release ("nparts" per releaseInterval hours per site)
            nparts, // Number of particles released per site (per hour in releaseScenario == 1
            recordsPerFile1, // Number of records per velocity file (allow two velocity files with different sizes)
            stepsPerStep, // Number of increments between each velocity record (also for time interpolations)
            parallelThreads, // Number of threads to use in parallel execution
            parallelThreadsHD, // Number of cores for reading HD files; Fails inexplicably on Salmon occasionally
            pstepsInterval, connectivityInterval,  // Interval in hours between recording element density summaries, connectivity
            vertDistrInterval, vertDistrMax; // Interval in hours for recording vertDistr, max depth for bins (1 m bins from 0 to vertDistrMax)

    double releaseTime, releaseTimeEnd, viabletime, // Time of particle release (if releaseScenario == "0") and end of particle release (if releaseScenario == 2), Time to attain settlement competency
            dt, // Time step (s) per record
            openBoundaryThresh, // distance threshold (m) to open boundary nodes for particles to be ejected
            D_h, // Horizontal diffusion parameter
            D_hVert, // Vertical diffusion parameter
            maxParticleAge, // Maximum age for particles. Set to <=0 to ignore.
            viableDegreeDays, maxDegreeDays, // Degree x days to use for settlement viability time and mortality
            lightThreshCopepodid, lightThreshNauplius,
            swimUpSpeedMean, swimUpSpeedStd,
            swimUpSpeedCopepodidMean, swimUpSpeedCopepodidStd,
            swimUpSpeedNaupliusMean, swimUpSpeedNaupliusStd,
            swimDownSpeedMean, swimDownSpeedStd, // Particle sinking distribution parameters
            swimDownSpeedCopepodidMean, swimDownSpeedCopepodidStd,
            swimDownSpeedNaupliusMean, swimDownSpeedNaupliusStd,
            salinityThreshold, // sink below threshold
            salinityThreshCopepodidMin, salinityThreshCopepodidMax, salinityThreshNaupliusMin, salinityThreshNaupliusMax, // Sandvik 2020 A3: linear increase in prSink from Max (none sink) to Min (all sink)
            passiveSinkingIntercept, passiveSinkingSlope,
            startDepth, // Particle initiation depth
            maxDepth, // maximum particle depth
            connectivityThresh, // Threshold distance for "settlement" (m)
            connectDepth1_max, // max depth for connectivity layer 1
            connectDepth1_min, // min depth for connectivity layer 1
            connectDepth2_max, // max depth for connectivity layer 2
            connectDepth2_min, // min depth for connectivity layer 2
            pstepsMaxDepth, // maximum depth for recording particle density in psteps output
            releaseInterval, // release frequency in hours
            restartParticlesCutoffDays; // when reading the specified restart particles file, cutoff in particle start date to apply (days before start date of run)
    List<Double> eggTemp_b, // temperature dependent egg production parameters; response = eggs per day per gravid female
            mortSal_b; // salinity dependent mortality parameters; response = mortality rate per hour for larvae

    public RunProperties(String filename) {
        System.out.println("Getting properties from " + filename);
        Properties properties = new Properties();
        try {
            properties.load(new FileInputStream(filename));
        } catch (IOException e) {
            System.err.println("--- Could not find properties file - check filename and working directory ---");
        }

        // Directories
        hfDir0 = properties.getProperty("hfDir0", "/home/sa04ts/hydro/WeStCOMS2/Archive/");
        hfDirPrefix0 = properties.getProperty("hfDirPrefix0", "netcdf_");
        hfDirSuffix0 = properties.getProperty("hfDirSuffix0", "");
        hfDir1 = properties.getProperty("hfDir1", "");
        hfDirPrefix1 = properties.getProperty("hfDirPrefix1", "");
        hfDirSuffix1 = properties.getProperty("hfDirSuffix1", "");
        hfDir2 = properties.getProperty("hfDir2", "");
        hfDirPrefix2 = properties.getProperty("hfDirPrefix2", "");
        hfDirSuffix2 = properties.getProperty("hfDirSuffix2", "");

        // Geography, hydrodynamic files, & mesh files
        coordOS = Boolean.parseBoolean(properties.getProperty("coordOS", "true"));
        recordsPerFile1 = Integer.parseInt(properties.getProperty("recordsPerFile1", "25"));
        mesh0 = properties.getProperty("mesh0", "/home/sa04ts/FVCOM_meshes/WeStCOMS2_mesh.nc");
        mesh1 = properties.getProperty("mesh1", "");
        meshType0 = properties.getProperty("meshType0", "");
        meshType1 = properties.getProperty("meshType1", "");
        hfFilePrefix0 = properties.getProperty("hfFilePrefix0", "westcoms2");
        hfFilePrefix1 = properties.getProperty("hfFilePrefix1", "westcoms2");
        hfFilePrefix2 = properties.getProperty("hfFilePrefix2", "swan");
        FVCOM = properties.getProperty("meshType0").equals("FVCOM");
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
        parallelThreadsHD = Integer.parseInt(properties.getProperty("parallelThreadsHD", "1"));
        readHydroVelocityOnly = Boolean.parseBoolean(properties.getProperty("readHydroVelocityOnly", "false"));
        duplicateLastDay = Boolean.parseBoolean(properties.getProperty("duplicateLastDay", "false"));
        dt = Double.parseDouble(properties.getProperty("dt", "3600"));
        fixDepth = Boolean.parseBoolean(properties.getProperty("fixDepth", "false"));
        maxDepth = Double.parseDouble(properties.getProperty("maxDepth", "10000"));

        // Release
        restartParticles = properties.getProperty("restartParticles", "");
        restartParticlesCutoffDays = Double.parseDouble(properties.getProperty("restartParticlesCutoffDays", "21"));
        siteDensityPath = properties.getProperty("siteDensityPath", "");
        startDepth = Double.parseDouble(properties.getProperty("startDepth", "1"));
        releaseScenario = Integer.parseInt(properties.getProperty("releaseScenario", "0"));
        releaseInterval = Double.parseDouble(properties.getProperty("releaseInterval", "1"));
        nparts = Integer.parseInt(properties.getProperty("nparts", "5"));
        releaseTime = Double.parseDouble(properties.getProperty("releaseTime", "0"));
        releaseTimeEnd = Double.parseDouble(properties.getProperty("releaseTimeEnd", "24"));

        // Arrival
        endOnArrival = Boolean.parseBoolean(properties.getProperty("endOnArrival", "false"));
        connectivityThresh = Double.parseDouble(properties.getProperty("connectivityThresh", "100"));

        // Advection & diffusion
        rk4 = Boolean.parseBoolean(properties.getProperty("rk4", "true"));
        stepsPerStep = Integer.parseInt(properties.getProperty("stepsPerStep", "30"));
        diffusion = Boolean.parseBoolean(properties.getProperty("diffusion", "true"));
        variableDh = Boolean.parseBoolean(properties.getProperty("variableDh", "false"));
        variableDhV = Boolean.parseBoolean(properties.getProperty("variableDhV", "true"));
        D_h = Double.parseDouble(properties.getProperty("D_h", "0.1"));
        D_hVert = Double.parseDouble(properties.getProperty("D_hVert", "0.001"));
        stokesDrift = Boolean.parseBoolean(properties.getProperty("stokesDrift", "false"));
        if (fixDepth) {
            D_hVert = 0;
            variableDhV = false;
            System.out.println("fixDepth = true; Setting D_hVert = 0 and variableDiffusion = false.");
        }

        // Behaviour
        species = properties.getProperty("species", "none");
        daylightPath = properties.getProperty("daylightPath", "");
        swimLightLevel = Boolean.parseBoolean(properties.getProperty("swimLightLevel", "false"));
        lightThreshCopepodid = Double.parseDouble(properties.getProperty("lightThreshCopepodid", "2.06e-5"));
        lightThreshNauplius = Double.parseDouble(properties.getProperty("lightThreshNauplius", "0.392"));
        swimUpSpeedMean = Double.parseDouble(properties.getProperty("swimUpSpeedMean", "0"));
        swimUpSpeedStd = Double.parseDouble(properties.getProperty("swimUpSpeedStd", "0"));
        swimUpSpeedCopepodidMean = Double.parseDouble(properties.getProperty("swimUpSpeedCopepodidMean", "" + swimUpSpeedMean));
        swimUpSpeedCopepodidStd = Double.parseDouble(properties.getProperty("swimUpSpeedCopepodidStd", "" + swimUpSpeedStd));
        swimUpSpeedNaupliusMean = Double.parseDouble(properties.getProperty("swimUpSpeedNaupliusMean", "" + swimUpSpeedMean/2));
        swimUpSpeedNaupliusStd = Double.parseDouble(properties.getProperty("swimUpSpeedNaupliusStd", "" + swimUpSpeedStd/2));
        swimDownSpeedMean = Double.parseDouble(properties.getProperty("swimDownSpeedMean", "0"));
        swimDownSpeedStd = Double.parseDouble(properties.getProperty("swimDownSpeedStd", "0"));
        swimDownSpeedCopepodidMean = Double.parseDouble(properties.getProperty("swimDownSpeedCopepodidMean", "" + swimDownSpeedMean));
        swimDownSpeedCopepodidStd = Double.parseDouble(properties.getProperty("swimDownSpeedCopepodidStd", "" + swimDownSpeedStd));
        swimDownSpeedNaupliusMean = Double.parseDouble(properties.getProperty("swimDownSpeedNaupliusMean", "" + swimDownSpeedMean/2));
        swimDownSpeedNaupliusStd = Double.parseDouble(properties.getProperty("swimDownSpeedNaupliusStd", "" + swimDownSpeedStd/2));
        salinityThreshold = Double.parseDouble(properties.getProperty("salinityThreshold", "20"));
        salinityThreshCopepodidMin = Double.parseDouble(properties.getProperty("salinityThreshCopepodidMin", "" + salinityThreshold));
        salinityThreshCopepodidMax = Double.parseDouble(properties.getProperty("salinityThreshCopepodidMax", "" + (salinityThreshold+0.001)));
        salinityThreshNaupliusMin = Double.parseDouble(properties.getProperty("salinityThreshNaupliusMin", "" + salinityThreshold));
        salinityThreshNaupliusMax = Double.parseDouble(properties.getProperty("salinityThreshNaupliusMax", "" + (salinityThreshold+0.001)));
        passiveSinkingIntercept = Double.parseDouble(properties.getProperty("passiveSinkingIntercept", "0.001527"));
        passiveSinkingSlope = Double.parseDouble(properties.getProperty("passiveSinkingSlope", "-0.0000168"));
        eggTemp_fn = properties.getProperty("eggTemp_fn", "constant");
        String eggTemp_defaults = switch (eggTemp_fn) {
            case "constant" -> "28.2";
            // Linear model comes from data in Kragesteen 2023 Table 1
            // eggs_d = b[0] + b[1]*temperature
            case "linear" -> "8.27,3.79";
            // Norwegian model uses Stien et al 2005
            // eggs_d = b[0] * (temperature + b[1])^2
            case "quadratic" -> "0.109,8.54";
            // Logistic model comes from data in Kragesteen 2023 Table 1
            // eggs_d = b[0]/(1 + exp(-b[1] * (temperature - b[2]))) + b[3]
            case "logistic" -> "65.5,0.550,8.49,10.51";
            default -> "28.2";
        };
        String eggTemp_b_Str = properties.getProperty("eggTemp_b", eggTemp_defaults);
        eggTemp_b = Arrays.stream(eggTemp_b_Str.split(",")).map(Double::parseDouble).toList();

        // Demographics
        mortSal_fn = properties.getProperty("mortSal_fn", "constant");
        String mortSal_defaults = switch (mortSal_fn) {
            // Stien et al 2005
            case "constant" -> "0.01";
            // estimated 2nd order polynomial fit to Bricknell et al. 2006 Table 1
            // mort_h = b[0] + b[1]*salinity + b[2]*salinity*salinity
            case "quadratic" -> "0.0011,1.1239,-0.07";
            // estimated logistic fit to Bricknell et al. 2006 Table 1
            // mort_h = b[0]/(1 + exp(-b[1] * (salinity - b[2]))) + b[3]
            case "logistic" -> "0.941,-0.594,13.5,0.0205";
            default -> "0.01";
        };
        String mortSal_b_Str = properties.getProperty("mortSal_b", mortSal_defaults);
        mortSal_b = Arrays.stream(mortSal_b_Str.split(",")).map(Double::parseDouble).toList();
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
        recordImmature = Boolean.parseBoolean(properties.getProperty("recordImmature", "false"));
        recordPsteps = Boolean.parseBoolean(properties.getProperty("recordPsteps", "true"));
        splitPsteps = Boolean.parseBoolean(properties.getProperty("splitPsteps", "true"));
        pstepsInterval = Integer.parseInt(properties.getProperty("pstepsInterval", "24"));
        pstepsMaxDepth = Double.parseDouble(properties.getProperty("pstepsMaxDepth", "10000"));
        recordConnectivity = Boolean.parseBoolean(properties.getProperty("recordConnectivity", "true"));
        connectImmature = Boolean.parseBoolean(properties.getProperty("connectImmature", "false"));
        connectivityInterval = Integer.parseInt(properties.getProperty("connectivityInterval", "24"));
        connectDepth1_max = Double.parseDouble(properties.getProperty("connectDepth1_max", "10000"));
        connectDepth1_min = Double.parseDouble(properties.getProperty("connectDepth1_min", "0"));
        connectDepth2_max = Double.parseDouble(properties.getProperty("connectDepth2_max", "10000"));
        connectDepth2_min = Double.parseDouble(properties.getProperty("connectDepth2_min", "10000"));
        recordLocations = Boolean.parseBoolean(properties.getProperty("recordLocations", "true"));
        recordArrivals = Boolean.parseBoolean(properties.getProperty("recordArrivals", "true"));
        recordMovement = Boolean.parseBoolean(properties.getProperty("recordMovement", "false"));
        recordActivity = Boolean.parseBoolean(properties.getProperty("recordElemActivity", "false"));
        recordVertDistr = Boolean.parseBoolean(properties.getProperty("recordVertDistr", "false"));
        vertDistrInterval = Integer.parseInt(properties.getProperty("vertDistrInterval", "1"));
        vertDistrMax = Integer.parseInt(properties.getProperty("vertDistrMax", "20"));

        // record connectivity between a second set of depths if set differently than defaults (i.e., realistic values)
        recordConnectivityDepth2 = connectDepth2_min < 10000 && connectDepth2_max < 10000;

        // hydrodynamic file requirements
        needS = (!fixDepth) || (!mortSal_fn.equals("constant"));
        needT = viableDegreeDays > -1;
        needZeta = false;
        needK = (!fixDepth) && variableDhV;
        needLight = (!fixDepth) && swimLightLevel;
        needVh = variableDh;
        needStokes = stokesDrift;

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
        if (!fixDepth && !meshType0.equals("FVCOM")) {
            conflict = true;
            conflictingProperties += "  verticalDynamics && !mesh1Type.equals(\"FVCOM\")\n";
        }
        if (needStokes && (hfDir2.isEmpty())) {
            conflict = true;
            conflictingProperties += "  Stokes drift on but no data directory\n";
        }

        if(conflict) {
            System.err.println(conflictingProperties);
            System.exit(1);
        }
    }

    public String printProperties() {
        return "\n------------ biotracker properties ------------\n" +
                "\n" +
                "-------- Run settings: \n" +
                "Start ymd: " + this.start_ymd + "\n" +
                "End ymd: " + this.end_ymd + "\n" +
                "Number of days: " + this.numberOfDays + "\n" +
                "Use EPSG:27700: " + this.coordOS + "\n" +
                "Run backwards: " + this.backwards + "\n" +
                "Number of cores: " + this.parallelThreads + "\n" +
                "Number of cores to read HD.nc: " + this.parallelThreadsHD + "\n" +
                "Duplicate last day of HD: " + this.duplicateLastDay + "\n" +
                "Check open boundaries: " + this.checkOpenBoundaries + "\n" +
                "Open boundary threshold (m): " + this.openBoundaryThresh + "\n" +
                "Particles per site per release: " + this.nparts + "\n" +
                "Release scenario: " + this.releaseScenario + "\n" +
                "Release time start: " + this.releaseTime + "\n" +
                "Release time end: " + this.releaseTimeEnd + "\n" +
                "Release interval (h): " + this.releaseInterval + "\n" +
                "Restart file: " + this.restartParticles + "\n" +
                "Restart particles cutoff days: " + this.restartParticlesCutoffDays + "\n" +
                "\n" +
                "-------- Meshes & Environment: \n" +
                "Records per HD file: " + this.recordsPerFile1 + "\n" +
                "dt (s per HD record): " + this.dt + "\n" +
                "Mesh 1: " + this.mesh0 + "\n" +
                "Mesh 1 type: " + this.meshType0 + "\n" +
                "HD 1 path: " + this.hfDir0 + this.hfDirPrefix0 + "YYYY" + this.hfDirSuffix0 + "/" + this.hfFilePrefix0 + "_YYYYMMDD.*nc" + "\n" +
                "Mesh 2: " + this.mesh1 + "\n" +
                "Mesh 2 type: " + this.meshType1 + "\n" +
                "HD 2 path: " + this.hfDir1 + this.hfDirPrefix1 + "YYYY" + this.hfDirSuffix1 + "/" + this.hfFilePrefix1 + "_YYYYMMDD.*nc" + "\n" +
                "Wave path: " + this.hfDir2 + this.hfDirPrefix2 + "YYYY" + this.hfDirSuffix2 + "/" + this.hfFilePrefix2 + "_YYYYMMDD.*nc" + "\n" +
                "Daily sunlight hours: " + this.daylightPath + "\n" +
                "Read salinity: " + this.needS + "\n" +
                "Read temperature: " + this.needT + "\n" +
                "Read shortwave: " + this.needLight + "\n" +
                "Read K: " + this.needK + "\n" +
                "Read viscovh: " + this.needVh + "\n" +
                "Read zeta: " + this.needZeta + "\n" +
                "Read wave data: " + this.needStokes + "\n" +
                "\n" +
                "-------- Sites: \n" +
                "Source site file: " + this.sitefile + "\n" +
                "Destination site file: " + this.sitefileEnd + "\n" +
                "Daily source lice counts: " + this.siteDensityPath + "\n" +
                "Site habitat (unused): " + this.habitat + "\n" +
                "Site suffix (unused): " + this.suffix + "\n" +
                "\n" +
                "-------- Dynamics: \n" +
                "RK4: " + this.rk4 + "\n" +
                "Interpolation step (s): " + (this.dt / this.stepsPerStep) + " (" + this.stepsPerStep + " / HD timestep)" + "\n" +
                "Diffusion: " + this.diffusion + "\n" +
                "Dh: " + (this.variableDh ? "variable" : this.D_h) + "\n" +
                "DhV: " + (this.variableDhV ? "variable" : this.D_hVert) + "\n" +
                "Stokes drift: " + this.stokesDrift + "\n" +
                "\n" +
                "-------- Biology: \n" +
                "End on arrival: " + this.endOnArrival + "\n" +
                "Fixed depth: " + this.fixDepth + "\n" +
                "Start depth (m): " + this.startDepth + "\n" +
                "Maximum depth (m): " + this.maxDepth + "\n" +
                "Viable time (h): " + this.viabletime + "\n" +
                "Max particle age (h): " + this.maxParticleAge + "\n" +
                "Viable degree days: " + this.viableDegreeDays + "\n" +
                "Max degree days: " + this.maxDegreeDays + "\n" +
                "Mortality rate function: " + this.mortSal_fn + "\n" +
                "Mortality rate parameters (/h): " + this.mortSal_b + "\n" +
                "Passive sinking intercept (m/s): " + this.passiveSinkingIntercept + "\n" +
                "Passive sinking slope (m/s/psu): " + this.passiveSinkingSlope + "\n" +
                "Copepodid P(swim down) = 0-1: " + this.salinityThreshCopepodidMax + " - " + salinityThreshCopepodidMin + "\n" +
                "Nauplius P(swim down) = 0-1: " + this.salinityThreshNaupliusMax + " - " + salinityThreshNaupliusMin + "\n" +
                "Copepodid swim down speed (m/s) ~ Norm(" + this.swimDownSpeedCopepodidMean + ", " + this.swimDownSpeedCopepodidStd + ")" + "\n" +
                "Nauplius swim down speed (m/s) ~ Norm(" + this.swimDownSpeedNaupliusMean + ", " + this.swimDownSpeedNaupliusStd + ")" + "\n" +
                "Swim up based on shortwave: " + this.swimLightLevel + "\n" +
                "Copepodid light threshold (umol/m2/s): " + this.lightThreshCopepodid + "\n" +
                "Nauplius light threshold (umol/m2/s): " + this.lightThreshNauplius + "\n" +
                "Copepodid swim up speed (m/s) ~ Norm(" + this.swimUpSpeedCopepodidMean + ", " + this.swimUpSpeedCopepodidStd + ")" + "\n" +
                "Nauplius swim up speed (m/s) ~ Norm(" + this.swimUpSpeedNaupliusMean + ", " + this.swimUpSpeedNaupliusStd + ")" + "\n" +
                "Egg production function: " + this.eggTemp_fn + "\n" +
                "Egg production parameters (/d): " + this.eggTemp_b + "\n" +
                "\n" +
                "-------- Recording: \n" +
                "Psteps (lice/element/time): \n" +
                "- Copepodid: " + this.recordPsteps + "\n" +
                "- Nauplius: " + this.recordImmature + "\n" +
                "- Split by source: " + this.splitPsteps + "\n" +
                "- Recording interval (h): " + this.pstepsInterval + "\n" +
                "- Max depth (m): " + this.pstepsMaxDepth + "\n" +
                "Connectivity:\n" +
                "- Copepodid: " + this.recordConnectivity + "\n" +
                "- Nauplius: " + this.connectImmature + "\n" +
                "- Threshold (m): " + this.connectivityThresh + "\n" +
                "- Recording interval (h): " + this.connectivityInterval + "\n" +
                "- Depth min 1 (m): " + this.connectDepth1_min + "\n" +
                "- Depth max 1 (m): " + this.connectDepth1_max + "\n" +
                (this.recordConnectivityDepth2 ? (
                        "- Depth min 2 (m): " + this.connectDepth2_min + "\n" +
                        "- Depth max 2 (m): " + this.connectDepth2_max + "\n"
                        ) : "") +
                "Vertical distributions:\n" +
                "- Copepodid: " + this.recordVertDistr + "\n" +
                "- Nauplius: " + this.recordImmature + "\n" +
                "- Recording interval (h): " + this.vertDistrInterval + "\n" +
                "- Max depth (m): " + this.vertDistrMax + "\n" +
                "Hourly particle locations: " + this.recordLocations + "\n" +
                "Particle arrivals: " + this.recordArrivals + "\n" +
                "Sub-hourly movement: " + this.recordMovement + "\n" +
                "Total sink-swim-passive by element: " + this.recordActivity + "\n" +
                "\n" +
                "-------------- end of properties --------------\n\n";
    }
}



