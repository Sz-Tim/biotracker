/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.io.*;
import java.util.Properties;
import static particle_track.Particle_track.dateIntParse;

/**
 *
 * @author sa01ta
 */
public class RunProperties {
    String basedir, sitedir, // Location of code and site data files
            datadir, datadirPrefix, datadirSuffix, // Location of default hydrodynamic data, with prefix and suffix for finding annual subdirectories
            datadir2, datadir2Prefix, datadir2Suffix, // Location of secondary (larger domain) hydrodynamic data, with prefix and suffix for finding annual subdirectories 
            mesh1, mesh2, // Full path to the mesh files used describing spatial structure of the hydrodynamic data (
            location, sitefile, habitat, suffix; // Descriptive strings
            
    boolean backwards, // run model backwards? Needs some work on loops to make this work correctly
            timeInterpolate, spatialInterpolate, // interpolate between hydro file data values?
            rk4, // use RK4 numerical integration (alternative is Euler; need about 10 times as many steps)
            cluster, // run on cluster? - POSSIBLY DEPRECATE DEPENDING ON USAGE
            parallel, // use multiple cores to speed up run?
            diffusion, variableDiff, // include random walk, use diffusion parameter from hydro output?
            salinityMort, // mortality calculated based on local salinity (sea lice - doesn't presently do anything)?
            endOnArrival, // stop at first suitable habitat site, or simply note arrival and move on?
            setDepth, // set particle depth at initiation?
            splitPsteps, // separate pSteps by source site? No longer produced (done in Matlab post processing) - POSSIBLY DEPRECATE
            pstepsIncMort, // Include mortality in pstep calculation (negative exponential, unless "calcMort" is true) - POSSIBLY DEPRECATE
            all_locs_out;
    
    int start_ymd, end_ymd, numberOfDays, // Start and end of run. If numberOfDays = 0, it is ignored and end_ymd is used instead
            releaseScenario, // 0 release all at "releaseTime", 1 continuous release ("nparts" per hour per site)
            nparts, // Number of particles released per site (per hour in releaseScenario == 1
            recordsPerFile, // Number of records per velocity file
            stepsPerStep, // Number of increments between each velocity record (also for time interpolations)
            depthLayers, // Number of depth layers in hydro output - IDEALLY DEPRECATE
            dumpInterval, // Interval in hours for printing particle locations and elements to file
            thresh, // Threshold distance for "settlement" (m)
            behaviour, // Particle behaviour - see Particle.java
            endlimit, // Maximum ID of startlocs to use as a destination (0 = use all)
            N, // Number of mesh element centroids - IDEALLY DEPRECATE
            M; // Number of mesh nodes - IDEALLY DEPRECATE

    double releaseTime, viabletime, // Time of particle release (if releaseScenario == "0"), Time to attain settlement competency
            dt, // Timestep (s) per record
            D_h, // Horizontal diffuision parameter
            diffusionMultiplier, // Additional diffusion transport distance linear multiplier
            mortalityRate, // Hourly mortality rate of particles
            maxParticleAge; // Maximum age for particles. Set to <=0 and it will be ignored.
  
    /**
     * 
     * @param filename 
     */
    public RunProperties(String filename)
    {
        System.out.println("GETTING PROPERTIES FROM "+filename);
        Properties properties = new Properties();
        try 
        {
            properties.load(new FileInputStream(filename));
        } 
        catch (IOException e) 
        {
            System.err.println("!!!Could not find properties file!!!");
        }
        
        for(String key : properties.stringPropertyNames()) {
            String value = properties.getProperty(key);
            System.out.println(key + " => " + value);
        }
        
        basedir = properties.getProperty("basedir","/home/sa01ta/particle_track/");
        sitedir = properties.getProperty("sitedir","/home/sa01ta/particle_track/minch_sites/");
        datadir = properties.getProperty("datadir","/home/sa01da/work/minch2/Archive/");
        datadirPrefix = properties.getProperty("datadirPrefix","netcdf_");
        datadirSuffix = properties.getProperty("datadirSuffix","");
        datadir2 = properties.getProperty("datadir2","/home/sa01da/data/oss/");
        datadir2Prefix = properties.getProperty("datadirPrefix","");
        datadir2Suffix = properties.getProperty("datadirSuffix","_OLD");
        
        mesh1 = properties.getProperty("mesh1","/home/sa01ta/particle_track/WestCOMS_mesh.nc");
        mesh2 = properties.getProperty("mesh2","");
    
        sitefile = properties.getProperty("sitefile","startlocations.dat");
        location = properties.getProperty("location","minch_continuous");
        habitat = properties.getProperty("habitat","");
        suffix = properties.getProperty("suffix","");     
        
        start_ymd = Integer.parseInt(properties.getProperty("start_ymd","20180101"));
        numberOfDays = Integer.parseInt(properties.getProperty("numberOfDays","0"));
        if (numberOfDays>0)
        {
            int[] startDate = dateIntParse(start_ymd);
            ISO_datestr tempIsoDate = new ISO_datestr(startDate[0], startDate[1], startDate[2]);
            for (int i = 1; i < numberOfDays; i++)
            {
                tempIsoDate.addDay();
            }
            end_ymd = Integer.parseInt(tempIsoDate.getDateStr());
        }
        else
        {
            end_ymd = Integer.parseInt(properties.getProperty("end_ymd","20180102"));
        }
        if (end_ymd<start_ymd)
        {
            System.err.println("***** End date before start date! *****");
            System.exit(1);
        }
        
        // Set variable values based on contents of the properties object
        backwards = Boolean.parseBoolean(properties.getProperty("backwards","false"));
        timeInterpolate = Boolean.parseBoolean(properties.getProperty("timeInterpolate","true"));
        spatialInterpolate = Boolean.parseBoolean(properties.getProperty("spatialInterpolate","true"));
    
        rk4 = Boolean.parseBoolean(properties.getProperty("rk4","true"));
        cluster = Boolean.parseBoolean(properties.getProperty("cluster","true"));
        
        diffusion = Boolean.parseBoolean(properties.getProperty("diffusion","true"));
        variableDiff = Boolean.parseBoolean(properties.getProperty("variableDiff","false"));
        salinityMort = Boolean.parseBoolean(properties.getProperty("salinityMort","false"));
        endOnArrival = Boolean.parseBoolean(properties.getProperty("endOnArrival","false"));
        // DEPRECATED
        //tidalRelease = Boolean.parseBoolean(properties.getProperty("tidalRelease","false"));
        setDepth = Boolean.parseBoolean(properties.getProperty("setDepth","false"));
        splitPsteps = Boolean.parseBoolean(properties.getProperty("splitPsteps","false"));
        
        pstepsIncMort = Boolean.parseBoolean(properties.getProperty("pstepsIncMort","true"));
        
        parallel = Boolean.parseBoolean(properties.getProperty("parallel","true"));
        // DEPRECATED
        //oldOutput = Boolean.parseBoolean(properties.getProperty("oldOutput"));
        
        all_locs_out = Boolean.parseBoolean(properties.getProperty("all_locs_out","true"));
    
        releaseScenario = Integer.parseInt(properties.getProperty("releaseScenario","0"));
        nparts = Integer.parseInt(properties.getProperty("nparts","5"));
        releaseTime = Double.parseDouble(properties.getProperty("releaseTime","0"));

        dt = Double.parseDouble(properties.getProperty("dt","3600"));
        
        // 21/11/2018 --- Ideally want to remove this. If using two models for hydro 
        // would ideally sense from the hydro file
        recordsPerFile = Integer.parseInt(properties.getProperty("recordsPerFile","7"));
        
        stepsPerStep = Integer.parseInt(properties.getProperty("stepsPerStep","25"));
        
        // 21/11/2018 --- possible to remove? If using two models for hydro 
        // would ideally sense from the hydro file
        depthLayers = Integer.parseInt(properties.getProperty("depthLayers","10"));
        
        dumpInterval = Integer.parseInt(properties.getProperty("dumpInterval","24"));

        thresh = Integer.parseInt(properties.getProperty("thresh","500"));
        viabletime = Double.parseDouble(properties.getProperty("viabletime","86"));
        behaviour = Integer.parseInt(properties.getProperty("behaviour","1"));
        
        D_h = Double.parseDouble(properties.getProperty("D_h","0.1"));
        diffusionMultiplier = Double.parseDouble(properties.getProperty("diffusionMultiplier","1"));
        mortalityRate = Double.parseDouble(properties.getProperty("mortalityRate","0.01"));
        maxParticleAge = Double.parseDouble(properties.getProperty("maxParticleAge","-1"));
        
        // 21/11/2018 --- possible to remove? If using two models for hydro 
        // would ideally sense from the hydro file
        N = Integer.parseInt(properties.getProperty("N","79244"));
        M = Integer.parseInt(properties.getProperty("M","46878"));
        
        endlimit = Integer.parseInt(properties.getProperty("endlimit","0"));
        
        properties.list(System.out);
        
    }

//    // default constructor for lack of property file
//    // NOT TESTED SINCE CODE REVISED.......
//    public RunProperties(String[] args)
//    {
//        if (args.length == 6 || args.length == 7 || args.length == 8)
//        {
//            System.out.println("Standard minch chapter args (5)");
//            // set values for many of the things that used to be read in
//            
//            habitat=args[0];
//            //firstday=Integer.parseInt(args[1]);
//            //lastday=Integer.parseInt(args[2]);
//            diffusionMultiplier=Double.parseDouble(args[3]);
//            D_h=0.1*diffusionMultiplier;
//            behaviour=Integer.parseInt(args[4]);
//            nparts=Integer.parseInt(args[5]);
//            if (args.length == 7)
//            {
//                stepsPerStep=Integer.parseInt(args[6]);
//            }
//            if (args.length == 8)
//            {
//                viabletime=Double.parseDouble(args[7]);
//            }
//        }
//        
//        if (args.length >= 12)
//        {
//            nparts=Integer.parseInt(args[0]);
//            //firstday=Integer.parseInt(args[1]);
//            //lastday=Integer.parseInt(args[2]);
//            recordsPerFile=Integer.parseInt(args[3]);
//            dt=Integer.parseInt(args[4]);
//            stepsPerStep=Integer.parseInt(args[5]);
//            thresh=Integer.parseInt(args[6]);
//            viabletime=Integer.parseInt(args[7]);
//            behaviour=Integer.parseInt(args[8]);
//            if (args[9].equalsIgnoreCase("true"))
//            {
//                timeInterpolate = true;
//            }
//            if (args[10].equalsIgnoreCase("true"))
//            {
//                spatialInterpolate = true;
//            }
//            D_h = Double.parseDouble(args[11]);
//        }
//        
//        if (args.length == 14)
//        {
//            suffix = args[12];
//            location = args[13];
//        }
//        
//        if (D_h==0)
//            {
//                diffusion = false;
//            }
//        // the case of using spatially variable diffusion as output by FVCOM
//        if (D_h<0)
//            {
//                variableDiff = true;
//            }
//    }
    
    
}



