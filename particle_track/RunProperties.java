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
    String basedir = "/home/sa01ta/particle_track/";
    String sitedir = "/home/sa01ta/particle_track/minch_sites";
    String datadir = "/home/sa01ta/particle_track/minch_continuous_overlap/";
    String datadir2 = "/home/sa01ta/particle_track/minch_continuous_overlap/";
    String location = "minch_continuous";
    String sitefile = "startlocations.dat";
    String habitat = "";
    String suffix = "";
    
    int start_ymd = 20150101;
    int end_ymd = 20150121;
    int numberOfDays = 0;
        
    boolean backwards = false;
    boolean timeInterpolate = true;
    boolean spatialInterpolate = true;
    
    boolean rk4 = true;
    boolean cluster = true;
    
    boolean diffusion = true;
    boolean variableDiff = false;
    // Particle properties ---------------------
    boolean calcMort = false;
    boolean endOnArrival = true; 
    boolean tidalRelease = true;
    boolean setDepth = false;
    
    boolean splitPsteps = true;
    // Include mortality in pstep calculation (negative exponential, unless "calcMort" is true)
    boolean pstepsIncMort = true;
    
    boolean parallel = false;
    boolean all_locs_out = true;
    
    int releaseScenario = 5; // Instantaneous/continuous/etc release - see switch statement in main method
    int nparts=100; // Number of particles released per site

    double dt=3600; // Timestep (s) per record
    //int firstday=1; // Start day for simulation
    //int lastday=21; // End day for simulation
    int recordsPerFile=7; // Number of records per velocity file
    int stepsPerStep=25; // Number of increments between each velocity record (also for time interpolations)
    int depthLayers=10;
    int dumpInterval=24; // Interval in hours for printing particle locations and elements to file 
    
    // The threshold distance, closer than which particles are deemed to have settled.
    int thresh = 500; // Threshold distance for "settlement" (m)
    double viabletime = 86; // Time to attain settlement competency
   
    int behaviour=1; // Particle behaviour - see Particle.java

    double D_h = 0.1;
        
    double diffusionMultiplier = 1;
    double mortalityRate = 0.01;
    double maxParticleAge = -1; // Maximum age for particles. Set to <=0 and it will be ignored.
    
    int N = 79244; // Number of mesh element centroids
    int M = 46878; // Number of mesh nodes
    
    int endlimit = 0;
    
    // default constructor for lack of property file
    // NOT TESTED SINCE CODE REVISED.......
    public RunProperties(String[] args)
    {
        if (args.length == 6 || args.length == 7 || args.length == 8)
        {
            System.out.println("Standard minch chapter args (5)");
            // set values for many of the things that used to be read in
            
            habitat=args[0];
            //firstday=Integer.parseInt(args[1]);
            //lastday=Integer.parseInt(args[2]);
            diffusionMultiplier=Double.parseDouble(args[3]);
            D_h=0.1*diffusionMultiplier;
            behaviour=Integer.parseInt(args[4]);
            nparts=Integer.parseInt(args[5]);
            if (args.length == 7)
            {
                stepsPerStep=Integer.parseInt(args[6]);
            }
            if (args.length == 8)
            {
                viabletime=Double.parseDouble(args[7]);
            }
        }
        
        if (args.length >= 12)
        {
            nparts=Integer.parseInt(args[0]);
            //firstday=Integer.parseInt(args[1]);
            //lastday=Integer.parseInt(args[2]);
            recordsPerFile=Integer.parseInt(args[3]);
            dt=Integer.parseInt(args[4]);
            stepsPerStep=Integer.parseInt(args[5]);
            thresh=Integer.parseInt(args[6]);
            viabletime=Integer.parseInt(args[7]);
            behaviour=Integer.parseInt(args[8]);
            if (args[9].equalsIgnoreCase("true"))
            {
                timeInterpolate = true;
            }
            if (args[10].equalsIgnoreCase("true"))
            {
                spatialInterpolate = true;
            }
            D_h = Double.parseDouble(args[11]);
        }
        
        if (args.length == 14)
        {
            suffix = args[12];
            location = args[13];
        }
        
        if (D_h==0)
            {
                diffusion = false;
            }
        // the case of using spatially variable diffusion as output by FVCOM
        if (D_h<0)
            {
                variableDiff = true;
            }
    }
    
    /**
     * Ultimately want to fix this so that if a property is not included in the .properties
     * file, the default value is uses automatically and a warning raised
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
        
        basedir = properties.getProperty("basedir");
        sitedir = properties.getProperty("sitedir");
        datadir = properties.getProperty("datadir");
        datadir2 = properties.getProperty("datadir2");
    
        sitefile = properties.getProperty("sitefile");
        location = properties.getProperty("location");
        habitat = properties.getProperty("habitat");
        suffix = properties.getProperty("suffix");     
        
        start_ymd = Integer.parseInt(properties.getProperty("start_ymd"));
        numberOfDays = Integer.parseInt(properties.getProperty("numberOfDays"));
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
            end_ymd = Integer.parseInt(properties.getProperty("end_ymd"));
        }
        if (end_ymd<start_ymd)
        {
            System.err.println("End date before start date!");
            System.exit(1);
        }
        
        // Set variable values based on contents of the properties object
        backwards = Boolean.parseBoolean(properties.getProperty("backwards"));
        timeInterpolate = Boolean.parseBoolean(properties.getProperty("timeInterpolate"));
        spatialInterpolate = Boolean.parseBoolean(properties.getProperty("spatialInterpolate"));
    
        rk4 = Boolean.parseBoolean(properties.getProperty("rk4"));
        cluster = Boolean.parseBoolean(properties.getProperty("cluster"));
        
        diffusion = Boolean.parseBoolean(properties.getProperty("diffusion"));
        variableDiff = Boolean.parseBoolean(properties.getProperty("variableDiff"));
        calcMort = Boolean.parseBoolean(properties.getProperty("calcMort"));
        endOnArrival = Boolean.parseBoolean(properties.getProperty("endOnArrival"));
        tidalRelease = Boolean.parseBoolean(properties.getProperty("tidalRelease"));
        setDepth = Boolean.parseBoolean(properties.getProperty("setDepth"));
        splitPsteps = Boolean.parseBoolean(properties.getProperty("splitPsteps"));
        
        pstepsIncMort = Boolean.parseBoolean(properties.getProperty("pstepsIncMort"));
        
        parallel = Boolean.parseBoolean(properties.getProperty("parallel"));
        
        all_locs_out = Boolean.parseBoolean(properties.getProperty("all_locs_out"));
    
        releaseScenario = Integer.parseInt(properties.getProperty("releaseScenario"));
        nparts = Integer.parseInt(properties.getProperty("nparts"));

        dt = Double.parseDouble(properties.getProperty("dt"));

        //firstday = Integer.parseInt(properties.getProperty("firstday"));
        //lastday = Integer.parseInt(properties.getProperty("lastday"));
        recordsPerFile = Integer.parseInt(properties.getProperty("recordsPerFile"));
        stepsPerStep = Integer.parseInt(properties.getProperty("stepsPerStep"));
        depthLayers = Integer.parseInt(properties.getProperty("depthLayers"));
        dumpInterval = Integer.parseInt(properties.getProperty("dumpInterval"));

        thresh = Integer.parseInt(properties.getProperty("thresh"));
        viabletime = Double.parseDouble(properties.getProperty("viabletime"));
        behaviour = Integer.parseInt(properties.getProperty("behaviour"));
        
        D_h = Double.parseDouble(properties.getProperty("D_h"));
        diffusionMultiplier = Double.parseDouble(properties.getProperty("diffusionMultiplier"));
        mortalityRate = Double.parseDouble(properties.getProperty("mortalityRate"));
        maxParticleAge = Double.parseDouble(properties.getProperty("maxParticleAge"));
        
        N = Integer.parseInt(properties.getProperty("N"));
        M = Integer.parseInt(properties.getProperty("M"));
        
        endlimit = Integer.parseInt(properties.getProperty("endlimit"));
    }
        
}



