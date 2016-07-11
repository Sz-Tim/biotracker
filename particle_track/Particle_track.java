/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.io.*;
import java.util.*;
import edu.cornell.lassp.houle.RngPack.*;

/**
 *
 * @author tomdude
 */
public class Particle_track {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        // TODO code application logic here
    
        System.out.println("Starting particle tracking program\n");
        Date date = new Date();
        // display time and date using toString()
        System.out.println(date.toString());
       
        //System.out.println(new Date().toString());
        long startTime = System.currentTimeMillis();
        
        RanMT ran = new RanMT(System.currentTimeMillis());
           
        System.out.println("Reading in data\n");
        
        /**
         * files are 
         * 
         * nodexy.dat, uvnode.dat, trinodes.dat, neighbours.dat
         * u_[3-31].dat
         * v_[3-31].dat
         */ 
        
        boolean backwards = false;
        boolean diffusion = true;
        boolean variableDiff = false;
        boolean calcMort = false;
        
        boolean cluster = true;
        
        // This actually isn't used, except in the case of defaulting on location
        String basedir = "C:\\Users\\sa01ta\\Documents\\particle_track\\";
        String scriptdir =  "C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\hydro_mesh_run\\minch2\\";
        if (cluster==true)
        {
            basedir = "/home/sa01ta/particle_track/";
        }
        
        // Set some basic default parameters
        int nparts=940;
        nparts=100;
        double dt=3600;
        int firstday=3;
        int lastday=31;
        int recordsPerFile=24;
        int stepsPerStep=200;
        
        // The threshold distance, closer than which particles are deemed to have settled.
        int thresh = 500;
        
        // Particles become viable when they are halway through their PLD. 
        // This is compared to individual particle age later.
        double viabletime = 0;
        String location = "lorn";
               
        int behaviour=1; // Particle behaviour - see Particle.java
        
        boolean timeInterpolate = true;
        boolean spatialInterpolate = true;
        boolean rk4 = true;
        if (rk4 == true)
        {
            System.out.println("Using RK4 algorithm");
        }
        else
        {
            System.out.println("Using Euler algorithm");
        }
        double D_h = 0.1;
        
        /**
         * SETUP CONFIGURATION FOR THE MINCH CHAPTER CASES
         */
        String habitat = "";
        String suffix = "";
        double diffusionMultiplier = 1;
        if (args.length == 6 || args.length == 7 || args.length == 8)
        {
            System.out.println("Standard minch chapter args (5)");
            // set values for many of the things that used to be read in
            nparts=Integer.parseInt(args[5]);
            recordsPerFile=6;
            dt=3600;
            thresh=500;
            
            location="minch_continuous";
            timeInterpolate = true;
            spatialInterpolate = true;
            
            habitat=args[0];
            firstday=Integer.parseInt(args[1]);
            lastday=Integer.parseInt(args[2]);
            diffusionMultiplier=Double.parseDouble(args[3]);
            D_h=0.1*diffusionMultiplier;
            if (D_h==0)
                {
                    diffusion = false;
                }
            // the case of using spatially variable diffusion as output by FVCOM
            if (D_h<0)
                {
                    variableDiff = true;
                }
            behaviour=Integer.parseInt(args[4]);
            stepsPerStep=25;
            if (args.length == 7)
            {
                stepsPerStep=Integer.parseInt(args[6]);
            }
            viabletime=9999;
            if (args.length == 8)
            {
                viabletime=Double.parseDouble(args[7]);
            }
        }
        
        if (args.length >= 12)
        {
            try 
            {
                nparts=Integer.parseInt(args[0]);
                firstday=Integer.parseInt(args[1]);
                lastday=Integer.parseInt(args[2]);
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
                if (D_h==0)
                {
                    diffusion = false;
                    System.out.println("No diffusion");
                }
                //System.out.println("nparts = "+nparts+"\nsettlement threshold distance = "+thresh+"\nviabletime = "+viabletime+" D_h = "+D_h);
            }
            catch (IndexOutOfBoundsException e)
            {
                System.err.println("Incorrect number of input parameters, found "+args.length);
            }
        }
        
        // Print all main arguments
        System.out.printf("-----------------------------------------------------------\n");
        System.out.printf("Location           = %s\n",location);
        System.out.printf("Habitat            = %s\n",habitat);
        System.out.printf("N_parts/site       = %d\n",nparts);
        System.out.printf("hydromod dt (s)    = %f\n",dt);
        System.out.printf("hydromod rec/file  = %d\n",recordsPerFile);
        System.out.printf("stepsperstep       = %d\n",stepsPerStep);
        System.out.printf("firstfile          = %d\n",firstday);
        System.out.printf("lastfile           = %d\n",lastday);
        System.out.printf("RK4                = %s\n",rk4);
        System.out.printf("Vertical behaviour = %d\n",behaviour);        
        System.out.printf("Viable time (h)    = %f\n",viabletime);
        System.out.printf("Viable time (d)    = %f\n",viabletime/24.0);
        System.out.printf("Threshold distance = %d\n",thresh);
        System.out.printf("Diffusion          = %f\n",D_h);
        System.out.printf("-----------------------------------------------------------\n");
        
        if ((location.equalsIgnoreCase("lorn") == true || location.equalsIgnoreCase("minch") == true) && viabletime==0)
        {
            // case when for viable at
            viabletime=(lastday-firstday+1)*recordsPerFile*dt/2;
            System.out.println("Ignoring input viable time, set to "+viabletime);
        }
        else if (viabletime==666)
        {
            // ensure that particles never become viable...
            viabletime=(lastday-firstday+1)*recordsPerFile*dt*2;
            System.out.println("Ignoring input viable time, set to "+viabletime);
        }
        else
        {
            viabletime=viabletime*3600;
            System.out.println("Using input viabletime: "+viabletime/3600.0+"hrs "+viabletime+"s");
        }
        System.out.printf("-----------------------------------------------------------\n");
        
        String datadir = "D:\\FVCOM_2011_jun_9d";
        // if only 12 arguments, just behave as previously (Firth of Lorn assumed, or set manually to Fyne)
        if (args.length == 14)
        {
            //datadir = args[13];
            suffix = args[12];
            location = args[13];
        }
        // default "Firth of Lorn" values
        int N = 25071;
        int M = 14000;
        if (location.equalsIgnoreCase("fyne"))
        {          
            N = 1593;
            M = 943;
        } 
        else if (location.equalsIgnoreCase("minch") || location.equalsIgnoreCase("minch_continuous"))
        {
            N = 79244;
            M = 46878;            
        }
        
        double[][] nodexy = new double[M][2];
        double[][] uvnode = new double[N][2];
        int[][] trinodes = new int[N][3];
        int[][] neighbours = new int[N][3];

        
        
        // --------------------------------------------------------------------------------------
        // File reading and domain configuration
        // --------------------------------------------------------------------------------------
        
        String datadir2 = "";
        if (location.equalsIgnoreCase("fyne"))
        {
            datadir = "D:\\fvcom_output\\"+datadir+"\\output\\ptrack_dat";
            datadir2 = "C:\\Users\\sa01ta\\Documents\\EFF\\mesh\\943dat2";
            System.out.println("DATA DIRECTORY = "+datadir);

        } 
        else if (location.equalsIgnoreCase("minch")) 
        {
            datadir = "D:\\dima_data\\minch1_results\\plots\\minch1_v05\\";
            datadir2 = "D:\\dima_data\\minch1_results\\plots\\minch1_v05\\";
            //datadir2 = "C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\hydro_mesh_run";
            if (cluster==true)
            {
                datadir=basedir+"minch1_v05/";
                datadir2=basedir+"minch1_v05/";
            }
            System.out.println("DATA DIRECTORY = "+datadir2);
        }
        else if (location.equalsIgnoreCase("minch_continuous"))
        {
            datadir = "D:\\dima_data\\minch_continuous\\";
            datadir2 = "C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\hydro_mesh_run\\minch2\\";
            if (cluster==true)
            {
                datadir=basedir+"minch_continuous/";
                datadir2=basedir+"minch_continuous/";
            }
            System.out.println("DATA DIRECTORY = "+datadir2);
        }
        else
        {
            System.out.println("No location......");
            datadir = "D:\\dima_data\\"+datadir;
            System.out.println("DATA DIRECTORY = "+datadir);
        }
        nodexy = IOUtils.readFileDoubleArray(datadir2+"nodexy.dat",M,2," ",true); // the original mesh nodes
        uvnode = IOUtils.readFileDoubleArray(datadir2+"uvnode.dat",N,2," ",true); // the centroids of the elements
        trinodes = IOUtils.readFileIntArray(datadir2+"trinodes.dat",N,3," ",true); // the corners of the elements
        neighbours = IOUtils.readFileIntArray(datadir2+"neighbours.dat",N,3," ",true);
        
        // reduce node/element IDs in files generated by matlab by one (loops start at zero, not one as in matlab)
        for (int i = 0; i < N; i++)
        {
            //System.out.println(bathymetry[i][0]);
            for (int j = 0; j < 3; j++)
            {
                trinodes[i][j]--;
                //System.out.printf("%d ",trinodes[i][j]);
                if (neighbours[i][j]>0)
                {
                    neighbours[i][j]--;
                }             
            }
            //System.out.printf("\n");         
            
        }
        int[] allelems = new int[uvnode.length];
        for (int j = 0; j < uvnode.length; j++)
        {
            allelems[j] = j;
        }

        
        
        int nrecords=(lastday-firstday+1)*recordsPerFile;
        int inviabletime=nrecords;
        
        dt=dt/(double)stepsPerStep;
        double dev_perstep = Math.pow(0.1,dt);
        System.out.println("Particle dt = "+dt+" dev_perstep = "+dev_perstep);

        System.out.println("behaviour = "+behaviour);
//        if (timeInterpolate==true)
//        {
//            System.out.println("WARNING: time interpolation not operational in last record of each file");
//        }

        // load array of start node IDs (as stored by matlab)
        double startlocs[][] = new double[10][3];
        double open_BC_locs[][] = new double[1][3];
        
        //String sitedir="C:\\Users\\sa01ta\\Documents\\westcoast_tracking_summary\\habitat_sites\\sites_to_use\\";
        String sitedir="C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\site_locations\\";
        if (cluster==true)
        {
            sitedir=basedir+"minch_sites/";
        }

        startlocs = IOUtils.setupStartLocs(location,habitat,basedir);
        open_BC_locs = IOUtils.setupOpenBCLocs(location,habitat,basedir); 
        
        if (habitat.equalsIgnoreCase("spill") || habitat.equalsIgnoreCase("test")){
            viabletime=(lastday+1)*recordsPerFile*dt*3600;
            System.out.println("ignoring input viable time, set to "+viabletime); 
        }
        
        // --------------------------------------------------------------------------------------
        // Creating initial particle array
        // --------------------------------------------------------------------------------------        
        int nparts_per_site=nparts;
        nparts=nparts*startlocs.length; 
        
        for (int i = 0; i < startlocs.length; i++)
        {
            startlocs[i][0]--;
            //System.out.println(startlocs[i][0]);
        }

        //nparts=startlocs.length;
        //int tot_psteps=nparts*nrecords;

        // an array to save the number of "particle-timesteps" in each cell
        double[][] pstepsImmature= new double[N][2];
        double[][] pstepsMature= new double[N][2];
        
   
        for (int i = 0; i < N; i++)
        {
            pstepsImmature[i][0] = i;
            pstepsMature[i][0] = i;
        }
        // array to save source, destination, and transport time for each particle
        int[][] particle_info= new int[nparts][3];
        double[][] settle_density = new double[nparts][1];
        
        double[] xstart = new double[nparts];
        double[] ystart = new double[nparts];
        int[] startElem = new int[nparts];
        int[] startid = new int[nparts];
        System.out.println("nparts = "+nparts);
        //System.out.println("hello");
        for (int i = 0; i < nparts; i++)
        {
            //startid[i]=(int)(startlocs.length*ran.raw());
            startid[i]=i%startlocs.length;
            xstart[i]=startlocs[startid[i]][1];
            ystart[i]=startlocs[startid[i]][2];

            // this boundary location is not actually in the mesh/an element, so set
            // new particle location to centre of nearest element.
            int closest=Particle.nearestCentroid(xstart[i],ystart[i],uvnode);
            startElem[i]=Particle.whichElement(xstart[i],ystart[i],allelems,nodexy,trinodes);
            if (startElem[i]<0)
            {
                xstart[i]=uvnode[closest][0];
                ystart[i]=uvnode[closest][1];
                startElem[i]=closest;
            }
            //System.out.printf("start location %d = %d %.4e %.4e %d %d\n",i,(int)startlocs[startid[i]][0],xstart[i],ystart[i],closest,startElem[i]);
        }

        // setup particles
        Particle[] particles = new Particle[nparts];
        System.out.println("particles.length = "+particles.length);
        for (int i = 0; i < particles.length; i++)
        {
            particles[i] = new Particle(xstart[i],ystart[i],i);
            particle_info[i][0]=startid[i];//(int)startlocs[startid[i]][0];
            particles[i].setElem(startElem[i]);
        }
        
        IOUtils.particleLocsToFile(particles,nparts,0,"particlelocations_start.out");

        // --------------------------------------------------------------------------------------
        // Final setup bits
        // --------------------------------------------------------------------------------------
        System.out.println("Starting time loop");

        int[] searchCounts= new int[5];

        double minDistTrav=10000000;
        double maxDistTrav=0;

        int stepcount=0;
        int calcCount=0;
        double time=0;
        
        double[][] particle1Velocity = new double[nrecords*stepsPerStep][3];
        double[][] particle1Location = new double[nrecords*stepsPerStep][2];
                
        int printCount=0;
        
        // Initial value for br set - this particular one is only used in case of "spill" habitat (identical to that file list)
        BufferedReader br = new BufferedReader(new FileReader(datadir2+"150408_spill_list1.dat"));
        if (location.equalsIgnoreCase("minch_continuous")) {
            if (habitat.equalsIgnoreCase("spill")){
                br = new BufferedReader(new FileReader("150408_spill_list1.dat"));
            }
            // Doing a minch_continuous run but using whatever set of files is appropriate for the desired period
            else
            {
                br = new BufferedReader(new FileReader("filelist.dat"));
            }
        }
        // Marking and reseting fails for some reason; use bash head command instead in submitting script
        //br.mark(0);
        //System.out.println("filelist top line: "+br.readLine());
        //br.reset();

        //int nfreeparts = 0;
        int nfreeparts = nparts;
        int nViable = 0;
        int nBoundary = 0;
        int nSettled = 0;
        
        String filenums = "";
        
        
        
        
        // --------------------------------------------------------------------------------------
        // Start time loop
        // --------------------------------------------------------------------------------------
        for (int fnum = firstday; fnum <= lastday; fnum++)
        // alternatively, run loop backwards
        //for (int day = lastday; day >= firstday; day--)
        {
            System.out.printf("\nfile %d - time %fsecs (%fhrs) \n",fnum,stepcount*dt*stepsPerStep,time);
            
            // clear any old data
            //clear FVCOM1
            // load the new data file. this puts variables straight into the
            // workspace
            int depthLayers = 1;       
            
            String ufile = "";
            String vfile = "";
            String elfile = "";
            String ufile1 = "";
            String vfile1 = "";
            String elfile1 = "";
            double[][] u = new double[recordsPerFile][N*depthLayers];
            double[][] v = new double[recordsPerFile][N*depthLayers];
            //double[][] el = new double[recordsPerFile][N*depthLayers];
            double[][] u1 = new double[recordsPerFile][N*depthLayers];
            double[][] v1 = new double[recordsPerFile][N*depthLayers];
            //double[][] el1 = new double[recordsPerFile][N*depthLayers];
            
            // Reading files on the first instance:
            // - Read two files to provide data for time interpolation during last time steps
            if (fnum == firstday)
            {
                System.out.println("Reading data day "+firstday);
                // maintain backward compatability with previous number reading
                filenums = ""+fnum;
                // replace filenums in its entirety               
                if (location.equalsIgnoreCase("minch_continuous"))
                {
                    filenums = br.readLine();
                } 
                System.out.println("t=0 Reading t: "+filenums);
                ufile = datadir+"u_"+filenums+".dat";
                vfile = datadir+"v_"+filenums+".dat";
                //String viscfile = datadir+"\\viscofm_"+fnum+".dat";
                //elfile = datadir+"el_"+filenums+".dat";
                System.out.println(ufile+" "+vfile+" "+elfile);
                u = IOUtils.readFileDoubleArray(ufile,recordsPerFile,N*depthLayers," ",true);
                v = IOUtils.readFileDoubleArray(vfile,recordsPerFile,N*depthLayers," ",true);
                //double[][] viscofm = readFileDoubleArray(viscfile,recordsPerFile,N*10," ",false);
                //el = readFileDoubleArray(elfile,recordsPerFile,M*depthLayers," ",false);
                //double[][] sal = readFileDoubleArray(sfile,recordsPerFile,M*10," ",false);
                
                filenums = ""+(fnum+1);
                if (location.equalsIgnoreCase("minch_continuous"))
                {
                    filenums = br.readLine();
                }
                System.out.println("t=0 Reading t+1: "+filenums);
                ufile1 = datadir+"u_"+filenums+".dat";
                vfile1 = datadir+"v_"+filenums+".dat";
                //elfile1 = datadir+"el_"+filenums+".dat";
                u1 = IOUtils.readFileDoubleArray(ufile1,recordsPerFile,N*depthLayers," ",true);
                v1 = IOUtils.readFileDoubleArray(vfile1,recordsPerFile,N*depthLayers," ",true);
                //el1 = readFileDoubleArray(elfile1,recordsPerFile,M*depthLayers," ",false);
//                double usum=0,u1sum=0;
//                for (int i = 0; i < recordsPerFile; i++)
//                {
//                    for (int j = 0; j < recordsPerFile; j++)
//                    {
//                        usum+=u[i][j];
//                        u1sum+=u1[i][j];
//                    }
//                }
//                System.out.println("U Array Sums: u="+usum+" u1="+u1sum);
            } 
            // Interim timesteps:
            // - switch "next file" to being "current file"
            // - read new "next file"
            else if (fnum > firstday && fnum < lastday)
            {
                System.out.println("**** Reading data day "+fnum);
                // At this point, "filenums" is the name prefix of the second set of files read
                // (as determined in second half of loop case above)
                ufile = datadir+"u_"+filenums+".dat";
                vfile = datadir+"v_"+filenums+".dat";
                //elfile = datadir+"el_"+filenums+".dat";
                u = IOUtils.readFileDoubleArray(ufile,recordsPerFile,N*depthLayers," ",true);
                v = IOUtils.readFileDoubleArray(vfile,recordsPerFile,N*depthLayers," ",true);
                //el = readFileDoubleArray(elfile,recordsPerFile,N*depthLayers," ",true);
                // Read in the next file, so that t2+1 is available for time interpolation
                filenums = ""+(fnum+1);
                if (location.equalsIgnoreCase("minch_continuous"))
                {
                    filenums = br.readLine();
                }
                System.out.println("t!=0 Reading t+1: "+filenums);
                ufile1 = datadir+"u_"+filenums+".dat";
                vfile1 = datadir+"v_"+filenums+".dat";
                //elfile1 = datadir+"el_"+filenums+".dat";
                u1 = IOUtils.readFileDoubleArray(ufile1,recordsPerFile,N*depthLayers," ",true);
                v1 = IOUtils.readFileDoubleArray(vfile1,recordsPerFile,N*depthLayers," ",true);
                //el1 = readFileDoubleArray(elfile1,recordsPerFile,N*depthLayers," ",true);
            } 
            // Last time step:
            // - switch "next file" to being "current file"
            else
            {
                System.out.println("t!=0 Reading t+1: "+filenums);
                ufile = datadir+"u_"+filenums+".dat";
                vfile = datadir+"v_"+filenums+".dat";
                //elfile = datadir+"el_"+filenums+".dat";
                u = IOUtils.readFileDoubleArray(ufile,recordsPerFile,N*depthLayers," ",true);
                v = IOUtils.readFileDoubleArray(vfile,recordsPerFile,N*depthLayers," ",true);
                //el = readFileDoubleArray(elfile,recordsPerFile,N*depthLayers," ",true);

            }

            
            int firsttime=0;
            int lasttime=recordsPerFile;
            
            // set an initial tide state
            String tideState = "flood";
            
            System.out.println("Free particles    = "+nfreeparts);
            System.out.println("Viable particles  = "+nViable);
            System.out.println("Arrived particles = "+nSettled);
            System.out.println("Boundary exits    = "+nBoundary);
            
            // default, run loop forwards
            //for (int tt = firsttime; tt <= 3; tt++)
            for (int tt = firsttime; tt <= recordsPerFile-1; tt++)
            // alternatively, run loop backwards
            //for (int tt = lasttime; tt >= firsttime; tt--)
            {
                //System.out.printf("--------- TIME %d ----------\n",tt);
                System.out.printf("%d ",tt);

                boolean debug = false;
                if (debug==true)
                    {
                        IOUtils.particleLocsToFile(particles,nparts,0,"particlelocations_t"+tt+".out");
                    }
                
                nfreeparts = (int)Math.min(nparts,nfreeparts+Math.floor((double)nparts/25.0));
                //nfreeparts = nparts;
                
                for (int st = 0; st < stepsPerStep; st++)
                {
                    //System.out.print(",");
                    //System.out.println("nfreeparts = "+nfreeparts);
                    for (int i = 0; i < nfreeparts; i++)
                    {
                        particles[i].increaseAge(dt/3600.0);
                        //System.out.printf("PARTICLE %d \n",i);
                        if (particles[i].getArrived()==false)
                        {
                            //System.out.println("particle able to move");
                            int elemPart = particles[i].getElem();
                            //System.out.printf("%d\n",elemPart);
                            // set and get the DEPTH layer for the particle based on tide state
//                            if (tt>0)
//                            {
//                                if (el[tt][trinodes[elemPart][0]]>el[tt-1][trinodes[elemPart][0]])
//                                {
//                                    particles[i].setDepthLayer(behaviour,"flood");
//                                } else {
//                                    particles[i].setDepthLayer(behaviour,"ebb");
//                                }
//                            }
                            int dep = particles[i].getDepthLayer();
                            //System.out.println("Depth ="+dep);
                            
                            // Find the salinity in the neighbourhood of the particle (used to compute instantaneous mortality rate).
                            // This is stored at NODES as opposed to ELEMENT CENTROIDS.
                            // So need to get the value from each of the corners and calculate 
                            // a value at the particle location (similar method to getting velocity from nearest centroids).
                            if (st == 0)
                            {
                                double salinity = 0;
                                double mort = 0;
//                                if (calcMort == true)
//                                {
//                                    salinity = particles[i].salinity(tt,sal,trinodes);
//                                    particles[i].setMortRate(salinity);
//                                }
                                particles[i].setDensity();
                            }
                              
                            double advectStep[] = new double[2];
//                            System.out.printf("ADVECT: Euler=[%.3e,%.3e] RK4=[%.3e,%.3e]\n",
//                                    advectStep[0],advectStep[1],advectStep2[0],advectStep2[1]);
                            if (rk4==true)
                            {
                                advectStep = particles[i].rk4Step(u, v, u1, v1, 
                                    neighbours, uvnode,  nodexy, trinodes, allelems,
                                    tt, st, dt, stepsPerStep, recordsPerFile, fnum, lastday, depthLayers);   
                            }
                            else
                            {
                                advectStep = particles[i].eulerStepOld(u, v, u1, v1, neighbours, uvnode, 
                                    tt, st, dt, stepsPerStep, recordsPerFile, fnum, lastday, depthLayers, 
                                    spatialInterpolate, timeInterpolate);
                            }
                            
                            // Reverse velocities if running backwards
                            if (backwards == true)
                            {
                                advectStep[0] = -advectStep[0];
                                advectStep[1] = -advectStep[1];
                            }
                            
                            // 3. Calculate diffusion (random walk step)
//                            if (variableDiff==true)
//                            {
//                                D_h=1000*viscofm[tt][elemPart*10+dep];
//                            }
                            // doing this (below) means that can set zero diffusion, even if read in diffusion values
                            if (diffusion==false)
                            {
                                D_h=0;
                            }
                            //rand('twister',sum(100*clock)); %resets it to a different state each time.
                            double diff_X = Math.sqrt(6*D_h*dt/(double)stepsPerStep);
                            double diff_Y = Math.sqrt(6*D_h*dt/(double)stepsPerStep);
                            double[] behave_uv = particles[i].behaveVelocity(behaviour);
                            // 4. update particle location
                            //double ran1 = BoxMueller.bmTransform(ran.raw(),ran.raw());
                            //double ran2 = BoxMueller.bmTransform(ran.raw(),ran.raw());
                            double ran1 = 2.0*ran.raw()-1.0;
                            double ran2 = 2.0*ran.raw()-1.0;
                            //System.out.println("D_h = "+D_h+" diff_X = "+diff_X+" diff_Y "+diff_Y+" ran1 = "+ran1+" ran2 = "+ran2);
                            //System.out.println("Distances travelled: X "+dt*water_U+" "+diff_X*ran1+" Y "+dt*water_U+" "+diff_Y*ran2);

                            double newlocx=particles[i].getLocation()[0]+advectStep[0]+dt*behave_uv[0]+diff_X*ran1; // simplest possible "Euler"
                            double newlocy=particles[i].getLocation()[1]+advectStep[1]+dt*behave_uv[1]+diff_Y*ran2;

                            //System.out.println("Old = ("+particles[i].getLocation()[0]+", "+particles[i].getLocation()[1]+") --- New = ("+newlocx+", "+newlocy+")");
                            
                            // find element containing particle and update seach counts for diagnosis
                            int[] c = particles[i].findContainingElement(newlocx, newlocy, elemPart, 
                                    nodexy, trinodes, neighbours, allelems);
                            int whereami = c[0];
                            for (int j = 0; j < 5; j++)
                            {
                                searchCounts[j] += c[j+1];
                            }

                            // if particle is within the mesh, update location normally and save the distance travelled
                            if (whereami != -1)
                            {
                                double distTrav = Math.sqrt((particles[i].getLocation()[0]-newlocx)*(particles[i].getLocation()[0]-newlocx)+
                                        (particles[i].getLocation()[1]-newlocy)*(particles[i].getLocation()[1]-newlocy));
                                if (distTrav>maxDistTrav)
                                {
                                    maxDistTrav=distTrav;
                                }
                                if (distTrav<minDistTrav)
                                {
                                    minDistTrav=distTrav;
                                }
                                particles[i].setLocation(newlocx,newlocy);
                                particles[i].setElem(whereami);
                                //System.out.printf("** MOVED **, new elem = %d (dist = %f)\n",particles[i].getElem(),Math.sqrt((newlocx-uvnode[particles[i].getElem()][0])*(newlocx-uvnode[particles[i].getElem()][0])+(newlocy-uvnode[particles[i].getElem()][1])*(newlocy-uvnode[particles[i].getElem()][1])));
                            }
                            
                            // if particle has skipped out of the model domain, place it at the nearest element centroid
                            if (whereami == -1)
                            {
                                int closest=Particle.nearestCentroid(particles[i].getLocation()[0],particles[i].getLocation()[1],uvnode);
                                //fprintf('x%d',closest);
                                particles[i].setLocation(uvnode[closest][0],uvnode[closest][1]);
                                particles[i].setElem(closest);
                            }

                            // set particle to become able to settle after a predefined time
                            //if ((day-firstday)*24+tt==viabletime)
                            if (particles[i].getAge()>viabletime/3600.0 && particles[i].getViable()==false)
                            //if (time>viabletime)
                            {
                                //System.out.println("Particle became viable");
                                //if (ran.raw()<dev_perstep)
                                //{
                                    particles[i].setViable(true);
                                    nViable++;
                                    // save location information
                            
                                
                                    //System.out.println(time+" "+particles[i].getAge());
                                    //System.out.println("I can settle");
                                //}
                                pstepsMature[elemPart][1]+=(dt/3600)*1.0/stepsPerStep;
                            } else {
                                pstepsImmature[elemPart][1]+=(dt/3600)*1.0/stepsPerStep;
                            }
                            
                            // check whether the particle has gone within a certain range of one of the boundary nodes
                            // (make it settle there, even if it is inviable)
                            for (int loc = 0; loc < open_BC_locs.length; loc++)
                                {
                                    double dist = Math.sqrt((particles[i].getLocation()[0]-open_BC_locs[loc][1])*(particles[i].getLocation()[0]-open_BC_locs[loc][1])+
                                            (particles[i].getLocation()[1]-open_BC_locs[loc][2])*(particles[i].getLocation()[1]-open_BC_locs[loc][2]));
                                    if (dist < 1500)
                                    {
                                        //System.out.printf("Boundary stop: %d at %d\n",i,loc);
                                        particles[i].setArrived(true);
                                        particle_info[i][1] = -loc;//(int)startlocs[loc][0];
                                        particle_info[i][2] = (int)particles[i].getAge();//((day-firstday)*24+tt);
                                        nBoundary++;
                                        break;
                                        
                                    }
                                }
                                
                            // if able to settle, is it close to a possible settlement
                            // location?
                            //System.out.println("Patricle age = "+particles[i].getAge()+" Viabletime/3600 = "+viabletime/3600.0+" viable = "+particles[i].getViable());
                            if (particles[i].getViable()==true)
                            {
                                //System.out.println(particles[i].getViable());
                                
                                for (int loc = 0; loc < startlocs.length; loc++)
                                {
                                    double dist = Math.sqrt((particles[i].getLocation()[0]-startlocs[loc][1])*(particles[i].getLocation()[0]-startlocs[loc][1])+
                                            (particles[i].getLocation()[1]-startlocs[loc][2])*(particles[i].getLocation()[1]-startlocs[loc][2]));
                                    if (dist < thresh)
                                    {
                                        //System.out.printf("settlement: %d at %d\n",i,loc);
                                        particles[i].setArrived(true);
                                        particle_info[i][1] = loc;//(int)startlocs[loc][0];
                                        particle_info[i][2] = (int)time;//((day-firstday)*24+tt);
                                        settle_density[i][0] = particles[i].getDensity();
                                        nSettled++;
                                        break;
                                    }
                                }    
                            }
                            
                            
                            
                            

                        }
                    }
                    time+=dt/3600.0;
//                    if ((day-firstday)*24+tt==viabletime && st==0)
//                    {
//                        System.out.println("print particle locations to file");
//                        particleLocsToFile(particles,"particlelocations"+thresh+"_"+viabletime+".out");
//                    }
                    //System.out.println("Saving particle locations");
//                    if (lochfyne != true)// && st%25==0)
//                    {   
//                        printCount++;
//                        //int[] testTrackSites = {28, 979, 78, 56, 1215, 505, 788, 359, 335, 857, 551};
//                        //for (int i = 0; i < testTrackSites.length; i++)
//                        //{
//                            particleLocsToFile(particles,printCount,"particlelocations_all.out");
//                            //particleLocsToFile2(time,particles,testTrackSites[i],"loc_"+i+".out");
//                        //}
//                    }
                    //System.out.println("-----------");
                    
                    // end of particle loop
                    calcCount++;
                }
                
                printCount++;
                //if(firstday==3)
                //{
                    IOUtils.particleLocsToFile(particles,nparts_per_site,printCount,"particlelocations_all"+suffix+".out");
                //}
//                for (int i = 0; i < Math.min(20,nparts); i++)
//                    {
//                        particleLocsToFile2(time,particles,i,"loc_"+i+".out");
//                    }
                // print particle locations here - once every hour
                stepcount++;
            }
            System.out.printf("\n");

        }  
        System.out.printf("\nelement search counts: %d %d %d %d %d\n",searchCounts[0],searchCounts[1],searchCounts[2],searchCounts[3],searchCounts[4]);
        System.out.printf("transport distances: min = %.4e, max = %.4e\n", minDistTrav, maxDistTrav);

        IOUtils.writeDoubleArrayToFile(pstepsImmature,"pstepsImmature"+suffix+".out");
        IOUtils.writeDoubleArrayToFile(pstepsMature,"pstepsMature"+suffix+".out");
        IOUtils.writeIntegerArrayToFile(particle_info,"particle_info"+suffix+".out");
        IOUtils.writeDoubleArrayToFile(settle_density,"settle_density"+suffix+".out");
        IOUtils.particleLocsToFile(particles,nparts,0,"particlelocations"+suffix+".out");
        IOUtils.writeDoubleArrayToFile(particle1Velocity,"particle1velocity"+suffix+".out");
        IOUtils.writeDoubleArrayToFile(particle1Location,"particle1location"+suffix+".out");
        IOUtils.particleLocsToFile1(particles,"particlelocations_end"+suffix+".out");
        
        
        //System.out.println(new Date().toString());
        long endTime = System.currentTimeMillis();
        System.out.println("Elapsed time = "+(endTime-startTime)/1000.0);
        
        // scale psteps
//        tot_psteps=tmax*nparts;
//        double[] psteps2 = new double[N];
//        for (int i = 0; i < N; i++)
//        {
//            psteps2[i]=(double)psteps[i][1]/tot_psteps;
//        }       
    }
    

    public double[][] makeCentroidXY()
    {
        double[][] out=new double[1][1];
        return out;
    }

    public void setupOutput()
    {

    }

    public void writeOutput()
    {
    
    }
    
     
    
    public double nc_varget(String filename, String variable)
    {
        double out=0;
        return out;
    }

}
