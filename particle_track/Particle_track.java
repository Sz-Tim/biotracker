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
        
        boolean cluster = false;
        
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
               
        // options for behaviour of particles - depth aspects are set at each
        // timestep by function
        /**
         * 1 - passive, stay on surface
         * 2 - passive, stay on bottom (layer 10)
         * 3 - passive, stay in mid layer (layer 5)
         * 6 - top during flood tides, mid during ebb (local)
         * 7 - mid during flood tides, bed during ebb (local)
         * 8 - top during flood tides, bed during ebb (local)
         */
        int behaviour=1;
        
        boolean timeInterpolate = false;
        boolean spatialInterpolate = false;
        double D_h = 0.1;
        
        /**
         * SETUP CONFIGURATION FOR THE MINCH CHAPTER CASES
         */
        String habitat = "";
        String suffix = "";
        double diffusionMultiplier = 1;
        if (args.length == 6 || args.length == 7)
        {
            System.out.println("Standard minch chapter args (5)");
            // set values for many of the things that used to be read in
            nparts=Integer.parseInt(args[5]);
            recordsPerFile=6;
            dt=3600;
            stepsPerStep=200;
            thresh=500;
            viabletime=87;
            //viabletime=1000000000;
            location="minch_continuous";
            timeInterpolate = false;
            spatialInterpolate = false;
            
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
            if (args.length == 7)
            {
                suffix=args[6];
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
                System.out.println("nparts = "+nparts+"\nsettlement threshold distance = "+thresh+"\nviabletime = "+viabletime+" D_h = "+D_h);
            }
            catch (IndexOutOfBoundsException e)
            {
                System.err.println("Incorrect number of input parameters, found "+args.length);
            }
        }
        
        
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
            System.out.println("DATA DIRECTORY = "+datadir2);
        }
        else
        {
            System.out.println("No location......");
            datadir = "D:\\dima_data\\"+datadir;
            System.out.println("DATA DIRECTORY = "+datadir);
        }
        nodexy = readFileDoubleArray(datadir2+"nodexy.dat",M,2," ",true); // the original mesh nodes
        uvnode = readFileDoubleArray(datadir2+"uvnode.dat",N,2," ",true); // the centroids of the elements
        trinodes = readFileIntArray(datadir2+"trinodes.dat",N,3," ",true); // the corners of the elements
        neighbours = readFileIntArray(datadir2+"neighbours.dat",N,3," ",true);
        
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

        if ((location.equalsIgnoreCase("lorn") == true || location.equalsIgnoreCase("minch") == true) && viabletime==0)
        {
            // case when for viable at
            viabletime=(lastday-firstday+1)*recordsPerFile*dt/2;
            System.out.println("ignoring input viable time, set to "+viabletime);
        }
        else if (viabletime==666)
        {
            // ensure that particles never become viable...
            viabletime=(lastday-firstday+1)*recordsPerFile*dt*2;
            System.out.println("ignoring input viable time, set to "+viabletime);
        }
        else
        {
            viabletime=viabletime*3600;
            System.out.println("Using input viabletime: "+viabletime/3600.0+"hrs "+viabletime+"s");
        }
        
        int nrecords=(lastday-firstday+1)*recordsPerFile;
        int inviabletime=nrecords;
        
        dt=dt/(double)stepsPerStep;
        double dev_perstep = Math.pow(0.1,dt);
        System.out.println("dt = "+dt+" dev_perstep = "+dev_perstep);

        System.out.println("behaviour = "+behaviour);
//        if (timeInterpolate==true)
//        {
//            System.out.println("WARNING: time interpolation not operational in last record of each file");
//        }

        // load array of start node IDs (as stored by matlab)
        double startlocs[][] = new double[10][3];
        double open_BC_locs[][] = new double[1][3];
        
        // FYNE configurations
        if (location.equalsIgnoreCase("fyne"))
        {
            System.out.println("****** Loch Fyne tracking ******");
            startlocs = setupStartLocationsFyneFarms();
        }
        // North Minch fish farm locations
        else if (location.equalsIgnoreCase("minch") || location.equalsIgnoreCase("minch_continuous"))
        {
            //String sitedir="C:\\Users\\sa01ta\\Documents\\westcoast_tracking_summary\\habitat_sites\\sites_to_use\\";
            String sitedir="C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\site_locations\\";
            if (cluster==true)
            {
                sitedir=basedir+"minch_sites/";
                
            }
            open_BC_locs = readFileDoubleArray(sitedir+"open_boundary_locs.dat",120,3," ",true);
            if (habitat.equalsIgnoreCase("coast"))
            {
                startlocs = readFileDoubleArray(sitedir+"minchcoast_1km_utm2.dat",3892,4," ",true);
            }
            else if (habitat.equalsIgnoreCase("fishfarm"))
            {
                //startlocs = readFileDoubleArray(sitedir+"140124_ms_site_utm_area.dat",122,4," ",true); 
                startlocs = readFileDoubleArray(sitedir+"150508_ms_site_utm_FMA.dat",122,4," ",true);
            }
            else if (habitat.equalsIgnoreCase("nephrops"))
            {
                startlocs = readFileDoubleArray(sitedir+"sublittoralmud_utm.dat",5896,4," ",true);
            }
            else if (habitat.equalsIgnoreCase("spill"))
            {
                sitedir="C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\site_locations\\";
                startlocs = readFileDoubleArray(sitedir+"150408_spill_utm_10vary.dat",16,4," ",true);
                //startlocs = setupStartLocationsSpill(uvnode, 36301);
                viabletime=(lastday+1)*recordsPerFile*dt*3600;
                System.out.println("ignoring input viable time, set to "+viabletime);
                
            }
            else
            {
                System.out.println("****** North Minch tracking ******");
                //startlocs = readFileDoubleArray("D:\\"+datadir+"\\site_xy.dat",17,3," ",true);
                //startlocs = readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\farms_utm2.dat",213,3," ",true);
                startlocs = readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\site_locations\\140124_ms_site_utm_area.dat",122,4," ",true);
                open_BC_locs = readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\hydro_mesh_run\\open_boundary_locs.dat",120,3," ",true);
            }
         
        }
        // Firth of Lorn boundary locations 
        // boundary_940pt in article submitted to Ecography (coastal sites c1km apart).
        // others (_with_interest or _with_search) include offshore renewable energy sites
        // (as in article submitted to Jornal of Applied Ecology).
        else
        {
            startlocs = readFileDoubleArray(basedir+"lorncoast\\boundary_940pt.dat",940,3," ",true);
        //double startlocs[][] = readFileDoubleArray(basedir+"lorncoast\\boundary_with_interest.dat",966,3," ",true);
            //startlocs = readFileDoubleArray(basedir+"lorncoast\\boundary_with_search.dat",1252,3," ",true);
            //startlocs = readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\lorn\\120903_renewableimpact\\131107_revision\\startlocs_array10.dat",971,3," ",true);
            //startlocs = readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\lorn\\120903_renewableimpact\\131107_revision\\startlocs_array25.dat",1018,3," ",true);
            //startlocs = readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\lorn\\120903_renewableimpact\\131107_revision\\startlocs_array50.dat",1095,3," ",true);
            //startlocs = readFileDoubleArray(basedir+"infrastructure\\portsmarinas_num.txt",41,3," ",true);
            //startlocs = readFileDoubleArray(basedir+"infrastructure\\portsmarinas_search.dat",353,3," ",true);
            open_BC_locs = readFileDoubleArray("C:\\Users\\sa01ta\\Documents\\lorn\\120903_renewableimpact\\131107_revision\\open_boundary_locs.dat",93,3," ",true);
        }
        
        
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
            particles[i] = new Particle(xstart[i],ystart[i]);
            particle_info[i][0]=startid[i];//(int)startlocs[startid[i]][0];
            particles[i].setElem(startElem[i]);
        }
        
        particleLocsToFile(particles,nparts,0,"particlelocations_start.out");

        // ------------------- loop 2 = timestep ----------------------------
        System.out.println("Starting time loop");


        int count0=0;
        int count1=0;
        int count2=0;
        int count3=0;
        int count4=0;

        double minDistTrav=10000000;
        double maxDistTrav=0;

        int stepcount=0;
        int calcCount=0;
        double time=0;
        
        double[][] particle1Velocity = new double[nrecords*stepsPerStep][3];
        double[][] particle1Location = new double[nrecords*stepsPerStep][2];
        
        System.out.printf("Viable time        = %f\n",viabletime);
        System.out.printf("Threshold distance = %d\n",thresh);
        System.out.printf("Diffusion          = %f\n",D_h);
        
        int printCount=0;
        
        // Initial value for br set - this particular one is only used in case of "spill" habitat (identical to that file list)
        BufferedReader br = new BufferedReader(new FileReader("C:\\Users\\sa01ta\\Documents\\Sealice_NorthMinch\\hydro_mesh_run\\150408_spill_list1.dat"));
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

        int nfreeparts = 0;
        int nViable = 0;
        int nBoundary = 0;
        int nSettled = 0;
        
        String filenums = "";
        
        // default, run loop forwards
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
            double[][] el = new double[recordsPerFile][N*depthLayers];
            double[][] u1 = new double[recordsPerFile][N*depthLayers];
            double[][] v1 = new double[recordsPerFile][N*depthLayers];
            double[][] el1 = new double[recordsPerFile][N*depthLayers];
            
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
                elfile = datadir+"el_"+filenums+".dat";
                System.out.println(ufile+" "+vfile+" "+elfile);
                u = readFileDoubleArray(ufile,recordsPerFile,N*depthLayers," ",true);
                v = readFileDoubleArray(vfile,recordsPerFile,N*depthLayers," ",true);
                //double[][] viscofm = readFileDoubleArray(viscfile,recordsPerFile,N*10," ",false);
                el = readFileDoubleArray(elfile,recordsPerFile,M*depthLayers," ",false);
                //double[][] sal = readFileDoubleArray(sfile,recordsPerFile,M*10," ",false);
                
                filenums = ""+(fnum+1);
                if (location.equalsIgnoreCase("minch_continuous"))
                {
                    filenums = br.readLine();
                }
                System.out.println("t=0 Reading t+1: "+filenums);
                ufile1 = datadir+"u_"+filenums+".dat";
                vfile1 = datadir+"v_"+filenums+".dat";
                elfile1 = datadir+"el_"+filenums+".dat";
                u1 = readFileDoubleArray(ufile1,recordsPerFile,N*depthLayers," ",true);
                v1 = readFileDoubleArray(vfile1,recordsPerFile,N*depthLayers," ",true);
                el1 = readFileDoubleArray(elfile1,recordsPerFile,M*depthLayers," ",false);
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
                // Cloning arrays doesn't seem to work: read data in again 
                // (what was file t+1 is now file t)
                // At this point, "filenums" is the name prefix of the second set of files read
                // (as determined in second half of loop case above)
                ufile = datadir+"u_"+filenums+".dat";
                vfile = datadir+"v_"+filenums+".dat";
                elfile = datadir+"el_"+filenums+".dat";
                u = readFileDoubleArray(ufile,recordsPerFile,N*depthLayers," ",true);
                v = readFileDoubleArray(vfile,recordsPerFile,N*depthLayers," ",true);
                el = readFileDoubleArray(elfile,recordsPerFile,N*depthLayers," ",true);
                //u = u1.clone();
                //v = v1.clone();
                //el = el1.clone();
                // Read in the next file, so that t2+1 is available for time interpolation
                filenums = ""+(fnum+1);
                if (location.equalsIgnoreCase("minch_continuous"))
                {
                    filenums = br.readLine();
                }
                System.out.println("t!=0 Reading t+1: "+filenums);
                ufile1 = datadir+"u_"+filenums+".dat";
                vfile1 = datadir+"v_"+filenums+".dat";
                elfile1 = datadir+"el_"+filenums+".dat";
                u1 = readFileDoubleArray(ufile1,recordsPerFile,N*depthLayers," ",true);
                v1 = readFileDoubleArray(vfile1,recordsPerFile,N*depthLayers," ",true);
                el1 = readFileDoubleArray(elfile1,recordsPerFile,N*depthLayers," ",true);
            } 
            // Last time step:
            // - switch "next file" to being "current file"
            else
            {
                System.out.println("t!=0 Reading t+1: "+filenums);
                ufile = datadir+"u_"+filenums+".dat";
                vfile = datadir+"v_"+filenums+".dat";
                elfile = datadir+"el_"+filenums+".dat";
                u = readFileDoubleArray(ufile,recordsPerFile,N*depthLayers," ",true);
                v = readFileDoubleArray(vfile,recordsPerFile,N*depthLayers," ",true);
                el = readFileDoubleArray(elfile,recordsPerFile,N*depthLayers," ",true);
                //u = u1.clone();
                //v = v1.clone();
                //el = el1.clone();
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
                        particleLocsToFile(particles,nparts,0,"particlelocations_t"+tt+".out");
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
                            // set the depth layer for the particle based on tide state
                            if (tt>0)
                            {
                                if (el[tt][trinodes[elemPart][0]]>el[tt-1][trinodes[elemPart][0]])
                                {
                                    particles[i].setDepthLayer(behaviour,"flood");
                                } else {
                                    particles[i].setDepthLayer(behaviour,"ebb");
                                }
                            }
                            
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
                            
                            // Set velocity array here first and fill with the right values. Need an array with velocity "now" 
                            // and one with velocity at "next" time step, to be populated before go to time interpolation bit 
                            // (which should happen last, after any possible spatial interpolation).
                            double water_U = 0;
                            double water_V = 0;
                            
                            if (spatialInterpolate==false)
                            {
                                // output from matlab is in rows for each time
                                // and then columns for each depth within each element
                                water_U = u[tt][elemPart*depthLayers+dep];
                                water_V = v[tt][elemPart*depthLayers+dep];
                                //System.out.printf("Vel(element %d) = %.4f %.4f\n",elemPart,water_U,water_V);
                                if (timeInterpolate == true)
                                {
                                    if (tt < recordsPerFile-1)
                                    {
                                        water_U = u[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(u[tt+1][elemPart*10+dep]-u[tt][elemPart*depthLayers+dep]);
                                        water_V = v[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(v[tt+1][elemPart*10+dep]-v[tt][elemPart*depthLayers+dep]);
                                    } 
                                    else
                                    {
                                        //water_U = u[tt][elemPart*depthLayers+dep];
                                        //water_V = v[tt][elemPart*depthLayers+dep];
                                        if (fnum < lastday)
                                        {
                                            
                                            water_U = u[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(u1[0][elemPart*10+dep]-u[tt][elemPart*depthLayers+dep] );
                                            water_V = v[tt][elemPart*depthLayers+dep] + ((double)st/(double)stepsPerStep)*(v1[0][elemPart*10+dep]-v[tt][elemPart*depthLayers+dep]);                             
                                        }
                                        // Very last timestep: just fix velocity to be same over entire last hour
                                        else
                                        {
                                            water_U = u[tt][elemPart*depthLayers+dep];
                                            water_V = v[tt][elemPart*depthLayers+dep];
                                        }
                                    }
                                } 
                                
                            }
                            else if (spatialInterpolate == true)
                            {
                                // Print element values for comparison
                                //System.out.printf("Vel(element %d) = %.4f %.4f\n",elemPart,u[tt][elemPart*10+dep],v[tt][elemPart*10+dep]);
//                                int nearest = Particle.nearestCentroid(particles[i].getLocation()[0], particles[i].getLocation()[1], uvnode);
//                                int which = Particle.whichElement(particles[i].getLocation()[0], particles[i].getLocation()[1], allelems, nodexy, trinodes);
                                
                                particles[i].setNrListToNeighbourCells(neighbours,uvnode);
//                                if (nearest != which)
//                                {
//                                    System.out.printf("Particle "+i+" Elem "+particles[i].getElem()+" nearest "+nearest+" whichElem "+which+"\n");
//                                    double d1 = Particle.distanceEuclid(particles[i].getLocation()[0],particles[i].getLocation()[1],uvnode[which][0],uvnode[which][1]);
//                                    double d2 = Particle.distanceEuclid(particles[i].getLocation()[0],particles[i].getLocation()[1],uvnode[nearest][0],uvnode[nearest][1]);
//                                    System.out.printf("dist(elemPart) = "+d1+" dist(nearest) = "+d2+"\n");
//                                    double[][] a = particles[i].getNrList();
//                                    System.out.printf("Neighbour cells: \n"+a[0][0]+" "+a[0][1]+"\n"+a[1][0]+" "+a[1][1]+"\n"+a[2][0]+" "+a[2][1]+"\n"+a[3][0]+" "+a[3][1]+"\n");
//                                }
                                
                                //particles[i].setNrList(uvnode);
                                double[] vel = particles[i].velocityFromNearestList(tt,u,v);
                                water_U = vel[0];
                                water_V = vel[1];
                                if (timeInterpolate == true)
                                {
                                    double[] velplus1 = new double[2];
                                    if (tt < recordsPerFile-1)
                                    {
                                        velplus1 = particles[i].velocityFromNearestList(tt+1,u,v);
                                        
                                    }
                                    else
                                    {
                                        //velplus1 = particles[i].velocityFromNearestList(0,u1,v1);
                                        if (fnum < lastday)
                                        {
                                            velplus1 = particles[i].velocityFromNearestList(0,u1,v1);
                                        }
                                        // Very last timestep: just fix velocity to be same over entire last hour
                                        else
                                        {
                                            velplus1 = particles[i].velocityFromNearestList(tt,u,v);
                                        }
                                    }
                                    water_U = vel[0] + ((double)st/(double)stepsPerStep)*(velplus1[0]-vel[0]);
                                    water_V = vel[1] + ((double)st/(double)stepsPerStep)*(velplus1[1]-vel[1]);
                                }
                            }
//                            if (st == 1)
//                            {
//                                System.out.printf("Vel (particle %d element %d) = %.4f %.4f\n",i,elemPart,water_U,water_V);
//                            }
                            
                            // below velocities used if running backwards
                            if (backwards == true)
                            {
                                water_U = -water_U;
                                water_V = -water_V;
                            }
                            
                            //[water_U,water_V]=calc_vel1(particles(i),[uchunk vchunk]);
                            // 3. Calculate diffusion 
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
                            double diff_X = Math.sqrt(6*D_h*dt/(double)stepsPerStep);// /(double)stepsPerStep);
                            double diff_Y = Math.sqrt(6*D_h*dt/(double)stepsPerStep);// /(double)stepsPerStep);    //+/- is random so direction doesn't matter
                            //diff_X=0; diff_Y=0;
                            double[] behave_uv = particles[i].behaveVelocity(behaviour);
                            // 4. update particle location
                            //double ran1 = BoxMueller.bmTransform(ran.raw(),ran.raw());
                            //double ran2 = BoxMueller.bmTransform(ran.raw(),ran.raw());
                            double ran1 = 2.0*ran.raw()-1.0;
                            double ran2 = 2.0*ran.raw()-1.0;
                            //System.out.println("D_h = "+D_h+" diff_X = "+diff_X+" diff_Y "+diff_Y+" ran1 = "+ran1+" ran2 = "+ran2);
                            //System.out.println("Distances travelled: X "+dt*water_U+" "+diff_X*ran1+" Y "+dt*water_U+" "+diff_Y*ran2);
                            double newlocx=particles[i].getLocation()[0]+dt*water_U+dt*behave_uv[0]+diff_X*ran1; // simplest possible "Euler"
                            double newlocy=particles[i].getLocation()[1]+dt*water_V+dt*behave_uv[1]+diff_Y*ran2;

                            //System.out.println("Old = ("+particles[i].getLocation()[0]+", "+particles[i].getLocation()[1]+") --- New = ("+newlocx+", "+newlocy+")");
                            
                            //double newlocx=particles[i].getLocation()[0]+diff_X*ran1; // simplest possible "Euler"
                            //double newlocy=particles[i].getLocation()[1]+diff_Y*ran2;

                            if (i == 4)
                            {                               
                                particle1Velocity[calcCount][0]=time;
                                particle1Velocity[calcCount][1]=Math.sqrt(water_U*water_U+water_V*water_V);
                                particle1Velocity[calcCount][2]=elemPart;
                                //System.out.println(particle1Velocity[calcCount][0]+" "+particle1Velocity[calcCount][1]);
                                
                                particle1Location[calcCount][0]=time;
                                particle1Location[calcCount][1]=Math.sqrt(Math.pow(xstart[i]-newlocx,2.0)+Math.pow(ystart[i]-newlocy,2.0));
                            }
                            
                            
                            // search progressively further from previous element for new particle location
                            //fprintf('+');
                            int[] elems = new int[1];
                            elems[0] = elemPart;
                            int whereami=Particle.whichElement(newlocx,newlocy,elems,nodexy,trinodes);
                            count0=count0+1;
                            if (whereami==-1)
                            {
                                int[] elems0 = neighbours[elemPart];
                                count1++;
                                whereami=Particle.whichElement(newlocx,newlocy,elems0,nodexy,trinodes);
                                // if fails, look in nearest 10 (id numerical)
                                if (whereami==-1)
                                {
                                    //fprintf('a');
                                    //checkfirst
                                    count2=count2+1;
                                
                                    int[] elems1 = new int[10];
                                    for (int j = 0; j < 10; j++)
                                    {
                                        elems1[j] = Math.min(Math.max(elemPart-5+j,0),N-1);
                                    }
                                    whereami=Particle.whichElement(newlocx,newlocy,elems1,nodexy,trinodes);
                                    // if fails, look in nearest 500 (id numerical)
                                    if (whereami==-1)
                                    {
                                        //fprintf('b');
                                        //checkfirst
                                        count3=count3+1;
                                        
                                        int[] elems2 = new int[500];
                                        for (int j = 0; j < 500; j++)
                                        {
                                            elems2[j] = Math.min(Math.max(elemPart-250+j,0),N-1);
                                        }
                                        whereami=Particle.whichElement(newlocx,newlocy,elems2,nodexy,trinodes);
                                        // if this fails, look in all elements
                                        if (whereami==-1)
                                        {
                                            //fprintf('c');
                                            count4=count4+1;
                                            whereami=Particle.whichElement(newlocx,newlocy,allelems,nodexy,trinodes);
                                            //elemPart(i)=nearest_centroid(particles(i).x,particles(i).y,centroids);
                                        }
                                    }
                                }
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
                    particleLocsToFile(particles,nparts_per_site,printCount,"particlelocations_all"+suffix+".out");
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
        System.out.printf("\nelement search counts: %d %d %d %d %d\n",count0,count1,count2,count3,count4);
        System.out.printf("transport distances: min = %.4e, max = %.4e\n", minDistTrav, maxDistTrav);

        writeDoubleArrayToFile(pstepsImmature,"pstepsImmature"+suffix+".out");
        writeDoubleArrayToFile(pstepsMature,"pstepsMature"+suffix+".out");
        writeIntegerArrayToFile(particle_info,"particle_info"+suffix+".out");
        writeDoubleArrayToFile(settle_density,"settle_density"+suffix+".out");
        particleLocsToFile(particles,nparts,0,"particlelocations"+suffix+".out");
        writeDoubleArrayToFile(particle1Velocity,"particle1velocity"+suffix+".out");
        writeDoubleArrayToFile(particle1Location,"particle1location"+suffix+".out");

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

    public static double[][] setupStartLocationsFyneFarms()
    {
        int nsites = 9;
        double startlocs[][]= new double[nsites][3];

        // S->N ordering (reordering sites from original nonsensical order based on:
        // order_SN=[9 7 4 6 5 8 3 1 2];)
        
        startlocs[0][0]=1; startlocs[0][1]=352401; startlocs[0][2]=6190933;
        startlocs[1][0]=2; startlocs[1][1]=354969; startlocs[1][2]=6193169;
        startlocs[2][0]=3; startlocs[2][1]=354246; startlocs[2][2]=6194759;
        startlocs[3][0]=4; startlocs[3][1]=348880; startlocs[3][2]=6199380;
        startlocs[4][0]=5; startlocs[4][1]=352745; startlocs[4][2]=6201735;
        startlocs[5][0]=6; startlocs[5][1]=348606; startlocs[5][2]=6204475;
        startlocs[6][0]=7; startlocs[6][1]=353078; startlocs[6][2]=6206339;
        startlocs[7][0]=8; startlocs[7][1]=357420; startlocs[7][2]=6217200; 
        startlocs[8][0]=9; startlocs[8][1]=361834; startlocs[8][2]=6223063;
        
        // "hypothetical" fishfarm used for SAMS newsletter
        //fishfarms(1,1)=351000;
        //fishfarms(1,2)=6195000;
        return startlocs;
    }
    
    public static double[][] setupStartLocationsSpill(double[][] uvnode, int startelem)
    {
        int nsites = 1;
        double startlocs[][]= new double[nsites][3];
        startlocs[0][0] = 1; startlocs[0][1] = uvnode[startelem][0]; startlocs[0][2] = uvnode[startelem][1];
        return startlocs;
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
    
    public static double[][] readFileDoubleArray(String filename, int rows, int cols, String sep, boolean note)
    {
        double[][] myDouble = new double[rows][cols];
        int x=0, y=0;
        boolean failed = false;
        double sum = 0;
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(filename));	//reading files in specified directory
 
            String line;
            
            outer: while ((line = in.readLine()) != null)	//file reading
            {
                y=0;
                if (x >= rows)
                {
                    System.out.println(filename+" has more rows than expected.");
                    failed = true;
                    break;
                }
                String[] values = line.split(" ");
                for (String str : values)
                {
                    if (y >= cols)
                    {
                        System.out.println(filename+" has more columns than expected: "+y);
                        failed = true;
                        break outer;
                    }
                    double str_double = Double.parseDouble(str);
                    myDouble[x][y]=str_double;
                    sum+=myDouble[x][y];
                    //System.out.print(myDouble[x][y] + " ");
                    y++;

                }

                //System.out.println("");
                x++;
            }
            in.close();
            
        } catch( IOException ioException ) {
            throw new RuntimeException(ioException);
            //System.err.println("******************* Cannot read from file "+filename+" ******************************");
            //failed = true;
        }
        if (note == true && failed == false)
        {
            System.out.printf("Created %dx%d array from file: %s\n",myDouble.length,myDouble[0].length,filename);
            System.out.println("Array sum at read time = "+sum);
        }
        else if (failed == true)
        {
            System.out.println("FAILED to read file "+filename);
            System.exit(1);
        }
        return myDouble;
    }
    
    public static int[][] readFileIntArray(String filename, int rows, int cols, String sep, boolean note)
    {
        int[][] myInt = new int[rows][cols];
        int x=0, y=0;
        boolean failed = false;
        try
        {
            BufferedReader in = new BufferedReader(new FileReader(filename));	//reading files in specified directory
            //System.out.println("in readFileIntArray");
            String line;
            while ((line = in.readLine()) != null)	//file reading
            {
                
                y=0;
                if (x >= rows)
                {
                    System.out.println(filename+" has more rows than expected.");
                    failed = true;
                    break;
                }
                String[] values = line.split(sep);
                for (String str : values)
                {
                    if (y >= cols)
                    {
                        System.out.println(filename+" has more columns than expected.");
                        failed = true;
                        break;
                    }
                    int str_int = Integer.parseInt(str);
                    myInt[x][y]=str_int;
                    //System.out.print(myDouble[x][y] + " ");
                    y++;

                }
                //System.out.println("");
                x++;
            }
            in.close();
            
        
        } catch( Exception e ) {
            System.err.println("******************* Cannot read from file "+filename+" ******************************");
            failed = true;
        }
        if (note == true && failed == false)
        {
            System.out.printf("Created %dx%d array from file: %s\n",myInt.length,myInt[0].length,filename);
        }
        else if (failed == true)
        {
            System.out.println("FAILED to read file "+filename);
            System.exit(1);
        }
        return myInt;
    }
    
    public static void writeDoubleArrayToFile(double[][] variable, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < variable.length; i++)
            {
                for (int j = 0; j < variable[0].length; j++)
                {
                    out.printf("%.4e ",variable[i][j]);
                }
                out.printf("\n");
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    public static void writeIntegerArrayToFile(int[][] variable, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < variable.length; i++)
            {
                for (int j = 0; j < variable[0].length; j++)
                {
                    out.printf("%d ",variable[i][j]);
                }
                out.printf("\n");
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    public static void particleLocsToFile(Particle[] particles, int nparts, int tt, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,true);
            PrintWriter out = new PrintWriter(fstream);
            for (int i = 0; i < particles.length/nparts; i++)
            {
                out.printf("%d %d %f %f %d\n",tt,i,particles[i].getLocation()[0],particles[i].getLocation()[1],particles[i].getElem());
            }
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }
    
    public static void particleLocsToFile2(double time, Particle[] particles, int i, String filename)
    {
        try
        {
            // Create file 
            FileWriter fstream = new FileWriter(filename,true);
            PrintWriter out = new PrintWriter(fstream);
            out.printf("%f %f %f %d\n",time,particles[i].getLocation()[0],particles[i].getLocation()[1],particles[i].getElem());
            //Close the output stream
            out.close();
        }catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    } 
    
    public double nc_varget(String filename, String variable)
    {
        double out=0;
        return out;
    }

}
