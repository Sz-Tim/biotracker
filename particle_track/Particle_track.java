/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

//import java.lang.Object;
////import ucar.nc2.dataset.*;
//import ucar.nc2.NCdumpW;
////import ucar.nc2.Attribute;
////import ucar.nc2.Netcdf;
//import ucar.nc2.NetcdfFile;
//import ucar.nc2.Variable;
//import ucar.nc2.iosp.IOServiceProvider;
//import ucar.ma2.*;
//
//import ucar.ma2.ArrayFloat;
//import ucar.ma2.InvalidRangeException;

import java.io.*;
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
        
        RanMT ran = new RanMT(System.currentTimeMillis());
        
        System.out.println("Reading in data\n");

        /**
         * files are 
         * meshdata.nc:
         * uv_[3-31].nc: u, v
         * 
         * or
         * 
         * nodexy.dat, uvnode.dat, trinodes.dat, neighbours.dat
         * u_[3-31].dat
         * v_[3-31].dat
         */ 
        int N = 25071;
        int M = 14000;
        
        //String meshfile = "D:\\sharedfolder\\May2010_mat\\meshdata.nc";
        double[][] nodexy = readFileDoubleArray("D:\\sharedfolder\\May2010_mat\\nodexy.dat",M,2," "); // the original mesh nodes
        double[][] uvnode = readFileDoubleArray("D:\\sharedfolder\\May2010_mat\\uvnode.dat",N,2," "); // the centroids of the elements
        int[][] trinodes = readFileIntArray("D:\\sharedfolder\\May2010_mat\\trinodes.dat",N,3," "); // the corners of the elements
        int[][] neighbours = readFileIntArray("D:\\sharedfolder\\May2010_mat\\neighbours.dat",N,3," "); // the neighbouring elements of each element
        int[] allelems = new int[uvnode.length];
        for (int j = 0; j < uvnode.length; j++)
        {
            allelems[j] = j;
        }
        
        // will this work directly? check netcdf file variables
        // TODO make function nc_varget -- same data as mesh.uvnodes
        // TODO actually need to interpolate using x,y and nv (list of corner nodes for each element)

        // setup farms and boundary
        int tstart=1000;
        int tmax=3000;
        int nparts=500;
        double dt=3600;
        
        int firstday=3;
        int lastday=4;
        int nrecords=(lastday-firstday+1)*24;
        // int nparts=20;
        // int tot_psteps=nparts*nrecords;

        dt=3600;

        //fprintf('----------------------------------------------------------------\n');
        int stepsPerStep=100;
        dt=dt/(double)stepsPerStep;

        // default start location
        int s=1;

        // three different options for computing velocity
        // 1 - velocity at nearest element centroid
        // 2 - velocity from nearest 5 centroids (inverse-distance weighted)
        // 3 - velocity from element that particle is located within
        int veltype=3;

        // options for behaviour of particles - depth aspects are set at each
        // timestep by function
        // 1 - passive, stay on surface
        // 2 - passive, stay on bottom (layer 10)
        // 3 - passive, stay in mid layer (layer 5)
        // 4 - vertical swimming: surface for hours 19-6, mid layer (5) hours 7-18
        // 5 - rapid drop (1->10) at hour 6, then gradually move back up
        // MORE...? "homing" ability
        int behaviour=1;

        // should we loops over start locations (=1) or over behaviour (=0)?
        //loopstart=0;

        int viabletime=336;
        int inviabletime=nrecords;

        // load array of start node IDs (as stored by matlab)
        double startlocs[][] = readFileDoubleArray("D:\\sharedfolder\\May2010_mat\\boundary_940pt.dat",940,3," ");

        nparts=startlocs.length;
        int tot_psteps=nparts*nrecords;

        // an array to save the number of "particle-timesteps" in each cell
        int[][] psteps= new int[N][2];
        for (int i = 0; i < N; i++)
        {
            psteps[i][0] = i;
        }
        // array to save source, destination, and transport time for each particle
        int[][] particle_info= new int[nparts][3];
        
        double[] xstart = new double[nparts];
        double[] ystart = new double[nparts];
        int[] startElem = new int[nparts];
        int[] startid = new int[nparts];
        System.out.println("nparts = "+nparts);
        for (int i = 0; i < nparts; i++)
        {
            //System.out.printf("PARTICLE %d\n",i);

            //startid=int16(length(boundary)*rand(1));
            startid[i]=i%startlocs.length;
            //System.out.printf("%d", startlocs.length);
            behaviour=1;

            xstart[i]=startlocs[startid[i]][1];
            ystart[i]=startlocs[startid[i]][2];
            //System.out.printf("start location %d = %d %.4e %.4e\n",i,startid[i],xstart[i],ystart[i]);

            // this boundary location is not actually in the mesh/an element, so set
            // new particle location to centre of nearest element.
            int closest=Particle.nearestCentroid(xstart[i],ystart[i],uvnode);

            xstart[i]=uvnode[closest][0];
            ystart[i]=uvnode[closest][1];
            
            startElem[i]=Particle.whichElement(xstart[i],ystart[i],allelems,nodexy,trinodes);
        }

        // setup particles
        Particle[] particles = new Particle[nparts];
        System.out.println("particles.length = "+particles.length);
        for (int i = 0; i < particles.length; i++)
        {
            particles[i] = new Particle(xstart[i],ystart[i]);
            particle_info[i][0]=startid[i];
            particles[i].setElem(startElem[i]);
        }


        // ------------------- loop 2 = timestep ----------------------------
        //fprintf('Starting time loop\n');


        int count0=0;
        int count1=0;
        int count2=0;
        int count3=0;
        int count4=0;

        double minDistTrav=10000000;
        double maxDistTrav=0;


        for (int day=firstday; day <= lastday; day++)
        {
            System.out.printf("\nday %d - ",day);
            // clear any old data
            //clear FVCOM1
            // load the new data file. this puts variables straight into the
            // workspace
            
            String ufile = "D:\\sharedfolder\\May2010_mat\\u_"+day+".dat";
            String vfile = "D:\\sharedfolder\\May2010_mat\\v_"+day+".dat";
            double[][] u = readFileDoubleArray(ufile,24,N*10," ");
            double[][] v = readFileDoubleArray(vfile,24,N*10," ");

            int lasttime=24;

            for (int tt=0; tt < lasttime; tt++)
            {
                System.out.printf(" %d",tt);
        //             if (mod(tt,100)==0) 
        //                 fprintf('%d ',tt); 
        //             end

                int dep=Particle.setParticleDepth(behaviour,tt);
                //uchunk=squeeze(u(tt,dep,:));
                //vchunk=squeeze(v(tt,dep,:));
                //velocities=[uchunk vchunk];

                for (int st = 0; st < stepsPerStep; st++)
                {
                    //fprintf(',\n');


                    for (int i = 0; i < nparts; i++)
                    {
                        //fprintf('PARTICLE %d\n',i);
                        if (particles[i].getArrived()==false)
                        {
                            //fprintf('%d',elemPart(i));
                            int elemPart = particles[i].getElem();
                            double water_U = u[tt][elemPart*10+dep];
                            double water_V = v[tt][elemPart*10+dep];
                            
                            // save location information "start of timestep"
                            psteps[elemPart][1]++;

                            //[water_U,water_V]=calc_vel1(particles(i),[uchunk vchunk]);
                            // 3. Calculate diffusion    
                            //rand('twister',sum(100*clock)); %resets it to a different state each time.
                            double diff_X = Math.sqrt(0.01*dt);
                            double diff_Y = Math.sqrt(0.01*dt);    //+/- is random so direction doesn't matter
                            double[] behave_uv = particles[i].behaveVelocity(behaviour);
                            // 4. update particle location
                            double newlocx=particles[i].getLocation()[0]+dt*water_U+dt*behave_uv[0]+diff_X*ran.raw(); // simplest possible "Euler"
                            double newlocy=particles[i].getLocation()[1]+dt*water_V+dt*behave_uv[1]+diff_Y*ran.raw();


                            //fprintf('+');
                            int[] elems = new int[1];
                            elems[0] = elemPart;
                            int whereami=Particle.whichElement(newlocx,newlocy,elems,nodexy,trinodes);
                            count0=count0+1;
                            if (whereami==0)
                            {
                                int[] elems0 = neighbours[elemPart];
                                count1++;
                                whereami=Particle.whichElement(newlocx,newlocy,elems0,nodexy,trinodes);
                                // if fails, look in nearest 10 (id numerical)
                                if (whereami==0)
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
                                    if (whereami==0)
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
                                        if (whereami==0)
                                        {
                                            //fprintf('c');
                                            count4=count4+1;
                                            whereami=Particle.whichElement(newlocx,newlocy,allelems,nodexy,trinodes);
                                            //elemPart(i)=nearest_centroid(particles(i).x,particles(i).y,centroids);
                                        }
                                    }
                                }
                            }


                            if (whereami != 0)
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
                            }
                            if (whereami == 0)
                            {
                                int closest=Particle.nearestCentroid(particles[i].getLocation()[0],particles[i].getLocation()[1],uvnode);
                                //fprintf('x%d',closest);
                                particles[i].setLocation(uvnode[closest][0],uvnode[closest][1]);
                                particles[i].setElem(closest);
                            }

                            // set particle to become able to settle after a predefined time
                            if ((day-1)*24+tt>viabletime)
                            {
                                particles[i].setViable(true);
                            }
                            // if able to settle, is it close to a possible settlement
                            // location?
                            if (particles[i].getViable()==true)
                            {
                                for (int loc = 0; loc < startlocs.length; loc++)
                                {
                                    double dist = Math.sqrt((particles[i].getLocation()[0]-startlocs[loc][0])*(particles[i].getLocation()[0]-startlocs[loc][0])+
                                            (particles[i].getLocation()[1]-startlocs[loc][1])*(particles[i].getLocation()[1]-startlocs[loc][1]));
                                    if (dist < 250)
                                    {
                                        particles[i].setArrived(true);
                                        particle_info[i][2] = loc;
                                        particle_info[i][3] = ((day-1)*24+tt);
                                        break;
                                    }
                                }    
                            }

                        }
                    }
                    // end of particle loop
                }
            }

        }  
        System.out.printf("element search counts: %d %d %d %d %d\n",count0,count1,count2,count3,count4);
        System.out.printf("transport distances: min = %.4e, max = %.4e\n", minDistTrav, maxDistTrav);

        writeIntegerArrayToFile(psteps,"psteps.out");
        writeIntegerArrayToFile(particle_info,"psteps.out");
       
        //profile viewer
        //p = profile("info");
        //profsave(p,"profile_results")
        // scale psteps
//        tot_psteps=tmax*nparts;
//        double[] psteps2 = new double[N];
//        for (int i = 0; i < N; i++)
//        {
//            psteps2[i]=(double)psteps[i][1]/tot_psteps;
//        }       
    }

    public static double[][] setupStartLocations()
    {
        int nsites = 9;
        double startlocs[][]= new double[nsites][2];
        
        startlocs[0][0]=357420; startlocs[0][1]=6217200; 
        startlocs[1][0]=361834; startlocs[1][1]=6223063;
        startlocs[2][0]=353078; startlocs[2][1]=6206339;
        startlocs[3][0]=354246; startlocs[3][1]=6194759; 
        startlocs[4][0]=352745; startlocs[4][1]=6201735;
        startlocs[5][0]=348880; startlocs[5][1]=6199380;
        startlocs[6][0]=354969; startlocs[6][1]=6193169;
        startlocs[7][0]=348606; startlocs[7][1]=6204475;
        startlocs[8][0]=352401; startlocs[8][1]=6190933;
        // non-fishfarms
//        double startlocs[][]=[354500 6188000; 
//            355200 6192000;
//            350000 6197000;
//            352800 6196000; 
//            348000 6209000;
//            354000 6204000;
//            357000 6213000;
//            360000 6219000;
//            370000 6227000];
        // "hypothetical" fishfarm used for SAMS newsletter
        //fishfarms(1,1)=351000;
        //fishfarms(1,2)=6195000;
        //"2" is Minard
        //"7" is Portavadie
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
    
    public static double[][] readFileDoubleArray(String filename, int rows, int cols, String sep)
    {
        double[][] myDouble = new double[rows][cols];
        int x=0, y=0;

        try
        {
            BufferedReader in = new BufferedReader(new FileReader(filename));	//reading files in specified directory
 
            String line;
            while ((line = in.readLine()) != null)	//file reading
            {
                y=0;
                if (x >= rows)
                {
                    System.out.println(filename+" has more rows than expected.");
                    break;
                }
                String[] values = line.split(" ");
                for (String str : values)
                {
                    if (y >= cols)
                    {
                        System.out.println(filename+" has more columns than expected.");
                        break;
                    }
                    double str_double = Double.parseDouble(str);
                    myDouble[x][y]=str_double;
                    //System.out.print(myDouble[x][y] + " ");
                    y++;

                }
                //System.out.println("");
                x++;
            }
            in.close();
        } catch( IOException ioException ) {}
        System.out.printf("Created %dx%d array from file: %s\n",myDouble.length,myDouble[0].length,filename);
        return myDouble;
    }
    
    public static int[][] readFileIntArray(String filename, int rows, int cols, String sep)
    {
        int[][] myInt = new int[rows][cols];
        int x=0, y=0;

        try
        {
            BufferedReader in = new BufferedReader(new FileReader(filename));	//reading files in specified directory
 
            String line;
            while ((line = in.readLine()) != null)	//file reading
            {
                y=0;
                if (x >= rows)
                {
                    System.out.println(filename+" has more rows than expected.");
                    break;
                }
                String[] values = line.split(sep);
                for (String str : values)
                {
                    if (y >= cols)
                    {
                        System.out.println(filename+" has more columns than expected.");
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
        } catch( IOException ioException ) {}
        System.out.printf("Created %dx%d array from file: %s\n",myInt.length,myInt[0].length,filename);
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
    
    public double nc_varget(String filename, String variable)
    {
        double out=0;
        return out;
    }

}
