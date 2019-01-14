/*
 * TODO LIST:
 * - Vertical migration/behaviour
 * - Parallelisation?
 *      + Threads using fork/join 
 *        http://www.oracle.com/technetwork/articles/java/fork-join-422606.html
 *        http://tutorials.jenkov.com/java-util-concurrent/java-fork-and-join-forkjoinpool.html
 *        (16 cores per cluster node)
 *      + MPI e.g. MPJ?
 */
package particle_track;

import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

//import java.util.Collection;
//import java.util.List;
//import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.CompletionService;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
//import java.util.concurrent.ThreadLocalRandom;
//import java.lang.InterruptedException;

//import edu.cornell.lassp.houle.RngPack.*;
//import static particle_track.IOUtils.countLines;

 
import java.io.IOException;
//import java.nio.file.Files;
//import java.nio.file.Path;
//import java.nio.file.Paths;
//import java.util.stream.Stream;
import org.apache.commons.io.FileUtils;
//import org.apache.commons.io.filefilter.TrueFileFilter;
import org.apache.commons.io.filefilter.WildcardFileFilter;
 



/**
 *
 * @author tomdude
 */
public class Particle_track {
   
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException {
        // TODO code application logic here

        System.out.println("Starting particle tracking program\n");
        Date date = new Date();
        // display time and date using toString()
        System.out.println(date.toString());

        long heapMaxSize = Runtime.getRuntime().maxMemory();
        System.out.println("Max heap " + heapMaxSize);

        //System.out.println(new Date().toString());
        long startTime = System.currentTimeMillis();

        //RanMT ran = new RanMT(System.currentTimeMillis());
        System.out.println("Reading in data\n");

        System.out.println(System.getProperty("user.dir"));

        //RunProperties rp = new RunProperties("model_setup.properties");
        RunProperties rp = new RunProperties(args[0]); // first (and only?) cmd line arg is properties filename e.g. model_setup.properties
        // Use this instead of previous to create runProps from CMD line args
        //RunProperties runProps = new RunProperties(args); 

        int[] startDate = dateIntParse(rp.start_ymd);
        ISO_datestr currentIsoDate = new ISO_datestr(startDate[0], startDate[1], startDate[2]);
        int[] endDate = dateIntParse(rp.end_ymd);
        ISO_datestr endIsoDate = new ISO_datestr(endDate[0], endDate[1], endDate[2]);

        int numberOfDays = endIsoDate.getDateNum() - currentIsoDate.getDateNum() + 1;

        // Print all main arguments
        System.out.printf("-----------------------------------------------------------\n");
        System.out.printf("Location           = %s\n", rp.location);
        System.out.printf("Habitat            = %s\n", rp.habitat);
        System.out.printf("N_parts/site       = %d\n", rp.nparts);
        System.out.printf("hydromod dt (s)    = %f\n", rp.dt);
        System.out.printf("hydromod rec/file  = %d, %d\n", rp.recordsPerFile1, rp.recordsPerFile2);
        System.out.printf("stepsperstep       = %d\n", rp.stepsPerStep);
        System.out.printf("firstfile          = %d\n", rp.start_ymd);
        System.out.printf("lastfile           = %d\n", rp.end_ymd);
        System.out.printf("Simulated dur. (d) = %f\n", (double) numberOfDays);
        System.out.printf("Simulated dur. (s) = %f\n", (double) numberOfDays * 86400);
        //System.out.printf("Simulated dur. (s) = %f\n",rp.dt*rp.recordsPerFile*(rp.lastday-rp.firstday+1));
        System.out.printf("RK4                = %s\n", rp.rk4);
        System.out.printf("Vertical behaviour = %d\n", rp.behaviour);
        System.out.printf("Viable time (h)    = %f\n", rp.viabletime);
        System.out.printf("Viable time (d)    = %f\n", rp.viabletime / 24.0);
        System.out.printf("Threshold distance = %d\n", rp.thresh);
        System.out.printf("Diffusion D_h      = %f (diffusion: %s)\n", rp.D_h, rp.diffusion);
        System.out.printf("-----------------------------------------------------------\n");

        // --------------------------------------------------------------------------------------
        // File reading and domain configuration
        // --------------------------------------------------------------------------------------
        // Set the directories required to run the model
        //String[] dirList = IOUtils.setDataDirectories(location, cluster);
//        String basedir = dirList[0];
//        String sitedir = dirList[2];
//        String datadir = dirList[3];
//        String datadir2 = dirList[4];

        
        List<Mesh> meshes = new ArrayList<>();
        meshes.add(new Mesh(rp.mesh1,rp.mesh1Type));
        if (rp.mesh2.equals("") != true)
        {
            meshes.add(new Mesh(rp.mesh2,rp.mesh2Type));
        }


        int[] allelems = IntStream.rangeClosed(0, meshes.get(0).getUvnode()[0].length-1).toArray();
        
        double subStepDt = rp.dt / (double) rp.stepsPerStep; // number of seconds per substep
        double dev_perstep = Math.pow(0.1, subStepDt);
        System.out.println("Particle subStepDt = " + subStepDt + " dev_perstep = " + dev_perstep);
        System.out.println("behaviour = " + rp.behaviour);

        //Need
        
        
        // --------------------------------------------------------------------------------------
        // Creating initial particle array
        // --------------------------------------------------------------------------------------
        // load array of start node IDs (as stored by matlab)
        double startlocs[][] = new double[10][3];       
        double endlocs[][] = new double[10][3];
        //double open_BC_locs[][] = new double[1][3];

        //startlocs = IOUtils.setupStartLocs(rp.sitefile, rp.sitedir, true);
        //startlocs = IOUtils.setupStartLocs(rp.location,rp.habitat,rp.basedir);
        //endlocs = IOUtils.setupEndLocs(rp.habitat, rp.sitedir, startlocs, rp.endlimit);
//        open_BC_locs = IOUtils.setupOpenBCLocs(rp.location, rp.datadir2);
        
        // A new way of creating habitat sites, allowing use of more information
        List<HabitatSite> habitat = new ArrayList<>();
        System.out.println("Creating start sites");
        habitat = IOUtils.createHabitatSites(rp.sitefile, rp.sitedir, 5, true);
        for (HabitatSite site : habitat)
        {
            System.out.println(site.toString());
        }
        // Need a list of end sites - have just used the same list for now
        List<HabitatSite> habitatEnd = new ArrayList<>();
        System.out.println("Creating end sites");
        habitatEnd = IOUtils.createHabitatSites(rp.sitefile, rp.sitedir, 5, true);
        

        int nparts_per_site = rp.nparts;
        int nTracksSavedPerSite = Math.min(1, nparts_per_site);
        int nparts = rp.nparts * startlocs.length;

        for (int i = 0; i < startlocs.length; i++) {
            startlocs[i][0]--;
            //System.out.println(startlocs[i][0]+" "+startlocs[i][1]+" "+startlocs[i][2]);
        }

        // --------------------------------------------------------------------------------------
        // Setup particles
        // --------------------------------------------------------------------------------------
        List<Particle> particles = new ArrayList<>(nparts);
        int numParticlesCreated = 0; // Counter to keep track of how many particles have been created
        boolean allowRelease = true; // boolean to be switched after a single release event
        

        // --------------------------------------------------------------------------------------
        // Particle data arrays for model output
        // --------------------------------------------------------------------------------------
        // an array to save the number of "particle-timesteps" in each cell
        // default case (non-split): two columns
        int pstepCols = 2; 
        // Alternatively, make a column for each source site
        if (rp.splitPsteps == true){
            pstepCols = startlocs.length + 1;            
        }
        
//        double[][] pstepsImmature = new double[rp.N][pstepCols];
//        double[][] pstepsMature = new double[rp.N][pstepCols];
//        for (int i = 0; i < rp.N; i++) {
//            pstepsImmature[i][0] = i;
//            pstepsMature[i][0] = i;
//        }

        // --------------------------------------------------------------------------------------
        // Set up times at which to print particle locations to file 
        // --------------------------------------------------------------------------------------
        int simLengthHours = numberOfDays * 24;
        System.out.println("simLengthHours " + simLengthHours);
        

        // --------------------------------------------------------------------------------------
        // Final setup bits
        // --------------------------------------------------------------------------------------
        System.out.println("Starting time loop");

        int[] searchCounts = new int[5];

        double minMaxDistTrav[] = new double[2];
        minMaxDistTrav[0] = 10000000;
        minMaxDistTrav[1] = 0;

        int stepcount = 0;
        int calcCount = 0;
        double time = 0;

        int printCount = 0;

        int[] freeViableSettleExit = new int[4];

        int numberOfExecutorThreads = Runtime.getRuntime().availableProcessors();
        if (rp.parallel == false) {
            numberOfExecutorThreads = 1;
        }
        //numberOfExecutorThreads = 1;
        System.out.println("Number of executor threads = " + numberOfExecutorThreads);
        ExecutorService executorService = Executors.newFixedThreadPool(numberOfExecutorThreads);
        CompletionService<List<Particle>> executorCompletionService = new ExecutorCompletionService<List<Particle>>(executorService);                           
        
        //final Collection<Callable<List<Particle>>> callables = new ArrayList<>();
        final Collection<Callable<List<Particle>>> callables = new ArrayList<Callable<List<Particle>>>();

        String locationHeader = "hour ID startDate age startLocation x y elem status density";
        String arrivalHeader = "ID startDate startTime startLocation endDate endTime endLocation age density";
        
        IOUtils.printFileHeader(arrivalHeader,"arrivals.out");
        
        HydroField hydroField1,hydroField2;
        // Set an initial values for the field. Subsequent value setting will produce arrays that have one more element in the time dimension
        List<File> files = (List<File>) FileUtils.listFiles(
                new File(rp.datadir+rp.datadirPrefix+currentIsoDate.getYear()+System.getProperty("file.separator")),
                new WildcardFileFilter("minch2_"+currentIsoDate.getYear()+String.format("%02d",currentIsoDate.getMonth())+String.format("%02d",currentIsoDate.getDay())+"*.nc"), 
                null);    
        String[] varNames = {"u","v","salinity","temp","zeta"};
        hydroField1 = new HydroField(files.get(0).getCanonicalPath(),varNames);
        
        try {
            // --------------------------------------------------------------------------------------
            // Start time loop
            // --------------------------------------------------------------------------------------
            long currTime = System.currentTimeMillis();
            //for (int fnum = rp.firstday; fnum <= rp.lastday; fnum++)
            for (int fnum = 0; fnum < numberOfDays; fnum++) {

                String today = currentIsoDate.getDateStr();
                System.out.println(today);
                IOUtils.printFileHeader(locationHeader,"locations_" + today + ".out");
                IOUtils.printFileHeader(arrivalHeader,"arrivals_" + today + ".out");
                
                long splitTime = System.currentTimeMillis();
                System.out.printf("\n------ Day %d (%s) - Stepcount %d (%f hrs) ------ \n",
                        fnum + 1, currentIsoDate.getDateStr(), stepcount, time);
                System.out.println("Elapsed time (s) = " + (splitTime - startTime) / 1000.0);
                System.out.println("Last 24hr time (s) = " + (splitTime - currTime) / 1000.0);
                currTime = System.currentTimeMillis();
                
                // set an initial tide state
                String tideState = "flood";

                // COUNT the number of particles in different states
                freeViableSettleExit = particleCounts(particles);
                System.out.println("Free particles    = " + freeViableSettleExit[0]);
                System.out.println("Viable particles  = " + freeViableSettleExit[1]);
                System.out.println("Arrival count     = " + freeViableSettleExit[2]);
                System.out.println("Boundary exits    = " + freeViableSettleExit[3]);

                // default, run loop forwards
                // ---- LOOP OVER ENTRIES IN THE HYDRO OUTPUT ------------------------
                for (int tt = 0; tt < 24; tt++) {
                    // alternatively, run loop backwards
                    //for (int tt = lasttime; tt >= firsttime; tt--)

                    System.out.printf("--------- HOUR %d ----------\n",tt);
                    // Calculate current time of the day (complete hours elapsed since midnight)
                    int currentHour = tt;
                    //System.out.printf("%d \n", tt + 1);
                    
                    // Read the PRIMARY hydrodynamic output file - ASSUME THAT THIS IS ONE FILE PER DAY
                    if (tt%rp.recordsPerFile1 == 0){
                        System.out.println("Reading file "+tt);
                        // Dima file naming format: minch2_20171229_0003.nc
                        List<File> files1 = (List<File>) FileUtils.listFiles(
                                new File(rp.datadir+rp.datadirPrefix+currentIsoDate.getYear()+System.getProperty("file.separator")),
                                new WildcardFileFilter("minch2_"+currentIsoDate.getYear()+String.format("%02d",currentIsoDate.getMonth())+String.format("%02d",currentIsoDate.getDay())+"*.nc"), 
                                null);    
                        String[] varNames1 = {"u","v","salinity","temp","zeta"};
                        //hydroField1 = new HydroField(files1.get(0).getCanonicalPath(),varNames1);
                        
                        ISO_datestr tomorrow = ISO_datestr.getTomorrow(currentIsoDate);
                        
                        // Check dates
//                        System.out.println("Today: "+currentIsoDate.getDateStr());
//                        System.out.println("Tomorrow: "+tomorrow.getDateStr());

                        List<File> files2 = (List<File>) FileUtils.listFiles(
                                new File(rp.datadir+rp.datadirPrefix+tomorrow.getYear()+System.getProperty("file.separator")),
                                new WildcardFileFilter("minch2_"+tomorrow.getYear()+String.format("%02d",tomorrow.getMonth())+String.format("%02d",tomorrow.getDay())+"*.nc"), 
                                null);                           
                        // Read both files and combine
                        hydroField1 = new HydroField(files1.get(0).getCanonicalPath(),files2.get(0).getCanonicalPath(),varNames1);
                    }   
                    
                    // Create new particles, if releases are scheduled hourly, or if release is scheduled for this
                    // exact hour
                    if (rp.releaseScenario==1 || (rp.releaseScenario==0 && time>rp.releaseTime && allowRelease==true))
                    {
                        System.out.printf("Release attempt: releaseScenario %d, releaseTime %f, allowRelease %s newParticlesCreated %d \n",
                            rp.releaseScenario,time,allowRelease,numParticlesCreated);
                        //System.out.printf("releaseScenario==1, releasing hourly (hour = %d)\n",currentHour);
                        List<Particle> newParts = createNewParticles(habitat,meshes.get(0),
                                rp,currentIsoDate,currentHour,numParticlesCreated);
                        particles.addAll(newParts);
                        numParticlesCreated = numParticlesCreated+(rp.nparts*habitat.size());
                        // If only one release to be made, prevent further releases
                        if (rp.releaseScenario==0)
                        {
                            allowRelease = false;
                        }
                    }

                    // ---- INTERPOLATE BETWEEN ENTRIES IN THE HYDRO OUTPUT ------------------------
                    for (int st = 0; st < rp.stepsPerStep; st++) {

                        // Update the element count arrays
                        //pstepUpdater(particles, rp, pstepsMature, pstepsImmature, subStepDt);

                        //System.out.print(",");
                        //System.out.println("nfreeparts = "+nfreeparts);
                        // MOVE the particles
                        if (rp.parallel == true) {
                            int particlesSize = particles.size();
                            int listStep = particlesSize / numberOfExecutorThreads;
                            for (int i = 0; i < numberOfExecutorThreads; i++) {
                                List<Particle> subList;
                                if(i==numberOfExecutorThreads-1){
                                    // Note: ArrayList.subList(a,b) is inclusive of a but exclusive of b => 
                                    subList = particles.subList(i * listStep, particlesSize);
                                    //System.out.println(listStep+" "+i+" "+(i*listStep)+" "+(particlesSize-1));
                                }else{
                                    subList = particles.subList(i * listStep, (i + 1) * listStep);
                                    //System.out.println(listStep+" "+i+" "+(i*listStep)+" "+((i + 1) * listStep - 1));
                                }
//                                callables.add(new ParallelParticleMover(subList, time, tt, st, subStepDt, rp,
//                                        hydroField1.getU(), hydroField1.getV(), mesh1.getNeighbours(), mesh1.getUvnode(), mesh1.getNodexy(), 
//                                        mesh1.getTrinodes(), allelems, mesh1.getDepthUvnode(), mesh1.getSiglay(),
//                                        habitatEnd, mesh1.getOpenBoundaryNodes(),
//                                        searchCounts,
//                                        minMaxDistTrav));
                                callables.add(new ParallelParticleMover(subList, time, tt, st, subStepDt, rp,
                                        hydroField1.getU(), hydroField1.getV(), meshes.get(0).getNeighbours(), meshes.get(0).getUvnode(), meshes.get(0).getNodexy(), 
                                        meshes.get(0).getTrinodes(), allelems, meshes.get(0).getDepthUvnode(), meshes.get(0).getSiglay(),
                                        habitatEnd, meshes.get(0).getOpenBoundaryNodes(),
                                        searchCounts,
                                        minMaxDistTrav));

                            }
                            for (Callable<List<Particle>> callable : callables) {
                                executorCompletionService.submit(callable);
                            }
                            for (Callable<List<Particle>> callable : callables) {
                                executorCompletionService.take().get();
                            }
                            callables.clear();
                        } else {
                            // Normal serial loop
                            for (Particle part : particles) {
                                // Can just use the move method that was shipped out to the ParallelParticleMover class
//                                ParallelParticleMover.move(part, time, tt, st, subStepDt, rp,
//                                        hydroField1.getU(), hydroField1.getV(), mesh1.getNeighbours(), mesh1.getUvnode(), mesh1.getNodexy(), mesh1.getTrinodes(), 
//                                        allelems, mesh1.getDepthUvnode(), mesh1.getSiglay(),
//                                        habitatEnd, mesh1.getOpenBoundaryNodes(),
//                                        searchCounts,
//                                        minMaxDistTrav);
                                ParallelParticleMover.move(part, time, tt, st, subStepDt, rp,
                                        hydroField1.getU(), hydroField1.getV(), meshes.get(0).getNeighbours(), meshes.get(0).getUvnode(), meshes.get(0).getNodexy(), meshes.get(0).getTrinodes(), 
                                        allelems, meshes.get(0).getDepthUvnode(), meshes.get(0).getSiglay(),
                                        habitatEnd, meshes.get(0).getOpenBoundaryNodes(),
                                        searchCounts,
                                        minMaxDistTrav);

                            }
                        }

                        // --------------- End of particle loop ---------------------
                        time += subStepDt / 3600.0;

                        // end of particle loop
                        calcCount++;
                    }

                    // New output: print ALL current particle location to a separate file, once each hour

                    //System.out.println("Print particle locations to file " + today + " " + currentHour);
                    //IOUtils.particleLocsToFile_full(particles, "locations_" + today + "_" + currentHour + ".out", true);

                    IOUtils.particleLocsToFile_full(particles,currentHour,"locations_" + today + ".out",true);
                                        
                    // It's the end of an hour, so if particles are allowed to infect more than once, reactivate them
                    for (Particle part : particles) {
                        if (part.getSettledThisHour()==true) // previously had clause oldOutput==false here
                        {
                            IOUtils.arrivalToFile(part, currentIsoDate, currentHour, "arrivals_" + today + ".out", true);
                            part.setSettledThisHour(false);
                        }
                    }

                    // Clean up "dead" (666) and "exited" (66) particles
                    List<Particle> particlesToRemove = new ArrayList<>(0);
                    for (Particle part : particles)
                    {
                        if (part.getStatus()==666 || part.getStatus()==66)
                        {
                            //System.out.printf("Removing particle %d, status %d\n",part.getID(),part.getStatus());
                            particlesToRemove.add(part);
                        }
                    }
                    particles.removeAll(particlesToRemove);
                    
                    System.out.println("Number of particles = "+particles.size());

                    printCount++;                    
                    stepcount++;
                }
                System.out.printf("\n");
                
                currentIsoDate.addDay();
            }
            System.out.printf("\nelement search counts: %d %d %d %d %d\n", searchCounts[0], searchCounts[1], searchCounts[2], searchCounts[3], searchCounts[4]);
            System.out.printf("transport distances: min = %.4e, max = %.4e\n", minMaxDistTrav[0], minMaxDistTrav[1]);
           
            executorService.shutdownNow();
        } finally {
            executorService.shutdownNow();
        }

        long endTime = System.currentTimeMillis();
        System.out.println("Elapsed time = " + (endTime - startTime) / 1000.0);
    }
    
    /**
     * Method to create new particles. These must be appended to the existing list
     * 
     * @param habitat
     * @param mesh
     * @param rp
     * @param currentDate
     * @param currentTime
     * @param numParticlesCreated
     * @return List of the new particles to be appended to existing list
     */
    public static List<Particle> createNewParticles(List<HabitatSite> habitat, Mesh mesh, 
            RunProperties rp, ISO_datestr currentDate, int currentTime, int numParticlesCreated)
    {
        //System.out.printf("In createNewParticles: nparts %d startlocsSize %d\n",rp.nparts,startlocs.length);
        List<Particle> newParts = new ArrayList<>(rp.nparts*habitat.size());
        for (int i = 0; i < rp.nparts*habitat.size(); i++)
            {
                int startid = i % habitat.size();
                double xstart = habitat.get(startid).getLocation()[0];
                double ystart = habitat.get(startid).getLocation()[1];
                
                //System.out.println(xstart+" "+ystart);
                
                //double xstart = startlocs[startid][1];
                //double ystart = startlocs[startid][2];
                
                // If start location is a boundary location it is not actually in the mesh/an element, so set
                // new particle location to centre of nearest element.
                int closest = Particle.nearestCentroid(xstart, ystart, mesh.getUvnode());
                int startElem = Particle.whichElement(new double[]{xstart,ystart}, 
                        IntStream.rangeClosed(0, mesh.getUvnode()[0].length-1).toArray(), 
                        mesh.getNodexy(), mesh.getTrinodes());
                if (startElem < 0) {
                    xstart = mesh.getUvnode()[0][closest];
                    ystart = mesh.getUvnode()[1][closest];
                    startElem = closest;
                }
            
                Particle p = new Particle(xstart, ystart, rp.startDepth, habitat.get(startid).getID(), numParticlesCreated+i, 
                        rp.mortalityRate, currentDate, currentTime, rp.coordRef);
                p.setElem(startElem);
                // No longer set releaseScenario for individual particles since they are created at the moment of release
                //p.setReleaseScenario(rp.releaseScenario, startlocs);
//                if (startlocs[startid].length > 4 && rp.setDepth == true) {
//                    p.setZ(startlocs[startid][4]);
//                }
                if (rp.setDepth == true) {
                    p.setZ(habitat.get(startid).getDepth());
                }
                newParts.add(p);
                //System.out.println(p.toString());
            }
        return newParts;
    }

    /**
     * Count the number of particles in different states (free, viable, settled,
     * exited domain)
     *
     * @param parts
     * @return
     */
    public static int[] particleCounts(List<Particle> parts) {
        int freeViableSettleExit[] = new int[4];
        // Add count 1 for each particle that satisfies this list of conditions
        // Lines below are equivalent to:
        //if (p.getFree()) {
        //    freeViableSettleExit[0] += 1;
        //} 
        for (Particle p : parts) {
            freeViableSettleExit[0] += p.getFree() ? 1 : 0;
            freeViableSettleExit[1] += p.getViable() ? 1 : 0;
            freeViableSettleExit[2] += p.getArrived() ? 1 : 0;
            freeViableSettleExit[3] += p.getBoundaryExit() ? 1 : 0;
        }
        return freeViableSettleExit;
    }
    
    // calculate a connectivity matrix detailing the 
    public static double[][] connectFromParticleArrivals(List<Particle> particles, int nStartLocs, int npartsPerSite)
    {
        double[][] connectMatrix = new double[nStartLocs][nStartLocs];
        for (Particle part : particles)
        {
            for (Arrival arrival: part.getArrivals())
            {
                connectMatrix[arrival.getSourceLocation()][arrival.getArrivalLocation()] += arrival.getArrivalDensity()/npartsPerSite;
            }
        }
        return connectMatrix;
    }

//    /**
//     * Make additions to the element presence counts (PSTEPS)
//     *
//     * @param particles
//     * @param rp
//     * @param pstepsMature
//     * @param pstepsImmature
//     * @param subStepDt
//     */
//    public static void pstepUpdater(List<Particle> particles, RunProperties rp,
//            double[][] pstepsMature, double[][] pstepsImmature, double subStepDt) {
//        for (Particle p : particles) {
//            double d = 1;
//            if (rp.pstepsIncMort == true) {
//                d = p.getDensity();
//            }
//            //System.out.println("density = "+d+" mortRate = "+p.getMortRate());
//            int elemPart = p.getElem();
//            // psteps arrays are updated by lots of threads
//            if (p.getViable() == true) {
//                if (rp.splitPsteps == false) {
//                    pstepsMature[elemPart][1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsMature[elemPart][p.getStartID() + 1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                }
//            } else if (p.getFree() == true) {
//                //System.out.println("Printing to pstepsImmature");
//                if (rp.splitPsteps == false) {
//                    pstepsImmature[elemPart][1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsImmature[elemPart][p.getStartID() + 1] += d * (subStepDt / 3600);//*1.0/rp.stepsPerStep;
//                }
//            }
//        }
//    }
    
    /**
     * Take a snapshot of the number of mature particles in each cell
     * @param particles
     * @param rp
     * @param nSourceSites
     * @return 
     */
//    public static double[][] pstepMatureSnapshot(List<Particle> particles, RunProperties rp,
//            int nSourceSites) {   
//        int pstepCols = 2; 
//        // Alternatively, make a column for each source site
//        if (rp.splitPsteps == true){
//            pstepCols = nSourceSites + 1;            
//        }
//        double[][] pstepsInstMature = new double[rp.N][pstepCols];
//        for (int i = 0; i < rp.N; i++) {
//            pstepsInstMature[i][0] = i;
//        }
//        
//        for (Particle p : particles) {
//            if (p.getViable() == true) {
//                double d = 1;
//                if (rp.pstepsIncMort == true) {
//                    d = p.getDensity();
//                }
//                //System.out.println("density = "+d+" mortRate = "+p.getMortRate());
//                int elemPart = p.getElem();
//                if (rp.splitPsteps == false) {
//                    pstepsInstMature[elemPart][1] += d;//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsInstMature[elemPart][p.getStartID() + 1] += d;//*1.0/rp.stepsPerStep;
//                }
//            }
//        }
//        return pstepsInstMature;
//    } 
    
//    /**
//     * Take a snapshot of the number of immature particles in each cell
//     * @param particles
//     * @param rp
//     * @param nSourceSites
//     * @return 
//     */
//    public static double[][] pstepImmatureSnapshot(List<Particle> particles, RunProperties rp,
//            int nSourceSites) {   
//        int pstepCols = 2; 
//        // Alternatively, make a column for each source site
//        if (rp.splitPsteps == true){
//            pstepCols = nSourceSites + 1;
//        }
//        double[][] pstepsInstImmature = new double[rp.N][pstepCols];
//        for (int i = 0; i < rp.N; i++) {
//            pstepsInstImmature[i][0] = i;
//        }
//        
//        for (Particle p : particles) {
//            if (p.getViable() == false && p.getFree() == true) {
//                double d = 1;
//                if (rp.pstepsIncMort == true) {
//                    d = p.getDensity();
//                }
//                //System.out.println("density = "+d+" mortRate = "+p.getMortRate());
//                int elemPart = p.getElem();
//                //System.out.println("Printing to pstepsImmature");
//                if (rp.splitPsteps == false) {
//                    pstepsInstImmature[elemPart][1] += d;//*1.0/rp.stepsPerStep;
//                } else {
//                    pstepsInstImmature[elemPart][p.getStartID() + 1] += d;//*1.0/rp.stepsPerStep;
//                }
//            }
//        }
//        return pstepsInstImmature;
//    }
    

    /**
     * work out the date from an integer in format YYYYMMDD - Mike Bedington
     * method
     *
     * @param ymd
     * @return
     */
    public static int[] dateIntParse(int ymd) {
        double start_ymd_mod = (double) (ymd);
        int startYear = (int) Math.floor(start_ymd_mod / 10000);
        start_ymd_mod = start_ymd_mod - startYear * 10000;
        int startMonth = (int) Math.floor(start_ymd_mod / 100);
        start_ymd_mod = start_ymd_mod - startMonth * 100;
        int startDay = (int) start_ymd_mod;

        int[] output = new int[]{startDay, startMonth, startYear};
        return output;
    }

    public static void memTest() {
        long heapSize = Runtime.getRuntime().totalMemory();
        System.out.println("Total heap memory " + heapSize);
        long heapFreeSize = Runtime.getRuntime().freeMemory();
        System.out.println("Free heap memory " + heapFreeSize);

    }

    public void setupOutput() {
    }

    public void writeOutput() {
    }
}
