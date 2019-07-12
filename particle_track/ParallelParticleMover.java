/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.util.Collection;
import java.util.List;
import java.util.ArrayList;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import java.lang.InterruptedException;
import java.util.concurrent.ThreadLocalRandom;
//import static particle_track.Particle_track.move;

/**
 *
 * @author sa01ta
 */
public class ParallelParticleMover implements Callable<List<Particle>> {
    
    private final List<Particle> particles;
    private final double time;
    private final int tt;
    private final int st;
    private final double subStepDt;
    private final RunProperties rp;
    
//    private final float[][][] u; 
//    private final float[][][] v;
//    private final int[][] neighbours;
//    private final float[][] uvnode;
//    private final float[][] nodexy;
//    private final int[][] trinodes;
    private final int[] allelems;
//    private final float[] depthUvnode;
//    private final float[] sigvec2;
//    private final int[] open_BC_locs;
    
    private final List<Mesh> meshes;
    private final List<HydroField> hydroFields;
    
    private final List<HabitatSite> habitatEnd;
    private final int[] searchCounts;
    private final double[] minMaxDistTrav;
    
    public ParallelParticleMover(List<Particle> particles, double time, int tt, int st, double subStepDt,
            RunProperties rp,
            List<Mesh> meshes,
            List<HydroField> hydroFields,
            List<HabitatSite> habitatEnd, 
            int[] allelems,
            int[] searchCounts,
            double[] minMaxDistTrav)
    {
        this.particles=particles;
        this.time=time;
        this.tt=tt;
        this.st=st;
        this.subStepDt=subStepDt;
        this.rp=rp;
        this.meshes=meshes;
        this.hydroFields=hydroFields;
        this.allelems=allelems;
        this.habitatEnd=habitatEnd;
        this.searchCounts=searchCounts;
        this.minMaxDistTrav=minMaxDistTrav;
    }
    
    @Override
    public ArrayList<Particle> call() throws Exception {
        for (Particle part : particles) {
            //System.out.println("Time "+time+" - Moving particle "+part.getID());
//            move(part, time, tt, st, subStepDt, rp, 
//                                    u, v, neighbours, uvnode, nodexy, trinodes, allelems, depthUvnode, sigvec2, 
//                                    habitatEnd, open_BC_locs,
//                                    searchCounts, 
//                                    minMaxDistTrav);
            move(part, time, tt, st, subStepDt, rp, meshes, hydroFields, habitatEnd, allelems, searchCounts, 
                                    minMaxDistTrav);
        }
        return new ArrayList<Particle>();
    }
    
    /**
     * Method to do all the actual particle movement stuff. This can be used either by "call" above to 
     * do particle movements split over the cores, or can just be used directly in a serial loop
     * (i.e. when parallel==false)
     * 
     * @param part
     * @param time
     * @param tt
     * @param st
     * @param subStepDt
     * @param rp
     * @param u
     * @param v
     * @param neighbours
     * @param uvnode
     * @param nodexy
     * @param trinodes
     * @param allelems
     * @param depthUvnode
     * @param siglay
     * @param habitatEnd
     * @param open_BC_locs
     * @param searchCounts
     * @param particle_info
     * @param settle_density
     * @param minMaxDistTrav 
     */
//    (Particle part, double time, int tt, int st, double subStepDt,
//            RunProperties rp, 
//            float u[][][], float v[][][],
//            int[][] neighbours, float[][] uvnode,  float[][] nodexy, 
//            int[][] trinodes, int[] allelems,
//            float[] depthUvnode, float[] siglay,
//            List<HabitatSite> habitatEnd, int[] open_BC_locs,
//            int[] searchCounts,
//            double[] minMaxDistTrav)
    public static void move(Particle part, double time, int tt, int st, double subStepDt,
            RunProperties rp, 
            List<Mesh> meshes,
            List<HydroField> hydroFields,
            List<HabitatSite> habitatEnd,
            int[] allelems,
            int[] searchCounts,
            double[] minMaxDistTrav)
    {
        //System.out.println("--- Moving particle "+part.getID()+" ---");
        //System.out.println("Particle location at start of ParallelParticleMover.move: "+part.printLocation());
        
        
        Mesh m = meshes.get(part.getMesh());
        HydroField hf = hydroFields.get(part.getMesh());
        int elemPart = part.getElem();
        
        //System.out.println("movepart");
        // Set particles free once they pass their defined release time (hours)
        if (part.getFree()==false)
        {
            if (time > part.getReleaseTime())
            {
                part.setFree(true);
                part.setStatus(1);
                //System.out.println("Particle "+i+" released from "+particles[i].getLocation()[0]+" "+particles[i].getLocation()[1]+" at t="+time);
                // This is updated by many threads - calculation now done outside of move
                //freeViableSettleExit[0]++;
            }
        }
        if (part.getFree()==true && part.getArrived()==false && part.getBoundaryExit()==false)
        {    
            // Increment in particle age
            part.incrementAge(subStepDt/3600.0); // particle age in hours
            // Increment in particle degree days accumulated (can be used for maturation or mortality in larvae)
            double temperature = hf.getT()[tt][part.getDepthLayer()][m.getTrinodes()[0][part.getElem()]];
            part.incrementDegreeDays(temperature,rp);
            
//            System.out.printf("PARTICLE %d : %s\n",part.getID(),part.printLocation());
            
            //System.out.println("particle able to move");
            
            //System.out.printf("%d\n",elemPart);
            // set and get the DEPTH layer for the particle based on tide state
//            if (tt>0)
//            {
//                if (el[tt][trinodes[elemPart][0]]>el[tt-1][trinodes[elemPart][0]])
//                {
//                    particles[i].setDepthLayer(behaviour,"flood");
//                } else {
//                    particles[i].setDepthLayer(behaviour,"ebb");
//                }
//            }

            if (meshes.get(part.getMesh()).getType().equalsIgnoreCase("FVCOM"))
            {
                part.setDepth(rp.D_hVert,rp.sinkingRateMean,rp.sinkingRateStd,subStepDt,m.getDepthUvnode()[elemPart]);

                // set depth layer based on depth in metres
                part.setLayerFromDepth(m.getDepthUvnode()[elemPart],m.getSiglay());
            }
            
            // Find the salinity in the neighbourhood of the particle (used to compute instantaneous mortality rate).
            // This is stored at NODES as opposed to ELEMENT CENTROIDS.
            // So need to get the value from each of the corners and calculate 
            // a value at the particle location (similar method to getting velocity from nearest centroids).
            if (st == 0)
            {
                double salinity = 0;
                double mort = 0;
//                                if (rp.salinityMort == true)
//                                {
//                                    salinity = particles[i].salinity(tt,sal,trinodes);
//                                    particles[i].setMortRate(salinity);
//                                }
                part.setDensity();
            }

            double advectStep[] = new double[2];

            if (rp.rk4==true)
            {
                advectStep = part.rk4Step(hydroFields, meshes,
                    tt, st, subStepDt, rp.stepsPerStep, rp.coordRef);   
            }
            else
            {
                advectStep = part.eulerStep(hydroFields, meshes,
                    tt, st, subStepDt, rp.stepsPerStep, rp.coordRef);
            }
//            System.out.printf("ADVECT: Euler=[%.3e,%.3e] RK4=[%.3e,%.3e]\n",
//                advectStep[0],advectStep[1],advectStep2[0],advectStep2[1]);

            // Reverse velocities if running backwards
            if (rp.backwards == true)
            {
                advectStep[0] = -advectStep[0];
                advectStep[1] = -advectStep[1];
            }

            // 3. Calculate diffusion (random walk step)
//                            if (variableDiff==true)
//                            {
//                                D_h=1000*viscofm[tt][elemPart*10+dep];
//                            }

            double diff_X = 0;
            double diff_Y = 0;
            if (rp.diffusion==true)
            {
                // Use in-built RNG that is intented for multithread concurrent use. Also saves importing anything.
                diff_X = ThreadLocalRandom.current().nextDouble(-1.0,1.0)*Math.sqrt(6*rp.D_h*subStepDt);///(double)rp.stepsPerStep);
                diff_Y = ThreadLocalRandom.current().nextDouble(-1.0,1.0)*Math.sqrt(6*rp.D_h*subStepDt);///(double)rp.stepsPerStep);
            }
            double[] behave_uv = part.behaveVelocity(rp.behaviour);

            //System.out.println("D_h = "+rp.D_h+" diff_X = "+diff_X+" diff_Y "+diff_Y+" ran1 = "+ran1+" ran2 = "+ran2);
//            System.out.println("Distances travelled: X "+advectStep[0]+" "+diff_X+" --- Y "+advectStep[1]+" "+diff_Y);

            double dx = advectStep[0]+subStepDt*behave_uv[0]+diff_X;
            double dy = advectStep[1]+subStepDt*behave_uv[1]+diff_Y;
            
            if (rp.coordRef.equalsIgnoreCase("WGS84"))
            {
                // These two methods always give the same first significant figure and generally the second, later sig. figs. differ in general
                //double[] dXY1 = distanceMetresToDegrees1(new double[]{dx,dy}, part.getLocation());
                double[] dXY2 = distanceMetresToDegrees2(new double[]{dx,dy}, part.getLocation());
                
//                System.out.printf("distance %.04f %.04f (lat: %.04f) --- Method1 x: %.04e y: %.04e --- Method2 x: %.04e y: %.04e\n",
//                        dx,dy,part.getLocation()[1],dXY1[0],dXY1[1],dXY2[0],dXY2[1]);
                
                dx = dXY2[0];
                dy = dXY2[1];
            }
            
            // 4. update particle location
            double newlocx=part.getLocation()[0] + dx; 
            double newlocy=part.getLocation()[1] + dy;
            //System.out.println("Old = ("+part.getLocation()[0]+", "+part.getLocation()[1]+") --- New = ("+newlocx+", "+newlocy+")");

            // find element containing particle and update seach counts for diagnosis
 
            int[] c = null;
          
            int[] el = new int[2];
            if (m.getType().equalsIgnoreCase("ROMS"))
            {
                el = part.getROMSnearestPointU();
                //System.out.println("NearestROMSPointU: "+el[0]+" "+el[1]);
            }
            else
            {
                el[0] = part.getElem();
            }
                       
            // Is the particle still in the present mesh, should it change mesh, and has it exited the overall domain?
            part.meshSelectOrExit(new double[]{newlocx,newlocy},meshes,rp);
           
            // ***************************** By this point, the particle has been allocated to a mesh and new locations set etc ***********************
            
            // set particle to become able to settle after a predefined time
            if (((part.getAge()>rp.viabletime) ||
                    (part.getDegreeDays() > rp.viableDegreeDays && rp.viableDegreeDays > 0)) && part.getViable()==false)
            {
                System.out.printf("Particle became viable, ID %d age %f degreeDays %f\n",part.getID(),part.getAge(),part.getDegreeDays());                                  
                part.setViable(true);
                part.setStatus(2);
                // This is updated by many threads - calculation now done outside of move
                //freeViableSettleExit[1]++;
            }

            // Stop particles in their tracks if they exceed a maximum age
            if ((part.getAge()>rp.maxParticleAge && rp.maxParticleAge > 0) ||
                    (part.getDegreeDays() > rp.maxDegreeDays && rp.maxDegreeDays > 0))
            {
                part.setFree(false);
                part.setStatus(666);
            }

            

            // **************** if able to settle, is it close to a possible settlement location? ******************************
            //System.out.println("Patricle age = "+particles.get(i).getAge()+" Viabletime/3600 = "+viabletime/3600.0+" viable = "+particles.get(i).getViable());
            if (part.getViable()==true)
            {
                //System.out.println(particles[i].getViable());

                for (HabitatSite site : habitatEnd)
                {
                    //System.out.println("In habitatEnd check "+part.getLocation()[0]+" "+part.getLocation()[1]);
//                    double dist = Math.sqrt((part.getLocation()[0]-site.getLocation()[0])*(part.getLocation()[0]-site.getLocation()[0])+
//                            (part.getLocation()[1]-site.getLocation()[1])*(part.getLocation()[1]-site.getLocation()[1]));
                    double dist = Particle.distanceEuclid2(part.getLocation()[0], part.getLocation()[1],
                            site.getLocation()[0], site.getLocation()[1], rp.coordRef);
                    
                    double distThresh = rp.thresh;
                    // Correct the threshold distance to a value in degrees, if required - REMOVED
//                    if (rp.coordRef.equalsIgnoreCase("WGS84"))
//                    {
//                        distThresh = distThresh / 111206;
//                    }
                    
                    if (dist < distThresh && part.getSettledThisHour()==false)
                    {
                        //IOUtils.arrivalToFile(part, currentDate, whereami, loc, filename, true);
                        
                        //System.out.printf("settlement: %d at %d\n",i,loc);
                        if (rp.endOnArrival==true)
                        {
                            part.setArrived(true);
                            part.setStatus(3);
                        }

                        part.setSettledThisHour(true);
                        part.setLastArrival(site.getID());
                        break;
                    }
                }    
            }            
        }
    }
    
    /**
     * Calculate a transport distance, provided in metres, in degrees.
     * This uses the equirectangular method (https://www.movable-type.co.uk/scripts/latlong.html),
     * which is not particularly accurate but should be OK for short distances.
     * 
     * @param distanceMetres    calculated transport distance [dx,dy]
     * @param location          particle location in degrees [lon,lat]
     * @return 
     */
    public static double[] distanceMetresToDegrees1(double[] distanceMetres, double[] location)
    {
        double[] distanceDegrees = new double[2];
        // Distance per degree of longitude changes with latitude
        distanceDegrees[0] = distanceMetres[0] / (111206*Math.cos(2*Math.PI*location[1]/360)); 
        // Distance per degree of latitude remains broadly constant
        distanceDegrees[1] = distanceMetres[1] / 111206; // 111206 = average radius of earth (6371000) * tan(1deg)
        
        return distanceDegrees;
    }
    
    /**
     * Calculate a transport distance, provided in metres, in degrees.
     * This should be more accurate than the equirectangular approximation.
     * (http://www.csgnetwork.com/degreelenllavcalc.html; view source).
     * 
     * @param distanceMetres    calculated transport distance [dx,dy]
     * @param location          particle location in degrees [lon,lat]
     * @return 
     */
    public static double[] distanceMetresToDegrees2(double[] distanceMetres, double[] location)
    {
        double[] distanceDegrees = new double[2];
        
        // Set up "Constants" for calculating distances
        double m1 = 111132.92;     // latitude calculation term 1
        double m2 = -559.82;       // latitude calculation term 2
        double m3 = 1.175;         // latitude calculation term 3
        double m4 = -0.0023;       // latitude calculation term 4
        double p1 = 111412.84;     // longitude calculation term 1
        double p2 = -93.5;         // longitude calculation term 2
        double p3 = 0.118;         // longitude calculation term 3

        double latRad = 2*Math.PI*location[1]/360.0;
        
        // Calculate the length of a degree of latitude and longitude in meters
        double latlen = m1 + (m2 * Math.cos(2 * latRad)) + (m3 * Math.cos(4 * latRad)) +
                (m4 * Math.cos(6 * latRad));
        double longlen = (p1 * Math.cos(latRad)) + (p2 * Math.cos(3 * latRad)) +
                    (p3 * Math.cos(5 * latRad));
        
        distanceDegrees[0] = distanceMetres[0]/longlen;
        distanceDegrees[1] = distanceMetres[1]/latlen;
        
        return distanceDegrees;
    }
    
    /**
     * Convert a difference in lat/longs to a distance in metres.
     * Essentially the inverse function of distanceMetresToDegrees2.
     * 
     * @param distanceDegrees
     * @param location
     * @return 
     */
    public static double[] distanceDegreesToMetres(double[] distanceDegrees, double[] location)
    {
        double[] distanceMetres = new double[2];
        
        // Set up "Constants" for calculating distances
        double m1 = 111132.92;     // latitude calculation term 1
        double m2 = -559.82;       // latitude calculation term 2
        double m3 = 1.175;         // latitude calculation term 3
        double m4 = -0.0023;       // latitude calculation term 4
        double p1 = 111412.84;     // longitude calculation term 1
        double p2 = -93.5;         // longitude calculation term 2
        double p3 = 0.118;         // longitude calculation term 3

        double latRad = 2*Math.PI*location[1]/360.0;
        
        // Calculate the length of a degree of latitude and longitude in meters
        double latlen = m1 + (m2 * Math.cos(2 * latRad)) + (m3 * Math.cos(4 * latRad)) +
                (m4 * Math.cos(6 * latRad));
        double longlen = (p1 * Math.cos(latRad)) + (p2 * Math.cos(3 * latRad)) +
                    (p3 * Math.cos(5 * latRad));
        
        distanceMetres[0] = distanceDegrees[0]*longlen;
        distanceMetres[1] = distanceDegrees[1]*latlen;
        
        return distanceMetres;
    }
    
    /**
     * 
     * @param x
     * @param y
     * @param m
     * @param rp
     * @return 
     */
    public static int openBoundaryCheck(float x, float y, Mesh m, RunProperties rp)
    {
        // check whether the particle has gone within a certain range of one of the boundary nodes
        // (make it settle there, even if it is inviable)
        int nBNode = 0;
        if (m.getType().equalsIgnoreCase("FVCOM"))
        {
            nBNode = m.getOpenBoundaryNodes().length;
        }
        else if (m.getType().equalsIgnoreCase("ROMS"))
        {
            nBNode = m.getConvexHull().length;
        }

        for (int loc = 0; loc < nBNode; loc++)
        {
            double dist = 3001;

            if (m.getType().equalsIgnoreCase("FVCOM"))
            {
                dist = Particle.distanceEuclid2(x, y,
                    m.getNodexy()[0][m.getOpenBoundaryNodes()[loc]], m.getNodexy()[1][m.getOpenBoundaryNodes()[loc]], rp.coordRef);
            }
            else if (m.getType().equalsIgnoreCase("ROMS"))
            {
                dist = Particle.distanceEuclid2(x, y,
                    m.getConvexHull()[loc][0], m.getConvexHull()[loc][1], rp.coordRef);

            }

            //System.out.println("dist to OBC loc = "+dist);

            double distThresh = 3000;

            if (dist < distThresh)
            {
                return loc;
            }
        }
        // If not close to any open boundary points, must be a land departure => return -1
        return -1;
    }

}
