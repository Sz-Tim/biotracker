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
    private final double[][] u; 
    private final double[][] v;
    private final int[][] neighbours;
    private final double[][] uvnode;
    private final double[][] nodexy;
    private final int[][] trinodes;
    private final int[] allelems;
    private final double[][] bathymetry;
    private final double[] sigvec2;
    private final double[][] startlocs;
    private final double[][] endlocs;
    private final double[][] open_BC_locs;
    private final int[] searchCounts;
    private final int[][] particle_info;
    private final double[][] settle_density;
    private final double[] minMaxDistTrav;
    
    public ParallelParticleMover(List<Particle> particles, double time, int tt, int st, double subStepDt,
            RunProperties rp, 
            double[][] u, double[][] v,
            int[][] neighbours, double[][] uvnode,  double[][] nodexy, 
            int[][] trinodes, int[] allelems,
            double[][] bathymetry, double[] sigvec2, double[][] startlocs,
            double[][] endlocs, double[][] open_BC_locs,
            int[] searchCounts,
            int[][] particle_info, double[][] settle_density,
            double[] minMaxDistTrav)
    {
        this.particles=particles;
        this.time=time;
        this.tt=tt;
        this.st=st;
        this.subStepDt=subStepDt;
        this.rp=rp;
        this.u=u; 
        this.v=v;
        this.neighbours=neighbours;
        this.uvnode=uvnode;
        this.nodexy=nodexy;
        this.trinodes=trinodes;
        this.allelems=allelems;
        this.bathymetry=bathymetry;
        this.sigvec2=sigvec2;
        this.startlocs=startlocs;
        this.endlocs=endlocs;
        this.open_BC_locs=open_BC_locs;
        this.searchCounts=searchCounts;
        this.particle_info=particle_info;
        this.settle_density=settle_density;
        this.minMaxDistTrav=minMaxDistTrav;
    }
    
    @Override
    public ArrayList<Particle> call() throws Exception {
        for (Particle part : particles) {
            //System.out.println("Time "+time+" - Moving particle "+part.getID());
            move(part, time, tt, st, subStepDt, rp, 
                                    u, v, neighbours, uvnode, nodexy, trinodes, allelems, bathymetry, sigvec2, 
                                    startlocs, endlocs, open_BC_locs,
                                    searchCounts, 
                                    particle_info, settle_density, minMaxDistTrav);
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
     * @param bathymetry
     * @param sigvec2
     * @param startlocs
     * @param endlocs
     * @param open_BC_locs
     * @param searchCounts
     * @param particle_info
     * @param settle_density
     * @param minMaxDistTrav 
     */
    public static void move(Particle part, double time, int tt, int st, double subStepDt,
            RunProperties rp, 
            double[][] u, double[][] v,
            int[][] neighbours, double[][] uvnode,  double[][] nodexy, 
            int[][] trinodes, int[] allelems,
            double[][] bathymetry, double[] sigvec2, double[][] startlocs,
            double[][] endlocs, double[][] open_BC_locs,
            int[] searchCounts,
            int[][] particle_info, double[][] settle_density,
            double[] minMaxDistTrav)
    {
        // Set particles free once they pass thier defined release time (hours)
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
            
            part.increaseAge(subStepDt/3600.0); // particle age in hours
            //System.out.printf("PARTICLE %d \n",i);
            
            //System.out.println("particle able to move");
            int elemPart = part.getElem();
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
            // set depth layer based on fixed depth in metres
            part.setLayerFromDepth(bathymetry[elemPart][0],sigvec2);

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
                part.setDensity();
            }

            double advectStep[] = new double[2];
//                            System.out.printf("ADVECT: Euler=[%.3e,%.3e] RK4=[%.3e,%.3e]\n",
//                                    advectStep[0],advectStep[1],advectStep2[0],advectStep2[1]);
            if (rp.rk4==true)
            {
                advectStep = part.rk4Step(u, v, 
                    neighbours, uvnode,  nodexy, trinodes, allelems,
                    tt, st, subStepDt, rp.stepsPerStep, rp.depthLayers);   
            }
            else
            {
                System.err.println("Euler step not calculated --- edit method to work with 7 row velocities");
//                                    advectStep = particles[i].eulerStepOld(u, v, u1, v1, neighbours, uvnode, 
//                                        tt, st, subStepDt, stepsPerStep, recordsPerFile, fnum, lastday, depthLayers, 
//                                        spatialInterpolate, timeInterpolate);
            }

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
                diff_X = ThreadLocalRandom.current().nextDouble(-1.0,1.0)*Math.sqrt(6*rp.D_h*subStepDt/(double)rp.stepsPerStep);
                diff_Y = ThreadLocalRandom.current().nextDouble(-1.0,1.0)*Math.sqrt(6*rp.D_h*subStepDt/(double)rp.stepsPerStep);
            }
            double[] behave_uv = part.behaveVelocity(rp.behaviour);

            //System.out.println("D_h = "+rp.D_h+" diff_X = "+diff_X+" diff_Y "+diff_Y+" ran1 = "+ran1+" ran2 = "+ran2);
            //System.out.println("Distances travelled: X "+subStepDt*water_U+" "+diff_X*ran1+" Y "+subStepDt*water_U+" "+diff_Y*ran2);

            // 4. update particle location
            double newlocx=part.getLocation()[0]+advectStep[0]+subStepDt*behave_uv[0]+diff_X; // simplest possible "Euler"
            double newlocy=part.getLocation()[1]+advectStep[1]+subStepDt*behave_uv[1]+diff_Y;
            //System.out.println("Old = ("+particles[i].getLocation()[0]+", "+particles[i].getLocation()[1]+") --- New = ("+newlocx+", "+newlocy+")");

            // find element containing particle and update seach counts for diagnosis
            int[] c = Particle.findContainingElement(newlocx, newlocy, elemPart, 
                    nodexy, trinodes, neighbours, allelems);
            int whereami = c[0];
            for (int j = 0; j < 4; j++)
            {
                //moveStats[j+4] = c[j];
                searchCounts[j] += c[j+1];
            }
            
            // if particle is within the mesh, update location normally and save the distance travelled
            double distTrav = 0;
            if (whereami != -1)
            {
                distTrav = Math.sqrt((part.getLocation()[0]-newlocx)*(part.getLocation()[0]-newlocx)+
                        (part.getLocation()[1]-newlocy)*(part.getLocation()[1]-newlocy));
                part.setLocation(newlocx,newlocy);
                part.setElem(whereami);
                //System.out.printf("** MOVED **, new elem = %d (dist = %f)\n",particles[i].getElem(),Math.sqrt((newlocx-uvnode[particles[i].getElem()][0])*(newlocx-uvnode[particles[i].getElem()][0])+(newlocy-uvnode[particles[i].getElem()][1])*(newlocy-uvnode[particles[i].getElem()][1])));
            }
            // What's done here is probably not thread-safe; but it's not critical to operation.
            // Maybe save particle own min/max travel distances and calculate overall at end?
            if (distTrav>minMaxDistTrav[1])
            {
                minMaxDistTrav[1]=distTrav;
            }
            if (distTrav<minMaxDistTrav[0])
            {
                minMaxDistTrav[0]=distTrav;
            }

            // if particle has skipped out of the model domain, place it at the nearest element centroid
            if (whereami == -1)
            {
                int closest=Particle.nearestCentroid(part.getLocation()[0],part.getLocation()[1],uvnode);
                //fprintf('x%d',closest);
                part.setLocation(uvnode[closest][0],uvnode[closest][1]);
                part.setElem(closest);
            }

            // set particle to become able to settle after a predefined time
            if (part.getAge()>rp.viabletime && part.getViable()==false)
            {
                //System.out.println("Particle became viable");                                  
                part.setViable(true);
                part.setStatus(2);
                // This is updated by many threads - calculation now done outside of move
                //freeViableSettleExit[1]++;
            }

            // Stop particles in their tracks if they exceed a maximum age
            if (part.getAge()>rp.maxParticleAge && rp.maxParticleAge > 0)
            {
                part.setFree(false);
                part.setStatus(666);
            }

            // check whether the particle has gone within a certain range of one of the boundary nodes
            // (make it settle there, even if it is inviable)
            for (int loc = 0; loc < open_BC_locs.length; loc++)
                {
                    double dist = Math.sqrt((part.getLocation()[0]-open_BC_locs[loc][1])*(part.getLocation()[0]-open_BC_locs[loc][1])+
                            (part.getLocation()[1]-open_BC_locs[loc][2])*(part.getLocation()[1]-open_BC_locs[loc][2]));
                    if (dist < 1500)
                    {
                        //System.out.printf("Boundary stop: %d at %d\n",i,loc);
                        //part.setArrived(true);
                        part.setBoundaryExit(true);
                        part.setStatus(66);
                        // elements of particle_info are only updated once by a single thread
                        //particle_info[part.getID()][1] = -loc;//(int)startlocs[loc][0];
                        //particle_info[part.getID()][2] = (int)part.getAge();//((day-firstday)*24+tt);
                        // This is updated by many threads - calculation now done outside of move
                        //freeViableSettleExit[3]++;
                        break;

                    }
                }

            // if able to settle, is it close to a possible settlement
            // location?
            //System.out.println("Patricle age = "+particles.get(i).getAge()+" Viabletime/3600 = "+viabletime/3600.0+" viable = "+particles.get(i).getViable());
            if (part.getViable()==true)
            {
                //System.out.println(particles[i].getViable());

                for (int loc = 0; loc < endlocs.length; loc++)
                {
                    double dist = Math.sqrt((part.getLocation()[0]-endlocs[loc][1])*(part.getLocation()[0]-endlocs[loc][1])+
                            (part.getLocation()[1]-endlocs[loc][2])*(part.getLocation()[1]-endlocs[loc][2]));
                    if (dist < rp.thresh && part.getSettledThisHour()==false)
                    {
                        //IOUtils.arrivalToFile(part, currentDate, whereami, loc, filename, true);
                        
                        //System.out.printf("settlement: %d at %d\n",i,loc);
                        if (rp.endOnArrival==true)
                        {
                            part.setArrived(true);
                            part.setStatus(3);
                        }
                        // elements of particle_info and settle_density are only updated once by a single thread
//                        particle_info[part.getID()][1] = loc;//(int)startlocs[loc][0];
//                        particle_info[part.getID()][2] = (int)time;//((day-firstday)*24+tt);
//                        settle_density[part.getID()][0] = part.getDensity();
//                        // This is updated by many threads - calculation now done outside of move                        
                        //freeViableSettleExit[2]++;
                        if (rp.oldOutput == true)
                        {
                            part.reportArrival(part.getStartID(), loc, time, part.getDensity());
                        }
                        
                        part.setSettledThisHour(true);
                        part.setLastArrival(loc);
                        break;
                    }
                }    
            }            
        }
    }
}
