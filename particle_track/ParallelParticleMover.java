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
 * @author sa01ta
 */
public class ParallelParticleMover implements Callable<List<Particle>> {

    private final List<Particle> particles;
    private final double time;
    private final int hour;
    private final int step;
    private final double subStepDt;
    private final RunProperties rp;

    private final int[] allElems;

    private final List<Mesh> meshes;
    private final List<HydroField> hydroFields;

    private final List<HabitatSite> habitatEnd;
    private final int[] searchCounts;
    private final double[] minMaxDistTrav;
    private final boolean isDaytime;

    public ParallelParticleMover(List<Particle> particles, double time, int hour, int step, double subStepDt,
                                 RunProperties rp,
                                 List<Mesh> meshes,
                                 List<HydroField> hydroFields,
                                 List<HabitatSite> habitatEnd,
                                 int[] allElems,
                                 int[] searchCounts,
                                 double[] minMaxDistTrav,
                                 boolean isDaytime) {
        this.particles = particles;
        this.time = time;
        this.hour = hour;
        this.step = step;
        this.subStepDt = subStepDt;
        this.rp = rp;
        this.meshes = meshes;
        this.hydroFields = hydroFields;
        this.allElems = allElems;
        this.habitatEnd = habitatEnd;
        this.searchCounts = searchCounts;
        this.minMaxDistTrav = minMaxDistTrav;
        this.isDaytime = isDaytime;
    }

    @Override
    public ArrayList<Particle> call() throws Exception {
        for (Particle part : particles) {
            move(part, time, hour, step, subStepDt, rp, meshes, hydroFields, habitatEnd, allElems, searchCounts,
                    minMaxDistTrav, isDaytime);
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
     * @param hour
     * @param step
     * @param subStepDt
     * @param rp
     * @param meshes
     * @param hydroFields
     * @param habitatEnd
     * @param allElems
     * @param searchCounts
     * @param minMaxDistTrav
     */
    public static void move(Particle part, double time, int hour, int step, double subStepDt,
                            RunProperties rp,
                            List<Mesh> meshes,
                            List<HydroField> hydroFields,
                            List<HabitatSite> habitatEnd,
                            int[] allElems,
                            int[] searchCounts,
                            double[] minMaxDistTrav,
                            boolean isDaytime) {

        Mesh m = meshes.get(part.getMesh());
        HydroField hf = hydroFields.get(part.getMesh());
        int elemPart = part.getElem();
        int nDims = rp.verticalDynamics ? 3 : 2;

        // Set particles free once they pass their defined release time (hours)
        if (!part.getFree()) {
            if (time > part.getReleaseTime()) {
                part.setFree(true);
                part.setStatus(1);
            }
        }
        if (part.getFree() && !part.getArrived() && !part.getBoundaryExit()) {

            // Three types of possible movement: advection, diffusion, and active swimming/sinking; [x,y,(z)]
            double[] advectStep = new double[nDims];
            double[] diffusion = new double[nDims];
            double[] activeMovement = new double[nDims];
            double[] displacement = new double[nDims];
            for (int i=0; i < nDims; i++) {
                advectStep[i] = 0;
                diffusion[i] = 0;
                activeMovement[i] = 0;
                displacement[i] = 0;
            }

            // Increment in particle age & degree days
            part.incrementAge(subStepDt / 3600.0); // particle age in hours
            if (!rp.readHydroVelocityOnly) {
                double temperature = hf.getT()[hour][part.getDepthLayer()][m.getTrinodes()[0][part.getElem()]];
                part.incrementDegreeDays(temperature, rp);
            }

            // Find particle depth layer
            double D_hVertDz = 0;
            if (rp.fixDepth) {
                part.setDepth(rp.startDepth, m.getDepthUvnode()[elemPart]);
                part.setLayerFromDepth(m.getDepthUvnode()[elemPart], m.getSiglay());
            } else if (meshes.get(part.getMesh()).getType().equalsIgnoreCase("FVCOM") || meshes.get(part.getMesh()).getType().equalsIgnoreCase("ROMS_TRI")) {
                float localSalinity = 35; // TODO: Why is this set to 35? 35 is used in subsequent if statements...
                // Calculate the gradient in vertical diffusion, if required
                if (rp.variableDiff || rp.salinityThreshold < 35) {
                    // Calculate the vertical diffusivity profile for the particle location
                    float dep = (float) part.getDepth();
                    float localDepth = m.getDepthUvnode()[elemPart]; // TODO: This ignores zeta -- use HydroField.getWaterDepthUvnode(), or just ignore
                    // Get the sigma depths at this location, and compare with particle depth
                    // For particles on surface or sea bed, set layerAbove = layerBelow = [0 | sigDepths.length-1]
                    float[] sigDepths = m.getSiglay();
                    int layerBelow = sigDepths.length;
                    for (int i = 0; i < sigDepths.length; i++) {
                        if (dep < sigDepths[i] * localDepth) {
                            layerBelow = i;
                            break;
                        }
                    }
                    int layerAbove = layerBelow - 1;
                    if (layerAbove < 0) {
                        layerAbove = 0;
                    }
                    if (layerBelow == sigDepths.length) {
                        layerBelow = layerAbove;
                    }

                    float sigLayerHeight = (sigDepths[layerBelow] - sigDepths[layerAbove]) * localDepth;  // height of sigma layer

                    if (rp.variableDiff) {
                        float diffAbove = hf.getDiffVert()[hour][layerAbove][m.getTrinodes()[0][part.getElem()]];  // TODO: vertical diffusion coefficient is not currently read in
                        float diffBelow = hf.getDiffVert()[hour][layerBelow][m.getTrinodes()[0][part.getElem()]];
                        float diffusionDifference = Math.abs(diffAbove - diffBelow);
                        // Calculate a gradient, as long as particle is not on the bed or the surface, else default = 0
                        if (sigLayerHeight != 0) {
                            D_hVertDz = diffusionDifference / sigLayerHeight;
                        }
                    }
                    if (rp.salinityThreshold < 35 && rp.species.equalsIgnoreCase("sealice")) {
                        // Check local salinity in order to induce lice sinking behaviour if too low
                        float salAbove = hf.getS()[hour][layerAbove][m.getTrinodes()[0][part.getElem()]];
                        float salBelow = hf.getS()[hour][layerBelow][m.getTrinodes()[0][part.getElem()]];
                        float salinityDifference = salBelow - salAbove;
                        float dzPosition = dep - localDepth * sigDepths[layerAbove];
                        if (sigLayerHeight != 0) {
                            localSalinity = salAbove + dzPosition * salinityDifference / sigLayerHeight;
                        } else {
                            localSalinity = salAbove;
                        }
                    }
                }
                if (part.getStatus()==2) {
                    if (localSalinity < rp.salinityThreshold) {
                        activeMovement[2] = part.sink(rp);
                    } else if (isDaytime) {
                        activeMovement[2] = part.swim(rp);
                    }
                }
            }

            // Implement mortality once per hour
            // TODO: This code is old and outdated. Also, salinity is calculated above if rp.salinityThreshold < 35
            if (step == 0) {
                double salinity = 0;
                double mort = 0;
//                if (rp.salinityMort == true) {
//                    salinity = part.salinity(tt,sal,trinodes);
//                    part.setMortRate(salinity);
//                }
                part.setDensity();
            }

            if (rp.rk4) {
                advectStep = part.rk4Step(hydroFields, meshes, hour, step, subStepDt, rp.stepsPerStep, rp.coordRef, rp.verticalDynamics);
            } else {
                advectStep = part.eulerStep(hydroFields, meshes, hour, step, subStepDt, rp.stepsPerStep, rp.coordRef, rp.verticalDynamics);
            }

            // Reverse velocities if running backwards
            if (rp.backwards) {
                for (int i = 0; i < advectStep.length; i++) {
                    advectStep[i] *= -1;
                }
            }

            if (rp.diffusion) {
                // Use in-built RNG that is intended for multithread concurrent use. Also saves importing anything.
                diffusion[0] = ThreadLocalRandom.current().nextDouble(-1.0, 1.0) * Math.sqrt(6 * rp.D_h * subStepDt);
                diffusion[1] = ThreadLocalRandom.current().nextDouble(-1.0, 1.0) * Math.sqrt(6 * rp.D_h * subStepDt);
                diffusion[2] = part.verticalDiffusion(rp, D_hVertDz, subStepDt);
            }
            // part.behaveVelocity was here for xy movement, but returned 0's with no plans to expand

            for (int i = 0; i < nDims; i++) {
                displacement[i] = advectStep[i] + subStepDt * activeMovement[i] + diffusion[i];
            }


            if (rp.coordRef.equalsIgnoreCase("WGS84")) {
                double[] dXY2 = distanceMetresToDegrees2(new double[]{displacement[0], displacement[1]}, part.getLocation());
                displacement[0] = dXY2[0];
                displacement[1] = dXY2[1];
            }

            // Update particle location, changing mesh or exit status if necessary
            double newlocx = part.getLocation()[0] + displacement[0];
            double newlocy = part.getLocation()[1] + displacement[1];
            part.addX(Math.abs(displacement[0]));
            part.addY(Math.abs(displacement[1]));

            part.meshSelectOrExit(new double[]{newlocx, newlocy}, meshes, rp);
            if (rp.verticalDynamics) {
                double newDepth = part.getDepth() + displacement[2];
                part.addZ(Math.abs(displacement[2]));
                part.setDepth(newDepth, m.getDepthUvnode()[elemPart]);
                part.setLayerFromDepth(m.getDepthUvnode()[elemPart], m.getSiglay());
            }

            if (step == 0) {
                IOUtils.writeMovements(part, hour, step, displacement, advectStep, activeMovement, diffusion, "movementFile.dat", true);
            }

            // ***************************** By this point, the particle has been allocated to a mesh and new locations set etc ***********************
            // set particle to become able to settle after a predefined time
            if (!part.getViable()) {
                if (part.canBecomeViable(rp)) {
                    part.setViable(true);
                    part.setStatus(2);
                }
            }
            // Stop particles in their tracks if they exceed a maximum age
            if (part.isTooOld(rp)) {
                part.setFree(false);
                part.setStatus(666);
            }

            // **************** if able to settle, is it close to a possible settlement location? ******************************
            if (part.getViable()) {
                for (HabitatSite site : habitatEnd) {
                    double dist = Particle.distanceEuclid2(part.getLocation()[0], part.getLocation()[1],
                            site.getLocation()[0], site.getLocation()[1], rp.coordRef);
                    if (dist < rp.thresh && !part.getSettledThisHour()) {
                        if (rp.endOnArrival) {
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
     * @param distanceMetres calculated transport distance [dx,dy]
     * @param location       particle location in degrees [lon,lat]
     * @return
     */
    public static double[] distanceMetresToDegrees1(double[] distanceMetres, double[] location) {
        double[] distanceDegrees = new double[2];
        // Distance per degree of longitude changes with latitude
        distanceDegrees[0] = distanceMetres[0] / (111206 * Math.cos(2 * Math.PI * location[1] / 360));
        // Distance per degree of latitude remains broadly constant
        distanceDegrees[1] = distanceMetres[1] / 111206; // 111206 = average radius of earth (6371000) * tan(1deg)

        return distanceDegrees;
    }

    /**
     * Calculate a transport distance, provided in metres, in degrees.
     * This should be more accurate than the equirectangular approximation.
     * (http://www.csgnetwork.com/degreelenllavcalc.html; view source).
     *
     * @param distanceMetres calculated transport distance [dx,dy]
     * @param location       particle location in degrees [lon,lat]
     * @return
     */
    public static double[] distanceMetresToDegrees2(double[] distanceMetres, double[] location) {
        double[] distanceDegrees = new double[2];

        // Set up "Constants" for calculating distances
        double m1 = 111132.92;     // latitude calculation term 1
        double m2 = -559.82;       // latitude calculation term 2
        double m3 = 1.175;         // latitude calculation term 3
        double m4 = -0.0023;       // latitude calculation term 4
        double p1 = 111412.84;     // longitude calculation term 1
        double p2 = -93.5;         // longitude calculation term 2
        double p3 = 0.118;         // longitude calculation term 3

        double latRad = 2 * Math.PI * location[1] / 360.0;

        // Calculate the length of a degree of latitude and longitude in meters
        double latlen = m1 + (m2 * Math.cos(2 * latRad)) + (m3 * Math.cos(4 * latRad)) +
                (m4 * Math.cos(6 * latRad));
        double longlen = (p1 * Math.cos(latRad)) + (p2 * Math.cos(3 * latRad)) +
                (p3 * Math.cos(5 * latRad));

        distanceDegrees[0] = distanceMetres[0] / longlen;
        distanceDegrees[1] = distanceMetres[1] / latlen;

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
    public static double[] distanceDegreesToMetres(double[] distanceDegrees, double[] location) {
        double[] distanceMetres = new double[2];

        // Set up "Constants" for calculating distances
        double m1 = 111132.92;     // latitude calculation term 1
        double m2 = -559.82;       // latitude calculation term 2
        double m3 = 1.175;         // latitude calculation term 3
        double m4 = -0.0023;       // latitude calculation term 4
        double p1 = 111412.84;     // longitude calculation term 1
        double p2 = -93.5;         // longitude calculation term 2
        double p3 = 0.118;         // longitude calculation term 3

        double latRad = 2 * Math.PI * location[1] / 360.0;

        // Calculate the length of a degree of latitude and longitude in meters
        double latlen = m1 + (m2 * Math.cos(2 * latRad)) + (m3 * Math.cos(4 * latRad)) +
                (m4 * Math.cos(6 * latRad));
        double longlen = (p1 * Math.cos(latRad)) + (p2 * Math.cos(3 * latRad)) +
                (p3 * Math.cos(5 * latRad));

        distanceMetres[0] = distanceDegrees[0] * longlen;
        distanceMetres[1] = distanceDegrees[1] * latlen;

        return distanceMetres;
    }

    /**
     * @param x
     * @param y
     * @param m
     * @param rp
     * @return
     */
    public static int openBoundaryCheck(float x, float y, Mesh m, RunProperties rp) {
        // check whether the particle has gone within a certain range of one of the boundary nodes
        // (make it settle there, even if it is inviable)
        int nBNode = 0;
        if (m.getType().equalsIgnoreCase("FVCOM") || m.getType().equalsIgnoreCase("ROMS_TRI")) {
            nBNode = m.getOpenBoundaryNodes().length;
        } else if (m.getType().equalsIgnoreCase("ROMS")) {
            nBNode = m.getConvexHull().length;
        }

        for (int loc = 0; loc < nBNode; loc++) {
            double dist = 6001;

            if (m.getType().equalsIgnoreCase("FVCOM") || m.getType().equalsIgnoreCase("ROMS_TRI")) {
                dist = Particle.distanceEuclid2(x, y,
                        m.getNodexy()[0][m.getOpenBoundaryNodes()[loc]], m.getNodexy()[1][m.getOpenBoundaryNodes()[loc]], rp.coordRef);
            } else if (m.getType().equalsIgnoreCase("ROMS")) {
                dist = Particle.distanceEuclid2(x, y,
                        m.getConvexHull()[loc][0], m.getConvexHull()[loc][1], rp.coordRef);

            }

            //System.out.println("dist to OBC loc = "+dist);

            double distThresh = 6000;

            if (dist < distThresh) {
                return loc;
            }
        }
        // If not close to any open boundary points, must be a land departure => return -1
        return -1;
    }

}
