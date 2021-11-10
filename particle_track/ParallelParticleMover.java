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
    private final double elapsedHours;
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
    private final String currentDate;

    public ParallelParticleMover(List<Particle> particles, double elapsedHours, int hour, int step, double subStepDt,
                                 RunProperties rp,
                                 List<Mesh> meshes,
                                 List<HydroField> hydroFields,
                                 List<HabitatSite> habitatEnd,
                                 int[] allElems,
                                 int[] searchCounts,
                                 double[] minMaxDistTrav,
                                 boolean isDaytime, String currentDate) {
        this.particles = particles;
        this.elapsedHours = elapsedHours;
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
        this.currentDate = currentDate;
    }

    @Override
    public ArrayList<Particle> call() throws Exception {
        for (Particle part : particles) {
            move(part, elapsedHours, hour, step, subStepDt, rp, meshes, hydroFields, habitatEnd, allElems, searchCounts,
                    minMaxDistTrav, isDaytime, currentDate);
        }
        return new ArrayList<Particle>();
    }

    /**
     * Method to do all the actual particle movement stuff. This can be used either by "call" above to
     * do particle movements split over the cores, or can just be used directly in a serial loop
     * (i.e. when parallel==false)
     *
     * @param part
     * @param elapsedHours
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
    public static void move(Particle part, double elapsedHours, int hour, int step, double subStepDt,
                            RunProperties rp,
                            List<Mesh> meshes,
                            List<HydroField> hydroFields,
                            List<HabitatSite> habitatEnd,
                            int[] allElems,
                            int[] searchCounts,
                            double[] minMaxDistTrav,
                            boolean isDaytime, String currentDate) {

        Mesh m = meshes.get(part.getMesh());
        HydroField hf = hydroFields.get(part.getMesh());
        int elemPart = part.getElem();
        int nDims = rp.verticalDynamics ? 3 : 2;

        // Set particles free once they pass their defined release time (hours)
        if (!part.getFree()) {
            if (elapsedHours > part.getReleaseTime()) {
                part.setFree(true);
                part.setStatus(1);
            }
        }
        if (part.getFree() && !part.getArrived() && !part.getBoundaryExit()) {

            // Three types of possible movement: advection, diffusion, and active swimming/sinking; [x,y,(z)]
            // Note: (z) stays 0 unless rp.verticalDynamics == true
            double[] advectStep = {0,0,0};
            double[] diffusion = {0,0,0};
            double[] activeMovement = {0,0,0};
            double[] displacement = {0,0,0};
            int sink = 0;
            int swim = 0;

            // sigma layers
            float dep = (float) part.getDepth();
            float localDepth = m.getDepthUvnode()[elemPart]; // TODO: This ignores zeta -- use HydroField.getWaterDepthUvnode(), or just ignore
            // Get the sigma depths at this location, and compare with particle depth
            // For particles on surface or sea bed, set surroundingLayers[0 = layerBelow = [0 | sigDepths.length-1]
            // returns [layerBelow, layerAbove][depthBelow, depthAbove]
            float[][] nearestLayers = Mesh.findNearestSigmas(part.getDepth(), m.getSiglay(), m.getDepthUvnode()[elemPart]);
            float[][] nearestLevels = Mesh.findNearestSigmas(part.getDepth(), m.getSiglev(), m.getDepthUvnode()[elemPart]);

            // Increment in particle age & degree days
            part.incrementAge(subStepDt / 3600.0); // particle age in hours
            if (!rp.readHydroVelocityOnly) {
                double temperature = 0;
                if (rp.fixDepth) {
                    temperature = hf.getAvgFromTrinodes(m, part.getLocation(), part.getDepthLayer(), elemPart, hour, "temp", rp);
                } else {
                    temperature = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "temp", rp, nearestLayers);
                }
                part.incrementDegreeDays(temperature, rp);
            }

            // Vertical diffusion and salinity
            double KmBelow;
            double KmAbove;
            double KmGradient = 0; // positive gradient = more turbulent below = particle gets pushed deeper
            double Km_z = rp.D_hVert;
            float[][] nearestLevelsAdj;
            double Km_zAdj = rp.D_hVert;
            double depthAdj;

            if (rp.fixDepth) {
                part.setDepth(rp.startDepth, m.getDepthUvnode()[elemPart]);
                part.setLayerFromDepth(m.getDepthUvnode()[elemPart], m.getSiglay());
            } else if (meshes.get(part.getMesh()).getType().equalsIgnoreCase("FVCOM") || meshes.get(part.getMesh()).getType().equalsIgnoreCase("ROMS_TRI")) {
                float localSalinity = 35;
                if (rp.variableDiffusion) {
                    KmBelow = hf.getAvgFromTrinodes(m, part.getLocation(), (int) nearestLevels[0][0], elemPart, hour, "km", rp);
                    KmAbove = hf.getAvgFromTrinodes(m, part.getLocation(), (int) nearestLevels[1][0], elemPart, hour, "km", rp);
                    if (nearestLevels[0][1] == nearestLevels[1][1]) {
                        KmGradient = 0;
                    } else {
                        KmGradient = (KmBelow - KmAbove) / (nearestLevels[0][1] - nearestLevels[1][1]); // positive gradient = more turbulent below = particle gets pushed deeper
                    }
                    Km_z = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "km", rp, nearestLevels);
                    depthAdj = part.getDepth() + subStepDt * KmGradient / 2; // = (z + D_hVertDz/2 * dt) following Visser 1997
                    nearestLevelsAdj = Mesh.findNearestSigmas(depthAdj, m.getSiglev(), m.getDepthUvnode()[elemPart]);
                    Km_zAdj = hf.getValueAtDepth(m, part, part.getLocation(), depthAdj, hour, "km", rp, nearestLevelsAdj);
                    if (Km_z > 0.1) { // following Johnsen et al 2016
                        Km_z = 0.1;
                    }
                    if (Km_zAdj > 0.1) {
                        Km_zAdj = 0.1;
                    }
                }

                double localSalinity = 35;
                if (rp.salinityThreshold < 35) {
                    localSalinity = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "salinity", rp, nearestLayers);
                }
                if (part.getStatus()<3) {
                    if (localSalinity < rp.salinityThreshold || Km_z > Math.abs(rp.vertSwimSpeedMean)) { // TODO: Is localKm the correct measure of turbulence?
                        activeMovement[2] = part.sink(rp);
                        sink++;
                    } else if (rp.swimLightLevel) {
                        activeMovement[2] = part.swim(rp, m, hf, hour);
                        if (Math.abs(activeMovement[2]) > 1e-10) {
                            swim++;
                        }
                    } else if (isDaytime) {
                        activeMovement[2] = part.swim(rp);
                        swim++;
                    }
                }
            }

            // Implement mortality once per hour
            if (step == 0) {
                part.setDensity();
            }

            // advection
            if (rp.rk4) {
                advectStep = part.rk4Step(hydroFields, meshes, hour, step, subStepDt, rp.stepsPerStep, rp.coordRef, rp.verticalDynamics);
            } else {
                advectStep = part.eulerStep(hydroFields, meshes, hour, step, subStepDt, rp.stepsPerStep, rp.coordRef, rp.verticalDynamics);
            }
            if (rp.backwards) {
                for (int i = 0; i < advectStep.length; i++) {
                    advectStep[i] *= -1;
                }
            }

            if (rp.diffusion) {
                diffusion = part.diffuse(rp, KmGradient, Km_zAdj, subStepDt, "uniform");
            }

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

            if (part.getID() % 2000 == 0) {
                IOUtils.writeMovements(part, isDaytime, elapsedHours, currentDate, hour, step, displacement, advectStep, activeMovement, diffusion, sink, swim, "movementFile.dat", true);
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
