/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ThreadLocalRandom;

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
    private final boolean isDaytime;
    private final String currentDate;
    private final int[][] elemActivity;
    private final int[][] hourActivity;

    public ParallelParticleMover(List<Particle> particles, double elapsedHours, int hour, int step, double subStepDt,
                                 RunProperties rp,
                                 List<Mesh> meshes,
                                 List<HydroField> hydroFields,
                                 List<HabitatSite> habitatEnd,
                                 int[] allElems,
                                 boolean isDaytime, String currentDate, int[][] elemActivity, int[][] hourActivity) {
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
        this.isDaytime = isDaytime;
        this.currentDate = currentDate;
        this.elemActivity = elemActivity;
        this.hourActivity = hourActivity;
    }

    @Override
    public ArrayList<Particle> call() throws Exception {
        for (Particle part : particles) {
            move(part, elapsedHours, hour, step, subStepDt, rp, meshes, hydroFields, habitatEnd, allElems,
                    isDaytime, currentDate, elemActivity, hourActivity);
        }
        return new ArrayList<>();
    }

    /**
     * Method to do all the actual particle movement stuff. This can be used either by "call" above to
     * do particle movements split over the cores, or can just be used directly in a serial loop
     * (i.e. when parallel==false)
     */
    public static void move(Particle part, double elapsedHours, int hour, int step, double subStepDt,
                            RunProperties rp,
                            List<Mesh> meshes,
                            List<HydroField> hydroFields,
                            List<HabitatSite> habitatEnd,
                            int[] allElems,
                            boolean isDaytime, String currentDate, int[][] elemActivity, int[][] hourActivity) {

        Mesh m = meshes.get(part.getMesh());
        HydroField hf = hydroFields.get(part.getMesh());

        // Set particles free once they pass their defined release time (hours)
        if (!part.isFree()) {
            if (elapsedHours > part.getReleaseTime()) {
                part.setFree(true);
                part.setStatus(1);
            }
        }
        if (part.isFree() && !part.hasArrived() && !part.hasExited()) {

            // Three types of possible movement: advection, diffusion, and active swimming/sinking; [x,y,(z)]
            // Note: (z) stays 0 if rp.fixDepth == true
            double[] advectStep;
            double[] diffusion = {0,0,0};
            double[] activeMovement = {0,0,0};
            double[] displacement = {0,0,0};
            int sink = 0;
            int swim = 0;
            double localTemperature = 12;
            double localSalinity = 35;
            double tempSurface = rp.needT ? hf.getAvgFromTrinodes(m, part.getLocation(), 1, part.getElem(), hour, "temp", rp) : localTemperature;

            // sigma layers
            float localDepth = m.getDepthUvnode()[part.getElem()]; // TODO: This ignores zeta -- use HydroField.getWaterDepthUvnode(), or just ignore
            // Get the sigma depths at this location, and compare with particle depth
            // For particles on surface or sea bed, set surroundingLayers[0 = layerBelow = [0 | sigDepths.length-1]
            // returns [layerBelow, layerAbove][depthBelow, depthAbove]
            float[][] nearestLayers = Mesh.findNearestSigmas(part.getDepth(), m.getSiglay(), m.getDepthUvnode()[part.getElem()]);
            float[][] nearestLevels = Mesh.findNearestSigmas(part.getDepth(), m.getSiglev(), m.getDepthUvnode()[part.getElem()]);

            // Get local salinity, temperature
            if (!rp.readHydroVelocityOnly) {
                if (rp.fixDepth) {
                    if (rp.needS) {
                        localSalinity = hf.getAvgFromTrinodes(m, part.getLocation(), part.getDepthLayer(), part.getElem(), hour, "salinity", rp);
                    }
                    if (rp.needT) {
                        localTemperature = hf.getAvgFromTrinodes(m, part.getLocation(), part.getDepthLayer(), part.getElem(), hour, "temp", rp);
                    }
                } else {
                    if (rp.needS) {
                        localSalinity = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "salinity", rp, nearestLayers);
                    }
                    if (rp.needT) {
                        localTemperature = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "temp", rp, nearestLayers);
                    }
                }
            }
            // Increment in particle age & degree days
            part.incrementAge(subStepDt / 3600.0); // particle age in hours
            part.incrementDegreeDays(localTemperature, rp);

            // Diffusion
            double D_h = rp.D_h;
            double K_below;
            double K_above;
            double K_gradient = 0; // positive gradient = more turbulent below = particle gets pushed deeper
            double K_z = rp.D_hVert;
            float[][] nearestLevelsAdj;
            double K_zAdj = rp.D_hVert;
            double depthAdj;

            if (rp.salinityMort) {
                part.setMortRateSalinity(localSalinity);
            }

            if (rp.variableDh) {
                D_h = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "vh", rp, nearestLayers);
            }

            if (rp.fixDepth) {
                part.setDepth(rp.startDepth, m.getDepthUvnode()[part.getElem()]);
                part.setLayerFromDepth(m.getDepthUvnode()[part.getElem()], m.getSiglay());
            } else if (rp.FVCOM || meshes.get(part.getMesh()).getType().equals("ROMS_TRI")) {
                if (rp.variableDhV) {
                    K_below = hf.getAvgFromTrinodes(m, part.getLocation(), (int) nearestLevels[0][0], part.getElem(), hour, "k", rp);
                    K_above = hf.getAvgFromTrinodes(m, part.getLocation(), (int) nearestLevels[1][0], part.getElem(), hour, "k", rp);
                    if (nearestLevels[0][1] == nearestLevels[1][1]) {
                        K_gradient = 0;
                    } else {
                        K_gradient = (K_below - K_above) / (nearestLevels[0][1] - nearestLevels[1][1]); // positive gradient = more turbulent below = particle gets pushed deeper
                    }
                    K_z = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "k", rp, nearestLevels);
                    depthAdj = part.getDepth() + subStepDt * K_gradient / 2; // = (z + K_gradient/2 * dt) following Visser 1997
                    nearestLevelsAdj = Mesh.findNearestSigmas(depthAdj, m.getSiglev(), m.getDepthUvnode()[part.getElem()]);
                    K_zAdj = hf.getValueAtDepth(m, part, part.getLocation(), depthAdj, hour, "k", rp, nearestLevelsAdj);
                    if (K_z > 0.1) { // following Johnsen et al 2016
                        K_z = 0.1;
                    }
                    if (K_zAdj > 0.1) {
                        K_zAdj = 0.1;
                    }
                }

                if (part.getStatus()<3) {
                    int activity = 2; // 0 = sink; 1 = swim; 2 = float
                    // following Sandvik et al 2020, citing on Crosbie 2019
                    double prSink = part.calcSwimDownProb(localSalinity, rp);
                    if(prSink > ThreadLocalRandom.current().nextDouble(0,1)  ||
                        (rp.variableDhV && Math.abs(K_z) > Math.abs(part.isInfectious() ? rp.swimUpSpeedCopepodidMean : rp.swimUpSpeedNaupliusMean))) {
                        activeMovement[2] = Math.max(part.sink(rp, subStepDt), part.passive(rp, localSalinity, subStepDt));
                        activity = 0;
                        sink++;
                    } else if (rp.swimLightLevel) {
                        activeMovement[2] = part.swim(rp, m, hf, hour, subStepDt);
                        if (Math.abs(activeMovement[2]) > 1e-10) {
                            activity = 1;
                            swim++;
                        }
                    } else if (isDaytime) {
                        activeMovement[2] = part.swim(rp, subStepDt);
                        activity = 1;
                        swim++;
                    } else {
                        activeMovement[2] = part.passive(rp, localSalinity, subStepDt);
                    }
                    if (rp.recordActivity && part.getMesh()==0) {
                        elemActivity[part.getElem()][activity]++;
                        hourActivity[Math.toIntExact(Math.round(elapsedHours))][activity]++;
                    }
                }
            }

            // Implement mortality once per hour
            boolean approxHour = elapsedHours % 1 < 0.01 || elapsedHours % 1 > 0.99;
            if (step == 0 && approxHour) {
                part.setDensity();
            }

            // advection
            if (rp.rk4) {
                advectStep = part.rk4Step(hydroFields, meshes, hour, step, subStepDt, rp.stepsPerStep, rp.coordOS, rp.FVCOM);
            } else {
                advectStep = part.eulerStep(hydroFields, meshes, hour, step, subStepDt, rp.stepsPerStep, rp.coordOS, rp.fixDepth);
            }
            if (rp.backwards) {
                for (int i = 0; i < advectStep.length; i++) {
                    advectStep[i] *= -1;
                }
            }

            if (rp.diffusion) {
                diffusion = part.diffuse(D_h, K_gradient, K_zAdj, subStepDt, "uniform");
            }

            for (int i = 0; i < 3; i++) {
                displacement[i] = advectStep[i] + activeMovement[i] + diffusion[i];
            }

            if (!rp.coordOS) {
                double[] dXY2 = distanceMetresToDegrees2(new double[]{displacement[0], displacement[1]}, part.getLocation());
                displacement[0] = dXY2[0];
                displacement[1] = dXY2[1];
            }

            // Calculate new particle location, changing mesh or exit status if necessary
            double[] posInit = {part.getLocation()[0], part.getLocation()[1], part.getDepth()};
            double[] dActual = new double[3];
            double newlocx = part.getLocation()[0] + displacement[0];
            double newlocy = part.getLocation()[1] + displacement[1];

            // Location actually updated here
            part.moveInMesh(new double[]{newlocx, newlocy}, meshes, rp);
            //part.meshSelectOrExit(new double[]{newlocx, newlocy}, meshes, rp);
            double newDepth = part.getDepth() + displacement[2];
            double maxAllowedDepth = meshes.get(part.getMesh()).getDepthUvnode()[part.getElem()] < rp.maxDepth ? meshes.get(part.getMesh()).getDepthUvnode()[part.getElem()] : rp.maxDepth;
            part.setDepth(newDepth, maxAllowedDepth);
            part.setLayerFromDepth(meshes.get(part.getMesh()).getDepthUvnode()[part.getElem()], meshes.get(part.getMesh()).getSiglay());

            dActual[0] = part.getLocation()[0] - posInit[0];
            dActual[1] = part.getLocation()[1] - posInit[1];
            dActual[2] = part.getDepth() - posInit[2];
            part.addX(Math.abs(dActual[0]));
            part.addY(Math.abs(dActual[1]));
            part.addXY(dActual[0], dActual[1]);
            part.addZ(Math.abs(dActual[2]));

            if (rp.recordMovement && (part.getID() % (rp.nparts * rp.numberOfDays / 3) == 0)) {  // * 10 = sample of ~485 particles
                IOUtils.writeMovements(part, currentDate, hour, step, sink, swim, localTemperature, tempSurface, localSalinity, dActual, "movementFile.csv", true);
            }

            // ***************************** By this point, the particle has been allocated to a mesh and new locations set etc ***********************
            // set particle to become able to settle after a predefined time
            if (!part.isInfectious()) {
                if (part.canBecomeInfectious(rp)) {
                    part.setInfectious(true);
                    part.setStatus(2);
                }
            }
            // Stop particles in their tracks if they exceed a maximum age
            if (part.isTooOld(rp)) {
                part.setFree(false);
                part.setStatus(666);
            }

            // **************** if able to settle, is it close to a possible settlement location? ******************************
            if (part.isInfectious() || rp.connectImmature) {
                for (HabitatSite site : habitatEnd) {
                    double dist = Particle.distanceEuclid2(part.getLocation()[0], part.getLocation()[1],
                            site.getLocation()[0], site.getLocation()[1], rp.coordOS);
                    if (dist < rp.connectivityThresh && !part.hasSettledThisHour()) {  // TODO: I think this doesn't make sense? Necessary for part.setLastArrival(site.getID()) for connectivity though.
                        if (rp.endOnArrival && part.isInfectious()) {
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
        double latlen = m1 + (m2 * Math.cos(2 * latRad)) + (m3 * Math.cos(4 * latRad)) + (m4 * Math.cos(6 * latRad));
        double longlen = (p1 * Math.cos(latRad)) + (p2 * Math.cos(3 * latRad)) + (p3 * Math.cos(5 * latRad));

        distanceDegrees[0] = distanceMetres[0] / longlen;
        distanceDegrees[1] = distanceMetres[1] / latlen;

        return distanceDegrees;
    }

    /**
     * Convert a difference in lat/longs to a distance in metres.
     * Essentially the inverse function of distanceMetresToDegrees2.
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


    public static int openBoundaryCheck(float x, float y, Mesh m, RunProperties rp) {
        // check whether the particle has gone within a certain range of one of the boundary nodes
        // note that this uses the mesh NODES (nodexy), not elements or centroids
        // (make it settle there, even if it is inviable)

        // most particles exit north or west; WeStCOMS indexes start south, so reverse for efficiency
        if(rp.FVCOM) {
            int nBNode = m.getOpenBoundaryNodes().length - 1;
            double dist;
            for (int loc = nBNode; loc >= 0; loc--) {
                dist = Particle.distanceEuclid2(x, y,
                        m.getNodexy()[0][m.getOpenBoundaryNodes()[loc]], m.getNodexy()[1][m.getOpenBoundaryNodes()[loc]], rp.coordOS);
                if (dist < rp.openBoundaryThresh) {
                    return loc;
                }
            }
        } else if (m.getType().equals("ROMS")) {
            int nBNode = m.getConvexHull().length;
            double dist;
            for (int loc = 0; loc < nBNode; loc++) {
               dist = Particle.distanceEuclid2(x, y,
                            m.getConvexHull()[loc][0], m.getConvexHull()[loc][1], rp.coordOS);
                if (dist < rp.openBoundaryThresh) {
                    return loc;
                }
            }
        }
        // If not close to any open boundary points, must be a land departure => return -1
        return -1;
    }

    public static int openBoundaryCheck(double[] loc, int elem, Mesh m, RunProperties rp) {
        // check whether the particle has entered into a boundary *element*
        // (make it settle there, even if it is inviable)
        int elemNew = Particle.findContainingElement(loc, elem, m, false, 10)[0];
        if (m.getOpenBoundaryElems().contains(elemNew)) {
            return elemNew;
        }
        return -1;
    }

    public static boolean openBoundaryCheck(int elem, ArrayList<Integer> openBoundaryElems, RunProperties rp) {
        if (rp.checkOpenBoundaries) {
            return openBoundaryElems.contains(elem);
        }
        return false;
    }
}
