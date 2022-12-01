/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

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

    public ParallelParticleMover(List<Particle> particles, double elapsedHours, int hour, int step, double subStepDt,
                                 RunProperties rp,
                                 List<Mesh> meshes,
                                 List<HydroField> hydroFields,
                                 List<HabitatSite> habitatEnd,
                                 int[] allElems,
                                 boolean isDaytime, String currentDate, int[][] elemActivity) {
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
    }

    @Override
    public ArrayList<Particle> call() throws Exception {
        for (Particle part : particles) {
            move(part, elapsedHours, hour, step, subStepDt, rp, meshes, hydroFields, habitatEnd, allElems,
                    isDaytime, currentDate, elemActivity);
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
                            boolean isDaytime, String currentDate, int[][] elemActivity) {

        Mesh m = meshes.get(part.getMesh());
        HydroField hf = hydroFields.get(part.getMesh());
        int elemPart = part.getElem();
        int nDims = rp.verticalDynamics ? 3 : 2;

        // Set particles free once they pass their defined release time (hours)
        if (!part.isFree()) {
            if (elapsedHours > part.getReleaseTime()) {
                part.setFree(true);
                part.setStatus(1);
            }
        }
        if (part.isFree() && !part.hasArrived() && !part.hasExited()) {

            // Three types of possible movement: advection, diffusion, and active swimming/sinking; [x,y,(z)]
            // Note: (z) stays 0 unless rp.verticalDynamics == true
            double[] advectStep;
            double[] diffusion = {0,0,0};
            double[] activeMovement = {0,0,0};
            double[] displacement = {0,0,0};
            int sink = 0;
            int swim = 0;
            double localTemperature = 12;
            double localSalinity = 35;
            double tempSurface = hf.getAvgFromTrinodes(m, part.getLocation(), 1, elemPart, hour, "temp", rp);

            // sigma layers
            float localDepth = m.getDepthUvnode()[elemPart]; // TODO: This ignores zeta -- use HydroField.getWaterDepthUvnode(), or just ignore
            // Get the sigma depths at this location, and compare with particle depth
            // For particles on surface or sea bed, set surroundingLayers[0 = layerBelow = [0 | sigDepths.length-1]
            // returns [layerBelow, layerAbove][depthBelow, depthAbove]
            float[][] nearestLayers = Mesh.findNearestSigmas(part.getDepth(), m.getSiglay(), m.getDepthUvnode()[elemPart]);
            float[][] nearestLevels = Mesh.findNearestSigmas(part.getDepth(), m.getSiglev(), m.getDepthUvnode()[elemPart]);

            // Get local salinity, temperature
            if (!rp.readHydroVelocityOnly) {
                if (rp.fixDepth) {
                    localSalinity = hf.getAvgFromTrinodes(m, part.getLocation(), part.getDepthLayer(), elemPart, hour, "salinity", rp);
                    localTemperature = hf.getAvgFromTrinodes(m, part.getLocation(), part.getDepthLayer(), elemPart, hour, "temp", rp);
                } else {
                    localSalinity = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "salinity", rp, nearestLayers);
                    localTemperature = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "temp", rp, nearestLayers);
                }
            }
            // Increment in particle age & degree days
            part.incrementAge(subStepDt / 3600.0); // particle age in hours
            part.incrementDegreeDays(localTemperature, rp);

            // Vertical diffusion and salinity
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

            if (rp.fixDepth) {
                part.setDepth(rp.startDepth, m.getDepthUvnode()[elemPart]);
                part.setLayerFromDepth(m.getDepthUvnode()[elemPart], m.getSiglay());
            } else if (meshes.get(part.getMesh()).getType().equalsIgnoreCase("FVCOM") || meshes.get(part.getMesh()).getType().equalsIgnoreCase("ROMS_TRI")) {
                if (rp.variableDiffusion) {
                    K_below = hf.getAvgFromTrinodes(m, part.getLocation(), (int) nearestLevels[0][0], elemPart, hour, "k", rp);
                    K_above = hf.getAvgFromTrinodes(m, part.getLocation(), (int) nearestLevels[1][0], elemPart, hour, "k", rp);
                    if (nearestLevels[0][1] == nearestLevels[1][1]) {
                        K_gradient = 0;
                    } else {
                        K_gradient = (K_below - K_above) / (nearestLevels[0][1] - nearestLevels[1][1]); // positive gradient = more turbulent below = particle gets pushed deeper
                    }
                    K_z = hf.getValueAtDepth(m, part, part.getLocation(), part.getDepth(), hour, "k", rp, nearestLevels);
                    depthAdj = part.getDepth() + subStepDt * K_gradient / 2; // = (z + K_gradient/2 * dt) following Visser 1997
                    nearestLevelsAdj = Mesh.findNearestSigmas(depthAdj, m.getSiglev(), m.getDepthUvnode()[elemPart]);
                    K_zAdj = hf.getValueAtDepth(m, part, part.getLocation(), depthAdj, hour, "k", rp, nearestLevelsAdj);
                    if (K_z > 0.1) { // following Johnsen et al 2016
                        K_z = 0.1;
                    }
                    if (K_zAdj > 0.1) {
                        K_zAdj = 0.1;
                    }
                    // TODO: Set mortality for HABS if turbulence / shear is to high
                }

                if (part.getStatus()<3) {
                    int activity = 2; // 0 = sink; 1 = swim; 2 = float
                    // following Sandvik et al 2020, citing on Crosbie 2019
                    double prSink = part.calcSinkProb(localSalinity, rp);
                    if(prSink > ThreadLocalRandom.current().nextDouble(0,1)  ||
                        (rp.variableDiffusion && Math.abs(K_z) > Math.abs(part.isViable() ? rp.vertSwimSpeedCopepodidMean : rp.vertSwimSpeedNaupliusMean))) {
                        activeMovement[2] = part.sink(rp);
                        activity = 0;
                        sink++;
                    } else if (rp.swimLightLevel) {
                        activeMovement[2] = part.swim(rp, m, hf, hour);
                        if (Math.abs(activeMovement[2]) > 1e-10) {
                            activity = 1;
                            swim++;
                        }
                    } else if (isDaytime) {
                        activeMovement[2] = part.swim(rp);
                        activity = 1;
                        swim++;
                    }
                    if (rp.recordElemActivity && part.getMesh()==0) {
                        elemActivity[elemPart][activity]++;
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
                diffusion = part.diffuse(rp, K_gradient, K_zAdj, subStepDt, "uniform");
            }

            if (rp.debug3D.contains("activity")) {
                activeMovement[2] = 0;
            }
            if (rp.debug3D.contains("currents")) {
                advectStep[2] = 0;
            }
            if (rp.debug3D.contains("diffusion")) {
                diffusion[2] = 0;
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
            double[] posInit = {part.getLocation()[0], part.getLocation()[1], part.getDepth()};
            double[] dActual = new double[3];
            double newlocx = part.getLocation()[0] + displacement[0];
            double newlocy = part.getLocation()[1] + displacement[1];

            part.meshSelectOrExit(new double[]{newlocx, newlocy}, meshes, rp);
            if (rp.verticalDynamics) {
                double newDepth = part.getDepth() + displacement[2];
                double maxAllowedDepth = m.getDepthUvnode()[elemPart] < rp.maxDepth ? m.getDepthUvnode()[elemPart] : rp.maxDepth;
                part.setDepth(newDepth, maxAllowedDepth);
                part.setLayerFromDepth(m.getDepthUvnode()[elemPart], m.getSiglay());
            }

            dActual[0] = part.getLocation()[0] - posInit[0];
            dActual[1] = part.getLocation()[1] - posInit[1];
            dActual[2] = part.getDepth() - posInit[2];
            part.addX(Math.abs(dActual[0]));
            part.addY(Math.abs(dActual[1]));
            part.addXY(dActual[0], dActual[1]);
            part.addZ(Math.abs(dActual[2]));

            if (rp.recordMovement && (part.getID() % (rp.nparts * rp.numberOfDays * 10) == 0)) {  // * 10 = sample of ~485 particles
                IOUtils.writeMovements(part, currentDate, hour, step, sink, swim, localTemperature, tempSurface, localSalinity, dActual, "movementFile.dat", true);
            }

            // ***************************** By this point, the particle has been allocated to a mesh and new locations set etc ***********************
            // set particle to become able to settle after a predefined time
            if (!part.isViable()) {
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
            if (part.isViable()) {
                for (HabitatSite site : habitatEnd) {
                    double dist = Particle.distanceEuclid2(part.getLocation()[0], part.getLocation()[1],
                            site.getLocation()[0], site.getLocation()[1], rp.coordRef);
                    if (dist < rp.thresh && !part.hasSettledThisHour()) {
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
        int nBNode = 0;
        if (m.getType().equalsIgnoreCase("FVCOM") || m.getType().equalsIgnoreCase("ROMS_TRI")) {
            nBNode = m.getOpenBoundaryNodes().length;
        } else if (m.getType().equalsIgnoreCase("ROMS")) {
            nBNode = m.getConvexHull().length;
        }

        for (int loc = 0; loc < nBNode; loc++) {
            double dist = 0;

            if (m.getType().equalsIgnoreCase("FVCOM") || m.getType().equalsIgnoreCase("ROMS_TRI")) {
                dist = Particle.distanceEuclid2(x, y,
                        m.getNodexy()[0][m.getOpenBoundaryNodes()[loc]], m.getNodexy()[1][m.getOpenBoundaryNodes()[loc]], rp.coordRef);
            } else if (m.getType().equalsIgnoreCase("ROMS")) {
                dist = Particle.distanceEuclid2(x, y,
                        m.getConvexHull()[loc][0], m.getConvexHull()[loc][1], rp.coordRef);
            }

            if (dist < rp.openBoundaryThresh) {
                return loc;
            }
        }
        // If not close to any open boundary points, must be a land departure => return -1
        return -1;
    }

    public static int openBoundaryCheck(Particle part, Mesh m, RunProperties rp) {
        // check whether the particle has entered into a boundary *element*
        // (make it settle there, even if it is inviable)
        int nBElems = m.getOpenBoundaryNodes().length;
        for (int elem = 0; elem < nBElems; elem++) {
            if (m.getOpenBoundaryNodes()[elem] == part.getElem()) {
                return part.getElem();
            }
        }
        return -1;
    }
}
