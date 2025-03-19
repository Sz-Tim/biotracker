/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

import java.awt.geom.Path2D;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;
import java.util.stream.IntStream;

//import org.apache.commons.lang.ArrayUtils;

/**
 * @author tomdude
 */
public class Particle {

    final private int id;
    // horizontal position
    private final double[] xy = new double[2];
    final private double[] startLoc = new double[2];
    private final String startSiteID;
    private int startSiteIndex = 0;

    private final boolean coordOS;

    private final ISO_datestr startDate;
    private double startTime = 0;

    private int mesh;
    private int elem;
    private int[] elemRomsU;
    private int[] elemRomsV;
    private int[] nearestROMSGridPointU;
    private int[] nearestROMSGridPointV;

    private double[][] nrList = new double[5][2];
    private final double[][] cornerList = new double[3][2];
    private int whichMesh;
    // release time
    private double releaseTime = 0;
    // vertical position
    private double depth;
    private int depLayer = 0;
    // settlement details
    private double age = 0;
    private double degreeDays = 0;
    private int status = 0;
    private double density = 1;
    private double mortRate = 0.01; // default hourly rate, based on results of Stein et al. 2005 (copepodid, nauplii rate marginally lower = 0.0078)
    private boolean arrived = false;
    private boolean viable = false;
    private boolean free = false;
    private boolean boundaryExit = false;

    private String species = "sealice";

    private double xTotal = 0;
    private double yTotal = 0;
    private double xyTotal = 0;
    private double zTotal = 0;

    // A list to store data on the arrivals made by each particle.
    // If rp.endOnArrival=true, this will contain a maximum of one element.
    private ArrayList<Arrival> arrivals = new ArrayList<>();

    public Particle(double xstart, double ystart, double startDepth, String startSiteID, int startIndex, int id, double mortalityRate,
                    ISO_datestr startDate, double startTime, boolean coordOS, String species) {
        this.id = id;
        this.xy[0] = xstart;
        this.xy[1] = ystart;
        this.startSiteID = startSiteID;
        this.startSiteIndex = startIndex;
        this.startDate = new ISO_datestr(startDate.getDateStr());
        this.startTime = startTime;
        this.startLoc[0] = xstart;
        this.startLoc[1] = ystart;
        this.mortRate = mortalityRate;
        this.depth = startDepth;
        this.coordOS = coordOS;
        this.species = species;
    }

    public Particle(double xstart, double ystart, double startDepth, String startSiteID, int startIndex, int id, double startDensity,
                    ISO_datestr startDate, double startTime, boolean coordOS) {
        this.id = id;
        this.xy[0] = xstart;
        this.xy[1] = ystart;
        this.startSiteID = startSiteID;
        this.startSiteIndex = startIndex;
        this.startDate = new ISO_datestr(startDate.getDateStr());
        this.startTime = startTime;
        this.startLoc[0] = xstart;
        this.startLoc[1] = ystart;
        this.density = startDensity;
        this.depth = startDepth;
        this.coordOS = coordOS;
    }

    /**
     * Create a new particle from the line of a location ASCII file.
     */
    public Particle(String locationFileLine, String species) {
        String[] values = locationFileLine.split(" ");
        // Don't use the 0 th  entry on the line (current time in locations) 
        this.id = Integer.parseInt(values[1]);
        this.startDate = new ISO_datestr(values[2]);
        this.age = Double.parseDouble(values[3]);
        this.startSiteID = values[4];
        this.xy[0] = Double.parseDouble(values[5]);
        this.xy[1] = Double.parseDouble(values[6]);
        this.elem = Integer.parseInt(values[7]);
        this.status = Integer.parseInt(values[8]);
        this.density = Double.parseDouble(values[9]);
        this.mesh = Integer.parseInt(values[10]);
        this.depth = Double.parseDouble(values[11]);
        this.degreeDays = Double.parseDouble(values[12]);
        this.species = species;
        this.coordOS = true;
    }

    @Override
    public String toString() {
        return this.getID() + " " + Arrays.toString(this.xy);
    }

    public ArrayList<Arrival> getArrivals() {
        return this.arrivals;
    }

    public void addArrival(Arrival arrival) {
        this.arrivals.add(arrival);
    }

    public void clearArrivals() {
        this.arrivals.clear();
    }

    public void setReleaseTime(double releaseTime) {
        this.releaseTime = releaseTime;
    }

    public double getReleaseTime() {
        return this.releaseTime;
    }

    public int getID() {
        return this.id;
    }

    public double[] getLocation() {
        return this.xy;
    }

    public double[] getStartLocation() {
        return this.startLoc;
    }

    public String getStartID() {
        return this.startSiteID;
    }

    public int getStartIndex() {
        return this.startSiteIndex;
    }

    public ISO_datestr getStartDate() {
        return this.startDate;
    }

    public double getStartTime() {
        return this.startTime;
    }

    public void setLocation(double x, double y) {
        this.xy[0] = x;
        this.xy[1] = y;
    }

    public void setElem(int elem) {
        this.elem = elem;
    }

    public int getElem() {
        return this.elem;
    }

    public void setROMSElemU(int[] elem) {
        this.elemRomsU = elem;
    }

    public void setROMSElemV(int[] elem) {
        this.elemRomsU = elem;
    }

    public int[] getROMSElemU() {
        return this.elemRomsU;
    }

    public int[] getROMSElemV() {
        return this.elemRomsV;
    }

    public void setROMSnearestPointU(int[] point) {
        this.nearestROMSGridPointU = point;
    }

    public void setROMSnearestPointV(int[] point) {
        this.nearestROMSGridPointV = point;
    }

    public int[] getROMSnearestPointU() {
        return this.nearestROMSGridPointU;
    }

    public int[] getROMSnearestPointV() {
        return this.nearestROMSGridPointV;
    }

    public void setMesh(int meshID) {
        this.mesh = meshID;
    }

    public int getMesh() {
        return this.mesh;
    }

    public String printLocation() {
        return "xy " + this.xy[0] + " " + this.xy[1] + " mesh " + this.mesh + " " + " elem " + this.elem;
    }

    public double[][] getNrList() {
        return this.nrList;
    }

    public void setCornerList(int elemPart, int[][] trinodes, double[][] nodexy) {
        this.cornerList[0][0] = nodexy[trinodes[elemPart][0]][0];
        this.cornerList[0][1] = nodexy[trinodes[elemPart][0]][1];
        this.cornerList[1][0] = nodexy[trinodes[elemPart][1]][0];
        this.cornerList[1][1] = nodexy[trinodes[elemPart][1]][1];
        this.cornerList[2][0] = nodexy[trinodes[elemPart][2]][0];
        this.cornerList[2][1] = nodexy[trinodes[elemPart][2]][1];
    }

    public double[][] getCornerList() {
        return this.cornerList;
    }

    public double getDepth() {
        return this.depth;
    }

    /**
     * Set depth of particle with check against maximum allowed depth
     */
    public void setDepth(double depth, double maxDepth) {
        if (depth > maxDepth) {
            depth = maxDepth;
        }
        if (depth < 0) {
            depth = 0;
        }
        this.depth = depth;
    }

    public int getDepthLayer() {
        return this.depLayer;
    }

    public void setDepthLayer(int depLayer) {
        this.depLayer = depLayer;
    }

    // see lower down for setMortRate (based on salinity)
    public double getMortRate() {
        return mortRate;
    }

    public void setDensity() {
        this.density = this.density * (1 - this.mortRate);
    }

    public double getDensity() {
        return this.density;
    }

    public void setStatus(int status) {
        // 0 = not free
        // 1 = free
        // 2 = viable (able to settle)
        // 3 = settled
        // 66 = boundary exit
        // 666 = dead (exceeded lifespan)
        this.status = status;
    }

    public int getStatus() {
        return this.status;
    }

    public void addX(double xTravel) {
        this.xTotal += xTravel;
    }

    public void addY(double yTravel) {
        this.yTotal += yTravel;
    }

    public void addXY(double xTravel, double yTravel) {
        this.xyTotal += Math.sqrt(xTravel*xTravel + yTravel*yTravel);
    }

    public void addZ(double zTravel) {
        this.zTotal += zTravel;
    }

    public double getxTotal() {
        return this.xTotal;
    }

    public double getyTotal() {
        return this.yTotal;
    }

    public double getxyTotal() {
        return this.xyTotal;
    }

    public double getzTotal() {
        return this.zTotal;
    }

    /**
     * put particle in the correct depth layer, based upon
     * its z position, and
     *
     * @param localDepth element bathymetry value
     */
    public void setLayerFromDepth(double localDepth, float[] layers) {
        int depNearest = 0;
        double dZmin = 1000;
        for (int i = 0; i < layers.length; i++) {
            if (Math.abs(this.depth - localDepth * layers[i]) < dZmin) {
                depNearest = i;
                dZmin = Math.abs(this.depth - localDepth * layers[i]);
            }
        }
        this.setDepthLayer(depNearest);
    }

    public int findLayerFromDepth(double depth, double localDepth, float[] layers) {
        int depNearest = 0;
        double dZmin = 1000;
        for (int i = 0; i < layers.length; i++) {
            if (Math.abs(depth - localDepth * layers[i]) < dZmin) {
                depNearest = i;
                dZmin = Math.abs(depth - localDepth * layers[i]);
            }
        }
        return depNearest;
    }

    public void setInfectious(boolean viable) {
        this.viable = viable;
    }

    public void setArrived(boolean arrived) {
        this.arrived = arrived;
    }

    public void setFree(boolean free) {
        this.free = free;
    }

    public void setBoundaryExit(boolean exit) {
        this.boundaryExit = exit;
    }

    public boolean isInfectious() {
        return this.viable;
    }

    public boolean hasArrived() {
        return this.arrived;
    }

    public boolean isFree() {
        return this.free;
    }

    public boolean hasExited() {
        return this.boundaryExit;
    }

    public boolean canBecomeInfectious(RunProperties rp) {
        boolean byAge = getAge() > rp.viabletime && rp.viabletime > 0;
        boolean byDegreeDays = getDegreeDays() > rp.viableDegreeDays && rp.viableDegreeDays > 0;
        boolean byPrevious = isInfectious();
        return byAge || byDegreeDays || byPrevious;
    }

    public boolean isTooOld(RunProperties rp) {
        boolean byAge = getAge() > rp.maxParticleAge && rp.maxParticleAge > 0;
        boolean byDegreeDays = getDegreeDays() > rp.maxDegreeDays && rp.maxDegreeDays > 0;
        boolean byPrevious = getStatus() == 666;
        boolean byDensity = getDensity() < 1e-15;
        return byAge || byDegreeDays || byPrevious || byDensity;
    }

    /**
     * Sets the mortality rate for the particle packet based upon local salinity
     */
    public void setMortRate(double salinity, RunProperties rp, double dt) {
        double mortRateHourly = switch (rp.mortSal_fn) {
            case "constant" -> rp.mortSal_b.get(0);
            case "quadratic" -> rp.mortSal_b.get(0) + rp.mortSal_b.get(1) * salinity + rp.mortSal_b.get(2) * salinity * salinity;
            // mort_h = b[0]/(1 + exp(-b[1] * (salinity - b[2]))) + b[3]
            case "logistic" -> rp.mortSal_b.get(0) / (1 + Math.exp(-rp.mortSal_b.get(1) * (salinity - rp.mortSal_b.get(2)))) + rp.mortSal_b.get(3);
            default -> rp.mortSal_b.get(0);
        };
        this.mortRate = 1 - Math.pow(1 - mortRateHourly, dt);
    }

    public void incrementAge(double increment) {
        this.age += increment;
    }

    public void incrementDegreeDays(double temp, RunProperties rp) {
        double inc = temp * (rp.dt / rp.stepsPerStep) / 86400;
        this.degreeDays += inc;
    }

    public double getAge() {
        return this.age;
    }

    public double getDegreeDays() {
        return this.degreeDays;
    }

    public double calcSwimDownProb(double salinity, RunProperties rp) {
        double prSwimDown = 0;
        double salThreshMin = 0;
        double salThreshMax = 0;
        if (this.isInfectious()) {
            salThreshMin = rp.salinityThreshCopepodidMin;
            salThreshMax = rp.salinityThreshCopepodidMax;
        } else {
            salThreshMin = rp.salinityThreshNaupliusMin;
            salThreshMax = rp.salinityThreshNaupliusMax;
        }
        if (salinity < salThreshMax & salinity > salThreshMin) {
            prSwimDown =  (salThreshMax - salinity) / (salThreshMax - salThreshMin);
        }
        if (salinity <= salThreshMin) {
            prSwimDown = 1;
        }
        return prSwimDown;
    }

    public double sink(RunProperties rp, double dt) {
        double swimDownSpeedMean = this.isInfectious() ? rp.swimDownSpeedCopepodidMean : rp.swimDownSpeedNaupliusMean;
        double swimDownSpeedStd = this.isInfectious() ? rp.swimDownSpeedCopepodidStd : rp.swimDownSpeedNaupliusStd;
        int forceDownward = rp.swimDownSpeedMean < 0 ? -1 : 1;
        return (forceDownward * swimDownSpeedMean + swimDownSpeedStd * ThreadLocalRandom.current().nextGaussian()) * dt;
    }

    public double swim(RunProperties rp, double dt) {
        double swimSpeedMean = this.isInfectious() ? rp.swimUpSpeedCopepodidMean : rp.swimUpSpeedNaupliusMean;
        double swimSpeedStd = this.isInfectious() ? rp.swimUpSpeedCopepodidStd : rp.swimUpSpeedNaupliusStd;
        int forceUpward = swimSpeedMean > 0 ? -1 : 1;
        return (forceUpward * swimSpeedMean + swimSpeedStd * ThreadLocalRandom.current().nextGaussian()) * dt;
    }

    public double swim(RunProperties rp, Mesh mesh, HydroField hydroField, int hour, double dt) {
        // following Johnsen et al. 2014, 2016, Myksvoll et al. 2018, Sandvik et al. 2020 for light attenuation and swimming thresholds
        // short_wave units = W m-2 ≈ 2.1 μmole m-2 s-1 (according to Tom Adams handover files)
        double lightAtSurface = 2.1 * hydroField.getAvgFromTrinodes(mesh, this.getLocation(), 0, this.elem, hour, "short_wave", rp);
        double lightAtDepth = lightAtSurface * Math.exp(-0.2 * this.depth);

        if ((this.isInfectious() && lightAtDepth > rp.lightThreshCopepodid) ||
                (!this.isInfectious() && lightAtDepth > rp.lightThreshNauplius)) {
            return swim(rp, dt);
        } else {
            return 0.0;
        }
    }

    public double passive(RunProperties rp, double localSalinity, double dt) {
        // based on Bricknell 2006 Fig. 3, assuming buoyancy is identical for copepodid & nauplius
        return (rp.passiveSinkingIntercept + rp.passiveSinkingSlope * localSalinity) * dt;
    }

    public double[] diffuse(double D_h, double K_gradient, double K_zAdj, double dt, String distribution) {
        double[] diffusion = {0,0,0};
        double[] randoms = {0,0,0};
        double r = 0.0;

        if (distribution.equals("uniform")) {
            for (int i = 0; i < 3; i++) {
                randoms[i] = ThreadLocalRandom.current().nextDouble(-1.0, 1.0);
            }
            r = 1.0 / 3.0;
        } else if (distribution.equals("gaussian")) {
            for (int i = 0; i < 3; i++) {
                randoms[i] = ThreadLocalRandom.current().nextGaussian();
            }
            r = 1.414;
        } else {
            System.err.println("Diffusion distribution must be uniform or gaussian.");
            System.exit(1);
        }
        diffusion[0] = randoms[0] * Math.sqrt(2 / r * D_h * dt);
        diffusion[1] = randoms[1] * Math.sqrt(2 / r * D_h * dt);

        // following Visser 1997 (eq. 6):
        // dZ = K_gradient * dt + Rand * sqrt(2 * dt/r * K_zAdj), where K_zAdj = K(z + K_gradient/2 * dt)
        // if !rp.variableDiffusion, K_gradient == 0, K_zAdj == Dh_vert
        diffusion[2] = K_gradient * dt + randoms[2] * Math.pow(2 / r * K_zAdj * dt, 0.5);

        return diffusion;
    }

    /**
     * Move the particle within the mesh, updating the particle location, mesh, and element.
     * In the case where meshes is length 2, particles are always preferentially in mesh 0. Particles can only exit
     * mesh 0 by passing through an open boundary, and are otherwise returned to their previous location if the previous
     * location was not near an open boundary. Returning from mesh 1 to mesh 0 does not require an open boundary. They
     * can similarly only exit mesh 1 (a full boundary exit) via an open boundary.
     *
     * @param newLoc Proposed particle location
     * @param meshes List of meshes
     * @param rp Run properties
     */
    public void moveInMesh(double[] newLoc, List<Mesh> meshes, RunProperties rp) {
        int oldMesh = this.getMesh(); // Current mesh index
        double[] oldLoc = this.getLocation(); // Current particle location
        int oldElem = this.getElem(); // Current element index
        int nMesh = meshes.size(); // Number of meshes

        if (oldMesh == 0) {
            // Handle movement within mesh 0
            handleMesh0(newLoc, meshes, rp, oldLoc, oldElem, nMesh);
        } else {
            // Handle movement within mesh 1
            handleMesh1(newLoc, meshes, rp, oldLoc, oldElem);
        }
    }

    /**
     * Handle particle movement within mesh 0.
     *
     * @param newLoc Proposed particle location
     * @param meshes List of meshes
     * @param rp Run properties
     * @param oldLoc Current particle location
     * @param oldElem Current element index
     * @param nMesh Number of meshes
     */
    private void handleMesh0(double[] newLoc, List<Mesh> meshes, RunProperties rp, double[] oldLoc, int oldElem, int nMesh) {
        int [] elemMesh0 = findContainingElement(newLoc, oldElem, meshes.get(0), false, 20);
        if (elemMesh0[0] > -1) {
            // New location is within mesh 0: update position
            updateLocation(newLoc, elemMesh0[0], meshes, 0, false);
        } else {
            // New location is outside mesh 0: check if near an open boundary
            boolean isNearBoundary = ParallelParticleMover.openBoundaryCheck(oldElem, meshes.get(0).getOpenBoundaryElems(), rp);

            if (nMesh == 1 && isNearBoundary) {
                // Near an open boundary and only 1 mesh: exit domain
                exitDomain();
            } else if (nMesh == 2 && isNearBoundary) {
                // Near an open boundary and 2 meshes: handle boundary transition from mesh 0 to 1
                handleMeshBoundary(newLoc, meshes);
            } else {
                // Not near an open boundary: bounce back to previous location
                bounceBack();
            }
        }
    }

    /**
     * Handle particle movement within mesh 1.
     *
     * @param newLoc Proposed particle displacement
     * @param meshes List of meshes
     * @param rp Run properties
     * @param oldLoc Current particle location
     * @param oldElem Current element index
     */
    private void handleMesh1(double[] newLoc, List<Mesh> meshes, RunProperties rp, double[] oldLoc, int oldElem) {
        int[] elemMesh0 = findContainingElement(newLoc, meshes.get(0).getAdjoiningElement(), meshes.get(0), false, 20);
        if (elemMesh0[0] > -1) {
            // New location is within mesh 0: update position
            updateLocation(newLoc, elemMesh0[0], meshes, 0, true);
        } else {
            int[] elemMesh1 = findContainingElement(newLoc, oldElem, meshes.get(1), false, 20);
            if (elemMesh1[0] > -1) {
                // New location is within mesh 1: update position
                updateLocation(newLoc, elemMesh1[0], meshes, 1, false);
            } else {
                // New location is outside both meshes: check if near an open boundary
                if (ParallelParticleMover.openBoundaryCheck(oldElem, meshes.get(1).getOpenBoundaryElems(), rp)) {
                    // Near an open boundary: exit domain
                    exitDomain();
                } else {
                    // Not near an open boundary: bounce back to previous location
                    bounceBack();
                }
            }
        }
    }

    /**
     * Update the particle location and place it in the specified mesh.
     *
     * @param newLoc New particle location
     * @param meshes List of meshes
     * @param meshIndex Index of the mesh to place the particle in
     * @param switchedMesh Flag indicating if the particle is moving from mesh 1
     */
    private void updateLocation(double[] newLoc, int newElem, List<Mesh> meshes, int meshIndex, boolean switchedMesh) {
        this.setLocation(newLoc[0], newLoc[1]);
        this.setMesh(meshIndex);
        this.setElem(newElem);
    }

    /**
     * Set the particle to exit the domain.
     */
    private void exitDomain() {
        this.setBoundaryExit(true);
        this.setStatus(66);
    }

    /**
     * Handle particle transition from mesh 0 to mesh 1.
     *
     * @param newLoc New particle location
     * @param meshes List of meshes
     */
    private void handleMeshBoundary(double[] newLoc, List<Mesh> meshes) {
        int[] elemMesh1 = findContainingElement(newLoc, meshes.get(1).getAdjoiningElement(), meshes.get(1), true, 40);
        if (elemMesh1[0] > -1) {
            // New location is within mesh 1: update position
            updateLocation(newLoc, elemMesh1[0], meshes, 1, true);
        } else {
            // Not within mesh 1: bounce back to previous location
            bounceBack();
        }
    }

    /**
     * Bounce the particle back to its previous location.
     *
     */
    private void bounceBack() {

    }


    /**
     * This method updates the particle mesh neighbourhood information, dealing with a change in mesh ID if required
     *
     * @param meshes            The mesh that has been shown to contain a particle
     * @param switchedMesh Logical, has the particle changed mesh?
     */
    public void placeInMesh(List<Mesh> meshes, int id, boolean switchedMesh) {
        Mesh m = meshes.get(id);
        this.setMesh(id);

        if (m.getType().equals("FVCOM") || m.getType().equals("ROMS_TRI")) {
            int el = 0;
            boolean checkAll = true;
            if (!switchedMesh) {
                el = this.getElem();
                //checkAll = false;
            }
            int[] c = findContainingElement(this.getLocation(), el, m,true, 40);
            // if particle is within the mesh, update location normally and save the distance travelled
            this.setElem(c[0]);

        } else if (m.getType().equals("ROMS")) {
            int[] searchCentreU = null, searchCentreV = null;
            if (!switchedMesh) {
                searchCentreU = this.getROMSnearestPointU();
                searchCentreV = this.getROMSnearestPointV();
            }
            int[] nearU = nearestROMSGridPoint((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonU(), m.getLatU(), searchCentreU);
            int[] nearV = nearestROMSGridPoint((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonV(), m.getLatV(), searchCentreV);

            // More to do here to turn the nearest grid point into the containing element
            int[] containingROMSElemU = whichROMSElement((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonU(), m.getLatU(),
                    nearU);
            int[] containingROMSElemV = whichROMSElement((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonV(), m.getLatV(),
                    nearV);

            // Need to save the nearest point, in order to read into
            this.setROMSnearestPointU(nearU);
            this.setROMSnearestPointV(nearV);
            this.setROMSElemU(containingROMSElemU);
            this.setROMSElemV(containingROMSElemV);
        }
    }


    /**
     * Find the nearest mesh element centroid
     */
    public static int nearestCentroid(double x, double y, float[][] uvnode) {
        int nearest = -1;
        double dist = 10000000;

        for (int i = 0; i < uvnode[0].length; i++) {
            double distnew = Math.sqrt((x - uvnode[0][i]) * (x - uvnode[0][i]) + (y - uvnode[1][i]) * (y - uvnode[1][i]));
            if (distnew < dist) {
                dist = distnew;
                nearest = i;
            }
        }
        //System.out.printf("In Particle.nearestCentroid "+nearest+"\n");
        return nearest;
    }

    /**
     * Make a list of the nearest mesh element centroids
     */
    public static double[][] nearestCentroidList(double x, double y, float[][] uvnode) {
        double[][] nearestList = new double[5][2];
        int nearest = -1;
        double dist = 10000000;

        for (int i = 0; i < uvnode[0].length; i++) {
            double distnew = Math.sqrt((x - uvnode[0][i]) * (x - uvnode[0][i]) + (y - uvnode[1][i]) * (y - uvnode[1][i]));
            if (distnew < dist) {
                dist = distnew;
                // Shift everything along one element
                nearestList[4][0] = nearestList[3][0];
                nearestList[4][1] = nearestList[3][1];
                nearestList[3][0] = nearestList[2][0];
                nearestList[3][1] = nearestList[2][1];
                nearestList[2][0] = nearestList[1][0];
                nearestList[2][1] = nearestList[1][1];
                nearestList[1][0] = nearestList[0][0];
                nearestList[1][1] = nearestList[0][1];

                nearestList[0][0] = i;
                nearestList[0][1] = dist;
            }
        }
        return nearestList;
    }

    /**
     * Look for the nearest point in a grid of points.
     * This method is somewhat naive in that it searches the entire grid. This
     * means that it can handle grids which are not oriented on a strict N/S
     * lattice, but is inefficient.
     */
    public static int[] nearestROMSGridPoint(float x, float y, float[][] xGrid, float[][] yGrid, int[] searchCentre) {
        int[] nearest = new int[]{-1, -1};
        float dist = 10000000;

        if (xGrid.length != yGrid.length || xGrid[0].length != yGrid[0].length) {
            System.err.println("Particle.nearestROMSGridPoint: ROMS x and y grids are not the same size");
        }

        int minX = 0;
        int maxX = xGrid.length;
        int minY = 0;
        int maxY = xGrid[0].length;

        if (searchCentre != null) {
            minX = Math.max(0, searchCentre[0] - 3);
            maxX = Math.min(xGrid.length, searchCentre[0] + 4);
            minY = Math.max(0, searchCentre[1] - 3);
            maxY = Math.min(xGrid[0].length, searchCentre[1] + 4);
        }

        for (int i = minX; i < maxX; i++) {
            for (int j = minY; j < maxY; j++) {
                float distnew = (float) distanceEuclid2(x, y, xGrid[i][j], yGrid[i][j], false);
                if (distnew < dist) {
                    dist = distnew;
                    nearest = new int[]{i, j};
                }
            }
        }
        return nearest;
    }

    /**
     * This makes a list of the nearest ROMS grid points - 4 of them.
     * Note that this might not actually be what we want, due to the form of the grid.
     * Instead, we probably want the 4 corners of the containing element - see nearestListROMS2.
     */
    @SuppressWarnings("Convert2Lambda")
    public static double[][] nearestListROMS(double[] xy, float[][] xGrid, float[][] yGrid, int[] nearestPoint) {
        double[][] allDists = new double[9][3];
        double[][] nearList = new double[5][3];

        if (nearestPoint == null) {
            nearestPoint = nearestROMSGridPoint((float) xy[0], (float) xy[1], xGrid, yGrid, null);
        }

        // We should now be able to define the 9 grid points closest to the location
        float dist = 10000000;
        int k = 0;
        for (int i = nearestPoint[0] - 1; i <= nearestPoint[0] + 1; i++) {
            for (int j = nearestPoint[1] - 1; j <= nearestPoint[1] + 1; j++) {
                float distnew = (float) distanceEuclid2(xy[0], xy[1], xGrid[i][j], yGrid[i][j], false);
                allDists[k][0] = i;
                allDists[k][1] = j;
                allDists[k][2] = distnew;
                k++;
            }
        }

        // noinspection Convert2Diamond,Convert2Lambda
        Arrays.sort(allDists, new Comparator<double[]>() {
            @Override
            public int compare(double[] o1, double[] o2) {
                return Double.compare(o1[2], o2[2]);
            }
        });

        int nPoints = 4;
        System.arraycopy(allDists, 0, nearList, 0, nPoints);
        return nearList;
    }


    public static double[][] nearestListROMS2(double[] xy, float[][] xGrid, float[][] yGrid, int[] nearestPoint) {
        double[][] nearList = new double[5][3];

        if (nearestPoint == null) {
            nearestPoint = nearestROMSGridPoint((float) xy[0], (float) xy[1], xGrid, yGrid, null);
        }

        // Get the index of the "top-left" corner of the containing element
        int[] whichElem = whichROMSElement((float) xy[0], (float) xy[1], xGrid, yGrid, nearestPoint);
        if (whichElem[0] == -1) {
            nearestPoint = nearestROMSGridPoint((float) xy[0], (float) xy[1], xGrid, yGrid, null);
            whichElem = whichROMSElement((float) xy[0], (float) xy[1], xGrid, yGrid, nearestPoint);
        }

        // Get the distances to the corners of that element
        nearList[0] = new double[]{whichElem[0], whichElem[1],
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0]][whichElem[1]], yGrid[whichElem[0]][whichElem[1]], false)};
        nearList[1] = new double[]{whichElem[0] + 1, whichElem[1],
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0] + 1][whichElem[1]], yGrid[whichElem[0] + 1][whichElem[1]], false)};
        nearList[2] = new double[]{whichElem[0], whichElem[1] + 1,
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0]][whichElem[1] + 1], yGrid[whichElem[0]][whichElem[1] + 1], false)};
        nearList[3] = new double[]{whichElem[0] + 1, whichElem[1] + 1,
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0] + 1][whichElem[1] + 1], yGrid[whichElem[0] + 1][whichElem[1] + 1], false)};
        nearList[4] = new double[]{0, 0, 10000000};

        return nearList;
    }


    /**
     * Identify the ROMS element (defined as four neighbouring points forming a
     * parallelogram) containing a particular coordinate.
     * Given the nearest UV point, the coordinate is in one of four parallelograms,
     * which share that point as a common corner.
     * <p>
     * The output index of which element is defined by the last xGrid,yGrid index
     * occurring
     *
     * @param nearP the nearest ROMS grid point
     */
    public static int[] whichROMSElement(float x, float y, float[][] lon, float[][] lat, int[] nearP) {
        int[] which = new int[]{-1, -1};

        // Look at each of the four possible boxes in turn.
        // In the Marine Institute ROMS grids, U/V points are arranged in a matrix starting
        // in the NW corner.

        // In arrays read in by MATLAB from file:
        // Increasing row (first) index means decreasing latitude.
        // Increasing column (second) index means increasing longitude.

        // In arrays read in by JAVA from file, arrays are transposed:
        // Increasing row (first) index means increasing longitude.
        // Increasing column (second) index means decreasing latitude.

        // "Top left" (in context of lon_u/lat_u index)
        if (nearP[0] != 0 && nearP[1] != 0) {
            float[][] box = new float[4][2];
            // create corners of box from top left corner, going anticlockwise
            box[0] = new float[]{lon[nearP[0] - 1][nearP[1] - 1], lat[nearP[0] - 1][nearP[1] - 1]};
            box[1] = new float[]{lon[nearP[0] - 1][nearP[1]], lat[nearP[0] - 1][nearP[1]]};
            box[2] = new float[]{lon[nearP[0]][nearP[1]], lat[nearP[0]][nearP[1]]};
            box[3] = new float[]{lon[nearP[0]][nearP[1] - 1], lat[nearP[0]][nearP[1] - 1]};
            Path2D.Float boxPath = Mesh.pointsToPath(box);
            // Return the top left corner
            if (boxPath.contains(x, y)) {
                which = new int[]{nearP[0] - 1, nearP[1] - 1};
            }
        }
        // "Top right"
        if (nearP[1] != 0) {
            float[][] box = new float[4][2];
            box[0] = new float[]{lon[nearP[0]][nearP[1] - 1], lat[nearP[0]][nearP[1] - 1]};
            box[1] = new float[]{lon[nearP[0]][nearP[1]], lat[nearP[0]][nearP[1]]};
            box[2] = new float[]{lon[nearP[0] + 1][nearP[1]], lat[nearP[0] + 1][nearP[1]]};
            box[3] = new float[]{lon[nearP[0] + 1][nearP[1] - 1], lat[nearP[0] + 1][nearP[1] - 1]};
            Path2D.Float boxPath = Mesh.pointsToPath(box);
            if (boxPath.contains(x, y)) {
                which = new int[]{nearP[0], nearP[1] - 1};
            }
        }
        // "Bottom right"
        if (nearP[0] != lon.length && nearP[1] != lon[0].length) {
            float[][] box = new float[4][2];
            box[0] = new float[]{lon[nearP[0]][nearP[1]], lat[nearP[0]][nearP[1]]};
            box[1] = new float[]{lon[nearP[0]][nearP[1] + 1], lat[nearP[0]][nearP[1] + 1]};
            box[2] = new float[]{lon[nearP[0] + 1][nearP[1] + 1], lat[nearP[0] + 1][nearP[1] + 1]};
            box[3] = new float[]{lon[nearP[0] + 1][nearP[1]], lat[nearP[0] + 1][nearP[1]]};
            Path2D.Float boxPath = Mesh.pointsToPath(box);
            if (boxPath.contains(x, y)) {
                which = new int[]{nearP[0], nearP[1]};
            }
        }
        // "Bottom left"
        if (nearP[0] != 0 && nearP[1] != 0) {
            float[][] box = new float[4][2];
            box[0] = new float[]{lon[nearP[0] - 1][nearP[1]], lat[nearP[0] - 1][nearP[1]]};
            box[1] = new float[]{lon[nearP[0] - 1][nearP[1] + 1], lat[nearP[0] - 1][nearP[1] + 1]};
            box[2] = new float[]{lon[nearP[0]][nearP[1] + 1], lat[nearP[0]][nearP[1] + 1]};
            box[3] = new float[]{lon[nearP[0]][nearP[1]], lat[nearP[0]][nearP[1]]};
            Path2D.Float boxPath = Mesh.pointsToPath(box);
            if (boxPath.contains(x, y)) {
                which = new int[]{nearP[0] - 1, nearP[1]};
            }
        }

        return which;
    }


    public static double[] velocityFromNearestList(double[][] nrList, float[][] uHour, float[][] vHour, float[][] wHour) {
        double[] velocity = {0,0,0};
        double usum = 0, vsum = 0, wsum = 0, sum = 0;
        for (int i = 0; i < nrList.length; i++) {
            int elem = (int) nrList[i][0];
            double distance = nrList[i][1];
            double weight = distance == 0 ? 1 : 1.0 / (distance * distance);
            int depLayer = (int) nrList[i][2];
            wsum += weight * wHour[depLayer][elem];
            usum += weight * uHour[depLayer][elem];
            vsum += weight * vHour[depLayer][elem];
            sum += weight;
        }
        velocity[0] = usum / sum;
        velocity[1] = vsum / sum;
        velocity[2] = - wsum / sum; // upward velocity; ww * -1 so that positive values = increased depth
        return velocity;
    }

    /**
     * A method to use the nearest grid point indexes and distances to interpolate a velocity.
     * This is done between the four the nearest values.
     * In the ROMS output, U and V are stored on different grids (Arakawa-C). Therefore,
     * two nrLists are read in, one for each grid.
     * <p>
     * This method of calculating velocity DOES NOT match the methods used in e.g. TRACMASS, LTRANS.
     * Those programs interpolate linearly between velocity values on opposing edges of model elements
     * ("horizontally" for U, "vertically" for V), adjusted according to latitude.
     */
    public static double[] velocityFromNearestListROMS(double[][] nrListU, double[][] nrListV, int hour, float[][][] u, float[][][] v) {
        double[] velocity = new double[2];
        double[] weightsU = new double[nrListU.length];
        double[] weightsV = new double[nrListV.length];
        double usum = 0, vsum = 0, sumWeightsU = 0, sumWeightsV = 0;
        // Just use the first four entries
        for (int i = 0; i < nrListU.length; i++) {
            if (nrListU[i][2] != 0) {
                weightsU[i] = 1.0 / (nrListU[i][2] * nrListU[i][2]);
            } else {
                weightsU[i] = 1;
            }
            if (nrListV[i][2] != 0) {
                weightsV[i] = 1.0 / (nrListV[i][2] * nrListV[i][2]);
            } else {
                weightsV[i] = 1;
            }
            usum = usum + weightsU[i] * u[hour][(int) nrListU[i][0]][(int) nrListU[i][1]];
            vsum = vsum + weightsV[i] * v[hour][(int) nrListV[i][0]][(int) nrListV[i][1]];
            sumWeightsU = sumWeightsU + weightsU[i];
            sumWeightsV = sumWeightsV + weightsV[i];
        }
        usum = usum / sumWeightsU;
        vsum = vsum / sumWeightsV;
        velocity[0] = usum;
        velocity[1] = vsum;
        return velocity;
    }


    /**
     * Euclidean distance between two points
     */
    public static double distanceEuclid(double x1, double y1, double x2, double y2) {
        return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
    }

    public static double distanceEuclid(double x1, double y1, double z1, double x2, double y2, double z2) {
        return Math.sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
    }

    /**
     * Compute the Euclidean distance (in metres) between two points.
     * If the supplied points are WGS84 coordinates, convert their x,y separation
     * to metres prior to calculating the Euclidean distance.
     */
    public static double distanceEuclid2(double x1, double y1, double x2, double y2, boolean coordOS) {
        double dx = x1 - x2;
        double dy = y1 - y2;

        if (!coordOS) {
            double[] distXY = new double[]{dx, dy};
            distXY = ParallelParticleMover.distanceDegreesToMetres(distXY, new double[]{x1, y2});
            return Math.sqrt(distXY[0] * distXY[0] + distXY[1] * distXY[1]);
        }

        return Math.sqrt(dx * dx + dy * dy);
    }

    public static double distanceEuclid2(double x1, double y1, double z1, double x2, double y2, double z2, boolean coordOS) {
        double dx = x1 - x2;
        double dy = y1 - y2;
        double dz = z1 - z2;

        if (!coordOS) {
            double[] distXY = new double[]{dx, dy};
            distXY = ParallelParticleMover.distanceDegreesToMetres(distXY, new double[]{x1, y2});
            return Math.sqrt(distXY[0] * distXY[0] + distXY[1] * distXY[1] + dz * dz);
        }

        return Math.sqrt(dx * dx + dy * dy + dz * dz);
    }


    /**
     * Find which element a particle resides within (edge checking).
     *
     * @param xy Coordinates of the particle
     * @param elems Elements to check
     * @param nodexy Coordinates of the nodes in the mesh
     * @param trinodes Nodes that make up each triangular element
     * @return Index of the element containing the particle, or -1 if not found
     */
    public static int whichElement(double[] xy, int[] elems, float[][] nodexy, int[][] trinodes) {
        for (int elem : elems) {
            double[] xt = new double[3];
            double[] yt = new double[3];

            for (int j = 0; j < 3; j++) {
                xt[j] = nodexy[0][trinodes[j][elem]];
                yt[j] = nodexy[1][trinodes[j][elem]];
            }

            if (isPointInBoundingBox(xy, xt, yt) && isPointInTriangle(xy, xt, yt)) {
                return elem;
            }
        }
        return -1;
    }

    /**
     * Check whether a point lies within the bounding box of a triangle.
     *
     * @param xy Coordinates of the point
     * @param xt X-coordinates of the triangle vertices
     * @param yt Y-coordinates of the triangle vertices
     * @return True if the point lies within the bounding box, otherwise false
     */
    private static boolean isPointInBoundingBox(double[] xy, double[] xt, double[] yt) {
        return xy[0] >= Math.min(Math.min(xt[0], xt[1]), xt[2]) &&
                xy[0] <= Math.max(Math.max(xt[0], xt[1]), xt[2]) &&
                xy[1] >= Math.min(Math.min(yt[0], yt[1]), yt[2]) &&
                xy[1] <= Math.max(Math.max(yt[0], yt[1]), yt[2]);
    }

    /**
     * Check whether a point lies within a triangle.
     *
     * @param xy Coordinates of the point
     * @param xt X-coordinates of the triangle vertices
     * @param yt Y-coordinates of the triangle vertices
     * @return True if the point lies within the triangle, otherwise false
     */
    private static boolean isPointInTriangle(double[] xy, double[] xt, double[] yt) {
        double f1 = (xy[1] - yt[0]) * (xt[1] - xt[0]) - (xy[0] - xt[0]) * (yt[1] - yt[0]);
        double f2 = (xy[1] - yt[2]) * (xt[0] - xt[2]) - (xy[0] - xt[2]) * (yt[0] - yt[2]);
        double f3 = (xy[1] - yt[1]) * (xt[2] - xt[1]) - (xy[0] - xt[1]) * (yt[2] - yt[1]);
        return f1 * f3 >= 0.0 && f3 * f2 >= 0.0;
    }

    /**
     * Search through lists of elements progressively further away from the last known
     * location in order to find the new location.
     * <p>
     * Returned array contains the element in which the particle is determined to be located,
     * plus the number of counts required at each scale (for diagnostics).
     *
     * @param xy Coordinates of the particle
     * @param elemPart Initial guess of the element containing the particle
     * @param m Mesh
     * @param numRounds Number of rounds of neighbours to search before giving up
     * @return Array containing the element index where the particle is located and diagnostic counts for each search scale
     */
    public static int[] findContainingElement(double[] xy, int elemPart, Mesh m, boolean checkAll, int numRounds) {
        int[] c = new int[numRounds+1];
        int[] elemInit = {elemPart};

        if (!m.isInMesh(xy)) {
            c[0] = -1;
            return c;
        }

        float[][] nodexy = m.getNodexy();
        int[][] trinodes = m.getTrinodes();
        int[][] neighbours = m.getNeighbours();

        int whereami = whichElement(xy, elemInit, nodexy, trinodes);
        c[1] = 1;

        int round = 2;
        int [] elemsCheck = elemInit;
        Set<Integer> elemsChecked = new HashSet<>();
        elemsChecked.add(elemPart);

        while (whereami == -1) {
            if (round <= numRounds) {
                c[round] = 1;
                elemsCheck = getNeighbours(elemsCheck, neighbours, elemsChecked);
                whereami = whichElement(xy, elemsCheck, nodexy, trinodes);
                for (int elem : elemsCheck) {
                    elemsChecked.add(elem);
                }
                round++;
            } else {
                if (checkAll) {
                    int[] allRemainingElems = IntStream.rangeClosed(0, trinodes[0].length - 1)
                            .filter(elem -> !elemsChecked.contains(elem))
                            .toArray();
                    whereami = whichElement(xy, allRemainingElems, nodexy, trinodes);
                }
                break;
            }
        }

        c[0] = whereami;
        return c;
    }

    /**
     * Get the immediate neighbours of an array of elements, excluding already checked elements.
     *
     * @param elems Array of elements
     * @param neighbours Neighboring elements for each element
     * @param elemsChecked Set of already checked elements
     * @return Array of neighbouring elements that have not been checked
     */
//    private static int[] getNeighbours(int[] elems, int[][] neighbours, Set<Integer> elemsChecked) {
//        return Arrays.stream(elems)
//                .flatMap(elem -> Arrays.stream(neighbours).mapToInt(neigh -> neigh[elem]))
//                .filter(neigh -> !elemsChecked.contains(neigh) && neigh != 0)
//                .distinct()
//                .toArray();
//    }

    private static int[] getNeighbours(int[] elems, int[][] neighbours, Set<Integer> elemsChecked) {
        Set<Integer> neighboursSet = new HashSet<>();

        for (int elem : elems) {
            for (int[] neighbourArray : neighbours) {
                int neigh = neighbourArray[elem];
                if (neigh != 0 && !elemsChecked.contains(neigh)) {
                    neighboursSet.add(neigh);
                }
            }
        }

        int[] result = new int[neighboursSet.size()];
        int index = 0;
        for (int neigh : neighboursSet) {
            result[index++] = neigh;
        }

        return result;
    }


    /**
     * Calculate a list of neighbour cells based on a specified (e.g. particle) location
     * (0: containing element, 1-3: neighbours of the containing element), and calculate
     * the Euclidean distances.
     */
    public static double[][] neighbourCellsList(double[] xy, int elemPart0, Mesh m, boolean coordOS) {
        double[][] nrList = new double[5][2];

        int[] elem = findContainingElement(xy, elemPart0, m, true, 10);
        // If particle is not within the mesh (value returned by findContainingElement = -1)
        // exit this method returning array of zeros.
        if (elem[0] == -1) {
            return nrList;
        }
        int thisElem = elem[0];
        nrList[0][0] = thisElem;
        nrList[0][1] = distanceEuclid2(xy[0], xy[1], m.getUvnode()[0][thisElem], m.getUvnode()[1][thisElem], coordOS);
        // distance to neighbouring elems
        nrList[1][0] = m.getNeighbours()[0][thisElem];
        nrList[1][1] = distanceEuclid2(xy[0], xy[1], m.getUvnode()[0][m.getNeighbours()[0][thisElem]], m.getUvnode()[1][m.getNeighbours()[0][thisElem]], coordOS);
        nrList[2][0] = m.getNeighbours()[1][thisElem];
        nrList[2][1] = distanceEuclid2(xy[0], xy[1], m.getUvnode()[0][m.getNeighbours()[1][thisElem]], m.getUvnode()[1][m.getNeighbours()[1][thisElem]], coordOS);
        nrList[3][0] = m.getNeighbours()[2][thisElem];
        nrList[3][1] = distanceEuclid2(xy[0], xy[1], m.getUvnode()[0][m.getNeighbours()[2][thisElem]], m.getUvnode()[1][m.getNeighbours()[2][thisElem]], coordOS);
        nrList[4][0] = 0;
        nrList[4][1] = 1000000;

        return nrList;
    }

    /**
     * Find the containing element plus 3 neighbors for the depth layer above and below the particle and calculate the
     * euclidean distances from the particle to the centroids
     *
     * @param xy particle location
     * @param elemPart particle element
     * @param depth particle depth
     * @param coordOS use OS coordinate system: EPSG 27700
     * @return double[8][3] with rows = elements and columns = [elementID, distance to particle (3D), depth layer]
     */
    public static double[][] neighbourCellsList3D(double[] xy, int elemPart, int depLayer, float[][] uvnode,
                                                  float[] depthUvnode, float[][] nodexy, int[][] trinodes,
                                                  float[] siglayers, int[][] neighbours, double depth, boolean coordOS) {

        double[][] nrList = new double[8][3]; // (3 neighbors per layer + containing element) * 2 layers = 8

        if (elemPart == -1) {
            return nrList;
        }

        int[][] nrLayerIndexes = new int[4][2];

        // layer 0 = below particle, layer 1 = above particle
        nrLayerIndexes[0] = Mesh.findNearestSigmaIndexes(depLayer, siglayers.length);
        for (int nbr = 0; nbr < 3; nbr++) {
            nrLayerIndexes[nbr+1] = Mesh.findNearestSigmaIndexes(depth, siglayers, depthUvnode[neighbours[nbr][elemPart]]);
        }

        // layer 0: below particle
        // containing element
        nrList[0][0] = elemPart;
        nrList[0][1] = distanceEuclid2(xy[0], xy[1], depth,
                uvnode[0][elemPart], uvnode[1][elemPart], depthUvnode[elemPart] * siglayers[nrLayerIndexes[0][0]], coordOS);
        nrList[0][2] = nrLayerIndexes[0][0];
        // 3 neighbours
        for (int nbr = 0; nbr < 3; nbr++) {
            int nbrElem = neighbours[nbr][elemPart];
            nrList[nbr+1][0] = nbrElem;
            nrList[nbr+1][1] = distanceEuclid2(xy[0], xy[1], depth,
                    uvnode[0][nbrElem], uvnode[1][nbrElem], depthUvnode[nbrElem] * siglayers[nrLayerIndexes[nbr+1][0]], coordOS);
            nrList[nbr+1][2] = nrLayerIndexes[nbr+1][0];
        }

        // layer 1: above particle
        // copy layer 0, then update if layer 0 != layer 1  --  faster than always re-calculating
        for (int i = 4; i < 8; i++) {
            nrList[i][0] = nrList[i-4][0];
            nrList[i][1] = nrList[i-4][1];
            nrList[i][2] = nrList[i-4][2];
        }
        // containing element
        if (nrLayerIndexes[0][0] != nrLayerIndexes[0][1]) {
            nrList[4][0] = elemPart;
            nrList[4][1] = distanceEuclid2(xy[0], xy[1], depth, uvnode[0][elemPart], uvnode[1][elemPart], depthUvnode[elemPart] * siglayers[nrLayerIndexes[0][1]], coordOS);
            nrList[4][2] = nrLayerIndexes[0][1];
        }
        // 3 neighbours
        for (int nbr = 0; nbr < 3; nbr++) {
            int nbrElem = neighbours[nbr][elemPart];
            if (nrLayerIndexes[nbr+1][0] != nrLayerIndexes[nbr+1][1]) {
                nrList[nbr+5][0] = nbrElem;
                nrList[nbr+5][1] = distanceEuclid2(xy[0], xy[1], depth,
                        uvnode[0][nbrElem], uvnode[1][nbrElem], depthUvnode[nbrElem] * siglayers[nrLayerIndexes[nbr+1][1]], coordOS);
                nrList[nbr+5][2] = nrLayerIndexes[nbr+1][1];
            }
        }
        return nrList;
    }

    /**
     * Update particle location using an RK4 integration step.
     */
    public double[] rk4Step(List<HydroField> hydroFields, // velocities
                            List<Mesh> meshes,      // other mesh info
                            int hour, int step, double dt,           // locate particle in space and time
                            int stepsPerStep,  // info on simulation length
                            boolean coordOS, boolean FVCOM) {
        int elemPart = this.getElem();
        int meshPart = this.getMesh();
        int depLayer = this.getDepthLayer();
        float[][] uvnode = meshes.get(meshPart).getUvnode();
        float[] depthUvnode = meshes.get(meshPart).getDepthUvnode();
        float[][] nodexy = meshes.get(meshPart).getNodexy();
        int[][] trinodes = meshes.get(meshPart).getTrinodes();
        float[] siglayers = meshes.get(meshPart).getSiglay();
        int[][] neighbours = meshes.get(meshPart).getNeighbours();
        float[][][] u = hydroFields.get(meshPart).getU();
        float[][][] v = hydroFields.get(meshPart).getV();
        float[][][] w = hydroFields.get(meshPart).getW();
        double[] advectStep = {0,0,0};

        // 1. Compute k_1 (initial spatial interpolation)
        double[] k1 = stepAhead(
                this.getLocation(), this.getDepth(),
                elemPart, depLayer, 0,
                uvnode, depthUvnode, nodexy, trinodes, siglayers, neighbours,
                u, v, w,
                hour, step, dt, stepsPerStep, coordOS, FVCOM
        );
        double[] k1Deg = new double[]{k1[0], k1[1]};
        if (!this.coordOS) {
            k1Deg = ParallelParticleMover.distanceMetresToDegrees2(k1Deg, this.getLocation());
        }
        double[] k1Loc = new double[]{this.getLocation()[0] + k1Deg[0] / 2.0, this.getLocation()[1] + k1Deg[1] / 2.0};
        double k1Depth = this.getDepth() + k1[2] / 2.0;  // k1[2] = 0 if fixDepth = true
        int k1ElemPart = findContainingElement(k1Loc, elemPart, meshes.get(meshPart), false, 10)[0];
        elemPart = k1ElemPart == -1 ? elemPart : k1ElemPart; // particles can't exit mesh during substeps
        depLayer = findLayerFromDepth(k1Depth, depthUvnode[elemPart], siglayers);

        // 2. Compute k_2 (spatial interpolation at half step, temporal interpolation at half step)
        // Estimated half-step location using Euler
        // NOTE that here, for the purposes of identifying the elements containing the part-step locations,
        // the steps in degrees are calculated. These are separate of the actual half-step values in metres
        // which are retained and summed at the end of the method to give a transport distance in metres
        double[] k2 = stepAhead(
                k1Loc, k1Depth,
                elemPart, depLayer, 1.0 / 2.0,
                uvnode, depthUvnode, nodexy, trinodes, siglayers, neighbours,
                u, v, w,
                hour, step, dt, stepsPerStep, coordOS, FVCOM);
        double[] k2Deg = new double[]{k2[0], k2[1]};
        if (!this.coordOS) {
            k2Deg = ParallelParticleMover.distanceMetresToDegrees2(k2Deg, this.getLocation());
        }
        double[] k2Loc = new double[]{this.getLocation()[0] + k2Deg[0] / 2.0, this.getLocation()[1] + k2Deg[1] / 2.0};
        double k2Depth = this.getDepth() + k2[2] / 2.0;  // k2[2] = 0 if fixDepth = true
        int k2ElemPart = findContainingElement(k1Loc, elemPart, meshes.get(meshPart), false, 10)[0];
        elemPart = k2ElemPart == -1 ? elemPart : k2ElemPart; // particles can't exit mesh during substeps
        depLayer = findLayerFromDepth(k2Depth, depthUvnode[elemPart], siglayers);

        // 3. Compute k_3 (spatial interpolation at half step, temporal interpolation at half step)
        double[] k3 = stepAhead(
                k2Loc, k2Depth,
                elemPart, depLayer, 1.0 / 2.0,
                uvnode, depthUvnode, nodexy, trinodes, siglayers, neighbours,
                u, v, w,
                hour, step, dt, stepsPerStep, coordOS, FVCOM);
        double[] k3Deg = new double[]{k3[0], k3[1]};
        if (!this.coordOS) {
            k3Deg = ParallelParticleMover.distanceMetresToDegrees2(k3Deg, this.getLocation());
        }
        double[] k3Loc = new double[]{this.getLocation()[0] + k3Deg[0], this.getLocation()[1] + k3Deg[1]};
        double k3Depth = this.getDepth() + k3[2];  // k3[2] = 0 if fixDepth = true
        int k3ElemPart = findContainingElement(k3Loc, elemPart, meshes.get(meshPart), false, 10)[0];
        elemPart = k3ElemPart == -1 ? elemPart : k3ElemPart; // particles can't exit mesh during substeps
        depLayer = findLayerFromDepth(k3Depth, depthUvnode[elemPart], siglayers);

        // 4. Compute k_4 (spatial interpolation at end step)
        double[] k4 = stepAhead(
                k3Loc, k3Depth,
                elemPart, depLayer, 1.0,
                uvnode, depthUvnode, nodexy, trinodes, siglayers, neighbours,
                u, v, w,
                hour, step, dt, stepsPerStep, coordOS, FVCOM);

        // 5. Add it all together
        if(k1[0] != 0 && k2[0] != 0 && k3[0] != 0 && k4[0] != 0) {
            for (int i = 0; i < 3; i++) {
                advectStep[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
            }
        }
        return advectStep;
    }

    /**
     * Calculate the correction steps required for the RK4 algorithm
     *
     * @param xy            Location of new velocity to be used (current location plus spatial step)
     * @param timeStepAhead Time step ahead i.e. "1.0/2.0" if half-step ahead, "1.0" if full step ahead
     */
    public static double[] stepAhead(double[] xy, double depth, int elemPart, int depLayer, double timeStepAhead,
                                     float[][] uvnode, float[] depthUvnode, float[][] nodexy, int[][] trinodes,
                                     float[] siglayers, int[][] neighbours,
                                     float[][][] u, float[][][] v, float[][][] w,
                                     int hour, int step, double dt,
                                     int stepsPerStep, boolean coordOS, boolean FVCOM) {
        double[] xyz_step = {0,0,0};
        double[] vel = {0,0,0};
        double[] velplus1 = {0,0,0};

        if (FVCOM) {
            double[][] xNrList;
            xNrList = neighbourCellsList3D(xy, elemPart, depLayer, uvnode, depthUvnode, nodexy, trinodes, siglayers, neighbours, depth, coordOS);

            // 2. Compute k_1 (spatial interpolation at start of step)
            // Velocity from start of time step
            vel = velocityFromNearestList(xNrList, u[hour], v[hour], w[hour]);
            // Velocity from end of this hour - will linearly interpolate to end of subTimeStep below and in stepAhead
            velplus1 = velocityFromNearestList(xNrList, u[hour+1], v[hour+1], w[hour+1]);

            // If predicted location of particle is outside mesh, return zero velocity
            if (xNrList[0][0] == 0) {
                return xyz_step;
            }

        } else {
//            if (verticalDynamics) {
//                System.err.println("Error in rk4Step: vertical dynamics not implemented for ROMS");
//                System.exit(1);
//            }
//            // nearestROMSGridPoint ---> whichROMSelement --->
//            double[][] nrListU = nearestListROMS(xy, meshes.get(meshPart).getLonU(), meshes.get(meshPart).getLatU(), null);
//            double[][] nrListV = nearestListROMS(xy, meshes.get(meshPart).getLonV(), meshes.get(meshPart).getLatV(), null);
//            vel = velocityFromNearestListROMS(nrListU, nrListV, hour, u, v);
//            velplus1 = velocityFromNearestListROMS(nrListU, nrListV, hour + 1, u, v);
//            if (nrListU[0][0] == 0) {
//                return xyz_step;
//            }
        }

        // Do the relevant temporal interpolation for this part of the step
        for (int i = 0; i < 3; i++) {
            xyz_step[i] = dt * (vel[i] + ((step + timeStepAhead) / (double) stepsPerStep) * (velplus1[i] - vel[i]));
        }

        return xyz_step;
    }

    /**
     * Compute an Euler integration step for particle movement
     */
    @SuppressWarnings("ConstantConditions")
    public double[] eulerStep(List<HydroField> hydroFields, // velocities
                              List<Mesh> meshes,     // other mesh info
                              int hour, int step, double dt,                                  // locate particle in space and time
                              int stepsPerStep, boolean coordOS, boolean fixDepth) {
        if(!fixDepth) {
            System.err.println("Error: vertical dynamics not implemented with Particle.eulerStep()");
            System.exit(1);
        }
        int elemPart = this.getElem();
        int meshPart = this.getMesh();
        int dep = this.getDepthLayer();

        float[][] uvnode = meshes.get(meshPart).getUvnode();
        float[][] nodexy = meshes.get(meshPart).getNodexy();
        int[][] trinodes = meshes.get(meshPart).getTrinodes();
        int[][] neighbours = meshes.get(meshPart).getNeighbours();
        float[][][] u = hydroFields.get(meshPart).getU();
        float[][][] v = hydroFields.get(meshPart).getV();
        float[][][] w = hydroFields.get(meshPart).getW();

        double[] thisLoc = this.getLocation();

        double[] advectStep = new double[2];
        // compute velocities at start and end of entire step, at the new location
        double[] vel = new double[2];
        double[] velplus1 = new double[2];

        if (meshes.get(meshPart).getType().equals("FVCOM") || meshes.get(meshPart).getType().equals("ROMS_TRI")) {
            double[][] xNrList = neighbourCellsList(thisLoc, elemPart, meshes.get(meshPart), coordOS);

            // 2. Compute k_1 (spatial interpolation at start of step)
            // Velocity from start of time step
            vel = velocityFromNearestList(xNrList, u[hour], v[hour], w[hour]);
            // Velocity from end of this hour - will linearly interpolate to end of subTimeStep below and in stepAhead
            velplus1 = velocityFromNearestList(xNrList, u[hour+1], v[hour+1], w[hour+1]);

            // If predicted location of particle is outside mesh, return zero velocity
            if (xNrList[0][0] == 0) {
                return advectStep;
            }
        }
        // Find the nearest nodes for interpolation
        else if (meshes.get(meshPart).getType().equals("ROMS")) {
            // Populate the list of points for interpolation
            double[][] nrListU = nearestListROMS2(thisLoc, meshes.get(meshPart).getLonU(), meshes.get(meshPart).getLatU(), this.nearestROMSGridPointU);
            double[][] nrListV = nearestListROMS2(thisLoc, meshes.get(meshPart).getLonV(), meshes.get(meshPart).getLatV(), this.nearestROMSGridPointV);
            // Interpolate the velocity
            vel = velocityFromNearestListROMS(nrListU, nrListV, hour, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV());
            velplus1 = velocityFromNearestListROMS(nrListU, nrListV, hour + 1, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV());
            if (nrListU[0][0] == 0) {
                return advectStep;
            }
        }
        // Compute the advection step based on temporally interpolated velocities
        advectStep[0] = dt * (vel[0] + ((double) (step) / (double) stepsPerStep) * (velplus1[0] - vel[0]));
        advectStep[1] = dt * (vel[1] + ((double) (step) / (double) stepsPerStep) * (velplus1[1] - vel[1]));
        return advectStep;
    }
}