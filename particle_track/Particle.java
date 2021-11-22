/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

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
    private String startSiteID = "0";
    private int startSiteIndex = 0;

    private final String coordRef;

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
    private double depth = 0;
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
    private boolean settledThisHour = false;
    private boolean boundaryExit = false;

    private String species = "none";

    private String lastArrival = "0";

    private double xTotal = 0;
    private double yTotal = 0;
    private double zTotal = 0;

    // A list to store data on the arrivals made by each particle.
    // If rp.endOnArrival=true, this will contain a maximum of one element.
    // --- Not presently used ---
    private List<Arrival> arrivals;

    public Particle(double xstart, double ystart, double startDepth, String startSiteID, int startIndex, int id, double mortalityRate,
                    ISO_datestr startDate, double startTime, String coordRef, String species) {
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
        this.coordRef = coordRef;
        this.species = species;
    }

    public Particle(double xstart, double ystart, double startDepth, String startSiteID, int startIndex, int id, double mortalityRate, double startDensity,
                    ISO_datestr startDate, double startTime, String coordRef, String species) {
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
        this.density = startDensity;
        this.depth = startDepth;
        this.coordRef = coordRef;
        this.species = species;
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
        this.coordRef = "OSGB1936";
    }

    @Override
    public String toString() {
        return this.getID() + " " + this.xy.toString();
    }

    // Not presently used
    public void reportArrival(int sourceLocation, int arrivalLocation, double arrivalTime, double arrivalDensity) {
        arrivals.add(new Arrival(sourceLocation, arrivalLocation, arrivalTime, arrivalDensity));
    }

    public List<Arrival> getArrivals() {
        return this.arrivals;
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

    public void setDepthLayer(int dep) {
        this.depLayer = dep;
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

    public void addZ(double zTravel) {
        this.zTotal += zTravel;
    }

    public double getxTotal() {
        return this.xTotal;
    }

    public double getyTotal() {
        return this.yTotal;
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

    public void setViable(boolean viable) {
        this.viable = viable;
    }

    public void setArrived(boolean arrived) {
        this.arrived = arrived;
    }

    public void setLastArrival(String loc) {
        this.lastArrival = loc;
    }

    public String getLastArrival() {
        return this.lastArrival;
    }

    public void setSettledThisHour(boolean settled) {
        this.settledThisHour = settled;
    }

    public void setFree(boolean free) {
        this.free = free;
    }

    public void setBoundaryExit(boolean exit) {
        this.boundaryExit = exit;
    }

    public boolean isViable() {
        return this.viable;
    }

    public boolean hasArrived() {
        return this.arrived;
    }

    public boolean hasSettledThisHour() {
        return this.settledThisHour;
    }

    public boolean isFree() {
        return this.free;
    }

    public boolean hasExited() {
        return this.boundaryExit;
    }

    public boolean canBecomeViable(RunProperties rp) {
        boolean byAge = getAge() > rp.viabletime && rp.viabletime > 0;
        boolean byDegreeDays = getDegreeDays() > rp.viableDegreeDays && rp.viableDegreeDays > 0;
        boolean byPrevious = isViable();
        return byAge || byDegreeDays || byPrevious;
    }

    public boolean isTooOld(RunProperties rp) {
        boolean byAge = getAge() > rp.maxParticleAge && rp.maxParticleAge > 0;
        boolean byDegreeDays = getDegreeDays() > rp.maxDegreeDays && rp.maxDegreeDays > 0;
        boolean byPrevious = getStatus() == 666;
        return byAge || byDegreeDays || byPrevious;
    }

    /**
     * Sets the mortality rate for the particle packet based upon local salinity
     */
    public void setMortRate(double salinity) {
        // estimated 2nd order polynomial fit to Bricknell et al.'s (2006) data
        this.mortRate = 0.0011 * salinity * salinity - 0.07 * salinity + 1.1439;
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

    public double sink(RunProperties rp) {
        int forceDownward = rp.sinkingRateMean < 0 ? -1 : 1;
        return forceDownward * rp.sinkingRateMean + rp.sinkingRateStd * ThreadLocalRandom.current().nextGaussian();
    }

    public double swim(RunProperties rp) {
        double swimSpeedMean = this.isViable() ? rp.vertSwimSpeedCopepodidMean : rp.vertSwimSpeedNaupliusMean;
        double swimSpeedStd = this.isViable() ? rp.vertSwimSpeedCopepodidStd : rp.vertSwimSpeedNaupliusStd;
        int forceUpward = swimSpeedMean > 0 ? -1 : 1;
        return forceUpward * swimSpeedMean + swimSpeedStd * ThreadLocalRandom.current().nextGaussian();
    }

    public double swim(RunProperties rp, Mesh mesh, HydroField hydroField, int hour) {
        // following Johnsen et al 2014, 2016, Myksvoll et al 2018, Sandvik et al 2020 for light attenuation and swimming thresholds
        // short_wave units = W m-2 ≈ 2.1 μmole m-2 s-1 (according to Tom Adams handover files)
        double lightAtSurface = 2.1 * hydroField.getAvgFromTrinodes(mesh, this.getLocation(), 0, this.elem, hour, "short_wave", rp);
        double lightAtDepth = lightAtSurface * Math.exp(-0.2 * this.depth);

        if ((this.isViable() && lightAtDepth > 2.06e-5) ||
                (!this.isViable() && lightAtDepth > 0.392)) {
            return swim(rp);
        } else {
            return 0.0;
        }
    }

    public double[] diffuse(RunProperties rp, double K_gradient, double K_zAdj, double dt, String distribution) {
        int nDims = rp.verticalDynamics ? 3 : 2;
        double[] diffusion = {0,0,0};
        double[] randoms = {0,0,0};
        double r = 0.0;

        if (distribution.equalsIgnoreCase("uniform")) {
            for (int i = 0; i < nDims; i++) {
                randoms[i] = ThreadLocalRandom.current().nextDouble(-1.0, 1.0);
            }
            r = 1.0 / 3.0;
        } else if (distribution.equalsIgnoreCase("gaussian")) {
            for (int i = 0; i < nDims; i++) {
                randoms[i] = ThreadLocalRandom.current().nextGaussian();
            }
            r = 1.414;
        } else {
            System.err.println("Diffusion distribution must be uniform or gaussian.");
            System.exit(1);
        }
        diffusion[0] = randoms[0] * Math.sqrt(2 / r * rp.D_h * dt);
        diffusion[1] = randoms[1] * Math.sqrt(2 / r * rp.D_h * dt);

        if (rp.verticalDynamics) {
            if (rp.variableDiffusion) {
                // following Visser 1997 (eq. 6):
                // dZ = K_gradient * dt + Rand * sqrt(2 * dt/r * K_zAdj), where K_zAdj = K(z + K_gradient/2 * dt)
                diffusion[2] = K_gradient * dt + randoms[2] * Math.pow(2 / r * K_zAdj * dt, 0.5);
            } else {
                diffusion[2] = randoms[2] * Math.sqrt(2 / r * rp.D_hVert * dt);
            }
        }
        return diffusion;
    }

    /**
     * Is the particle still in the present mesh, should it change mesh, and has it exited the overall domain?
     */
    public void meshSelectOrExit(double[] newLoc, List<Mesh> meshes, RunProperties rp) {
        //  i) Is particle in present mesh?
        //      i.i) If YES and presently in mesh n > 0, is it also within mesh < n?
        //          If YES, place in that mesh
        //          If NO, stay in same mesh
        //      i.ii) If YES and in mesh 0, is it within certain range of boundary?
        //          i.ii.i) If YES, is it in other mesh > 0?
        //              If YES, switch to other mesh
        //              If NO, BOUNDARY EXIT!! (nearest bnode)
        //          If NO: CARRY ON!!!
        //      i.iii) If NO, is particle within other mesh > n?
        //          If YES, switch mesh
        //          If NO, BOUNDARY EXIT!! (nearest bnode)

        // Get current mesh and element
        int meshID = this.getMesh();
        double[] oldLoc = this.getLocation();
        int[] el = new int[2];
        if (meshes.get(meshID).getType().equalsIgnoreCase("ROMS")) {
            el = this.getROMSnearestPointU();
        } else if (meshes.get(meshID).getType().equalsIgnoreCase("ROMS_TRI")) {
            el[0] = this.getElem();
        } else {
            el[0] = this.getElem();
        }

        // i)
        Mesh m = meshes.get(meshID);
        if (m.isInMesh(newLoc, true, false, el)) {
            // i.i) in mesh > 0
            if (meshID > 0) {
                m = meshes.get(0);
                if (m.isInMesh(newLoc, true, true, null)) {
                    // switch to mesh 0
                    this.setLocation(newLoc[0], newLoc[1]);
                    this.placeInMesh(meshes, 0, true);
                } else {
                    // stay in current mesh
                    this.setLocation(newLoc[0], newLoc[1]);
                    this.placeInMesh(meshes, 1, false);
                }
            }
            // i.ii) in mesh 0
            else {
                // boundary check
                int bnode = -1;
                if (rp.checkOpenBoundaries) {
                    bnode = ParallelParticleMover.openBoundaryCheck((float) newLoc[0], (float) newLoc[1], m, rp);
                }
                if (bnode != -1) {
                    // Close to boundary
                    if (meshes.size() == 2) {
                        m = meshes.get(1);
                        if (m.isInMesh(newLoc, true, true, null)) {
                            // switch to other mesh
                            this.setLocation(newLoc[0], newLoc[1]);
                            this.placeInMesh(meshes, 1, true);
                        }
                    } else {
                        this.setBoundaryExit(true);
                        this.setStatus(66);
                    }
                } else {
                    this.setLocation(newLoc[0], newLoc[1]);
                    this.placeInMesh(meshes, 0, false);
                    // stay in same mesh and keep going!
                }
            }
        } else {
            if (meshes.size() == 2) {
                // Not in original mesh, so check the other one
                int otherMesh = 1;
                if (meshID == 1) {
                    otherMesh = 0;
                }
                m = meshes.get(otherMesh);
                if (m.isInMesh(newLoc, true, true, null)) {
                    // switch to other mesh
                    this.setLocation(newLoc[0], newLoc[1]);
                    this.placeInMesh(meshes, otherMesh, true);
                } else {
                    // boundary exit
                    this.setBoundaryExit(true);
                    this.setStatus(66);
                }
            }
        }
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

        if (m.getType().equalsIgnoreCase("FVCOM") || m.getType().equalsIgnoreCase("ROMS_TRI")) {
            int el = 0;
            boolean checkAll = true;
            if (!switchedMesh) {
                el = this.getElem();
                checkAll = false;
            }
            int[] c = findContainingElement(this.getLocation(), el,
                    m.getNodexy(), m.getTrinodes(), m.getNeighbours(), checkAll);
            // if particle is within the mesh, update location normally and save the distance travelled
            this.setElem(c[0]);

        } else if (m.getType().equalsIgnoreCase("ROMS")) {
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

            // Need to save nearest point, in order to read into
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
     * Make a list of nearest mesh element centroids
     */
    public static double[][] nearestCentroidList(double x, double y, float[][] uvnode) {
        double[][] nearestList = new double[5][2];
        int nearest = -1;
        double dist = 10000000;

        for (int i = 0; i < uvnode[0].length; i++) {
            double distnew = Math.sqrt((x - uvnode[0][i]) * (x - uvnode[0][i]) + (y - uvnode[1][i]) * (y - uvnode[1][i]));
            if (distnew < dist) {
                dist = distnew;
                nearest = i;
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
                float distnew = (float) distanceEuclid2(x, y, xGrid[i][j], yGrid[i][j], "WGS84");
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
     * Note that this might no actually be what we want, due to the form of the grid.
     * Instead, we probably want the 4 corners of the containing element - see nearestListROMS2.
     */
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
                float distnew = (float) distanceEuclid2(xy[0], xy[1], xGrid[i][j], yGrid[i][j], "WGS84");
                allDists[k][0] = i;
                allDists[k][1] = j;
                allDists[k][2] = distnew;
                k++;
            }
        }

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
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0]][whichElem[1]], yGrid[whichElem[0]][whichElem[1]], "WGS84")};
        nearList[1] = new double[]{whichElem[0] + 1, whichElem[1],
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0] + 1][whichElem[1]], yGrid[whichElem[0] + 1][whichElem[1]], "WGS84")};
        nearList[2] = new double[]{whichElem[0], whichElem[1] + 1,
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0]][whichElem[1] + 1], yGrid[whichElem[0]][whichElem[1] + 1], "WGS84")};
        nearList[3] = new double[]{whichElem[0] + 1, whichElem[1] + 1,
                (float) distanceEuclid2(xy[0], xy[1], xGrid[whichElem[0] + 1][whichElem[1] + 1], yGrid[whichElem[0] + 1][whichElem[1] + 1], "WGS84")};
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


    public static double[] velocityFromNearestList(double[][] nrList, int hour, float[][][] u, float[][][] v, float[][][] w, int depLayer, boolean verticalDynamics) {
        int nDims = verticalDynamics ? 3 : 2;
        double[] velocity = {0,0,0};
        double[] weights = new double[nrList.length];
        double usum = 0, vsum = 0, wsum = 0, sum = 0;
        for (int i = 0; i < nrList.length; i++) {
            int elem = (int) nrList[i][0];
            double distance = nrList[i][1];
            if (nrList[i][1] != 0) {
                weights[i] = 1.0 / (distance * distance);
            } else {
                weights[i] = 1;
            }
            if (verticalDynamics) {
                depLayer = (int) nrList[i][2];
                wsum = wsum - weights[i] * w[hour][depLayer][elem];  // upward velocity; ww * -1 so that positive values = increased depth
            }
            usum = usum + weights[i] * u[hour][depLayer][(int) nrList[i][0]];
            vsum = vsum + weights[i] * v[hour][depLayer][(int) nrList[i][0]];
            sum = sum + weights[i];
        }
        velocity[0] = usum / sum;
        velocity[1] = vsum / sum;
        if (verticalDynamics) {
            velocity[2] = wsum / sum;
        }
        return velocity;
    }

    /**
     * A method to use the nearest grid point indexes and distances to interpolate a velocity.
     * This is done between the four nearest values.
     * In the ROMS output, U and V are stored on different grids (Arakawa-C). Therefore
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
    public static double distanceEuclid2(double x1, double y1, double x2, double y2, String coordRef) {
        double dx = x1 - x2;
        double dy = y1 - y2;

        double[] distXY = new double[]{dx, dy};

        if (coordRef.equalsIgnoreCase("WGS84")) {
            distXY = ParallelParticleMover.distanceDegreesToMetres(distXY, new double[]{x1, y2});
        }

        return Math.sqrt(distXY[0] * distXY[0] + distXY[1] * distXY[1]);
    }

    public static double distanceEuclid2(double x1, double y1, double z1, double x2, double y2, double z2, String coordRef) {
        double dx = x1 - x2;
        double dy = y1 - y2;

        double[] distXY = new double[]{dx, dy};

        if (coordRef.equalsIgnoreCase("WGS84")) {
            distXY = ParallelParticleMover.distanceDegreesToMetres(distXY, new double[]{x1, y2});
        }

        return Math.sqrt(distXY[0] * distXY[0] + distXY[1] * distXY[1] + (z1 - z2) * (z1 - z2));
    }


    /**
     * Search through lists of elements progressively further away from the last known
     * location in order to find the new location.
     * <p>
     * Returned array contains the element in which the particle is determined to be located,
     * plus the number of counts required at each scale (for diagnostics).
     */
    public static int[] findContainingElement(double[] xy, int elemPart,
                                              float[][] nodexy, int[][] trinodes, int[][] neighbours, boolean checkAll) {
        int[] c = new int[6];
        int[] elems = new int[1];
        elems[0] = elemPart;

        int whereami = whichElement(xy, elems, nodexy, trinodes);
        c[1] = 1;

        // If this fails, look in my neighbours
        if (whereami == -1) {
            c[2] = 1;
            int[] elems0 = new int[]{neighbours[0][elemPart], neighbours[1][elemPart], neighbours[2][elemPart]};
            whereami = whichElement(xy, elems0, nodexy, trinodes);

            // if fails, look in the neighbours of my neighbours
            if (whereami == -1) {
                // Find the non-0 neighbours out of the ones just searched
                int[] nonZeroNeighbours = Arrays.stream(elems0).filter(num -> num != 0).toArray();
                // List THEIR neighbours
                int[] elems1 = new int[nonZeroNeighbours.length * 3];
                for (int i = 0; i < nonZeroNeighbours.length; i++) {
                    for (int j = 0; j < 3; j++) {
                        elems1[i * 3 + j] = neighbours[j][nonZeroNeighbours[i]];
                    }
                }
                c[3] = 1;
                whereami = whichElement(xy, elems1, nodexy, trinodes);

                // If fails, look in neighbours' neighbours' neighbours
                // Shouldn't really need to get this far? Maybe you should.
                if (whereami == -1) {
                    c[4] = 1;
                    // Find the non-0 neighbours out of the ones just searched
                    int[] nonZeroNeighbours1 = Arrays.stream(elems1).filter(num -> num != 0).toArray();
                    // List THEIR neighbours
                    int[] elems2 = new int[nonZeroNeighbours1.length * 3];
                    for (int i = 0; i < nonZeroNeighbours1.length; i++) {
                        for (int j = 0; j < 3; j++) {
                            elems2[i * 3 + j] = neighbours[j][nonZeroNeighbours1[i]];
                        }
                    }
                    whereami = whichElement(xy, elems2, nodexy, trinodes);

                    // if this fails, look in all elements (desperation) - unless movement is attempting to place particle outside mesh)
                    if (whereami == -1 && checkAll) {
                        c[5] = 1;
                        int[] allelems = IntStream.rangeClosed(0, trinodes[0].length - 1).toArray();
                        whereami = whichElement(xy, allelems, nodexy, trinodes);
                    }
                }
            }
        }
        c[0] = whereami;
        return c;
    }

    /**
     * Find which element a particle resides within (edge checking)
     */
    public static int whichElement(double[] xy, int[] elems, float[][] nodexy, int[][] trinodes) {
        int which = -1;
        int res = 0;
        for (int i = 0; i < elems.length; i++) {
            double[] xt = new double[3];
            double[] yt = new double[3];

            for (int j = 0; j < 3; j++) {
                try {
                    xt[j] = nodexy[0][trinodes[j][elems[i]]];
                    yt[j] = nodexy[1][trinodes[j][elems[i]]];
                } catch (Exception e) {
                    assert trinodes != null;
                    assert trinodes[j] != null;
                    assert nodexy != null;
                    assert nodexy[0] != null;
                    System.err.println(i + " " + j + " " + elems[i] + " " + trinodes[j][elems[i]] + " " + nodexy[0][trinodes[j][elems[i]]]);
                }
            }
            // check whether (x,y) lies within this
            double f1 = (xy[1] - yt[0]) * (xt[1] - xt[0]) - (xy[0] - xt[0]) * (yt[1] - yt[0]);
            double f2 = (xy[1] - yt[2]) * (xt[0] - xt[2]) - (xy[0] - xt[2]) * (yt[0] - yt[2]);
            double f3 = (xy[1] - yt[1]) * (xt[2] - xt[1]) - (xy[0] - xt[1]) * (yt[2] - yt[1]);
            if (f1 * f3 >= 0.0 && f3 * f2 >= 0.0) {
                res = 1;
            }
            if (res == 1) {
                which = elems[i];
                break;
            }
        }
        return which;
    }

    /**
     * Calculate a list of neighbour cells based on a specified (e.g. particle) location
     * (0: containing element, 1-3: neighbours of the containing element), and calculate
     * the Euclidean distances.
     */
    public static double[][] neighbourCellsList(double[] xy, int elemPart0, int meshPart, List<Mesh> meshes, String coordRef) {
        double[][] nrList = new double[5][2];

        float[][] nodexy = meshes.get(meshPart).getNodexy();
        float[][] uvnode = meshes.get(meshPart).getUvnode();
        int[][] trinodes = meshes.get(meshPart).getTrinodes();
        int[][] neighbours = meshes.get(meshPart).getNeighbours();

        int[] elem = findContainingElement(xy, elemPart0, nodexy, trinodes, neighbours, false);
        // If particle is not within the mesh (value returned by findContainingElement = -1)
        // exit this method returning array of zeros.
        if (elem[0] == -1) {
            return nrList;
        }
        int thisElem = elem[0];
        nrList[0][0] = thisElem;
        nrList[0][1] = distanceEuclid2(xy[0], xy[1], uvnode[0][thisElem], uvnode[1][thisElem], coordRef);
        // distance to neighbouring elems
        nrList[1][0] = neighbours[0][thisElem];
        nrList[1][1] = distanceEuclid2(xy[0], xy[1], uvnode[0][neighbours[0][thisElem]], uvnode[1][neighbours[0][thisElem]], coordRef);
        nrList[2][0] = neighbours[1][thisElem];
        nrList[2][1] = distanceEuclid2(xy[0], xy[1], uvnode[0][neighbours[1][thisElem]], uvnode[1][neighbours[1][thisElem]], coordRef);
        nrList[3][0] = neighbours[2][thisElem];
        nrList[3][1] = distanceEuclid2(xy[0], xy[1], uvnode[0][neighbours[2][thisElem]], uvnode[1][neighbours[2][thisElem]], coordRef);
        nrList[4][0] = 0;
        nrList[4][1] = 1000000;

        return nrList;
    }

    /**
     * Find the containing element plus 3 neighbors for the depth layer above and below the particle and calculate the
     * euclidean distances from the particle to the centroids
     *
     * @param xy particle location
     * @param elemPart0 particle element
     * @param meshPart particle mesh
     * @param meshes meshes
     * @param depth particle depth
     * @param coordRef coordinate system
     * @return double[8][3] with rows = elements and columns = [elementID, distance to particle (3D), depth layer]
     */
    public static double[][] neighbourCellsList3D(double[] xy, int elemPart0, int meshPart, List<Mesh> meshes, double depth, String coordRef) {

        float[][] uvnode = meshes.get(meshPart).getUvnode();
        float[] depthUvnode = meshes.get(meshPart).getDepthUvnode();
        int[][] neighbours = meshes.get(meshPart).getNeighbours();
        double[][] nrList = new double[8][3]; // (3 neighbors per layer + containing element) * 2 layers = 8
        float[][] nearestLayers;

        int thisElem = findContainingElement(xy, elemPart0, meshes.get(meshPart).getNodexy(), meshes.get(meshPart).getTrinodes(), neighbours, false)[0];
        if (thisElem == -1) {
            return nrList;
        }

        int nrRow = 0;
        for (int layer = 0; layer < 2; layer++) {
            // containing element
            nearestLayers = Mesh.findNearestSigmas(depth, meshes.get(meshPart).getSiglay(), depthUvnode[thisElem]);
            nrList[nrRow][0] = thisElem;
            nrList[nrRow][1] = distanceEuclid2(xy[0], xy[1], depth, uvnode[0][thisElem], uvnode[1][thisElem], (int) nearestLayers[layer][0], coordRef);
            nrList[nrRow][2] = (int) nearestLayers[layer][0];
            nrRow++;
            // neighboring elements
            for (int nbr = 0; nbr < 3; nbr++) {
                nearestLayers = Mesh.findNearestSigmas(depth, meshes.get(meshPart).getSiglay(), depthUvnode[neighbours[nbr][thisElem]]);
                nrList[nrRow][0] = neighbours[nbr][thisElem];
                nrList[nrRow][1] = distanceEuclid2(xy[0], xy[1], depth, uvnode[0][neighbours[nbr][thisElem]], uvnode[1][neighbours[nbr][thisElem]], (int) nearestLayers[layer][0], coordRef);
                nrList[nrRow][2] = (int) nearestLayers[layer][0];
                nrRow++;
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
                            int stepsPerStep, String coordRef, boolean verticalDynamics) {  // info on simulation length
        int elemPart = this.getElem();
        int meshPart = this.getMesh();
        int nDims = verticalDynamics ? 3 : 2;
        int depLayer = this.getDepthLayer();

        double[] advectStep = {0,0,0};
        double[] vel = {0,0,0};
        double[] velplus1 = {0,0,0};

        if (meshes.get(meshPart).getType().equalsIgnoreCase("FVCOM") || meshes.get(meshPart).getType().equalsIgnoreCase("ROMS_TRI")) {
            if (verticalDynamics) {
                this.nrList = neighbourCellsList3D(this.getLocation(), elemPart, meshPart, meshes, this.depth, coordRef);
            } else {
                this.nrList = neighbourCellsList(this.getLocation(), elemPart, meshPart, meshes, coordRef);
            }

            // 2. Compute k_1 (spatial interpolation at start of step)
            // Velocity from start of timestep
            vel = velocityFromNearestList(this.getNrList(), hour, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV(), hydroFields.get(meshPart).getW(), depLayer, verticalDynamics);
            // Velocity from end of this hour - will linearly interpolate to end of subTimeStep below and in stepAhead
            velplus1 = velocityFromNearestList(this.getNrList(), hour + 1, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV(), hydroFields.get(meshPart).getW(), depLayer, verticalDynamics);

        } else if (meshes.get(meshPart).getType().equalsIgnoreCase("ROMS")) {
            if (verticalDynamics) {
                System.err.println("Error: vertical dynamics not implemented for ROMS");
                System.exit(1);
            }
            double[][] nrListROMSU = nearestListROMS(this.getLocation(), meshes.get(meshPart).getLonU(), meshes.get(meshPart).getLatU(), null);
            double[][] nrListROMSV = nearestListROMS(this.getLocation(), meshes.get(meshPart).getLonV(), meshes.get(meshPart).getLatV(), null);
            vel = velocityFromNearestListROMS(nrListROMSU, nrListROMSV, hour, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV());
            velplus1 = velocityFromNearestListROMS(nrListROMSU, nrListROMSV, hour + 1, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV());
        }
        double[] k1 = {0,0,0};
        for (int i = 0; i < nDims; i++) {
            k1[i]  = dt * (vel[i] + ((double) step / (double) stepsPerStep) * (velplus1[i] - vel[i]));
        }

        // 3. Compute k_2 (spatial interpolation at half step, temporal interp at half step)
        // Estimated half-step location using Euler
        // NOTE that here, for the purposes of identifying the elements containing the part-step locations,
        // the steps in degrees are calculated. These are separate of the actual half-step values in metres
        // which are retained and summed at the end of the method to give a transport distance in metres
        double[] k1Deg = new double[]{k1[0], k1[1]};
        if (this.coordRef.equalsIgnoreCase("WGS84")) {
            k1Deg = ParallelParticleMover.distanceMetresToDegrees2(k1Deg, this.getLocation());
        }
        double k1Depth = this.getDepth() + k1[2] / 2.0;  // k1[2] = 0 if verticalDynamics = false
        depLayer = findLayerFromDepth(k1Depth, meshes.get(meshPart).getDepthUvnode()[elemPart], meshes.get(meshPart).getSiglay());
        double[] k2 = stepAhead(
                new double[]{this.getLocation()[0] + k1Deg[0] / 2.0, this.getLocation()[1] + k1Deg[1] / 2.0}, k1Depth,
                elemPart, depLayer, 1.0 / 2.0,
                meshes, meshPart, hydroFields,
                hour, step, dt, stepsPerStep, coordRef, verticalDynamics);

        // 4. Compute k_3 (spatial interpolation at half step, temporal interp at half step)
        double[] k2Deg = new double[]{k2[0], k2[1]};
        if (this.coordRef.equalsIgnoreCase("WGS84")) {
            k2Deg = ParallelParticleMover.distanceMetresToDegrees2(k2Deg, this.getLocation());
        }
        double k2Depth = this.getDepth() + k2[2] / 2.0;  // k2[2] = 0 if verticalDynamics = false
        depLayer = findLayerFromDepth(k2Depth, meshes.get(meshPart).getDepthUvnode()[elemPart], meshes.get(meshPart).getSiglay());
        double[] k3 = stepAhead(
                new double[]{this.getLocation()[0] + k2Deg[0] / 2.0, this.getLocation()[1] + k2Deg[1] / 2.0}, k2Depth,
                elemPart, depLayer, 1.0 / 2.0,
                meshes, meshPart, hydroFields,
                hour, step, dt, stepsPerStep, coordRef, verticalDynamics);

        // 5. Compute k_4 (spatial interpolation at end step)
        double[] k3Deg = new double[]{k3[0], k3[1]};
        if (this.coordRef.equalsIgnoreCase("WGS84")) {
            k3Deg = ParallelParticleMover.distanceMetresToDegrees2(k3Deg, this.getLocation());
        }
        double k3Depth = this.getDepth() + k3[2];  // k3[2] = 0 if verticalDynamics = false
        depLayer = findLayerFromDepth(k3Depth, meshes.get(meshPart).getDepthUvnode()[elemPart], meshes.get(meshPart).getSiglay());
        double[] k4 = stepAhead(
                new double[]{this.getLocation()[0] + k3Deg[0], this.getLocation()[1] + k3Deg[1]}, k3Depth,
                elemPart, depLayer, 1.0,
                meshes, meshPart, hydroFields,
                hour, step, dt, stepsPerStep, coordRef, verticalDynamics);

        // 6. Add it all together
        if(k1[0] != 0 && k2[0] != 0 && k3[0] != 0 && k4[0] != 0) {
            for (int i = 0; i < nDims; i++) {
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
                                     List<Mesh> meshes, int meshPart,
                                     List<HydroField> hydroFields,
                                     int hour, int step, double dt,
                                     int stepsPerStep, String coordRef, boolean verticalDynamics) {
        int nDims = verticalDynamics ? 3 : 2;
        double[] xyz_step = {0,0,0};
        double[] vel = {0,0,0};
        double[] velplus1 = {0,0,0};

        if (meshes.get(meshPart).getType().equalsIgnoreCase("FVCOM") || meshes.get(meshPart).getType().equalsIgnoreCase("ROMS_TRI")) {
            double[][] xNrList;
            if (verticalDynamics) {
                xNrList = neighbourCellsList3D(xy, elemPart, meshPart, meshes, depth, coordRef);
            } else {
                xNrList = neighbourCellsList(xy, elemPart, meshPart, meshes, coordRef);
            }

            // 2. Compute k_1 (spatial interpolation at start of step)
            // Velocity from start of timestep
            vel = velocityFromNearestList(xNrList, hour, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV(), hydroFields.get(meshPart).getW(), depLayer, verticalDynamics);
            // Velocity from end of this hour - will linearly interpolate to end of subTimeStep below and in stepAhead
            velplus1 = velocityFromNearestList(xNrList, hour + 1, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV(), hydroFields.get(meshPart).getW(), depLayer, verticalDynamics);

            // If predicted location of particle is outside mesh, return zero velocity
            if (xNrList[0][0] == 0) {
                return xyz_step;
            }

        } else if (meshes.get(meshPart).getType().equalsIgnoreCase("ROMS")) {
            if (verticalDynamics) {
                System.err.println("Error in rk4Step: vertical dynamics not implemented for ROMS");
                System.exit(1);
            }
            // nearestROMSGridPoint ---> whichROMSelement --->
            double[][] nrListU = nearestListROMS(xy, meshes.get(meshPart).getLonU(), meshes.get(meshPart).getLatU(), null);
            double[][] nrListV = nearestListROMS(xy, meshes.get(meshPart).getLonV(), meshes.get(meshPart).getLatV(), null);
            vel = velocityFromNearestListROMS(nrListU, nrListV, hour, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV());
            velplus1 = velocityFromNearestListROMS(nrListU, nrListV, hour + 1, hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV());
            if (nrListU[0][0] == 0) {
                return xyz_step;
            }
        }

        // Do the relevant temporal interpolation for this part of the step
        for (int i = 0; i < nDims; i++) {
            xyz_step[i] = dt * (vel[i] + ((step + timeStepAhead) / (double) stepsPerStep) * (velplus1[i] - vel[i]));
        }

        return xyz_step;
    }

    /**
     * Compute an Euler integration step for particle movement
     */
    public double[] eulerStep(List<HydroField> hydroFields, // velocities
                              List<Mesh> meshes,     // other mesh info
                              int hour, int step, double dt,                                  // locate particle in space and time
                              int stepsPerStep, String coordRef, boolean verticalDynamics) {
        if(verticalDynamics) {
            System.err.println("Error: vertical dynamics not implemented with Particle.eulerStep()");
            System.exit(1);
        }
        int elemPart = this.getElem();
        int meshPart = this.getMesh();
        int dep = this.getDepthLayer();

        double[] thisLoc = this.getLocation();

        double[] advectStep = new double[2];
        // compute velocities at start and end of entire step, at the new location
        double[] vel = new double[2];
        double[] velplus1 = new double[2];

        if (meshes.get(meshPart).getType().equalsIgnoreCase("FVCOM") || meshes.get(meshPart).getType().equalsIgnoreCase("ROMS_TRI")) {
            double[][] xNrList = neighbourCellsList(thisLoc, elemPart, meshPart, meshes, coordRef);

            // 2. Compute k_1 (spatial interpolation at start of step)
            // Velocity from start of timestep
            vel = velocityFromNearestList(xNrList, hour,
                    hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV(), hydroFields.get(meshPart).getW(), dep, verticalDynamics);
            // Velocity from end of this hour - will linearly interpolate to end of subTimeStep below and in stepAhead
            velplus1 = velocityFromNearestList(xNrList, hour + 1,
                    hydroFields.get(meshPart).getU(), hydroFields.get(meshPart).getV(), hydroFields.get(meshPart).getW(), dep, verticalDynamics);

            // If predicted location of particle is outside mesh, return zero velocity
            if (xNrList[0][0] == 0) {
                return advectStep;
            }
        }
        // Find nearest nodes for interpolation
        else if (meshes.get(meshPart).getType().equalsIgnoreCase("ROMS")) {
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