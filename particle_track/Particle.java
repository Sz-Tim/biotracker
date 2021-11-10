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
    private double[] xy = new double[2];
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
     *
     * @param locationFileLine
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
     * Set depth of particle with check against local depth
     *
     * @param depth          (this is a negative value)
     * @param localDepth (this is a positive value)
     */
    public void setDepth(double depth, double localDepth) {
        if (depth > localDepth) {
            depth = localDepth;
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
     * @param layers
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

    public boolean getViable() {
        return this.viable;
    }

    public boolean getArrived() {
        return this.arrived;
    }

    public boolean getSettledThisHour() {
        return this.settledThisHour;
    }

    public boolean getFree() {
        return this.free;
    }

    public boolean getBoundaryExit() {
        return this.boundaryExit;
    }

    public boolean canBecomeViable(RunProperties rp) {
        boolean byAge = getAge() > rp.viabletime && rp.viabletime > 0;
        boolean byDegreeDays = getDegreeDays() > rp.viableDegreeDays && rp.viableDegreeDays > 0;
        boolean byPrevious = getViable();
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
     *
     * @param salinity
     * @return
     */
    public void setMortRate(double salinity) {
        // estimated 2nd order polynomial fit to Bricknell et al.'s (2006) data
        this.mortRate = 0.0011 * salinity * salinity - 0.07 * salinity + 1.1439;
        //System.out.println("salinity = "+salinity+" mortrate = "+this.mortRate);
    }

    public void incrementAge(double increment) {
        this.age += increment;
    }

    public void incrementDegreeDays(double temp, RunProperties rp) {
        double inc = temp * (rp.dt / rp.stepsPerStep) / 86400;
        //System.out.println("DD increment: T="+temp+" dt="+rp.dt+" stepsPerStep="+rp.stepsPerStep+" dt/step="+(rp.dt/rp.stepsPerStep)+" inc="+inc);
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
        int forceUpward = rp.vertSwimSpeedMean > 0 ? -1 : 1;
        return forceUpward * rp.vertSwimSpeedMean + rp.vertSwimSpeedStd * ThreadLocalRandom.current().nextGaussian();
    }

    public double swim(RunProperties rp, Mesh mesh, HydroField hydroField, int hour) {
        // following Johnsen et al 2014, 2016, Myksvoll et al 2018, Sandvik et al 2020 for light attenuation and swimming thresholds
        // short_wave units = W m-2 ≈ 2.1 μmole m-2 s-1 (according to Tom Adams handover files)
        double lightAtSurface = 2.1 * hydroField.getAvgFromTrinodes(mesh, this.getLocation(), 0, this.elem, hour, "short_wave", rp);
        double lightAtDepth = lightAtSurface * Math.exp(-0.2 * this.depth);

        if ((this.getStatus() == 1 && lightAtDepth > 2.06e-5) || (this.getStatus() == 2 && lightAtDepth > 0.392)) {
            return swim(rp);
        } else {
            return 0.0;
        }
    }

    public double[] diffuse(RunProperties rp, double D_hVertDz, double dt, String distribution) {
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
            System.out.println("Diffusion distribution must be uniform or gaussian.");
            System.exit(0);
        }
        diffusion[0] = randoms[0] * Math.sqrt(2 / r * rp.D_h * dt);
        diffusion[1] = randoms[1] * Math.sqrt(2 / r * rp.D_h * dt);

        if (rp.verticalDynamics) {
            if (rp.variableDiffusion) {
                // folowing Visser 1997 (eq. 6):
                // dZ = K_grad * dt + Rand * sqrt(2 * 1/r * K * (z_n + 1/2 * K_grad * dt) * dt)
                // dZ = D_hVertDz * dt + Rand * Math.pow(2 * dt/r * D_hVert * (this.depth + dt/2 * D_hVertDz) , 0.5)
                diffusion[2] = D_hVertDz * dt +
                        randoms[2] * Math.pow(2 * dt/r * rp.D_hVert * (this.depth + dt/2 * D_hVertDz) , 0.5);
            } else {
                diffusion[2] = randoms[2] * Math.sqrt(2 / r * rp.D_hVert * dt);
            }
        }
        return diffusion;
    }

    /**
     * Is the particle still in the present mesh, should it change mesh, and has it exited the overall domain?
     *
     * @param newLoc
     * @param meshes
     * @param rp
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
        int move = 0;
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

        // Only move the particle at the end of this method, IFF it is allocated to move
        //this.setLocation(newLoc[0],newLoc[1]);

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
                    move = 1;
                } else {
                    // stay in current mesh
                    //System.out.println("In mesh "+meshID+", which is of type "+meshes.get(meshID).getType());
                    this.setLocation(newLoc[0], newLoc[1]);
                    this.placeInMesh(meshes, 1, false);
                    move = 1;
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
                            move = 1;
                        }
                    } else {
                        this.setBoundaryExit(true);
                        this.setStatus(66);
                        move = 0;
                    }
                } else {
                    this.setLocation(newLoc[0], newLoc[1]);
                    this.placeInMesh(meshes, 0, false);
                    move = 1;
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
                    move = 1;
                } else {
                    // boundary exit
                    this.setBoundaryExit(true);
                    this.setStatus(66);
                    move = 0;
                }
            } else {
                move = 0;
            }
        }
    }


    /**
     * This method updates the particle mesh neighbourhood information, dealing with a change in mesh ID if required
     *
     * @param meshes            The mesh that has been shown to contain a particle
     * @param id
     * @param switchedMesh Logical, has the particle changed mesh?
     */
    public void placeInMesh(List<Mesh> meshes, int id, boolean switchedMesh) {
        //System.out.println("In placeInMesh "+m.getType());
        Mesh m = meshes.get(id);
        this.setMesh(id);

        if (m.getType().equalsIgnoreCase("FVCOM") || m.getType().equalsIgnoreCase("ROMS_TRI")) {
            //System.out.println("mesh verified as FVCOM");
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
            //System.out.println("mesh verified as ROMS");
            int[] searchCentreU = null, searchCentreV = null;
            if (!switchedMesh) {
                searchCentreU = this.getROMSnearestPointU();
                searchCentreV = this.getROMSnearestPointV();
            }
            int[] nearU = nearestROMSGridPoint((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonU(), m.getLatU(), searchCentreU);
            int[] nearV = nearestROMSGridPoint((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonV(), m.getLatV(), searchCentreV);

            //System.out.println("found nearest point");
            // More to do here to turn the nearest grid point into the containing element
            int[] containingROMSElemU = whichROMSElement((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonU(), m.getLatU(),
                    nearU);
            int[] containingROMSElemV = whichROMSElement((float) this.getLocation()[0], (float) this.getLocation()[1],
                    m.getLonV(), m.getLatV(),
                    nearV);
            //System.out.println("found which element");

            // Need to save nearest point, in order to read into
            this.setROMSnearestPointU(nearU);
            this.setROMSnearestPointV(nearV);
            this.setROMSElemU(containingROMSElemU);
            this.setROMSElemV(containingROMSElemV);
        }
    }


    /**
     * Find the nearest mesh element centroid
     *
     * @param x
     * @param y
     * @param uvnode
     * @return
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
     *
     * @param x
     * @param y
     * @param uvnode
     * @return
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
     *
     * @param x
     * @param y
     * @param xGrid
     * @param yGrid
     * @return
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
            //System.out.println("Particle.nearestROMSGridPoint searchCentre: "+searchCentre[0]+","+searchCentre[1]);
            minX = Math.max(0, searchCentre[0] - 3);
            maxX = Math.min(xGrid.length, searchCentre[0] + 4);
            minY = Math.max(0, searchCentre[1] - 3);
            maxY = Math.min(xGrid[0].length, searchCentre[1] + 4);
        } else {
            //System.out.println("Particle.nearestROMSGridPoint searchCentre: NULL");
        }


//        System.out.println("Particle.nearestROMSGridPoint limits: X=["+minX+","+maxX+"] Y=["+minY+","+maxY+"]");

        for (int i = minX; i < maxX; i++) {
            for (int j = minY; j < maxY; j++) {
                float distnew = (float) distanceEuclid2(x, y, xGrid[i][j], yGrid[i][j], "WGS84");
                //float distnew = (float)Math.sqrt((x-xGrid[i][j])*(x-xGrid[i][j])+(y-yGrid[i][j])*(y-yGrid[i][j]));
//                System.out.printf("In --- Particle.nearestROMSGridPoint "+i+" "+j+" "+distnew+" ---\n");
                if (distnew < dist) {
                    dist = distnew;
                    nearest = new int[]{i, j};
                }
            }
        }
//        System.out.printf("In Particle.nearestROMSGridPoint; RESULT = "+nearest[0]+" "+nearest[1]+" "+dist+"\n");
        return nearest;
    }

    /**
     * This makes a list of the nearest ROMS grid points - 4 of them.
     * Note that this might no actually be what we want, due to the form of the grid.
     * Instead, we probably want the 4 corners of the containing element - see nearestListROMS2.
     *
     * @param xy
     * @param xGrid
     * @param yGrid
     * @param nearestPoint
     * @return
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
//            System.out.println("i "+i);
            for (int j = nearestPoint[1] - 1; j <= nearestPoint[1] + 1; j++) {
//                System.out.println("j "+j);
                //float distnew = (float)Math.sqrt((xy[0]-xGrid[i][j])*(xy[0]-xGrid[i][j])+(xy[1]-yGrid[i][j])*(xy[1]-yGrid[i][j]));
                float distnew = (float) distanceEuclid2(xy[0], xy[1], xGrid[i][j], yGrid[i][j], "WGS84");
//                System.out.println("distnew "+distnew);

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

//        for (double[] d : allDists)
//            System.out.println(Arrays.toString(d));

        //int nPoints = nearList.length;
        int nPoints = 4;
        System.arraycopy(allDists, 0, nearList, 0, nPoints);

        //System.out.printf("In Particle.nearestCentroid "+nearest+"\n");
        return nearList;
    }

    /**
     * @param xy
     * @param xGrid
     * @param yGrid
     * @param nearestPoint
     * @return
     */
    public static double[][] nearestListROMS2(double[] xy, float[][] xGrid, float[][] yGrid, int[] nearestPoint) {
        //double[][] allDists = new double[9][3];
        double[][] nearList = new double[5][3];

        if (nearestPoint == null) {
            //System.out.println("nearestListROMS2: null nearestPoint supplied");
            nearestPoint = nearestROMSGridPoint((float) xy[0], (float) xy[1], xGrid, yGrid, null);
        }

        // Get the index of the "top-left" corner of the containing element
        int[] whichElem = whichROMSElement((float) xy[0], (float) xy[1], xGrid, yGrid, nearestPoint);
        //System.out.println("nearestListROMS2: nearestPoint (supplied or calc'd) "+nearestPoint[0]+" "+nearestPoint[1]+" whichElem (calc'd) "+whichElem[0]+" "+whichElem[1]);
        if (whichElem[0] == -1) {
            //System.out.println("Distance to nearestPoint: "+distanceEuclid2(xy[0], xy[1], xGrid[nearestPoint[0]][nearestPoint[1]], yGrid[nearestPoint[0]][nearestPoint[1]], "WGS84"));
            nearestPoint = nearestROMSGridPoint((float) xy[0], (float) xy[1], xGrid, yGrid, null);
            //System.out.println("Calc'd a new nearest point: "+nearestPoint[0]+" "+nearestPoint[1]);
            whichElem = whichROMSElement((float) xy[0], (float) xy[1], xGrid, yGrid, nearestPoint);
            //System.out.println("nearestListROMS2: whichElem (calc'd a second time) "+whichElem[0]+" "+whichElem[1]);
            //System.exit(1);
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
     * @param x
     * @param y
     * @param lon
     * @param lat
     * @param nearP the nearest ROMS grid point
     * @return
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
//        if (which[0]==-1)
//        {
//            System.err.println("Particle.whichROMSElement: Location not found in any of 4 elements surrounding nearest point");
//        }

        return which;
    }


    /**
     * @param nrList
     * @param hour
     * @param u
     * @param v
     * @param depLayer
     * @return
     */
    public static double[] velocityFromNearestList(double[][] nrList, int hour, float u[][][], float v[][][], float w[][][], int depLayer, boolean verticalDynamics) {
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
     *
     * @param nrListU
     * @param nrListV
     * @param hour
     * @param u
     * @param v
     * @return
     */
    public static double[] velocityFromNearestListROMS(double[][] nrListU, double[][] nrListV, int hour, float u[][][], float v[][][]) {
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

            //System.out.printf("tt %d i %d",tt,i);
//            System.out.printf(" --U-- elem %d %d dist %.4f weight %.4e --- vel = %.4f\n",
//                (int)nrListU[i][0],(int)nrListU[i][1],nrListU[i][2],weightsU[i],u[tt][(int)nrListU[i][0]][(int)nrListU[i][1]]);
//            System.out.printf(" --V-- elem %d %d dist %.4f weight %.4e --- vel = %.4f\n",
//                (int)nrListV[i][0],(int)nrListV[i][1],nrListV[i][2],weightsV[i],v[tt][(int)nrListV[i][0]][(int)nrListV[i][1]]);


            usum = usum + weightsU[i] * u[hour][(int) nrListU[i][0]][(int) nrListU[i][1]];
            vsum = vsum + weightsV[i] * v[hour][(int) nrListV[i][0]][(int) nrListV[i][1]];

            sumWeightsU = sumWeightsU + weightsU[i];
            sumWeightsV = sumWeightsV + weightsV[i];
        }
        usum = usum / sumWeightsU;
        vsum = vsum / sumWeightsV;
        velocity[0] = usum;
        velocity[1] = vsum;
//        System.out.printf("Interpolated Velocity = %.4f %.4f\n",velocity[0],velocity[1]);
        return velocity;
    }


    /**
     * Euclidean distance between two points
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @return
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
     *
     * @param x1
     * @param y1
     * @param x2
     * @param y2
     * @param coordRef
     * @return
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
     *
     * @param xy
     * @param elemPart
     * @param nodexy
     * @param trinodes
     * @param neighbours
     * @return
     */
    public static int[] findContainingElement(double[] xy, int elemPart,
                                              float[][] nodexy, int[][] trinodes, int[][] neighbours, boolean checkAll) {
        int[] c = new int[6];
        int[] elems = new int[1];
        elems[0] = elemPart;

        //System.out.println("findContainingElement: elems[0]="+elems[0]);
        int whereami = whichElement(xy, elems, nodexy, trinodes);
        //System.out.println("findContainingElement: whereami="+whereami);
        c[1] = 1;

        // If this fails, look in my neighbours
        if (whereami == -1) {
            //int[] elems0 = neighbours[elemPart];
            c[2] = 1;
            //System.out.println("Looking in neighbours");
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
                //System.out.println("Looking in neighbours' neighbours");
                c[3] = 1;
                whereami = whichElement(xy, elems1, nodexy, trinodes);

                // If fails, look in neighbours' neighbours' neighbours
                // Shouldn't really need to get this far? Maybe you should.
                if (whereami == -1) {
                    //System.out.println("WARNING: Looking in neighbours' neighbours' neighbours!!!");
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
                        //System.out.printf("Looking in nearest 80000.... ");
                        c[5] = 1;
                        int[] allelems = IntStream.rangeClosed(0, trinodes[0].length - 1).toArray();
                        whereami = whichElement(xy, allelems, nodexy, trinodes);
//                        if (whereami == -1)
//                        {
//                            System.out.printf("Movement placed particle outside mesh\n");
//                        }
                    }
                }
            }
        }
        c[0] = whereami;
//        System.out.printf("%d %d %d %d %d %d %d\n",elemPart,c[0],c[1],c[2],c[3],c[4],c[5]);
//        if (c[0] == 0 || c[0] == -1)
//        {
//            System.out.println("Element out of bounds, fixing location");
//            System.out.printf("whereami=0 --- %.6e %.6e %d\n",newlocx,newlocy,elemPart);
//        }
        return c;
    }

    /**
     * Find which element a particle resides within (edge checking)
     *
     * @param xy
     * @param elems
     * @param nodexy
     * @param trinodes
     * @return
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
                    System.err.println(i + " " + j + " " + elems[i] + " " + trinodes[j][elems[i]] + " " + nodexy[0][trinodes[j][elems[i]]]);
                }
            }
            // check whether (x,y) lies within this
            //fprintf('check %d\n', possibleElems(i));

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
        //System.out.printf("whichElement: particle in %d\n", which);
        return which;
    }

    /**
     * Calculate a list of neighbour cells based on a specified (e.g. particle) location
     * (0: containing element, 1-3: neighbours of the containing element), and calculate
     * the Euclidean distances.
     *
     * @param xy
     * @param elemPart0
     * @param meshes
     * @param meshPart
     * @param coordRef
     * @return
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
        float[][] nearestLayers = new float[2][2];

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
       /* System.out.println("Neighborhood:");
        for (int i = 0; i < nrRow; i++) {
            System.out.printf("%d: %.1f, %d   ", (int) nrList[i][0], nrList[i][1], (int) nrList[i][2]);
        }
        System.out.println("");*/

        return nrList;
    }

    /**
     * Update particle location using an RK4 integration step.
     *
     * @param hydroFields
     * @param meshes
     * @param hour
     * @param step
     * @param dt
     * @param stepsPerStep
     * @param coordRef
     * @return
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
                System.out.println("Error: vertical dynamics not implemented for ROMS");
                System.exit(0);
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
        //System.out.println("Half step (k1 -> k2)");
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
        //System.out.println("Half step (k2 -> k3)");
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
        //System.out.println("End step (k3 -> k4)");
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
        if (k1[0] == 0 || k2[0] == 0 || k3[0] == 0 || k4[0] == 0) {
//            System.out.printf("RK4 attempt to step out of mesh: Location = [%.6e,%.6e], Element = %d\n",
//                    this.getLocation()[0],this.getLocation()[1],elemPart);
        } else {
            for (int i = 0; i < nDims; i++) {
                advectStep[i] = (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6.0;
            }
        }

        // Check that things approximately tie up with corresponding Euler step distance
        // Generally they do, but sometimes they really don't (calculating Euler steps with
        // similar time step to RK4 is not advisable as sensitive to complex current features).
        // Further the approximation for printing below doesn't include the time interpolation.
//        System.out.printf("RK4 calc (%d): advectStep=[%.3e,%.3e] vel0=[%.3e,%.3e] dt=%.3e vel0*dt=[%.3e,%.3e]\n",
//                step,advectStep[0],advectStep[1],vel[0],vel[1],dt,vel[0]*dt,vel[1]*dt);
        return advectStep;
    }

    /**
     * Calculate the correction steps required for the RK4 algorithm
     *
     * @param xy            Location of new velocity to be used (current location plus spatial step)
     * @param elemPart
     * @param depLayer
     * @param timeStepAhead Time step ahead i.e. "1.0/2.0" if half-step ahead, "1.0" if full step ahead
     * @param meshes
     * @param meshPart
     * @param hydroFields
     * @param hour
     * @param step
     * @param dt
     * @param stepsPerStep
     * @param coordRef
     * @return
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
                System.out.println("Error in rk4Step: vertical dynamics not implemented for ROMS");
                System.exit(0);
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
            xyz_step[i] = dt * (vel[i] + ((double) (step + timeStepAhead) / (double) stepsPerStep) * (velplus1[i] - vel[i]));
        }

        return xyz_step;
    }

    /**
     * Compute an Euler integration step for particle movement
     *
     * @param hydroFields
     * @param meshes
     * @param hour
     * @param step
     * @param dt
     * @param stepsPerStep
     * @param coordRef
     * @return
     */
    public double[] eulerStep(List<HydroField> hydroFields, // velocities
                              List<Mesh> meshes,     // other mesh info
                              int hour, int step, double dt,                                  // locate particle in space and time
                              int stepsPerStep, String coordRef, boolean verticalDynamics) {
        if(verticalDynamics) {
            System.out.println("Error: vertical dynamics not implemented with Particle.eulerStep()");
            System.exit(0);
        }
        int elemPart = this.getElem();
        int meshPart = this.getMesh();
        int dep = this.getDepthLayer();

        double[] thisLoc = this.getLocation();

        //System.out.printf("RK4Step: Location = [%.6e,%.6e], Element = %d\n",this.getLocation()[0],this.getLocation()[1],elemPart);
        double[] advectStep = new double[2];
        // compute velocities at start and end of entire step, at the new location
        double[] vel = new double[2];
        double[] velplus1 = new double[2];

//        this.nrList = neighbourCellsList(this.getLocation(), elemPart, 
//            meshPart, meshes, allelems,coordRef);

        if (meshes.get(meshPart).getType().equalsIgnoreCase("FVCOM") || meshes.get(meshPart).getType().equalsIgnoreCase("ROMS_TRI")) {
            double[][] xNrList = neighbourCellsList(thisLoc, elemPart, meshPart, meshes, coordRef);
            //this.setNrListToNeighbourCells(neighbours,uvnode);

            // 2. Compute k_1 (spatial interpolation at start of step)
            //System.out.println("Start step");
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
//            System.out.println("in eulerStep, about to find nearest points");
            // nearestROMSGridPoint ---> whichROMSelement --->

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

        // Sense check calculated velocity and advection step (should be better when step close to 0)
//        System.out.printf("Euler calc (%d): advectStep=[%.3e,%.3e] vel0=[%.3e,%.3e] dt=%.3e vel0*dt=[%.3e,%.3e]\n",
//            step,advectStep[0],advectStep[1],vel[0],vel[1],dt,vel[0]*dt,vel[1]*dt);

        return advectStep;
    }

    /**
     * Calculate velocity at a particle's location, given known
     * @param tt
     * @param u
     * @param v
     * @return
     */
//    public double[] velPart(int tt, float u[][][], float v[][][])
//    {
//        double[] velocity = new double[2];
//        velocity = velocityFromNearestList(this.nrList,tt,u,v,1);
//        return velocity;
//    }

    /**
     * Method to compute velocity at a certain time and location, directly from an element of a velocity field.
     *
     * This does the same thing as "velPart", except it generates a list of cells
     * to spatially interpolate from.
     *
     * @param tt
     * @param xy
     * @param elemPart0
     * @param u
     * @param v
     * @param neighbours
     * @param uvnode
     * @param nodexy
     * @param trinodes
     * @param allelems
     * @return
     */
//    public static double[] velocityInterpAtLoc(int tt, double xy[], int elemPart0, float u[][][], float v[][][],
//            int[][] neighbours, float[][] uvnode, float[][] nodexy, int[][] trinodes, int[] allelems, int depLayer, String coordRef)
//    {
//        double[] velocity = new double[2];
//        // elemPart0 is the starting search element. Set = 1 if outside range.
//        if (elemPart0 < 1 || elemPart0 > allelems.length)
//        {
//            elemPart0 = 1;
//            System.out.println("velocityInterpAtLoc picking default start element = 1");
//        }
//        double[][] cellList = neighbourCellsList(xy, elemPart0, neighbours, uvnode, nodexy, trinodes, allelems, coordRef);
//        velocity = velocityFromNearestList(cellList,tt,u,v,1);
//        return velocity;
//    }

    /**
     * Set the particle's nearestList (NrList) to be identical to the list of elements which 
     * neighbour the element that the particle is actually in
     * @param neighbours
     * @param uvnode
     */
//    public void setNrListToNeighbourCells(int[][] neighbours, float[][] uvnode)
//    {      
//        // distance to elem
//        //int elem = nearestCentroid(this.xy[0],this.xy[1],uvnode);
//        int elem = this.elem;
//        this.nrList[0][0] = elem;
//        this.nrList[0][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[0][elem],uvnode[1][elem]);
//        // distance to neighbouring elems
//        this.nrList[1][0] = neighbours[0][elem];
//        this.nrList[1][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[0][neighbours[0][elem]],uvnode[1][neighbours[0][elem]]);
//        this.nrList[2][0] = neighbours[1][elem];
//        this.nrList[2][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[0][neighbours[1][elem]],uvnode[1][neighbours[1][elem]]);
//        this.nrList[3][0] = neighbours[2][elem];
//        this.nrList[3][1] = distanceEuclid(this.xy[0],this.xy[1],uvnode[0][neighbours[2][elem]],uvnode[1][neighbours[2][elem]]);   
//        this.nrList[4][0] = 0;
//        this.nrList[4][1] = 1000000;     
//    }

}