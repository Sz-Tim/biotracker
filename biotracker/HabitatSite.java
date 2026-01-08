/*
 * Class to represent a site at which particles can be released/settle in particle tracking
 */
package biotracker;

//import java.util.Arrays;

import java.util.*;
import java.util.stream.IntStream;

import static biotracker.Particle.distanceEuclid2;

/**
 * @author SA01TA
 */
public class HabitatSite {
    private final String ID; // A unique identifier for a site (string e.g. SEPA site IDs for fish farms
    private float[] xy; // Coordinate for site
    private float depth;
    private double scale; // A scale factor to be applied if required (e.g. site biomass for fish farms)
    private int containingMesh; // the mesh containing the habitat (possibly multiple in case of polygon?)
    private String containingMeshType;
    private int nearestFVCOMCentroid; // the mesh u,v point that is closest to this location
    private int containingFVCOMElem; // the element containing the habitat (possibly multiple in case of polygon?)
    private int[] nearestROMSGridPointU;
    private int[] containingROMSElemU;
    private int[] nearestROMSGridPointV;
    private int[] containingROMSElemV;
    private double[] envConditionDay;
    private double[] envConditionSum;
    private int[] envConditionCount;
    private int[] envConditionCountDay;

    public HabitatSite(String ID, float x, float y, float depth, float scale, List<Mesh> meshes, RunProperties rp) {
        this.ID = ID;
        this.xy = new float[]{x, y};
        this.scale = scale;
        // [u, v, w, uv, salinity, temp, light, Vh, km, stokesU0, stokesV0, stokesUV0, Hsig, waveDir, waveT]
        this.envConditionDay = new double[15];
        this.envConditionSum = new double[15];
        this.envConditionCount = new int[15];
        this.envConditionCountDay = new int[15];

        double[] xy2 = new double[]{this.xy[0], this.xy[1]};

        // Find out which mesh the particle is in, based on bounding box.
        // We make the assumption that meshes appear in the list in their
        // order of precedence, and allocate particles to the first mesh that
        // contains their location
        // Assume as a default that the habitat site is in the first mesh, if not change it to something else
        this.containingMesh = 0;
        boolean insideMesh = false;
        // Use this to force checking of FVCOM mesh (single mesh scenario)
        boolean exitIfNotInMesh = true;

        for (int m = 0; m < meshes.size(); m++) {
            if (meshes.get(m).isInMesh(xy2, true, true, null)) {
                this.containingMesh = m;
                insideMesh = true;
                break;
            }
        }

        this.containingMeshType = "FVCOM";
        this.nearestFVCOMCentroid = -1;
        this.containingFVCOMElem = -1;

        //noinspection ConstantConditions
        if (!insideMesh && exitIfNotInMesh) {
            if(rp.verboseSetUp) {
                System.out.println(ID + " " + x + " " + y + " Habitat site not within any provided mesh --- : Not creating site");
            }
            this.containingMeshType = "NONE";
        } else {
            if (meshes.get(this.containingMesh).getType().equalsIgnoreCase("FVCOM") || meshes.get(this.containingMesh).getType().equalsIgnoreCase("ROMS_TRI")) {
                this.containingMeshType = meshes.get(this.containingMesh).getType();
                this.nearestFVCOMCentroid = Particle.nearestCentroid(xy[0], xy[1], meshes.get(this.containingMesh).getUvnode());
                this.depth = meshes.get(this.containingMesh).getDepthUvnode()[this.nearestFVCOMCentroid];

                //noinspection ConstantConditions
                if (!insideMesh) {
                    double d1 = distanceEuclid2(xy[0], xy[1], meshes.get(this.containingMesh).getUvnode()[0][this.nearestFVCOMCentroid], meshes.get(this.containingMesh).getUvnode()[1][this.nearestFVCOMCentroid], rp.coordOS);
                    if(rp.verboseSetUp) {
                        System.out.println("Habitat site (" + xy[0] + "," + xy[1] + ") outside mesh. Nearest centroid: " + this.nearestFVCOMCentroid
                                + " (" + meshes.get(this.containingMesh).getUvnode()[0][this.nearestFVCOMCentroid] + "," + meshes.get(this.containingMesh).getUvnode()[1][this.nearestFVCOMCentroid] + ")"
                                + " distance = " + d1);
                    }
                    if (d1 < 5000) {
                        if(rp.verboseSetUp) {
                            System.out.println("Within 5000m from mesh edge; moving to nearest element centroid");
                        }
                        this.xy = new float[]{meshes.get(this.containingMesh).getUvnode()[0][this.nearestFVCOMCentroid], meshes.get(this.containingMesh).getUvnode()[1][this.nearestFVCOMCentroid]};
                        xy2 = new double[]{this.xy[0], this.xy[1]};
                    } else {
                        System.err.println("More than 5000m from mesh edge; cannot create habitat site");
                    }
                }
                this.containingFVCOMElem = Particle.whichElement(xy2,
                        IntStream.rangeClosed(0, meshes.get(this.containingMesh).getUvnode()[0].length - 1).toArray(),
                        meshes.get(this.containingMesh).getNodexy(), meshes.get(this.containingMesh).getTrinodes());
                int a = this.containingFVCOMElem;
                if(rp.verboseSetUp) {
                    System.out.println(ID + " " + x + " " + y + " " + a);
                }
            } else if (meshes.get(this.containingMesh).getType().equalsIgnoreCase("ROMS")) {
                this.containingMeshType = "ROMS";
                System.out.println("habitat site in ROMS mesh");

                this.nearestROMSGridPointU = Particle.nearestROMSGridPoint(xy[0], xy[1],
                        meshes.get(this.containingMesh).getLonU(), meshes.get(this.containingMesh).getLatU(), null);
                this.nearestROMSGridPointV = Particle.nearestROMSGridPoint(xy[0], xy[1],
                        meshes.get(this.containingMesh).getLonV(), meshes.get(this.containingMesh).getLatV(), null);

                System.out.println("found nearest point");
                // More to do here to turn the nearest grid point into the containing element
                this.containingROMSElemU = Particle.whichROMSElement(xy[0], xy[1],
                        meshes.get(this.containingMesh).getLonU(), meshes.get(this.containingMesh).getLatU(),
                        this.nearestROMSGridPointU);
                this.containingROMSElemV = Particle.whichROMSElement(xy[0], xy[1],
                        meshes.get(this.containingMesh).getLonV(), meshes.get(this.containingMesh).getLatV(),
                        this.nearestROMSGridPointV);
                System.out.println("found which element");

                System.out.println("Habitat site --- nearestROMSpointU=[" + this.nearestROMSGridPointU[0] + "," + this.nearestROMSGridPointU[1] + "]("
                        + meshes.get(this.containingMesh).getLonU()[this.nearestROMSGridPointU[0]][this.nearestROMSGridPointU[1]] + ","
                        + meshes.get(this.containingMesh).getLatU()[this.nearestROMSGridPointU[0]][this.nearestROMSGridPointU[1]] + ")"
                        + " --- containingROMSElemU=[" + this.containingROMSElemU[0] + "," + this.containingROMSElemU[1] + "]");
                System.out.println("Habitat site --- nearestROMSpointV=[" + this.nearestROMSGridPointV[0] + "," + this.nearestROMSGridPointV[1] + "]("
                        + meshes.get(this.containingMesh).getLonV()[this.nearestROMSGridPointV[0]][this.nearestROMSGridPointV[1]] + ","
                        + meshes.get(this.containingMesh).getLatV()[this.nearestROMSGridPointV[0]][this.nearestROMSGridPointV[1]] + ")"
                        + " --- containingROMSElemV=[" + this.containingROMSElemV[0] + "," + this.containingROMSElemV[1] + "]");
            }
        }

    }

    @Override
    public String toString() {
        String details = this.ID + "," + this.xy[0] + "," + this.xy[1] + "," + this.depth + "," + this.containingMesh;
        if (this.containingMeshType.equalsIgnoreCase("FVCOM") || this.containingMeshType.equalsIgnoreCase("ROMS_TRI")) {
            return details + "," + this.nearestFVCOMCentroid + "," + this.containingFVCOMElem + "," + this.containingMeshType;
        } else if (this.containingMeshType.equalsIgnoreCase("ROMS")) {
            return details + "," +
                    "U:" + this.nearestROMSGridPointU[0] + ";" + this.nearestROMSGridPointU[1] + ";" +
                    "V:" + this.nearestROMSGridPointV[0] + ";" + this.nearestROMSGridPointV[1] + "," +
                    "U:" + this.containingROMSElemU[0] + ";" + this.containingROMSElemU[1] + ";" +
                    "V:" + this.containingROMSElemV[0] + ";" + this.containingROMSElemV[1] + "," +
                    "ROMS";
        }
        return details;
    }

    public String avgIntervalToString() {
        StringBuilder intervalAvg = new StringBuilder();
        double[] avgEnv = this.getAvgEnvConditionDay();
        for (int i = 0; i < avgEnv.length; i++) {
            if (avgEnv[i] > -9999) {
                intervalAvg.append(avgEnv[i]);
                if (i < avgEnv.length - 1) {
                    intervalAvg.append(",");
                }
            }
        }
        return intervalAvg.toString();
    }

    public String avgOverallToString() {
        StringBuilder overallAvg = new StringBuilder();
        double[] avgEnv = this.getAvgEnvCondition();
        for (int i = 0; i < avgEnv.length; i++) {
            if(avgEnv[i] > -9999) {
                overallAvg.append(avgEnv[i]);
                if (i < avgEnv.length - 1) {
                    overallAvg.append(",");
                }
            }
        }
        return overallAvg.toString();
    }

    public String getID() {
        return this.ID;
    }

    public float[] getLocation() {
        return this.xy;
    }

    public double getDepth() {
        return this.depth;
    }

    public double getScale() {
        return this.scale;
    }

    public void setScale(double scale) {
        this.scale = scale;
    }

    public double[] getEnvConditionSum() {
        return this.envConditionSum;
    }

    public int[] getEnvConditionCount() {
        return this.envConditionCount;
    }

    public double[] getEnvConditionDay() {
        return this.envConditionDay;
    }

    public int[] getEnvConditionCountDay() { return this.envConditionCountDay; }

    public void setEnvConditionDay(double[] envConditionDay) {
        this.envConditionDay = envConditionDay;
    }

    public void setEnvConditionCountDay(int[] envConditionCountDay) {
        this.envConditionCountDay = envConditionCountDay;
    }

    public void addEnvCondition(double[] currentConditions) {
        for (int i=0; i < this.envConditionSum.length; i++) {
            this.envConditionSum[i] += currentConditions[i];
            this.envConditionDay[i] += currentConditions[i];
            this.envConditionCount[i] ++;
            this.envConditionCountDay[i] ++;
        }
    }

    public double[] getAvgEnvCondition() {
        double[] avgCondition = new double[this.envConditionSum.length];
        for (int i=0; i < this.envConditionSum.length; i++) {
            avgCondition[i] = this.getEnvConditionSum()[i] / this.envConditionCount[i];
        }
        return avgCondition;
    }

    public double[] getAvgEnvConditionDay() {
        double[] avgCondition = new double[this.envConditionDay.length];
        for (int i=0; i < this.envConditionDay.length; i++) {
            avgCondition[i] = this.getEnvConditionDay()[i] / this.envConditionCountDay[i];
        }
        return avgCondition;
    }

    public int getContainingMesh() {
        return this.containingMesh;
    }

    public String getContainingMeshType() {
        return this.containingMeshType;
    }

    public int getContainingFVCOMElem() {
        return this.containingFVCOMElem;
    }

    public int[] getContainingROMSElemU() {
        return this.containingROMSElemU;
    }

    public int[] getContainingROMSElemV() {
        return this.containingROMSElemV;
    }

    public int[] getNearestROMSPointU() {
        return this.nearestROMSGridPointU;
    }

    public int[] getNearestROMSPointV() {
        return this.nearestROMSGridPointV;
    }

}
