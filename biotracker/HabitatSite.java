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
    private float scale; // A scale factor to be applied if required (e.g. site biomass for fish farms)
    private int containingMesh; // the mesh containing the habitat (possibly multiple in case of polygon?)
    private String containingMeshType;
    private int nearestFVCOMCentroid; // the mesh u,v point that is closest to this location
    private int containingFVCOMElem; // the element containing the habitat (possibly multiple in case of polygon?)
    private int[] nearestROMSGridPointU;
    private int[] containingROMSElemU;
    private int[] nearestROMSGridPointV;
    private int[] containingROMSElemV;
    private double[] envConditionSum;
    private int[] envConditionCount;

    public HabitatSite(String ID, float x, float y, float depth, float scale, List<Mesh> meshes, RunProperties rp) {
        this.ID = ID;
        this.xy = new float[]{x, y};
        this.scale = scale;
        this.envConditionSum = new double[4]; // [u, v, w, salinity]
        this.envConditionCount = new int[4]; // [u, v, w, salinity]

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
                    double d1 = distanceEuclid2(xy[0], xy[1], meshes.get(this.containingMesh).getUvnode()[0][this.nearestFVCOMCentroid], meshes.get(this.containingMesh).getUvnode()[1][this.nearestFVCOMCentroid], rp.coordRef);
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
        String details = this.ID + "\t" + this.xy[0] + "\t" + this.xy[1] + "\t" + this.containingMesh;
        if (this.containingMeshType.equalsIgnoreCase("FVCOM") || this.containingMeshType.equalsIgnoreCase("ROMS_TRI")) {
            details = this.ID + "\t" + this.xy[0] + "\t" + this.xy[1] + "\t" + this.depth + "\t"
                    + this.containingMesh + "\t" + this.nearestFVCOMCentroid + "\t" + this.containingFVCOMElem
                    + "\t" + this.containingMeshType
                    + "\t" + this.getAvgEnvCondition()[0] + "\t" + this.getAvgEnvCondition()[1]
                    + "\t" + this.getAvgEnvCondition()[2] + "\t" + this.getAvgEnvCondition()[3] ;
        } else if (this.containingMeshType.equalsIgnoreCase("ROMS")) {
            details = this.ID + "\t" + this.xy[0] + "\t" + this.xy[1] + "\t" + this.containingMesh + " U_grid: ("
                    + this.nearestROMSGridPointU[0] + "," + this.nearestROMSGridPointU[1] + ") ("
                    + this.containingROMSElemU[0] + "," + this.nearestROMSGridPointU[1] + ")";
        }
        return details;
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

    public void setScale(float scale) {
        this.scale = scale;
    }

    public double[] getEnvConditionSum() {
        return this.envConditionSum;
    }

    public int[] getEnvConditionCount() {
        return this.envConditionCount;
    }

    public void addEnvCondition(double[] currentConditions) {
        for (int i=0; i < this.envConditionSum.length; i++) {
            this.envConditionSum[i] += currentConditions[i];
            this.envConditionCount[i] ++;
        }
    }

    public double[] getAvgEnvCondition() {
        double[] avgCondition = new double[this.envConditionSum.length];
        for (int i=0; i < this.envConditionSum.length; i++) {
            avgCondition[i] = this.getEnvConditionSum()[i] / this.envConditionCount[i];
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
