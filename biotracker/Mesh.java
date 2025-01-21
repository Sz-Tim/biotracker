/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

import java.awt.geom.Path2D;
import java.lang.reflect.Array;
import java.util.ArrayList;

//import extUtils.ConcaveHull;

/**
 * @author SA01TA
 */
public class Mesh {

    // Define the mesh type
    private final String meshType;

    // Variables which are used in case of FVCOM triangular mesh
    private float[][] uvnode; // element centroid locations
    private float[][] nodexy; // element node locations
    private float[] depthUvnode; // depth below geoid at each element centroid
    private float[] depthNodexy; // depth below geoid at each mesh (xy)node; positive values
    private int[][] trinodes; // the nodes at the corners of each element
    private int[][] neighbours;
    private int[] openBoundaryNodes;
    private ArrayList<Integer> openBoundaryElems;
    private int[] boundaryNodes;
    private int adjoiningElement; // ONE element index at the boundary of the nested mesh (element of mesh 1 at boundary of mesh 0)

    // Variables used in the case of ROMS grid
    private float[][] lon_u;
    private float[][] lat_u;
    private float[][] lon_v;
    private float[][] lat_v;
    private float[][] lon_rho;
    private float[][] lat_rho;
    private float[][] depth_rho;

    private final int[][] rangeUV = new int[2][2];

    // Variables used in all cases
    private float[] siglay;
    private float[] siglev;
    private float[][] convexHull;
    private Path2D.Float convexHullPath;

    /**
     * Create a Mesh from a NetCDF file, and supplementary text files detailing
     * the open and closed boundary nodes.
     */
    public Mesh(String meshFilename, String type, boolean coordOS) {
        System.out.println("Reading mesh file: " + meshFilename);

        meshType = type;
        if (type.equals("ROMS_TRI") || type.equals("FVCOM")) {
            if (!coordOS) {
                uvnode = IOUtils.readNetcdfFloat2D(meshFilename, "uvnode", null, null);
                nodexy = IOUtils.readNetcdfFloat2D(meshFilename, "nodexy", null, null);
                if (uvnode == null) {
                    float[] lonc = IOUtils.readNetcdfFloat1D(meshFilename, "lonc");
                    if (lonc[0] > 180) {
                        for (int i = 0; i < lonc.length; i++) {
                            lonc[i] = lonc[i] - 360;
                        }
                    }
                    float[] latc = IOUtils.readNetcdfFloat1D(meshFilename, "latc");
                    uvnode = new float[2][lonc.length];
                    for (int i = 0; i < lonc.length; i++) {
                        uvnode[0][i] = lonc[i];
                        uvnode[1][i] = latc[i];
                    }

                    float[] lon = IOUtils.readNetcdfFloat1D(meshFilename, "lon");
                    if (lon[0] > 180) {
                        for (int i = 0; i < lon.length; i++) {
                            lon[i] = lon[i] - 360;
                        }
                    }
                    float[] lat = IOUtils.readNetcdfFloat1D(meshFilename, "lat");
                    nodexy = new float[2][lon.length];
                    for (int i = 0; i < lon.length; i++) {
                        nodexy[0][i] = lon[i];
                        nodexy[1][i] = lat[i];
                    }
                }
            } else {
                uvnode = IOUtils.readNetcdfFloat2D(meshFilename, "uvnode_os", null, null);
                nodexy = IOUtils.readNetcdfFloat2D(meshFilename, "nodexy_os", null, null);
            }

            depthUvnode = IOUtils.readNetcdfFloat1D(meshFilename, "depthUvnode");
            if (depthUvnode == null) {
                depthUvnode = IOUtils.readNetcdfFloat1D(meshFilename, "h_center");
            }
            depthNodexy = IOUtils.readNetcdfFloat1D(meshFilename, "depthNodexy");
            if (depthNodexy == null) {
                depthNodexy = IOUtils.readNetcdfFloat1D(meshFilename, "h");
            }
            trinodes = IOUtils.readNetcdfInteger2D(meshFilename, "trinodes");
            if (trinodes == null) {
                trinodes = IOUtils.readNetcdfInteger2D(meshFilename, "nv");
            }
            neighbours = IOUtils.readNetcdfInteger2D(meshFilename, "nbe");
            siglay = IOUtils.readNetcdfFloat1D(meshFilename, "siglay");
            if (siglay == null) {
                System.out.println("  Siglay is defined on mesh, read and extract first line");
                float[][] siglayGrid = IOUtils.readNetcdfFloat2D(meshFilename, "siglay", null, null);
                siglay = new float[siglayGrid.length];
                for (int i = 0; i < siglayGrid.length; i++) {
                    siglay[i] = siglayGrid[i][0];
                }
            }
            // Ensure values of siglay & siglev are positive
            for (int i = 0; i < siglay.length; i++) {
                if (siglay[i] < 0) {
                    siglay[i] *= -1;
                }
            }
            siglev = IOUtils.readNetcdfFloat1D(meshFilename, "siglev");
            for (int i = 0; i < siglev.length; i++) {
                if (siglev[i] < 0) {
                    siglev[i] *= -1;
                }
            }

            openBoundaryNodes = IOUtils.readNetcdfInteger1D(meshFilename, "boundaryNodesOpen");
            boundaryNodes = IOUtils.readNetcdfInteger1D(meshFilename, "boundaryNodesAll");
            // reduce node/element IDs in files generated by matlab by one (loops start at zero, not one as in matlab)
            for (int i = 0; i < trinodes.length; i++) {
                for (int j = 0; j < trinodes[1].length; j++) {
                    trinodes[i][j]--;
                    if (neighbours[i][j] > 0) {
                        neighbours[i][j]--;
                    }
                }
            }
            if (openBoundaryNodes != null) {
                for (int i = 0; i < openBoundaryNodes.length; i++) {
                    openBoundaryNodes[i]--;
                }
                // find elements with at least one trinode as an open boundary
                System.out.println("Finding boundary elements...");
                openBoundaryElems = findOpenBoundaryElements();
            } else {
                System.err.println("No openBoundaryNode information in mesh file");
            }
            if (boundaryNodes != null) {
                for (int i = 0; i < boundaryNodes.length; i++) {
                    boundaryNodes[i]--;
                }
            } else {
                System.err.println("No boundaryNode information in mesh file");
            }

            System.out.println("neighbours: " + neighbours.length + " " + neighbours[1].length);

            float[][] uvnodeT = new float[uvnode[1].length][2];
            for (int i = 0; i < uvnode[1].length; i++) {
                uvnodeT[i][0] = uvnode[0][i];
                uvnodeT[i][1] = uvnode[1][i];
            }
            convexHull = ConvexHull.convexHull(uvnodeT);
            IOUtils.writeFloatArrayToFile(convexHull, "convexHull_" + type + ".dat", false, false);
        } else if (type.equals("ROMS")) {
            rangeUV[0][0] = 300;
            rangeUV[0][1] = 748;
            rangeUV[1][0] = 250;
            rangeUV[1][1] = 900;

            int[] origin = new int[]{rangeUV[0][0], rangeUV[1][0]};
            int[] shape = new int[]{rangeUV[0][1] - rangeUV[0][0], rangeUV[1][1] - rangeUV[1][0]};

            lon_u = IOUtils.readNetcdfFloat2D(meshFilename, "lon_u", origin, shape);
            lat_u = IOUtils.readNetcdfFloat2D(meshFilename, "lat_u", origin, shape);
            lon_v = IOUtils.readNetcdfFloat2D(meshFilename, "lon_v", origin, shape);
            lat_v = IOUtils.readNetcdfFloat2D(meshFilename, "lat_v", origin, shape);
            lon_rho = IOUtils.readNetcdfFloat2D(meshFilename, "lon_rho", origin, shape);
            lat_rho = IOUtils.readNetcdfFloat2D(meshFilename, "lat_rho", origin, shape);
            siglay = IOUtils.readNetcdfFloat1D(meshFilename, "s_rho");

            float[][] mask_u = IOUtils.readNetcdfFloat2D(meshFilename, "mask_u", origin, shape);
            float[][] mask_v = IOUtils.readNetcdfFloat2D(meshFilename, "mask_v", origin, shape);

            // Set the index limits. These are the ones to be used in the case of reading the full ROMS domain
            float[][] ux1 = IOUtils.reshapeFloat(lon_u, lon_u.length * lon_u[0].length, 1);
            float[][] uy1 = IOUtils.reshapeFloat(lat_u, lat_u.length * lat_u[0].length, 1);
            float[][] uXy = new float[lat_u.length * lat_u[0].length][2];
            float[][] vx1 = IOUtils.reshapeFloat(lon_v, lon_v.length * lon_v[0].length, 1);
            float[][] vy1 = IOUtils.reshapeFloat(lat_v, lat_v.length * lat_v[0].length, 1);
            float[][] vXy = new float[lat_v.length * lat_v[0].length][2];
            for (int i = 0; i < lat_v.length * lat_v[0].length; i++) {
                uXy[i] = new float[]{ux1[i][0], uy1[i][0]};
                vXy[i] = new float[]{vx1[i][0], vy1[i][0]};
            }
            convexHull = ConvexHull.convexHull(uXy);
            IOUtils.writeFloatArrayToFile(convexHull, "convexHull_ROMS.dat", false, false);

            // List the corners to verify which way around the netCDF array is read in
            System.out.println("ROMS grid size: " + lon_u.length + " " + lon_u[0].length);
            System.out.println("Corner 1: " + lon_u[0][0] + " " + lat_u[0][0]);
            System.out.println("Corner 2: " + lon_u[lon_u.length - 1][0] + " " + lat_u[lon_u.length - 1][0]);
            System.out.println("Corner 3: " + lon_u[lon_u.length - 1][lon_u[0].length - 1] + " " + lat_u[lon_u.length - 1][lon_u[0].length - 1]);
            System.out.println("Corner 4: " + lon_u[0][lon_u[0].length - 1] + " " + lat_u[0][lon_u[0].length - 1]);
            // increasing first index => increasing longitude
            // increasing second index => decreasing latitude
        }
        convexHullPath = Mesh.pointsToPath(this.getConvexHull());
    }

    private ArrayList<Integer> findOpenBoundaryElements() {
        ArrayList<Integer> obElems = new ArrayList<>();
        for (int openBoundaryNode : openBoundaryNodes) {
            for (int j = 0; j < trinodes[0].length; j++) {
                if (!obElems.contains(j) &&
                        (trinodes[0][j] == openBoundaryNode ||
                                trinodes[1][j] == openBoundaryNode ||
                                trinodes[2][j] == openBoundaryNode)) {
                    obElems.add(j);
                }
            }
        }
        return obElems;
    }

    /**
     * Find sigma layers or levels above and below particle
     *
     * @param particleDepth depth of particle: Particle.getDepth()
     * @param sigmas sigma layers or levels: Mesh.getSigLay() or Mesh.getSigLev()
     * @param localDepth water depth at particle location: Mesh.getDepthUvnode(element)
     * @return float[0=below,1=above][0=index,1=depth]
     */
    public static float[][] findNearestSigmas(double particleDepth, float[] sigmas, float localDepth) {
        int layerBelow = sigmas.length;
        float[] sigDepths = new float[sigmas.length];
        for (int i = 0; i < sigmas.length; i++) {
            sigDepths[i] = sigmas[i] * localDepth;
            if (particleDepth < sigDepths[i]) {
                layerBelow = i;
                break;
            }
        }
        int layerAbove = layerBelow - 1;
        if (layerAbove < 0) {
            layerAbove = 0;
        }
        if (layerBelow == sigmas.length) {
            layerBelow = layerAbove;
        }
        float[][] sigmaInfo = new float[2][2];
        sigmaInfo[0][0] = layerBelow;
        sigmaInfo[0][1] = sigDepths[layerBelow];
        sigmaInfo[1][0] = layerAbove;
        sigmaInfo[1][1] = sigDepths[layerAbove];
        return sigmaInfo;
    }

    public static int[] findNearestSigmaIndexes(double particleDepth, float[] sigmas, float localDepth) {
        int[] layerBelowAbove = new int[2];
        layerBelowAbove[0] = sigmas.length;
        for (int i = 0; i < sigmas.length; i++) {
            if (particleDepth < sigmas[i] * localDepth) {
                layerBelowAbove[0] = i;
                break;
            }
        }
        layerBelowAbove[1] = layerBelowAbove[0] - 1;
        if (layerBelowAbove[1] < 0) {
            layerBelowAbove[1] = 0;
        }
        if (layerBelowAbove[0] == sigmas.length) {
            layerBelowAbove[0] = layerBelowAbove[1];
        }
        return layerBelowAbove;
    }

    public static int[] findNearestSigmaIndexes(int particleDepthLayer, int nSigmas) {
        int[] layerBelowAbove = new int[]{particleDepthLayer+1, particleDepthLayer};
        if (layerBelowAbove[1] < 0) {
            layerBelowAbove[1] = 0;
        }
        if (layerBelowAbove[0] == nSigmas) {
            layerBelowAbove[0] = layerBelowAbove[1];
        }
        return layerBelowAbove;
    }


    public float[][] getNodexy() {
        return nodexy;
    }

    public float[][] getUvnode() {
        return uvnode;
    }

    public float[] getDepthUvnode() {
        return depthUvnode;
    }

    public float[] getDepthNodexy() {
        return depthNodexy;
    }

    public int[][] getTrinodes() {
        return trinodes;
    }

    public int[][] getNeighbours() {
        return neighbours;
    }

    public float[][] getLonU() {
        return lon_u;
    }

    public float[][] getLatU() {
        return lat_u;
    }

    public float[][] getLonV() {
        return lon_v;
    }

    public float[][] getLatV() {
        return lat_v;
    }

    public float[][] getLatRho() {
        return lat_rho;
    }

    public float[][] getLonRho() {
        return lon_rho;
    }

    public int[][] getRange() {
        return rangeUV;
    }

    public float[] getSiglay() {
        return siglay;
    }

    public float[] getSiglev() {
        return siglev;
    }

    public int[] getOpenBoundaryNodes() {
        return openBoundaryNodes;
    }

    public ArrayList<Integer> getOpenBoundaryElems() {
        return openBoundaryElems;
    }

    public int[] getBoundaryNodes() {
        return boundaryNodes;
    }

    public void setAdjoiningElement(int adjoiningElement) { this.adjoiningElement = adjoiningElement; }

    public int getAdjoiningElement() { return adjoiningElement; }

    public float[][] getConvexHull() {
        return convexHull;
    }

    public Path2D.Float getConvexHullPath() {
        return convexHullPath;
    }

    public String getType() {
        return meshType;
    }

    public int getNElems() {
        return depthUvnode.length;
    }

    public int getNNodes() {
        return depthNodexy.length;
    }

    /**
     * Create a Path2D object from a list of x,y points
     */
    public static Path2D.Float pointsToPath(float[][] points) {
        if (points[0].length != 2) {
            System.err.print("Array provided to pointsToPath must have second dimension = 2");
        }
        Path2D.Float path = new Path2D.Float();
        path.moveTo(points[0][0], points[0][1]);
        for (int i = 1; i < points.length; i++) {
            path.lineTo(points[i][0], points[i][1]);
        }
        path.closePath();

        return path;
    }

    /**
     * Is a specific location within the mesh convex hull?
     * Note that this does not guarantee the point is within a mesh element, so we
     * should normally find the nearest/containing node or element in addition to the
     * mesh ID.
     */
    public boolean isInMesh(double[] xy)//, boolean checkElements, int elemLoc)
    {
        Path2D.Float cHull = this.getConvexHullPath();
        return cHull.contains(xy[0], xy[1]);
    }

    public boolean isInMesh(double[] xy, boolean checkElements, boolean checkAllElems, int[] elemLoc) {
        Path2D.Float cHull = this.getConvexHullPath();
        boolean inMesh = cHull.contains(xy[0], xy[1]);
        // Do the element check, if required. This will identify when a point is NOT in the mesh
        // when it has fallen inside the convex hull but is outside any element
        if (checkElements && inMesh) {
            int c = -1;
            if (this.getType().equals("FVCOM") || this.getType().equals("ROMS_TRI")) {
                int eL = 0;
                if (elemLoc != null) {
                    eL = elemLoc[0];
                }
                c = Particle.findContainingElement(xy, eL, this, checkAllElems, 1)[0];
            } else if (this.getType().equals("ROMS")) {
                // Use the U grid to determine whether within the mesh or not.
                // This leads to some particles that are outside the V grid being identified as in the mesh.
                // The same would happen with rho
                // If we can make it fast, check U and V grids, and chuck out if not in one of them?
                // Change isInMesh references in ParallelParticleMover, too
                int[] nearestPointU = Particle.nearestROMSGridPoint((float) xy[0], (float) xy[1], this.getLonU(), this.getLatU(), elemLoc);
                int[] cRomsU = Particle.whichROMSElement((float) xy[0], (float) xy[1], this.getLonU(), this.getLatU(), nearestPointU);
                int[] nearestPointV = Particle.nearestROMSGridPoint((float) xy[0], (float) xy[1], this.getLonV(), this.getLatV(), elemLoc);
                int[] cRomsV = Particle.whichROMSElement((float) xy[0], (float) xy[1], this.getLonV(), this.getLatV(), nearestPointV);
                c = Math.min(cRomsU[0], cRomsV[0]);
            }
            if (c == -1) {
                inMesh = false;
            }
        }
        return inMesh;
    }
}
