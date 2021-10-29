/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.io.IOException;

import ucar.ma2.InvalidRangeException;

import java.io.File;
import java.util.List;
import java.util.stream.IntStream;

import org.apache.commons.io.FileUtils;
import org.apache.commons.io.filefilter.TrueFileFilter;
import org.apache.commons.io.filefilter.WildcardFileFilter;

/**
 * @author SA01TA
 */
public class HydroField {
    private float[][][] u;  // eastward velocity [time][sigLay][elem]
    private float[][][] v;  // northward velocity [time][sigLay][elem]
    private float[][][] w;  // upward velocity [time][sigLay][elem]
    private float[][][] s;  // salinity [time][sigLay][node]
    private float[][][] t;  // temperature [time][sigLay][node]
    private float[][] el;  // zeta = temporally varying sea surface height above geoid [time][node]
    private float[] h;  // h = bathymetric sea surface depth below geoid [node]
    private float[][][] diffVert;  // vertical diffusion coefficient [time][sigLay][node]

    /**
     * Default constructor for a single day's data
     *
     * @param filename
     * @param varNames
     * @param origin   the first element indices contained in a vector
     * @param shape    shape vector for U/V arrays
     * @param shapeST  an additional array shape vector for T/S arrays
     * @param type     type of hydro field (FVCOM/ROMS)
     */
    public HydroField(String filename, String[] varNames, int[] origin, int[] shape, int[] shapeST, String type, boolean readHydroVelocityOnly) {
        boolean readUpwardVelocity = false;
        for (String name: varNames) {
            if (name.equals("ww")) {
                readUpwardVelocity = true;
            }
        }

        System.out.println("Reading hydro file: " + filename);

        // Create additional shape matrices for the 2D variables (elevation)
        int[] origin2 = null;
        if (origin != null) {
            origin2 = new int[origin.length - 1];
            origin2[0] = origin[0];
            origin2[1] = origin[2];
        }
        int[] shape2 = null;
        int[] shapeST2 = null;
        if (shape != null) {
            shape2 = new int[shape.length - 1];
            shape2[0] = shape[0];
            shape2[1] = shape[2];
        }
        if (shapeST != null) {
            shapeST2 = new int[shapeST.length - 1];
            shapeST2[0] = shapeST[0];
            shapeST2[1] = shapeST[2];
        }


        if (type.equalsIgnoreCase("FVCOM") || type.equalsIgnoreCase("ROMS_TRI")) {
            u = IOUtils.readNetcdfFloat3D(filename, varNames[0], origin, shape);
            v = IOUtils.readNetcdfFloat3D(filename, varNames[1], origin, shape);
            if (!readHydroVelocityOnly) {
                s = IOUtils.readNetcdfFloat3D(filename, varNames[2], origin, shapeST);
                t = IOUtils.readNetcdfFloat3D(filename, varNames[3], origin, shapeST);
                el = IOUtils.readNetcdfFloat2D(filename, varNames[4], origin2, shapeST2);
                //diffVert = IOUtils.readNetcdfFloat3D(filename,varNames[5],origin2,shapeST2);
            }

            if (readUpwardVelocity) {
                w = IOUtils.readNetcdfFloat3D(filename, "ww", origin, shape);
            }

        } else if (type.equalsIgnoreCase("ROMS")) {
            u = IOUtils.readNetcdfFloat3D(filename, varNames[0], origin, shape);
            v = IOUtils.readNetcdfFloat3D(filename, varNames[1], origin, shape);
            if (!readHydroVelocityOnly) {
                float[][][] elTmp = IOUtils.readNetcdfFloat3D(filename, varNames[4], origin, shape);
                System.out.println("elTmp (" + elTmp.length + "," + elTmp[0].length + "," + elTmp[0][0].length + ")");
                el = new float[elTmp[0].length][elTmp[0][0].length];
                for (int xInd = 0; xInd < elTmp[0].length; xInd++) {
                    for (int yInd = 0; yInd < elTmp[0][0].length; yInd++) {
                        el[xInd][yInd] = elTmp[0][xInd][yInd];
                    }
                }
            }

//                   // Can add extra stuff here for reading 4D velocity fields 
//        // (necessary for reading ROMS output even though first dimension is on size 1)
////        if (origin != null)
////        {
////            if (origin.length == 4)
////            {
////                float uTmp[][][][] = IOUtils.readNetcdfFloat4D(filename,varNames[0],origin,shape);
////            }
////        }
////        else
////        {
//            u = IOUtils.readNetcdfFloat3D(filename,varNames[0],origin,shape);
//            v = IOUtils.readNetcdfFloat3D(filename,varNames[1],origin,shape);
//            el = IOUtils.readNetcdfFloat2D(filename,varNames[4],origin2,shape2);
        }

    }

    /**
     * Second constructor to read two files at once and combine into a single field
     *
     * @param filename1
     * @param filename2
     * @param varNames
     * @param origin
     * @param shape
     */
    public HydroField(String filename1, String filename2, String[] varNames, int[] origin, int[] shape, int shapeST[], String type, boolean readHydroVelocityOnly) {
        boolean readUpwardVelocity = false;
        for (String name: varNames) {
            if (name.equals("ww")) {
                readUpwardVelocity = true;
            }
        }
        if ((readUpwardVelocity && varNames.length != 6) || (!readUpwardVelocity && varNames.length != 5)) {
            System.err.println("\nIncorrect number of variable names for hydro data extraction with verticalDynamics = " + readUpwardVelocity);
        }

        // Create additional shape matrices for the 2D variables (elevation)
        // Check whether the requested dimensions are non-null and remove the
        // middle dimension (depth) if so.
        int[] origin2 = null;
        if (origin != null) {
            origin2 = new int[origin.length - 1];
            origin2[0] = origin[0];
            origin2[1] = origin[2];
        }
        int[] shape2 = null;
        int[] shapeST2 = null;
        if (shape != null) {
            shape2 = new int[shape.length - 1];
            shape2[0] = shape[0];
            shape2[1] = shape[2];
        }
        if (shapeST != null) {
            shapeST2 = new int[shapeST.length - 1];
            shapeST2[0] = shapeST[0];
            shapeST2[1] = shapeST[2];
        }

        System.out.println("Reading two hydro files and combining");
        float[][][] u1 = null, v1 = null, s1 = null, t1 = null, u2 = null, v2 = null, s2 = null, t2 = null, diffVert1 = null, diffVert2 = null, w1 = null, w2 = null;
        float[][] el1 = null, el2 = null;

        if (type.equalsIgnoreCase("FVCOM") || type.equalsIgnoreCase("ROMS_TRI")) {
            System.out.println("Reading hydro file: " + filename1);
            u1 = IOUtils.readNetcdfFloat3D(filename1, varNames[0], origin, shape);
            v1 = IOUtils.readNetcdfFloat3D(filename1, varNames[1], origin, shape);
            if (!readHydroVelocityOnly) {
                s1 = IOUtils.readNetcdfFloat3D(filename1, varNames[2], origin, shapeST);
                t1 = IOUtils.readNetcdfFloat3D(filename1, varNames[3], origin, shapeST);
                el1 = IOUtils.readNetcdfFloat2D(filename1, varNames[4], origin2, shapeST2);  // origin and shape need to lose a dimension (middle one) here
                //diffVert1 = IOUtils.readNetcdfFloat3D(filename1,varNames[5],origin2,shapeST2);
            }
            if (readUpwardVelocity) {
                w1 = IOUtils.readNetcdfFloat3D(filename1, "ww", origin, shape);
            }


            // When reading two files, we ALWAYS want to start from t=0 for the second one
            if (origin != null) {
                origin[0] = 0;
                origin2[0] = 0;
            }
            System.out.println("Reading hydro file: " + filename2);
            u2 = IOUtils.readNetcdfFloat3D(filename2, varNames[0], origin, shape);
            v2 = IOUtils.readNetcdfFloat3D(filename2, varNames[1], origin, shape);
            if (!readHydroVelocityOnly) {
                s2 = IOUtils.readNetcdfFloat3D(filename2, varNames[2], origin, shapeST);
                t2 = IOUtils.readNetcdfFloat3D(filename2, varNames[3], origin, shapeST);
                el2 = IOUtils.readNetcdfFloat2D(filename2, varNames[4], origin2, shapeST2);  // origin and shape need to lose a dimension here (depth)
                //diffVert2 = IOUtils.readNetcdfFloat3D(filename2,varNames[5],origin2,shapeST2);
            }
            if (readUpwardVelocity) {
                w2 = IOUtils.readNetcdfFloat3D(filename2, "ww", origin, shape);
            }
        } else if (type.equalsIgnoreCase("ROMS")) {
            System.out.println("Reading hydro file: " + filename1);
            u1 = IOUtils.readNetcdfFloat3D(filename1, varNames[0], origin, shape);
            v1 = IOUtils.readNetcdfFloat3D(filename1, varNames[1], origin, shape);
            float[][][] elTmp1 = IOUtils.readNetcdfFloat3D(filename1, varNames[4], origin, shape);
            el1 = new float[elTmp1[0].length][elTmp1[0][0].length];
            for (int xInd = 0; xInd < elTmp1[0].length; xInd++) {
                for (int yInd = 0; yInd < elTmp1[0][0].length; yInd++) {
                    el1[xInd][yInd] = elTmp1[0][xInd][yInd];
                }
            }

            System.out.println("Reading hydro file: " + filename2);
            u2 = IOUtils.readNetcdfFloat3D(filename2, varNames[0], origin, shape);
            v2 = IOUtils.readNetcdfFloat3D(filename2, varNames[1], origin, shape);
            float[][][] elTmp2 = IOUtils.readNetcdfFloat3D(filename2, varNames[4], origin, shape);
            el2 = new float[elTmp2[0].length][elTmp2[0][0].length];
            for (int xInd = 0; xInd < elTmp2[0].length; xInd++) {
                for (int yInd = 0; yInd < elTmp2[0][0].length; yInd++) {
                    el2[xInd][yInd] = elTmp2[0][xInd][yInd];
                }
            }
        } else {
            System.err.println("Could not read hydro file data - incorrect type");
        }

        //System.out.println("u1 length ("+u1.length+","+u1[0].length+","+u1[0][0].length+")");

        // These two are recorded at element centroids in FVCOM
        u = new float[u1.length + 1][u1[0].length][u1[0][0].length];
        v = new float[u1.length + 1][u1[0].length][u1[0][0].length];
        if (readUpwardVelocity) {
            w = new float[u1.length + 1][u1[0].length][u1[0][0].length];
        }
        // Next three quantities are recorded at nodes in FVCOM
        if (!readHydroVelocityOnly) {
            if (s1 != null) {
                s = new float[s1.length + 1][s1[0].length][s1[0][1].length];
            }
            if (t1 != null) {
                t = new float[t1.length + 1][t1[0].length][t1[0][0].length];
            }
            if (diffVert1 != null) {
                diffVert = new float[diffVert1.length + 1][diffVert1[0].length][diffVert1[0][0].length];
            }
            el = new float[el1.length + 1][el1[0].length];
        }

        double sumU = 0;

        // Use this to test zero velocity/diffusion cases
        boolean createTest = false;
        float testU = 0;
        float testV = (float) 0.1;
        float testW = (float) 0.1;
        if (createTest) {
            for (int tt = 0; tt < u.length; tt++) {
                for (int dep = 0; dep < u[0].length; dep++) {
                    for (int elem = 0; elem < u[0][0].length; elem++) {
                        u[tt][dep][elem] = testU;
                        v[tt][dep][elem] = testV;
                        if (readUpwardVelocity) {
                            w[tt][dep][elem] = testW;
                        }
                        sumU += u[tt][dep][elem];
                    }
                    if (!readHydroVelocityOnly) {
                        for (int node = 0; node < el1[0].length; node++) {

                            if (s2 != null) {
                                s[u.length - 1][dep][node] = 0;
                            }
                            if (t2 != null) {
                                t[u.length - 1][dep][node] = 0;
                            }
                            if (diffVert2 != null) {
                                diffVert[u.length - 1][dep][node] = 0;
                            }
                            if (dep == 0) {
                                el[u.length - 1][node] = 0;
                            }
                        }
                    }
                }
            }
        } else {

            float[] sumUDep = new float[u1[0].length];
            for (int tt = 0; tt < u1.length; tt++) {
                for (int dep = 0; dep < u1[0].length; dep++) {
                    for (int elem = 0; elem < u1[0][0].length; elem++) {
                        u[tt][dep][elem] = u1[tt][dep][elem];
                        v[tt][dep][elem] = v1[tt][dep][elem];
                        if (readUpwardVelocity) {
                            w[tt][dep][elem] = w1[tt][dep][elem];
                        }
                        sumU += u[tt][dep][elem];
                        sumUDep[dep] += u[tt][dep][elem];
                    }
                    if (!readHydroVelocityOnly) {
                        for (int node = 0; node < el1[0].length; node++) {
                            if (s1 != null) {
                                s[tt][dep][node] = s1[tt][dep][node];
                            }
                            if (t1 != null) {
                                t[tt][dep][node] = t1[tt][dep][node];
                            }
                            if (diffVert1 != null) {
                                diffVert[tt][dep][node] = t1[tt][dep][node];
                            }
                            if (dep == 0) {
                                el[tt][node] = el1[tt][node];
                            }
                        }
                    }
                }
            }
//                for (int dep = 0; dep < u1[0].length; dep++) {
//                    System.out.println("sumUDep["+dep+"] = "+sumUDep[dep]);
//                }
            for (int dep = 0; dep < u1[0].length; dep++) {
                for (int elem = 0; elem < u1[0][0].length; elem++) {
                    u[u.length - 1][dep][elem] = u2[0][dep][elem];
                    v[u.length - 1][dep][elem] = v2[0][dep][elem];
                    if (readUpwardVelocity) {
                        w[u.length - 1][dep][elem] = w2[0][dep][elem];
                    }
                    sumU += u[u.length - 1][dep][elem];
                }
                if (!readHydroVelocityOnly) {
                    for (int node = 0; node < el1[0].length; node++) {
                        if (s2 != null) {
                            s[u.length - 1][dep][node] = s2[0][dep][node];
                        }
                        if (t2 != null) {
                            t[u.length - 1][dep][node] = t2[0][dep][node];
                        }
                        if (diffVert2 != null) {
                            diffVert[u.length - 1][dep][node] = t2[0][dep][node];
                        }
                        if (dep == 0) {
                            el[u.length - 1][node] = el2[0][node];
                        }
                    }
                }
            }
        }

        System.out.println("Combined files to single arrays (e.g. velocity dimensions " + u.length + " " + u[1].length + " " + u[0][1].length + "; sum = " + sumU + ")");
//            sumIndex(u,0,0,true);
//            sumIndex(u,0,3,true);
//            sumIndex(u,0,6,true);
//            sumIndex(u,0,9,true);
//            sumIndex(u,0,12,true);
//            sumIndex(u,1,0,true);
//            sumIndex(u,1,5,true);
//            sumIndex(u,1,9,true);
//            sumIndex(u,2,0,true);
//            sumIndex(u,2,79243,true);
    }

    /**
     * Public getter methods for each internal field
     */
    public float[][][] getU() {
        return u;
    }

    public float[][][] getV() {
        return v;
    }

    public float[][][] getW() {
        return w;
    }

    public float[][][] getS() {
        return s;
    }

    public float[][][] getT() {
        return t;
    }

    public float[][] getEl() {
        return el;
    }

    public float[][][] getDiffVert() {
        return diffVert;
    }

    public float[][] getWaterDepthNodexy(Mesh mesh) {
        float[][] waterDepth = getEl();
        float[] bathDepth = mesh.getDepthNodexy();
        for (int hour=0; hour < waterDepth.length; hour++) {
            for (int node=0; node < waterDepth[hour].length; node++) {
                waterDepth[hour][node] += bathDepth[node];
            }
        }
        return waterDepth;
    }

    // TODO: Maybe. No measure of sea surface height above geoid at Uvnode
//    public float[][] getWaterDepthUvnode(Mesh mesh) {
//        float[][] waterDepth = getEl();
//        float[] bathDepth = mesh.getDepthUvnode();
//        for (int hour=0; hour < waterDepth.length; hour++) {
//            for (int node=0; node < waterDepth[hour].length; node++) {
//                waterDepth[hour][node] += bathDepth[node];
//            }
//        }
//        return waterDepth;
//    }

    /**
     * Compute a sum of an array over a particular dimension, for a particular
     * index of that dimension - to assist in validation of data read in
     *
     * @param var
     * @param dimension
     * @param index
     * @param print
     * @return
     */
    public float sumIndex(float[][][] var, int dimension, int index, boolean print) {
        float s = 0;

        int[] d0Range, d1Range, d2Range;

        // Set ranges and indices for sum
        if (dimension == 0) {
            d0Range = new int[1];
            d0Range[0] = index;
        } else {
            d0Range = IntStream.rangeClosed(0, var.length - 1).toArray();
        }
        if (dimension == 1) {
            d1Range = new int[1];
            d1Range[0] = index;
        } else {
            d1Range = IntStream.rangeClosed(0, var[1].length - 1).toArray();
        }
        if (dimension == 2) {
            d2Range = new int[1];
            d2Range[0] = index;
        } else {
            d2Range = IntStream.rangeClosed(0, var[0][1].length - 1).toArray();
        }

        for (int d0 = 0; d0 < d0Range.length; d0++) {
            for (int d1 = 0; d1 < d1Range.length; d1++) {
                for (int d2 = 0; d2 < d2Range.length; d2++) {
                    s += var[d0Range[d0]][d1Range[d1]][d2Range[d2]];
                }
            }
        }

        if (print) {
            System.out.println("Array sum (dimension:" + dimension + " index:" + index + ") = " + s);
        }

        return s;
    }

}
