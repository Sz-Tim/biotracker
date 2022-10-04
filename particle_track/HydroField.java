/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package particle_track;

import java.util.stream.IntStream;

/**
 * @author SA01TA
 */
public class HydroField {
    // names based on WeStCOMS naming
    private float[][][] u;  // eastward velocity [time][sigLay][elem]
    private float[][][] v;  // northward velocity [time][sigLay][elem]
    private float[][][] w;  // upward velocity [time][sigLay][elem]  --  Note: depth is positive, with 0 = surface; w is multiplied by -1 in Particle.velocityFromNearestList() so that positive velocity = downward
    private float[][][] s;  // salinity [time][sigLay][node]
    private float[][][] t;  // temperature [time][sigLay][node]
    private float[][] zeta;  // zeta = temporally varying sea surface height above geoid [time][node]
    private float[][][] k;  // turbulent eddy viscosity for momentum [time][sigLev][node] -- Km in FVCOM
    private float[][] light;  // shortwave radiation at surface [time][node]

    /**
     * Default constructor for a single day's data
     *
     * @param origin  the first element indices contained in a vector
     * @param shape   shape vector for U/V arrays
     * @param shapeST an additional array shape vector for T/S arrays
     * @param type    type of hydro field (FVCOM/ROMS)
     */
    public HydroField(String filename, String[] varNames, int[] origin, int[] shape, int[] shapeST, String type, RunProperties rp) {

        System.out.println("Reading hydro file: " + filename);

        // Create additional shape matrices for the 2D variables (elevation)
        int[] origin2 = null;
        if (origin != null) {
            origin2 = new int[origin.length - 1];
            origin2[0] = origin[0];
            origin2[1] = origin[2];
        }
        int[] shapeST2 = null;
        if (shapeST != null) {
            shapeST2 = new int[shapeST.length - 1];
            shapeST2[0] = shapeST[0];
            shapeST2[1] = shapeST[2];
        }

        if (type.equalsIgnoreCase("FVCOM") || type.equalsIgnoreCase("ROMS_TRI")) {
            u = IOUtils.readNetcdfFloat3D(filename, varNames[0], origin, shape);
            v = IOUtils.readNetcdfFloat3D(filename, varNames[1], origin, shape);
            w = IOUtils.readNetcdfFloat3D(filename, varNames[2], origin, shape);
            if (!rp.readHydroVelocityOnly) {
                s = IOUtils.readNetcdfFloat3D(filename, varNames[3], origin, shapeST);
                t = IOUtils.readNetcdfFloat3D(filename, varNames[4], origin, shapeST);
                light = IOUtils.readNetcdfFloat2D(filename, varNames[7], origin2, shapeST2);
                zeta = IOUtils.readNetcdfFloat2D(filename, varNames[5], origin2, shapeST2);
                if (rp.variableDiffusion) {
                    k = IOUtils.readNetcdfFloat3D(filename, varNames[6], origin2, shapeST2);
                }
            }
        } else if (type.equalsIgnoreCase("ROMS")) {
            u = IOUtils.readNetcdfFloat3D(filename, varNames[0], origin, shape);
            v = IOUtils.readNetcdfFloat3D(filename, varNames[1], origin, shape);
            if (!rp.readHydroVelocityOnly) {
                float[][][] elTmp = IOUtils.readNetcdfFloat3D(filename, varNames[5], origin, shape);
                System.out.println("elTmp (" + elTmp.length + "," + elTmp[0].length + "," + elTmp[0][0].length + ")");
                zeta = new float[elTmp[0].length][elTmp[0][0].length];
                for (int xInd = 0; xInd < elTmp[0].length; xInd++) {
                    System.arraycopy(elTmp[0][xInd], 0, zeta[xInd], 0, elTmp[0][0].length);
                }
            }
        }

    }

    /**
     * Second constructor to read two files at once and combine into a single field
     */
    public HydroField(String filename1, String filename2, String[] varNames, int[] origin, int[] shape, int[] shapeST, String type, RunProperties rp) {

        // Create additional shape matrices for the 2D variables (elevation)
        // Check whether the requested dimensions are non-null and remove the
        // middle dimension (depth) if so.
        int[] origin2 = null;
        if (origin != null) {
            origin2 = new int[origin.length - 1];
            origin2[0] = origin[0];
            origin2[1] = origin[2];
        }
        int[] shapeST2 = null;
        if (shapeST != null) {
            shapeST2 = new int[shapeST.length - 1];
            shapeST2[0] = shapeST[0];
            shapeST2[1] = shapeST[2];
        }

        System.out.println("Reading two hydro files and combining");
        float[][][] u1 = null, v1 = null, s1 = null, t1 = null, u2 = null, v2 = null, s2 = null, t2 = null, k1 = null, k2 = null, w1 = null, w2 = null;
        float[][] zeta1 = null, zeta2 = null, light1 = null, light2 = null;

        if (type.equalsIgnoreCase("FVCOM") || type.equalsIgnoreCase("ROMS_TRI")) {
            System.out.println("Reading hydro file: " + filename1);
            u1 = IOUtils.readNetcdfFloat3D(filename1, varNames[0], origin, shape);
            v1 = IOUtils.readNetcdfFloat3D(filename1, varNames[1], origin, shape);
            w1 = IOUtils.readNetcdfFloat3D(filename1, varNames[2], origin, shape);
            s1 = IOUtils.readNetcdfFloat3D(filename1, varNames[3], origin, shapeST);
            t1 = IOUtils.readNetcdfFloat3D(filename1, varNames[4], origin, shapeST);
            zeta1 = IOUtils.readNetcdfFloat2D(filename1, varNames[5], origin2, shapeST2);  // origin and shape need to lose a dimension (middle one) here
            light1 = IOUtils.readNetcdfFloat2D(filename1, varNames[7], origin2, shapeST2);
            if (rp.variableDiffusion) {
                k1 = IOUtils.readNetcdfFloat3D(filename1, varNames[6], origin, shapeST);
            }

            // When reading two files, we ALWAYS want to start from t=0 for the second one
            if (origin != null) {
                origin[0] = 0;
                origin2[0] = 0;
            }
            System.out.println("Reading hydro file: " + filename2);
            u2 = IOUtils.readNetcdfFloat3D(filename2, varNames[0], origin, shape);
            v2 = IOUtils.readNetcdfFloat3D(filename2, varNames[1], origin, shape);
            w2 = IOUtils.readNetcdfFloat3D(filename2, varNames[2], origin, shape);
            s2 = IOUtils.readNetcdfFloat3D(filename2, varNames[3], origin, shapeST);
            t2 = IOUtils.readNetcdfFloat3D(filename2, varNames[4], origin, shapeST);
            zeta2 = IOUtils.readNetcdfFloat2D(filename2, varNames[5], origin2, shapeST2);  // origin and shape need to lose a dimension here (depth)
            light2 = IOUtils.readNetcdfFloat2D(filename2, varNames[7], origin2, shapeST2);
            if (rp.variableDiffusion) {
                k2 = IOUtils.readNetcdfFloat3D(filename2, varNames[6], origin, shapeST);
            }
        } else if (type.equalsIgnoreCase("ROMS")) {
            System.out.println("Reading hydro file: " + filename1);
            u1 = IOUtils.readNetcdfFloat3D(filename1, varNames[0], origin, shape);
            v1 = IOUtils.readNetcdfFloat3D(filename1, varNames[1], origin, shape);
            float[][][] elTmp1 = IOUtils.readNetcdfFloat3D(filename1, varNames[5], origin, shape);
            zeta1 = new float[elTmp1[0].length][elTmp1[0][0].length];
            for (int xInd = 0; xInd < elTmp1[0].length; xInd++) {
                System.arraycopy(elTmp1[0][xInd], 0, zeta1[xInd], 0, elTmp1[0][0].length);
            }
            System.out.println("Reading hydro file: " + filename2);
            u2 = IOUtils.readNetcdfFloat3D(filename2, varNames[0], origin, shape);
            v2 = IOUtils.readNetcdfFloat3D(filename2, varNames[1], origin, shape);
            float[][][] elTmp2 = IOUtils.readNetcdfFloat3D(filename2, varNames[4], origin, shape);
            zeta2 = new float[elTmp2[0].length][elTmp2[0][0].length];
            for (int xInd = 0; xInd < elTmp2[0].length; xInd++) {
                System.arraycopy(elTmp2[0][xInd], 0, zeta2[xInd], 0, elTmp2[0][0].length);
            }
        } else {
            System.err.println("Could not read hydro file data - incorrect type");
        }

        // These are recorded at element centroids in FVCOM
        assert u1 != null;
        u = new float[u1.length + 1][u1[0].length][u1[0][0].length];
        v = new float[u1.length + 1][u1[0].length][u1[0][0].length];
        w = new float[u1.length + 1][u1[0].length][u1[0][0].length];
        // Next quantities are recorded at nodes in FVCOM
        if (s1 != null) {
            s = new float[s1.length + 1][s1[0].length][s1[0][1].length];
        }
        if (t1 != null) {
            t = new float[t1.length + 1][t1[0].length][t1[0][0].length];
        }
        zeta = new float[zeta1.length + 1][zeta1[0].length];
        assert light1 != null;
        light = new float[light1.length + 1][light1[0].length];
        if (rp.variableDiffusion) {
            if (k1 != null) {
                k = new float[k1.length + 1][k1[0].length][k1[0][0].length];
            }
        }

        double sumU = 0;

        // Use this to test zero velocity/diffusion cases
        boolean createTest = false;
        float testU = 0;
        float testV = (float) 0.1;
        float testW = (float) 0.1;
        //noinspection ConstantConditions
        if (createTest) {
            for (int hour = 0; hour < u.length; hour++) {
                for (int dep = 0; dep < u[0].length; dep++) {
                    for (int elem = 0; elem < u[0][0].length; elem++) {
                        u[hour][dep][elem] = testU;
                        v[hour][dep][elem] = testV;
                        w[hour][dep][elem] = testW;
                        sumU += u[hour][dep][elem];
                    }
                    if (!rp.readHydroVelocityOnly) {
                        for (int node = 0; node < zeta1[0].length; node++) {

                            if (s2 != null) {
                                s[u.length - 1][dep][node] = 0;
                            }
                            if (t2 != null) {
                                t[u.length - 1][dep][node] = 0;
                            }
                            if (k2 != null) {
                                k[u.length - 1][dep][node] = 0;
                                // easiest way to adjust for extra sigma level
                                if (dep == u[0].length - 1) {
                                    k[u.length - 1][dep + 1][node] = 0;
                                }
                            }
                            if (dep == 0) {
                                zeta[u.length - 1][node] = 0;
                                light[u.length - 1][node] = 0;
                            }
                        }
                    }
                }
            }
        } else {
            for (int hour = 0; hour < u1.length; hour++) {
                for (int dep = 0; dep < u1[0].length; dep++) {
                    for (int elem = 0; elem < u1[0][0].length; elem++) {
                        u[hour][dep][elem] = u1[hour][dep][elem];
                        v[hour][dep][elem] = v1[hour][dep][elem];
                        assert w1 != null;
                        w[hour][dep][elem] = w1[hour][dep][elem];
                        sumU += u[hour][dep][elem];
                    }
                    for (int node = 0; node < zeta1[0].length; node++) {
                        if (s1 != null) {
                            s[hour][dep][node] = s1[hour][dep][node];
                        }
                        if (t1 != null) {
                            t[hour][dep][node] = t1[hour][dep][node];
                        }
                        if (dep == 0) {
                            zeta[hour][node] = zeta1[hour][node];
                            light[hour][node] = light1[hour][node];
                        }
                        if (rp.variableDiffusion) {
                            if (k1 != null) {
                                k[hour][dep][node] = k1[hour][dep][node];
                                // easiest way to adjust for extra sigma level
                                if (dep == u1[0].length - 1) {
                                    k[hour][dep + 1][node] = k1[hour][dep + 1][node];
                                }
                            }
                        }
                    }
                }
            }
            for (int dep = 0; dep < u1[0].length; dep++) {
                for (int elem = 0; elem < u1[0][0].length; elem++) {
                    u[u.length - 1][dep][elem] = u2[0][dep][elem];
                    v[u.length - 1][dep][elem] = v2[0][dep][elem];
                    w[u.length - 1][dep][elem] = w2[0][dep][elem];
                    sumU += u[u.length - 1][dep][elem];
                }
                for (int node = 0; node < zeta1[0].length; node++) {
                    if (s2 != null) {
                        s[u.length - 1][dep][node] = s2[0][dep][node];
                    }
                    if (t2 != null) {
                        t[u.length - 1][dep][node] = t2[0][dep][node];
                    }
                    if (dep == 0) {
                        zeta[u.length - 1][node] = zeta2[0][node];
                        light[u.length - 1][node] = light2[0][node];
                    }
                    if (rp.variableDiffusion) {
                        if (k2 != null) {
                            k[u.length - 1][dep][node] = k2[0][dep][node];
                            // easiest way to adjust for extra sigma level
                            if (dep == u1[0].length - 1) {
                                k[u.length - 1][dep + 1][node] = k2[0][dep + 1][node];
                            }
                        }
                    }
                }
            }
        }

//        System.out.printf("Combined arrays: %d hours, %d layers, %d levels, %d elements, %d nodes\n",
//                u.length, u[1].length, k[1].length, u[0][1].length, k[0][1].length);
    }

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

    public float[][] getZeta() {
        return zeta;
    }

    public float[][][] getK() {
        return k;
    }

    public float[][] getLight() {
        return light;
    }

    public float[][] getWaterDepthNodexy(Mesh mesh) {
        float[][] waterDepth = getZeta();
        float[] bathDepth = mesh.getDepthNodexy();
        for (int hour = 0; hour < waterDepth.length; hour++) {
            for (int node = 0; node < waterDepth[hour].length; node++) {
                waterDepth[hour][node] += bathDepth[node];
            }
        }
        return waterDepth;
    }

    public double getWaterDepthUvnode(Mesh mesh, Particle part, int hour, RunProperties rp) {
        return getAvgFromTrinodes(mesh, part.getLocation(), part.getDepthLayer(), part.getElem(), hour, "zeta", rp) + mesh.getDepthUvnode()[part.getElem()];
    }

    public double getAvgFromTrinodes(Mesh mesh, double[] location, int depthLayer, int elem, int hour, String varName, RunProperties rp) {
        double sum = 0;
        double weightSum = 0;
        int[] particleNodes = new int[3];
        double[] dist = new double[3];
        double[] weights = new double[3];
        for (int i = 0; i < 3; i++) {
            particleNodes[i] = mesh.getTrinodes()[i][elem];
            dist[i] = Particle.distanceEuclid2(location[0], location[1], mesh.getNodexy()[0][particleNodes[i]], mesh.getNodexy()[1][particleNodes[i]], rp.coordRef);
            weights[i] = 1.0 / (dist[i] * dist[i]);
            weightSum += weights[i];
            if (varName.equalsIgnoreCase("temp")) {
                sum += getT()[hour][depthLayer][particleNodes[i]] * weights[i];
            } else if (varName.equalsIgnoreCase("salinity")) {
                sum += getS()[hour][depthLayer][particleNodes[i]] * weights[i];
            } else if (varName.equalsIgnoreCase("k")) {
                sum += getK()[hour][depthLayer][particleNodes[i]] * weights[i];
            } else if (varName.equalsIgnoreCase("zeta")) {
                sum += getZeta()[hour][particleNodes[i]] * weights[i];
            } else if (varName.equalsIgnoreCase("short_wave")) {
                sum += getLight()[hour][particleNodes[i]] * weights[i];
            }
        }

        return sum / weightSum;
    }

    public double getValueAtDepth(Mesh m, Particle part, double[] location, double depth, int hour, String varName, RunProperties rp, float[][] nearestSigmas) {

        float localDepth = m.getDepthUvnode()[part.getElem()]; // TODO: This ignores zeta -- use HydroField.getWaterDepthUvnode(), or just ignore
        float sigmaHeight = nearestSigmas[0][1] - nearestSigmas[1][1];
        double dzPartVsAbove = depth - nearestSigmas[1][1];

        double varBelow = getAvgFromTrinodes(m, location, (int) nearestSigmas[0][0], part.getElem(), hour, varName, rp);
        double varAbove = getAvgFromTrinodes(m, location, (int) nearestSigmas[1][0], part.getElem(), hour, varName, rp);
        if (sigmaHeight != 0) {
            return varAbove + dzPartVsAbove * (varBelow - varAbove) / sigmaHeight;
        } else {
            return varAbove;
        }
    }

    /**
     * Compute a sum of an array over a particular dimension, for a particular
     * index of that dimension - to assist in validation of data read in
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

        for (int i : d0Range) {
            for (int j : d1Range) {
                for (int k : d2Range) {
                    s += var[i][j][k];
                }
            }
        }

        if (print) {
            System.out.println("Array sum (dimension:" + dimension + " index:" + index + ") = " + s);
        }

        return s;
    }

}
