/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package biotracker;

import java.util.stream.IntStream;
import java.util.concurrent.*;


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
    private float[][][] k;  // turbulent eddy viscosity for scalars [time][sigLev][node] -- kh in WeStCOMS
    private float[][][] vh; // horizontal turbulent eddy viscosity for scalars [time][sigLay][node] -- viscofh in WeStCOMS
    private float[][] light;  // shortwave radiation at surface [time][node]

    public HydroField(float[][][] u, float[][][] v, float[][][] w, float[][][] s, float[][][] t, float[][] zeta, float[][][] k, float[][][] vh, float[][] light) {
        this.u = u;
        this.v = v;
        this.w = w;
        this.s = s;
        this.t = t;
        this.zeta = zeta;
        this.k = k;
        this.vh = vh;
        this.light = light;
    }


    /**
     * Default constructor for a single day's data
     *
     * @param origin  the first element indices contained in a vector
     * @param shape   shape vector for U/V arrays
     * @param shapeST an additional array shape vector for T/S arrays
     * @param type    type of hydro field (FVCOM/ROMS)
     */
    public HydroField(String filename, String[] varNames, int[] origin, int[] shape, int[] shapeST, String type, RunProperties rp) {
        System.out.println("Reading " + filename);
        final int MAX_RETRIES = 3;
        final long RETRY_DELAY_MS = 5000;

        // Create additional shape matrices for the 2D variables (elevation)
        final int[] origin2 = origin != null ? new int[]{origin[0], origin[2]} : null;
        final int[] shapeST2 = shapeST != null ? new int[]{shapeST[0], shapeST[2]} : null;

        if (rp.parallelThreadsHD > 1) {
            ExecutorService executor = Executors.newFixedThreadPool(rp.parallelThreadsHD);
            try {
                Future<float[][][]> uFuture = executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat3D(filename, varNames[0], origin, shape), MAX_RETRIES, RETRY_DELAY_MS));
                Future<float[][][]> vFuture = executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat3D(filename, varNames[1], origin, shape), MAX_RETRIES, RETRY_DELAY_MS));
                Future<float[][][]> wFuture = rp.fixDepth ? null : executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat3D(filename, varNames[2], origin, shape), MAX_RETRIES, RETRY_DELAY_MS));
                Future<float[][][]> sFuture = rp.needS ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat3D(filename, varNames[3], origin, shapeST), MAX_RETRIES, RETRY_DELAY_MS)) : null;
                Future<float[][][]> tFuture = rp.needT ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat3D(filename, varNames[4], origin, shapeST), MAX_RETRIES, RETRY_DELAY_MS)) : null;
                Future<float[][]> zetaFuture = rp.needZeta ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat2D(filename, varNames[5], origin2, shapeST2), MAX_RETRIES, RETRY_DELAY_MS)) : null;
                Future<float[][][]> kFuture = rp.needK ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat3D(filename, varNames[6], origin2, shapeST2), MAX_RETRIES, RETRY_DELAY_MS)) : null;
                Future<float[][][]> vhFuture = rp.needVh ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat3D(filename, varNames[7], origin, shapeST), MAX_RETRIES, RETRY_DELAY_MS)) : null;
                Future<float[][]> lightFuture = rp.needLight ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat2D(filename, varNames[8], origin2, shapeST2), MAX_RETRIES, RETRY_DELAY_MS)) : null;

                u = uFuture.get();
                v = vFuture.get();
                w = rp.fixDepth ? new float[u.length][u[0].length][u[0][0].length] : wFuture.get();
                s = rp.needS ? sFuture.get() : null;
                t = rp.needT ? tFuture.get() : null;
                zeta = rp.needZeta ? zetaFuture.get() : null;
                k = rp.needK ? kFuture.get() : null;
                vh = rp.needVh ? vhFuture.get() : null;
                light = rp.needLight ? lightFuture.get() : null;
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                executor.shutdown();
            }
        } else {
            try {
                u = IOUtils.readNetcdfFloat3D(filename, varNames[0], origin, shape);
                v = IOUtils.readNetcdfFloat3D(filename, varNames[1], origin, shape);
                w = rp.fixDepth ? new float[u.length][u[0].length][u[0][0].length] : IOUtils.readNetcdfFloat3D(filename, varNames[2], origin, shape);
                s = rp.needS ? IOUtils.readNetcdfFloat3D(filename, varNames[3], origin, shapeST) : null;
                t = rp.needT ? IOUtils.readNetcdfFloat3D(filename, varNames[4], origin, shapeST) : null;
                zeta = rp.needZeta ? IOUtils.readNetcdfFloat2D(filename, varNames[5], origin2, shapeST2) : null;
                k = rp.needK ? IOUtils.readNetcdfFloat3D(filename, varNames[6], origin2, shapeST2) : null;
                vh = rp.needVh ? IOUtils.readNetcdfFloat3D(filename, varNames[7], origin, shapeST) : null;
                light = rp.needLight ? IOUtils.readNetcdfFloat2D(filename, varNames[8], origin2, shapeST2) : null;
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    private static <T> T retry(Callable<T> task, int MAX_RETRIES, long RETRY_DELAY_MS) throws Exception {
        int attempts = 0;
        while (true) {
            try {
                return task.call();
            } catch (Exception e) {
                if (++attempts >= MAX_RETRIES) {
                    throw e;
                }
                System.err.println("Attempt " + attempts + " failed, retrying in " + RETRY_DELAY_MS + "ms...");
                Thread.sleep(RETRY_DELAY_MS);
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

        int nHour = 0, nDep = 0, nElem = 0, nNode = 0;

        System.out.println("Reading two hydro files and combining");
        float[][][] u1 = null, v1 = null, w1 = null, s1 = null, t1 = null, k1 = null, vh1 = null, u2 = null, v2 = null, w2 = null, s2 = null, t2 = null, k2 = null, vh2 = null;
        float[][] zeta1 = null, zeta2 = null, light1 = null, light2 = null;

        if (type.equalsIgnoreCase("FVCOM") || type.equalsIgnoreCase("ROMS_TRI")) {
            System.out.println("Reading hydro file: " + filename1);
            u1 = IOUtils.readNetcdfFloat3D(filename1, varNames[0], origin, shape);
            v1 = IOUtils.readNetcdfFloat3D(filename1, varNames[1], origin, shape);
            w1 = IOUtils.readNetcdfFloat3D(filename1, varNames[2], origin, shape);
            nHour = u1.length;
            nDep = u1[0].length;
            nElem = u1[0][0].length;
            if (rp.needS) {
                s1 = IOUtils.readNetcdfFloat3D(filename1, varNames[3], origin, shapeST);
                nNode = s1[0][0].length;
            }
            if (rp.needT) {
                t1 = IOUtils.readNetcdfFloat3D(filename1, varNames[4], origin, shapeST);
                nNode = t1[0][0].length;
            }
            if (rp.needZeta) {
                zeta1 = IOUtils.readNetcdfFloat2D(filename1, varNames[5], origin2, shapeST2);  // origin and shape rp.need to lose a dimension (middle one) here
                nNode = zeta1[0].length;
            }
            if (rp.needK) {
                k1 = IOUtils.readNetcdfFloat3D(filename1, varNames[6], origin, shapeST);
            }
            if (rp.needVh) {
                vh1 = IOUtils.readNetcdfFloat3D(filename1, varNames[7], origin, shapeST);
            }
            if (rp.needLight) {
                light1 = IOUtils.readNetcdfFloat2D(filename1, varNames[8], origin2, shapeST2);
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
            if (rp.needS) {
                s2 = IOUtils.readNetcdfFloat3D(filename2, varNames[3], origin, shapeST);
            }
            if (rp.needT) {
                t2 = IOUtils.readNetcdfFloat3D(filename2, varNames[4], origin, shapeST);
            }
            if (rp.needZeta) {
                zeta2 = IOUtils.readNetcdfFloat2D(filename2, varNames[5], origin2, shapeST2);  // origin and shape rp.need to lose a dimension here (depth)
            }
            if (rp.needK) {
                k2 = IOUtils.readNetcdfFloat3D(filename2, varNames[6], origin, shapeST);
            }
            if (rp.needVh) {
                vh2 = IOUtils.readNetcdfFloat3D(filename2, varNames[7], origin, shapeST);
            }
            if (rp.needLight) {
                light2 = IOUtils.readNetcdfFloat2D(filename2, varNames[8], origin2, shapeST2);
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
        u = new float[nHour + 1][nDep][nElem];
        v = new float[nHour + 1][nDep][nElem];
        w = new float[nHour + 1][nDep][nElem];
        // Next quantities are recorded at nodes in FVCOM
        if (s1 != null) {
            s = new float[nHour + 1][nDep][nNode];
        }
        if (t1 != null) {
            t = new float[nHour + 1][nDep][nNode];
        }
        if (zeta1 != null) {
            zeta = new float[nHour + 1][nNode];
        }
        if (k1 != null) {
            k = new float[nHour + 1][nDep+1][nNode];
        }
        if (vh1 != null) {
            vh = new float[nHour + 1][nDep][nNode];
        }
        if (light1 != null) {
            light = new float[nHour + 1][nNode];
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
                        for (int node = 0; node < nNode; node++) {

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
                            if (vh2 != null) {
                                vh[u.length - 1][dep][node] = 0;
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
//            for (int hour = 0; hour < nHour; hour++) {
//                for (int dep = 0; dep < nDep; dep++) {
//                    for (int elem = 0; elem < nElem; elem++) {
//                        u[hour][dep][elem] = u1[hour][dep][elem];
//                        v[hour][dep][elem] = v1[hour][dep][elem];
//                        w[hour][dep][elem] = w1[hour][dep][elem];
//                        sumU += u[hour][dep][elem];
//                    }
//                    for (int node = 0; node < nNode; node++) {
//                        if (s1 != null) {
//                            s[hour][dep][node] = s1[hour][dep][node];
//                        }
//                        if (t1 != null) {
//                            t[hour][dep][node] = t1[hour][dep][node];
//                        }
//                        if (dep == 0) {
//                            if (zeta != null) {
//                                zeta[hour][node] = zeta1[hour][node];
//                            }
//                            if (light != null) {
//                                light[hour][node] = light1[hour][node];
//                            }
//                        }
//                        if (k1 != null) {
//                            k[hour][dep][node] = k1[hour][dep][node];
//                            // easiest way to adjust for extra sigma level
//                            if (dep == nDep - 1) {
//                                k[hour][dep + 1][node] = k1[hour][dep + 1][node];
//                            }
//                        }
//                        if (vh1 != null) {
//                            vh[hour][dep][node] = vh1[hour][dep][node];
//                        }
//                    }
//                }
//            }
            for (int hour = 0; hour < nHour; hour++) {
                for (int dep = 0; dep < nDep; dep++) {
                    System.arraycopy(u1[hour][dep], 0, u1[hour][dep], 0, nElem);
                    System.arraycopy(v1[hour][dep], 0, v1[hour][dep], 0, nElem);
                    System.arraycopy(w1[hour][dep], 0, w1[hour][dep], 0, nElem);
                    for (int elem = 0; elem < nElem; elem++) {
                        sumU += u[hour][dep][elem];
                    }
                }
            }
            if (s1 != null) {
                for (int hour = 0; hour < nHour; hour++) {
                    for (int dep = 0; dep < nDep; dep++) {
                        System.arraycopy(s1[hour][dep], 0, s[hour][dep], 0, nNode);
                    }
                }
            }
            if (t1 != null) {
                for (int hour = 0; hour < nHour; hour++) {
                    for (int dep = 0; dep < nDep; dep++) {
                        System.arraycopy(t1[hour][dep], 0, t[hour][dep], 0, nNode);
                    }
                }
            }
            if (zeta != null) {
                for (int hour = 0; hour < nHour; hour++) {
                    System.arraycopy(zeta1[hour], 0, zeta[hour], 0, nNode);
                }
            }
            if (light != null) {
                for (int hour = 0; hour < nHour; hour++) {
                    System.arraycopy(light1[hour], 0, light[hour], 0, nNode);
                }
            }
            if (k1 != null) {
                for (int hour = 0; hour < nHour; hour++) {
                    for (int dep = 0; dep < nDep; dep++) {
                        System.arraycopy(k1[hour][dep], 0, k[hour][dep], 0, nNode);
                        if (dep == nDep - 1) {
                            System.arraycopy(k1[hour][dep + 1], 0, k[hour][dep + 1], 0, nNode);
                        }
                    }
                }
            }
            if (vh1 != null) {
                for (int hour = 0; hour < nHour; hour++) {
                    for (int dep = 0; dep < nDep; dep++) {
                        System.arraycopy(vh1[hour][dep], 0, vh[hour][dep], 0, nNode);
                    }
                }
            }
            for (int dep = 0; dep < nDep; dep++) {
                for (int elem = 0; elem < nElem; elem++) {
                    u[u.length - 1][dep][elem] = u2[0][dep][elem];
                    v[v.length - 1][dep][elem] = v2[0][dep][elem];
                    w[w.length - 1][dep][elem] = w2[0][dep][elem];
                    sumU += u[u.length - 1][dep][elem];
                }
                for (int node = 0; node < nNode; node++) {
                    if (s2 != null) {
                        s[s.length - 1][dep][node] = s2[0][dep][node];
                    }
                    if (t2 != null) {
                        t[t.length - 1][dep][node] = t2[0][dep][node];
                    }
                    if (dep == 0) {
                        zeta[zeta.length - 1][node] = zeta2[0][node];
                        if (light != null) {
                            light[light.length - 1][node] = light2[0][node];
                        }
                    }
                    if (rp.variableDhV) {
                        if (k2 != null) {
                            k[k.length - 1][dep][node] = k2[0][dep][node];
                            // easiest way to adjust for extra sigma level
                            if (dep == nDep - 1) {
                                k[k.length - 1][dep + 1][node] = k2[0][dep + 1][node];
                            }
                        }
                    }
                    if (vh2 != null) {
                        vh[vh.length - 1][dep][node] = vh2[0][dep][node];
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
    public float[][][] getVh() {
        return vh;
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
            dist[i] = Particle.distanceEuclid2(location[0], location[1], mesh.getNodexy()[0][particleNodes[i]], mesh.getNodexy()[1][particleNodes[i]], rp.coordOS);
            weights[i] = 1.0 / (dist[i] * dist[i]);
            weightSum += weights[i];
            switch (varName) {
                case "temp" -> sum += getT()[hour][depthLayer][particleNodes[i]] * weights[i];
                case "salinity" -> sum += getS()[hour][depthLayer][particleNodes[i]] * weights[i];
                case "short_wave" -> sum += getLight()[hour][particleNodes[i]] * weights[i];
                case "k" -> sum += getK()[hour][depthLayer][particleNodes[i]] * weights[i];
                case "vh" -> sum += getVh()[hour][depthLayer][particleNodes[i]] * weights[i];
                case "zeta" -> sum += getZeta()[hour][particleNodes[i]] * weights[i];
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


    public double getValueAtDepth(Mesh m, int elem, double[] location, double depth, int hour, String varName, RunProperties rp, float[][] nearestSigmas) {

        float localDepth = m.getDepthUvnode()[elem]; // TODO: This ignores zeta -- use HydroField.getWaterDepthUvnode(), or just ignore
        float sigmaHeight = nearestSigmas[0][1] - nearestSigmas[1][1];
        double dzPartVsAbove = depth - nearestSigmas[1][1];

        double varBelow = getAvgFromTrinodes(m, location, (int) nearestSigmas[0][0], elem, hour, varName, rp);
        double varAbove = getAvgFromTrinodes(m, location, (int) nearestSigmas[1][0], elem, hour, varName, rp);
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
