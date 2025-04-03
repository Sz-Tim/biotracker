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
    private float[][] Hsig;  // significant wave height [time][node]
    private float[][] Dir;  // mean wave direction [time][node]
    private float[][] Tm;  // mean absolute wave period [time][node]

    public HydroField(float[][][] u, float[][][] v, float[][][] w, float[][][] s, float[][][] t, float[][] zeta, float[][][] k, float[][][] vh, float[][] light, float[][] Hsig, float[][] Dir, float[][] Tm) {
        this.u = u;
        this.v = v;
        this.w = w;
        this.s = s;
        this.t = t;
        this.zeta = zeta;
        this.k = k;
        this.vh = vh;
        this.light = light;
        this.Hsig = Hsig;
        this.Dir = Dir;
        this.Tm = Tm;
    }


    /**
     * Default constructor for a single day's data
     *
     * @param origin  the first element indices contained in a vector
     * @param shape   shape vector for U/V arrays
     * @param shapeST an additional array shape vector for T/S arrays
     */
    public HydroField(String filename, String[] varNames, int[] origin, int[] shape, int[] shapeST, RunProperties rp) {
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

    public HydroField(String filename, String[] varNames, int[] origin, int[] shape, int[] shapeST, Boolean hfStokes, RunProperties rp) {
        System.out.println("Reading " + filename);
        final int MAX_RETRIES = 3;
        final long RETRY_DELAY_MS = 5000;

        // Create additional shape matrices for the 2D variables (elevation)
        final int[] origin2 = origin != null ? new int[]{origin[0], origin[2]} : null;
        final int[] shapeST2 = shapeST != null ? new int[]{shapeST[0], shapeST[2]} : null;

        if (rp.parallelThreadsHD > 1) {
            ExecutorService executor = Executors.newFixedThreadPool(rp.parallelThreadsHD);
            try {
                Future<float[][]> HsigFuture = rp.needStokes ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat2D(filename, varNames[9], origin2, shapeST2), MAX_RETRIES, RETRY_DELAY_MS)) : null;
                Future<float[][]> DirFuture = rp.needStokes ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat2D(filename, varNames[10], origin2, shapeST2), MAX_RETRIES, RETRY_DELAY_MS)) : null;
                Future<float[][]> TmFuture = rp.needStokes ? executor.submit(() -> retry(() -> IOUtils.readNetcdfFloat2D(filename, varNames[11], origin2, shapeST2), MAX_RETRIES, RETRY_DELAY_MS)) : null;

                Hsig = rp.needStokes ? HsigFuture.get() : null;
                Dir = rp.needStokes ? DirFuture.get() : null;
                Tm = rp.needStokes ? TmFuture.get() : null;
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                executor.shutdown();
            }
        } else {
            try {
                Hsig = rp.needStokes ? IOUtils.readNetcdfFloat2D(filename, varNames[9], origin2, shapeST2) : null;
                Dir = rp.needStokes ? IOUtils.readNetcdfFloat2D(filename, varNames[10], origin2, shapeST2) : null;
                Tm = rp.needStokes ? IOUtils.readNetcdfFloat2D(filename, varNames[11], origin2, shapeST2) : null;
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

    public float[][] getHsig() {
        return Hsig;
    }

    public float[][] getDir() {
        return Dir;
    }

    public float[][] getTm() { return Tm; }

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


    // Returns m/s
    public static double[] calcStokesDepthExponential(float Hsig, double DirRadians, float Tm, double depth) {
        double k = (2 * Math.PI) / (9.81 * Math.pow(Tm, 2));
        double stokesTransportMonochromatic = (2 * Math.PI / Tm) * Math.pow(Hsig, 2) / 16;
        double[] stokesSurfaceUV = new double[2];
        stokesSurfaceUV[0] = stokesTransportMonochromatic * Math.cos(DirRadians);
        stokesSurfaceUV[1] = stokesTransportMonochromatic * Math.sin(DirRadians);
        double stokesSurfaceSpeed = Math.sqrt(Math.pow(stokesSurfaceUV[0], 2) + Math.pow(stokesSurfaceUV[1], 2));
        double ke = stokesSurfaceSpeed / (2 * stokesTransportMonochromatic) / 3;
        double stokesSpeed = stokesSurfaceSpeed * Math.exp(-2*ke*depth) / (1+8*ke*depth); // check: depth needs to be negative!
        double[] stokesDepthUV = new double[2];
        stokesDepthUV[0] = stokesSpeed * stokesSurfaceUV[0] / stokesSpeed;
        stokesDepthUV[1] = stokesSpeed * stokesSurfaceUV[1] / stokesSpeed;
        return stokesDepthUV;
    }

    public double[] getAvgFromTrinodes(int elem, int[][] trinodes, double[] nodeDist, int hour, double depth) {
        double weightSum = 0;
        double [] stokesUV = {0,0};
        double stokesSumU = 0;
        double stokesSumV = 0;
        double[] weights = new double[3];
        for (int i=0; i<3; i++) {
            weights[i] = 1.0 / (nodeDist[i] * nodeDist[i]);
            weightSum += weights[i];
            stokesUV = calcStokesDepthExponential(getHsig()[hour][trinodes[i][elem]], Math.toRadians(getDir()[hour][trinodes[i][elem]]), getTm()[hour][trinodes[i][elem]], depth);
            stokesSumU += stokesUV[0] * weights[i];
            stokesSumV += stokesUV[1] * weights[i];
        }
        stokesUV[0] = stokesSumU / weightSum;
        stokesUV[1] = stokesSumV / weightSum;
        return stokesUV;
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
                case "Hsig" -> sum += getHsig()[hour][particleNodes[i]] * weights[i];
                case "Dir" -> sum += getDir()[hour][particleNodes[i]] * weights[i];
                case "Tm" -> sum += getTm()[hour][particleNodes[i]] * weights[i];
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
