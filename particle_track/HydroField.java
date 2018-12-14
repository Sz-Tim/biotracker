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
 *
 * @author SA01TA
 */
public class HydroField {
    private float[][][] u;
    private float[][][] v;
    private float[][][] s;
    private float[][][] t;
    private float[][] el;
    
    /**
     * Default constructor for a single day's data
     * 
     * @param filename 
     * @param varNames
     */
    public HydroField(String filename, String[] varNames)
    {
        System.out.println("Reading hydro file: "+filename);

        u = IOUtils.readNetcdfFloat3D(filename,varNames[0]);
        v = IOUtils.readNetcdfFloat3D(filename,varNames[1]);
        s = IOUtils.readNetcdfFloat3D(filename,varNames[2]);
        t = IOUtils.readNetcdfFloat3D(filename,varNames[3]);
        el = IOUtils.readNetcdfFloat2D(filename,varNames[4]); 
    }
    
    /**
     * Second constructor to read two files at once and combine into a single field
     * 
     * @param filename1
     * @param filename2 
     * @param varNames
     */
    public HydroField(String filename1, String filename2, String[] varNames)
    {
            if (varNames.length != 5)
            {
                System.err.println("Incorrect number of variable names for hydro data extraction");
            }
            System.out.println("Reading two hydro files and combining");
            System.out.println("Reading hydro file: "+filename1);
            float[][][] u1 = IOUtils.readNetcdfFloat3D(filename1,varNames[0]);
            float[][][] v1 = IOUtils.readNetcdfFloat3D(filename1,varNames[1]);
            float[][][] s1 = IOUtils.readNetcdfFloat3D(filename1,varNames[2]);
            float[][][] t1 = IOUtils.readNetcdfFloat3D(filename1,varNames[3]);
            float[][] el1 = IOUtils.readNetcdfFloat2D(filename1,varNames[4]); 
        
            System.out.println("Reading hydro file: "+filename2);
            float[][][] u2 = IOUtils.readNetcdfFloat3D(filename2,varNames[0]);
            float[][][] v2 = IOUtils.readNetcdfFloat3D(filename2,varNames[1]);
            float[][][] s2 = IOUtils.readNetcdfFloat3D(filename2,varNames[2]);
            float[][][] t2 = IOUtils.readNetcdfFloat3D(filename2,varNames[3]);
            float[][] el2 = IOUtils.readNetcdfFloat2D(filename2,varNames[4]);
            
            // These two are recorded at element centroids in FVCOM
            u = new float[u1.length+1][u1[1].length][u1[0][1].length];
            v = new float[u1.length+1][u1[1].length][u1[0][1].length];
            // Next three quantities are recorded at nodes in FVCOM
            s = new float[s1.length+1][s1[1].length][s1[0][1].length];
            t = new float[s1.length+1][s1[1].length][s1[0][1].length];
            el = new float[s1.length+1][s1[0][1].length];
            
            double sumU = 0;
            for (int tt = 0; tt < u1.length; tt++) {
                for (int dep = 0; dep < u1[1].length; dep++) {
                    for (int elem = 0; elem < u1[0][1].length; elem++) {
                        u[tt][dep][elem] = u1[tt][dep][elem];
                        v[tt][dep][elem] = v1[tt][dep][elem];
                        sumU += u[tt][dep][elem];
                    }
                    for (int node = 0; node < s1[0][1].length; node++) {
                        s[tt][dep][node] = s1[tt][dep][node];
                        t[tt][dep][node] = t1[tt][dep][node];
                        if (dep == 0)
                        {
                            el[tt][node] = el1[tt][node];
                        }
                    }
                }
            }
            for (int dep = 0; dep < u1[1].length; dep++) {
                for (int elem = 0; elem < u1[0][1].length; elem++) {
                    u[u.length-1][dep][elem] = u2[0][dep][elem];
                    v[u.length-1][dep][elem] = v2[0][dep][elem];
                    sumU+=u[u.length-1][dep][elem];
                }
                for (int node = 0; node < s1[0][1].length; node++) {
                    s[u.length-1][dep][node] = s2[0][dep][node];
                    t[u.length-1][dep][node] = t2[0][dep][node];
                    if (dep == 0)
                    {
                        el[u.length-1][node] = el2[0][node];
                    }
                }
            }
            
            System.out.println("Combined files to single arrays (e.g. velocity dimensions "+u.length+" "+u[1].length+" "+u[0][1].length+"; sum = "+sumU+")");
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
    public float[][][] getU()
    {
        return u;
    }
    public float[][][] getV()
    {
        return v;
    }
    public float[][][] getS()
    {
        return s;
    }
    public float[][][] getT()
    {
        return t;
    }
    public float[][] getEl()
    {
        return el;
    }
    
    /**
     * Compute a sum of an array over a particular dimension, for a particular 
     * index of that dimension
     * 
     * @param var
     * @param dimension
     * @param index
     * @param print
     * @return 
     */
    public float sumIndex(float[][][] var, int dimension, int index, boolean print)
    {
        float s = 0;
        
        int[] d0Range, d1Range, d2Range;
        
        // Set ranges and indices for sum
        if (dimension==0){
            d0Range = new int[1];
            d0Range[0] = index;
        } else {
            d0Range = IntStream.rangeClosed(0, var.length-1).toArray();
        }
        if (dimension==1){
            d1Range = new int[1];
            d1Range[0] = index;  
        } else {
            d1Range = IntStream.rangeClosed(0, var[1].length-1).toArray();
        }
        if (dimension==2){
            d2Range = new int[1];
            d2Range[0] = index;  
        } else {
            d2Range = IntStream.rangeClosed(0, var[0][1].length-1).toArray();
        }
        
        for (int d0 = 0; d0 < d0Range.length; d0++) {
            for (int d1 = 0; d1 < d1Range.length; d1++) {
                for (int d2 = 0; d2 < d2Range.length; d2++) {
                    s += var[d0Range[d0]][d1Range[d1]][d2Range[d2]];
                }
            }
        }
        
        if (print == true)
        {
            System.out.println("Array sum (dimension:"+dimension+" index:"+index+") = "+s);
        }
        
        return s;
    }
    
}
